"""Independent conduction checks and real-output failure injection."""

import copy
import json
import os
from pathlib import Path
import shutil
import tempfile
import unittest

from open_cae_evidence.adapters.calculix import CalculiXAdapter
from open_cae_evidence.evidence import manifest, write_json
from open_cae_evidence.results import validate_contract
from open_cae_evidence.runner import replay, run_benchmark
from open_cae_evidence.specification import load_problem
from open_cae_evidence.thermal_input import make_deck, audit_deck
from open_cae_evidence.thermal_problem import analytical, validate
from open_cae_evidence.verification import ThermalVerifier

CONFIG = Path(__file__).resolve().parents[1] / "benchmarks/slab_steady_heat.json"


class ThermalInputTests(unittest.TestCase):
    def setUp(self):
        self.config = json.loads(CONFIG.read_text())

    def test_invalid_conductivity_units_and_values(self):
        for value in (0, -1, float("nan"), float("inf"), True):
            config = copy.deepcopy(self.config)
            config["material"]["conductivity"] = value
            with self.subTest(value=value), self.assertRaises(ValueError):
                validate(config)
        for value in ("W/K", "W/(mm*K)", "W/m^2"):
            config = copy.deepcopy(self.config)
            config["units"]["conductivity"] = value
            with self.subTest(unit=value), self.assertRaises(ValueError):
                validate(config)

    def test_wrong_or_missing_physical_temperatures(self):
        for changes in (
            {"x_min_temperature": 300, "x_max_temperature": 400},
            {"x_min_temperature": 300},
            {"x_max_temperature": 0},
            {"x_min_temperature": float("nan")},
            {"other_faces": "convection"},
        ):
            config = copy.deepcopy(self.config)
            config["boundary"].update(changes)
            with self.subTest(changes=changes), self.assertRaises(ValueError):
                validate(config)
        config = copy.deepcopy(self.config)
        del config["boundary"]["x_max_temperature"]
        with self.assertRaises(ValueError):
            validate(config)
        config = copy.deepcopy(self.config)
        config["internal_heat_generation"] = 1
        with self.assertRaises(ValueError):
            validate(config)

    def test_serialized_material_boundary_and_mesh_corruption(self):
        deck, _, _, _ = make_deck(self.config, [2, 2, 2])
        corruptions = [
            deck.replace("*CONDUCTIVITY\n10", "*CONDUCTIVITY\n-10"),
            deck.replace("*CONDUCTIVITY\n10", "*CONDUCTIVITY\n10000"),
            deck.replace("HOT, 11, 11, 400", "HOT, 11, 11, 300"),
            deck.replace("COLD, 11, 11, 300\n", ""),
            deck.replace("HOT, 11, 11, 400", "ALL, 11, 11, 400"),
            deck.replace("*HEAT TRANSFER, STEADY STATE", "*HEAT TRANSFER"),
            deck.replace("*END STEP", "*CFLUX\n1, 11, 1\n*END STEP"),
        ]
        rows = deck.splitlines()
        start = rows.index("*ELEMENT, TYPE=DC3D8, ELSET=SLAB") + 1
        row = rows[start].split(",")
        row[1], row[2] = row[2], row[1]
        rows[start] = ",".join(row)
        corruptions.append("\n".join(rows))
        for corrupt in corruptions:
            with self.subTest(length=len(corrupt)), self.assertRaises(ValueError):
                audit_deck(corrupt, self.config, [2, 2, 2])

    def test_fourier_reference_dimensions_and_scaling(self):
        ref = analytical(self.config, 0.05)
        self.assertEqual(ref["temperature_K"], 350)
        self.assertEqual(ref["gradient_K_m"], [-1000, 0, 0])
        self.assertEqual(ref["heat_flux_W_m2"], [10000, 0, 0])
        self.assertEqual(ref["heat_flow_W"], 600)
        config = copy.deepcopy(self.config)
        config["material"]["conductivity"] *= 2
        self.assertEqual(analytical(config, 0.05)["temperature_K"], 350)
        self.assertEqual(analytical(config, 0.05)["heat_flow_W"], 1200)


class RealThermalTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.config = json.loads(CONFIG.read_text())
        cls.problem = load_problem(cls.config)
        cls.solver = shutil.which(os.environ.get("CCX", "ccx"))
        if not cls.solver:
            raise RuntimeError("Real CalculiX required for all thermal tests")
        root = Path(os.environ.get("TEST_EVIDENCE_DIR", "evidence/tests"))
        root.mkdir(parents=True, exist_ok=True)
        cls.path = Path(tempfile.mkdtemp(prefix="thermal-", dir=root)).resolve()
        cls.archive = cls.path / "archive"
        if run_benchmark(cls.config, cls.archive, cls.solver) != 0:
            raise RuntimeError("Thermal vertical slice failed")
        cls.contract = json.loads((cls.archive / "normalized-results.json").read_text())
        cls.normalized = cls.contract["runs"][0]
        cls.adapter = CalculiXAdapter()
        cls.verifier = ThermalVerifier()
        cls.prepared = cls.adapter.prepare(cls.problem, [2, 2, 2])
        cls.version = cls.normalized["solver"]["version"]

    def test_three_meshes_and_conservation(self):
        validate_contract(self.contract)
        self.assertEqual(len(self.contract["runs"]), 3)
        for run in self.contract["runs"]:
            result = self.verifier.verify(self.problem, run)
            self.assertEqual(result["status"], "passed")
            self.assertAlmostEqual(
                result["metrics"]["surface_outward_heat_flow_W"]["x_max"], 600, places=7
            )
            self.assertLess(result["metrics"]["normalized_energy_imbalance"], 2e-6)
        study = self.contract["verification"]["mesh_study"]
        self.assertEqual(study["kind"], "mesh_consistency")
        self.assertIsNone(study["demonstrated_convergence_order"])
        replay(self.archive)

    def test_real_solver_wrong_temperature_assignment_fails_physics(self):
        path = self.path / "swapped-boundaries"
        path.mkdir()
        deck = (
            self.prepared.files["slab.inp"]
            .replace("HOT, 11, 11, 400", "HOT, 11, 11, 300")
            .replace("COLD, 11, 11, 300", "COLD, 11, 11, 400")
        )
        with self.assertRaises(ValueError):
            audit_deck(deck, self.config, [2, 2, 2])
        (path / "slab.inp").write_text(deck)
        version = self.adapter.execute(self.prepared, path, self.solver)
        normalized = self.adapter.extract(self.prepared, path, version)
        result = self.verifier.verify(self.problem, normalized)
        write_json(path / "negative-verification.json", result)
        self.assertEqual(result["status"], "failed")
        self.assertFalse(result["checks"]["nodal_temperature"])
        self.assertFalse(result["checks"]["solver_flux_direction_magnitude"])
        # Conservation alone can hold for a physically wrong boundary assignment.
        self.assertTrue(result["checks"]["energy_conservation"])

    def test_corrupted_real_temperature_and_flux_sign_cannot_pass(self):
        for field in ("temperature", "heat_flux", "nodal_heat_rate"):
            changed = copy.deepcopy(self.normalized)
            if field == "temperature":
                changed["quantities"][field]["values"][0]["value"] += 50
            elif field == "heat_flux":
                for r in changed["quantities"][field]["values"]:
                    r["value"] = [-v for v in r["value"]]
            else:
                changed["quantities"][field]["values"][0]["value"] *= 2
            with self.subTest(field=field):
                result = self.verifier.verify(self.problem, changed)
                write_json(self.path / f"negative-{field}.json", result)
                self.assertEqual(result["status"], "failed")

    def test_incomplete_or_wrong_unit_normalized_output_is_rejected(self):
        for field in ("temperature", "heat_flux", "nodal_heat_rate"):
            changed = copy.deepcopy(self.normalized)
            changed["quantities"][field]["values"].pop()
            with self.subTest(field=field), self.assertRaises(ValueError):
                self.verifier.verify(self.problem, changed)
        changed = copy.deepcopy(self.normalized)
        changed["quantities"]["heat_flux"]["unit"] = "W"
        with self.assertRaises(ValueError):
            self.verifier.verify(self.problem, changed)

    def test_corrupted_dat_and_frd_and_missing_boundary_output(self):
        origin = self.archive / "mesh-2x2x2"
        path = self.path / "corrupt-raw"
        shutil.copytree(origin, path)
        dat = path / "slab.dat"
        original = dat.read_text()
        dat.write_text(original.replace("4.000000E+02", "9.000000E+02", 1))
        with self.assertRaises(ValueError):
            self.adapter.extract(self.prepared, path, self.version)
        dat.write_text(original[: len(original) // 2])
        with self.assertRaises(ValueError):
            self.adapter.extract(self.prepared, path, self.version)
        dat.write_text(original)
        frd = path / "slab.frd"
        frd.write_text(frd.read_text().replace(" -5  T ", " -5  U ", 1))
        with self.assertRaises(ValueError):
            self.adapter.extract(self.prepared, path, self.version)

    def test_missing_and_rehashed_inconsistent_evidence(self):
        path = self.path / "incomplete-archive"
        shutil.copytree(self.archive, path)
        (path / "mesh-2x2x2/slab.dat").unlink()
        with self.assertRaises(ValueError):
            replay(path)
        write_json(path / "manifest.json", manifest(path))
        with self.assertRaises((ValueError, OSError)):
            replay(path)
        path = self.path / "rehashed-inconsistent"
        shutil.copytree(self.archive, path)
        normalized = json.loads((path / "normalized-results.json").read_text())
        normalized["runs"][0]["quantities"]["temperature"]["values"][0]["value"] += 50
        write_json(path / "normalized-results.json", normalized)
        write_json(path / "manifest.json", manifest(path))
        with self.assertRaisesRegex(ValueError, "differs from raw evidence"):
            replay(path)

    def test_real_conductivity_scaling(self):
        config = copy.deepcopy(self.config)
        config["material"]["conductivity"] *= 2
        path = self.path / "double-conductivity"
        from open_cae_evidence.runner import run_mesh

        run = run_mesh(config, [2, 2, 2], path, self.solver)
        self.assertTrue(run["passed"])
        self.assertAlmostEqual(
            run["thermal"]["surface_outward_heat_flow_W"]["x_max"], 1200, places=7
        )
        self.assertLess(run["thermal"]["temperature_max_error_K"], 1e-3)

    def test_contract_version_units_and_false_outcomes_fail(self):
        for name in ("version", "units", "false_pass", "empty_mesh_study"):
            changed = copy.deepcopy(self.contract)
            if name == "version":
                changed["schema_version"] = "99.0.0"
            elif name == "units":
                changed["runs"][0]["quantities"]["temperature"]["unit"] = "degC"
            elif name == "false_pass":
                changed["runs"][0]["verification"]["checks"][
                    "energy_conservation"
                ] = False
            else:
                changed["runs"] = []
            with self.subTest(name=name), self.assertRaises(ValueError):
                validate_contract(changed)
