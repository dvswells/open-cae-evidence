"""Independent analytic checks and adversarial tests against real sourced solves."""

import copy
import json
import math
import os
from pathlib import Path
import shutil
import tempfile
import unittest
import xml.etree.ElementTree as ET

from open_cae_evidence import q1, source_reference
from open_cae_evidence.adapters.calculix import CalculiXAdapter
from open_cae_evidence.adapters.elmer import ElmerAdapter
from open_cae_evidence.adapters import elmer_input
from open_cae_evidence.evidence import manifest, write_json
from open_cae_evidence.results import validate_contract
from open_cae_evidence.runner import run_benchmark, replay
from open_cae_evidence.source_comparison import main as compare_main, compare_contracts
from open_cae_evidence.source_verification import study
from open_cae_evidence.specification import load_problem
from open_cae_evidence.thermal_input import audit_deck
from open_cae_evidence.verification import for_problem

CONFIG = Path("benchmarks/slab_volumetric_source.json")


class ReferenceTests(unittest.TestCase):
    def setUp(self):
        self.config = json.loads(CONFIG.read_text())

    def test_exact_temperature_flux_and_global_balance(self):
        for x, T, q in (
            (0.0, 400.0, 5000.0),
            (0.05, 362.5, 10000.0),
            (0.1, 300.0, 15000.0),
        ):
            a = source_reference.analytical(self.config, x)
            self.assertAlmostEqual(a["temperature_K"], T)
            self.assertAlmostEqual(a["heat_flux_W_m2"][0], q)
        b = source_reference.balance(self.config)
        self.assertAlmostEqual(b["outward_heat_flow_W"]["x_min"], -300)
        self.assertAlmostEqual(b["outward_heat_flow_W"]["x_max"], 900)
        self.assertAlmostEqual(
            sum(b["outward_heat_flow_W"].values()), b["generation_W"]
        )
        self.assertAlmostEqual(b["generation_W"], 600)

    def test_governing_equation_by_independent_finite_differences(self):
        x, dx = 0.037, 1e-4
        T = lambda x: source_reference.analytical(self.config, x)["temperature_K"]
        q = lambda x: source_reference.analytical(self.config, x)["heat_flux_W_m2"][0]
        second = (T(x + dx) - 2 * T(x) + T(x - dx)) / (dx * dx)
        self.assertAlmostEqual(-10 * second, 100000, delta=1e-3)
        self.assertAlmostEqual((q(x + dx) - q(x - dx)) / (2 * dx), 100000, delta=1e-6)
        self.assertAlmostEqual(
            -10 * (T(x + dx) - T(x - dx)) / (2 * dx), q(x), delta=1e-6
        )

    def test_closed_form_interpolant_norms_independent_of_solver(self):
        a = source_reference.interpolation_norms(self.config, 0.05)
        self.assertAlmostEqual(a["temperature_L2_K_m3_2"], math.sqrt(0.03125))
        self.assertAlmostEqual(a["gradient_H1_K_m1_2"], math.sqrt(125))
        self.assertAlmostEqual(a["raw_flux_L2_W_m_minus1_2"], math.sqrt(12500))
        self.assertAlmostEqual(sum(w for _, w in q1.G3), 8)

    def test_invalid_units_sources_and_regions(self):
        for value in (0, -1, float("nan"), True):
            c = copy.deepcopy(self.config)
            c["internal_heat_generation"] = value
            with self.assertRaises(ValueError):
                load_problem(c)
        for unit in ("W/m^2", "W/kg", "W/mm^3"):
            c = copy.deepcopy(self.config)
            c["units"]["volumetric_heat_generation"] = unit
            with self.assertRaises(ValueError):
                load_problem(c)
        c = copy.deepcopy(self.config)
        c["source_region"] = "x_max"
        with self.assertRaises(ValueError):
            load_problem(c)

    def test_source_free_definition_is_unchanged_and_rejects_sources(self):
        c = json.loads(Path("benchmarks/slab_steady_heat.json").read_text())
        self.assertAlmostEqual(
            source_reference.analytical(self.config, 0.05)["temperature_K"] - 350, 12.5
        )
        load_problem(c)
        c["internal_heat_generation"] = 100000
        with self.assertRaises(ValueError):
            load_problem(c)

    def test_precision_limited_orders_are_not_invented(self):
        rows = []
        for h in (0.05, 0.025, 0.0125):
            m = {"h_x_m": h, "precision_floors": {}}
            for key in (
                "temperature_L2_K_m3_2",
                "gradient_H1_K_m1_2",
                "raw_flux_L2_W_m_minus1_2",
                "recovered_flux_L2_W_m_minus1_2",
                "raw_boundary_heat_flow_error_W",
            ):
                m[key] = 1e-15
                m["precision_floors"][key] = 1e-8
            rows.append({"passed": True, "thermal": m})
        r = study(rows)
        self.assertFalse(r["passed"])
        self.assertIsNone(r["observed_rates"][0]["orders"]["temperature_L2_K_m3_2"])


class SourceSolverTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temp = tempfile.TemporaryDirectory()
        cls.root = Path(os.environ.get("CAE_SOURCE_TEST_EVIDENCE", cls.temp.name))
        cls.root.mkdir(parents=True, exist_ok=True)
        cls.config = json.loads(CONFIG.read_text())
        cls.config["meshes"] = cls.config["meshes"][:3]
        cls.problem = load_problem(cls.config)
        cls.verifier = for_problem(cls.problem)
        cls.adapters = {"CalculiX": CalculiXAdapter(), "Elmer": ElmerAdapter()}
        cls.exes = {
            "CalculiX": shutil.which(os.environ.get("CCX", "ccx")),
            "Elmer": shutil.which(os.environ.get("ELMER_SOLVER", "ElmerSolver")),
        }
        cls.archives = {}
        cls.contracts = {}
        cls.prepared = {}
        for name, a in cls.adapters.items():
            if cls.exes[name] is None:
                raise RuntimeError("Real " + name + " is required")
            path = cls.root / name
            cls.archives[name] = path
            if run_benchmark(cls.config, path, cls.exes[name], a):
                raise RuntimeError("Real sourced " + name + " benchmark failed")
            cls.contracts[name] = json.loads(
                (path / "normalized-results.json").read_text()
            )
            cls.prepared[name] = a.prepare(cls.problem, [2, 2, 2])

    @classmethod
    def tearDownClass(cls):
        cls.temp.cleanup()

    def directory(self):
        d = self.root / self._testMethodName
        d.mkdir()
        return d

    def test_both_real_solvers_nonzero_nonnodal_errors_and_convergence(self):
        for c in self.contracts.values():
            validate_contract(c)
            self.assertEqual(c["schema_version"], "1.0.0")
            m = c["runs"][0]["verification"]["metrics"]
            self.assertEqual(m["temperature_max_error_K"], 0)
            self.assertGreater(m["temperature_L2_K_m3_2"], 0.17)
            self.assertGreater(m["gradient_H1_K_m1_2"], 11)
            self.assertTrue(c["verification"]["mesh_study"]["passed"])
            self.assertAlmostEqual(
                c["verification"]["mesh_study"]["observed_rates"][0]["orders"][
                    "temperature_L2_K_m3_2"
                ],
                2,
            )

    def test_center_only_gradient_error_would_falsely_be_zero(self):
        from open_cae_evidence.results import mesh_maps, thermal_fields

        run = self.contracts["Elmer"]["runs"][0]
        nodes, elements = mesh_maps(run)
        data = thermal_fields(run)
        center_error = 0.0
        for conn in elements.values():
            start, span = q1.cell(nodes, conn)
            _, grad = q1.evaluate(
                tuple(data["temperature"][n] for n in conn), span, (0.0, 0.0, 0.0)
            )
            exact = source_reference.analytical(self.config, start[0] + span[0] / 2)[
                "gradient_K_m"
            ]
            center_error += sum((x - y) ** 2 for x, y in zip(grad, exact)) * math.prod(
                span
            )
        self.assertLess(math.sqrt(center_error), 1e-10)
        self.assertGreater(run["verification"]["metrics"]["gradient_H1_K_m1_2"], 11)

    def test_boundary_reactions_and_raw_flux_estimates_are_distinct(self):
        for c in self.contracts.values():
            m = c["runs"][0]["verification"]["metrics"]
            self.assertAlmostEqual(m["actual_source_integral_W"], 600)
            self.assertAlmostEqual(m["reaction_net_outward_W"], 600)
            self.assertAlmostEqual(m["raw_net_outward_W"], 300)
            self.assertAlmostEqual(m["surface_outward_heat_flow_W"]["x_max"], 750)
            self.assertAlmostEqual(
                m["boundary_reaction_outward_heat_flow_W"]["x_max"], 900
            )

    def test_calculix_native_RFL_is_not_a_boundary_reaction(self):
        r = self.contracts["CalculiX"]["runs"][0]
        self.assertAlmostEqual(
            sum(
                v["value"]
                for v in r["quantities"]["native_total_nodal_heat_rate"]["values"]
            ),
            0,
        )
        self.assertAlmostEqual(
            sum(v["value"] for v in r["quantities"]["nodal_heat_rate"]["values"]), -600
        )
        self.assertEqual(
            r["quantities"]["nodal_heat_rate"]["source"]["kind"],
            "postprocessed_from_solver",
        )

    def test_actual_inputs_reject_missing_wrong_sign_magnitude_region(self):
        for name, prep in self.prepared.items():
            if name == "CalculiX":
                deck = prep.files["slab.inp"]
                changed = [deck.replace("*DFLUX\nSLAB, BF, 100000\n", "")]
                changed += [
                    deck.replace("SLAB, BF, 100000", s)
                    for s in (
                        "SLAB, BF, 50000",
                        "SLAB, BF, -100000",
                        "SLAB, S1, 100000",
                        "HOT, BF, 100000",
                    )
                ]
                for text in changed:
                    with self.assertRaises(ValueError):
                        audit_deck(text, self.config, [2, 2, 2])
            else:
                for old, new in (
                    (
                        "Volumetric Heat Source = Real 100000",
                        "Volumetric Heat Source = Real 50000",
                    ),
                    (
                        "Volumetric Heat Source = Real 100000",
                        "Heat Source = Real 100000",
                    ),
                    (
                        "Volumetric Heat Source = Real 100000",
                        "Volumetric Heat Source = Real -100000",
                    ),
                    ("  Body Force = 1\n", ""),
                    ("Target Bodies(1) = 1", "Target Bodies(1) = 2"),
                ):
                    files = dict(prep.files)
                    files["slab.sif"] = files["slab.sif"].replace(old, new)
                    with self.assertRaises(ValueError):
                        elmer_input.audit(files, self.config, [2, 2, 2])

    def test_missing_or_wrong_normalized_source_coverage_and_units(self):
        for name in self.contracts:
            for change in ("missing", "coverage", "unit"):
                r = copy.deepcopy(self.contracts[name]["runs"][0])
                if change == "missing":
                    del r["raw_output_verification"]["volumetric_source"]
                elif change == "coverage":
                    r["raw_output_verification"]["volumetric_source"][
                        "element_values"
                    ].pop()
                else:
                    r["raw_output_verification"]["volumetric_source"]["unit"] = "W/m^2"
                with self.assertRaises(ValueError):
                    self.verifier.verify(self.problem, r)

    def test_corrupt_actual_native_temperature_output(self):
        prep = self.prepared["Elmer"]
        d = self.directory()
        shutil.copytree(self.archives["Elmer"] / "mesh-2x2x2", d / "mesh")
        d = d / "mesh"
        p = next(d.glob("slab*.vtu"))
        xml = ET.fromstring(p.read_text())
        a = xml.find('.//DataArray[@Name="temperature"]')
        a.text = " ".join(str(float(v) + 5) for v in a.text.split())
        p.write_text(ET.tostring(xml, encoding="unicode"))
        r = self.adapters["Elmer"].extract(prep, d, "9.0")
        result = self.verifier.verify(self.problem, r)
        write_json(d / "failure-verification.json", result)
        self.assertEqual(result["status"], "failed")
        self.assertFalse(result["checks"]["nodal_temperature"])

    def test_reversed_real_native_recovered_flux(self):
        prep = self.prepared["Elmer"]
        d = self.directory()
        shutil.copytree(self.archives["Elmer"] / "mesh-2x2x2", d / "mesh")
        d = d / "mesh"
        p = next(d.glob("slab*.vtu"))
        xml = ET.fromstring(p.read_text())
        a = xml.find('.//DataArray[@Name="temperature flux"]')
        a.text = " ".join(str(-float(v)) for v in a.text.split())
        p.write_text(ET.tostring(xml, encoding="unicode"))
        r = self.adapters["Elmer"].extract(prep, d, "9.0")
        result = self.verifier.verify(self.problem, r)
        write_json(d / "failure-verification.json", result)
        self.assertFalse(result["checks"]["native_flux_consistency"])
        self.assertEqual(result["status"], "failed")

    def test_incomplete_native_output(self):
        prep = self.prepared["Elmer"]
        d = self.directory()
        shutil.copytree(self.archives["Elmer"] / "mesh-2x2x2", d / "mesh")
        d = d / "mesh"
        p = next(d.glob("slab*.vtu"))
        xml = ET.fromstring(p.read_text())
        a = xml.find('.//DataArray[@Name="temperature"]')
        a.text = " ".join(a.text.split()[:-1])
        p.write_text(ET.tostring(xml, encoding="unicode"))
        with self.assertRaises(ValueError):
            self.adapters["Elmer"].extract(prep, d, "9.0")

    def test_wrong_coordinate_transforms_in_actual_inputs(self):
        prep = self.prepared["Elmer"]
        files = dict(prep.files)
        files["slab.sif"] = files["slab.sif"].replace(
            "Coordinate Mapping(3) = 1 2 3", "Coordinate Mapping(3) = 3 2 1"
        )
        with self.assertRaises(ValueError):
            elmer_input.audit(files, self.config, [2, 2, 2])

    def test_fabricated_balanced_reactions_do_not_pass_physics(self):
        r = copy.deepcopy(self.contracts["Elmer"]["runs"][0])
        nodes = {n["id"]: n["coordinates"] for n in r["mesh"]["nodes"]}
        hot = next(
            v
            for v in r["quantities"]["nodal_heat_rate"]["values"]
            if nodes[v["node"]][0] == 0
        )
        cold = next(
            v
            for v in r["quantities"]["nodal_heat_rate"]["values"]
            if nodes[v["node"]][0] == 0.1
        )
        hot["value"] += 50
        cold["value"] -= 50
        result = self.verifier.verify(self.problem, r)
        self.assertTrue(result["checks"]["reaction_source_energy_conservation"])
        self.assertFalse(result["checks"]["boundary_reaction_reference"])
        self.assertEqual(result["status"], "failed")

    def test_real_wrong_source_converges_conserves_but_fails_reference(self):
        d = self.directory()
        wrong = copy.deepcopy(self.config)
        wrong["internal_heat_generation"] = 50000
        actual = load_problem(wrong)
        a = self.adapters["Elmer"]
        prep = a.prepare(actual, [2, 2, 2])
        for name, text in prep.files.items():
            (d / name).write_text(text)
        version = a.execute(prep, d, self.exes["Elmer"])
        r = a.extract(prep, d, version)
        result = self.verifier.verify(self.problem, r)
        write_json(d / "actual-physical-spec.json", wrong)
        write_json(d / "expected-physical-spec.json", self.config)
        write_json(d / "failure-verification.json", result)
        self.assertTrue(result["checks"]["reaction_source_energy_conservation"])
        self.assertFalse(result["checks"]["audited_source_matches_specification"])
        self.assertFalse(result["checks"]["nodal_temperature"])
        self.assertEqual(result["status"], "failed")
        with self.assertRaises(ValueError):
            elmer_input.audit(prep.files, self.config, [2, 2, 2])

    def test_native_source_archives_replay_without_solver_execution(self):
        for name, a in self.adapters.items():
            replay(self.archives[name], a)

    def test_rehashed_native_temperature_corruption_is_detected(self):
        d = self.root / self._testMethodName
        shutil.copytree(self.archives["Elmer"], d)
        p = next((d / "mesh-2x2x2").glob("slab*.vtu"))
        xml = ET.fromstring(p.read_text())
        a = xml.find('.//DataArray[@Name="temperature"]')
        a.text = " ".join(str(float(v) + 5) for v in a.text.split())
        p.write_text(ET.tostring(xml, encoding="unicode"))
        write_json(d / "manifest.json", manifest(d))
        with self.assertRaises(ValueError):
            replay(d, self.adapters["Elmer"])

    def test_rehashed_convergence_and_energy_metrics_are_detected(self):
        for field, value in (
            ("temperature_L2_K_m3_2", 0),
            ("actual_source_integral_W", 123),
        ):
            d = self.root / (self._testMethodName + "-" + field)
            shutil.copytree(self.archives["Elmer"], d)
            p = d / "normalized-results.json"
            r = json.loads(p.read_text())
            r["runs"][0]["verification"]["metrics"][field] = value
            write_json(p, r)
            write_json(d / "manifest.json", manifest(d))
            with self.assertRaises(ValueError):
                replay(d, self.adapters["Elmer"])

    def test_rehashed_nonnodal_sample_evidence_is_detected(self):
        d = self.root / self._testMethodName
        shutil.copytree(self.archives["Elmer"], d)
        p = d / "mesh-2x2x2" / "nonnodal-samples.csv"
        text = p.read_text()
        p.write_text(text.replace("362.5", "363.5", 1))
        if p.read_text() == text:
            p.write_text(text + "corrupt extra row\n")
        write_json(d / "manifest.json", manifest(d))
        with self.assertRaises(ValueError):
            replay(d, self.adapters["Elmer"])

    def test_common_raw_flux_comparison_and_altered_metrics(self):
        a, b = self.contracts["CalculiX"], self.contracts["Elmer"]
        result = compare_contracts(a, b)
        self.assertEqual(result["status"], "passed")
        changed = copy.deepcopy(b)
        changed["runs"][0]["verification"]["metrics"]["temperature_L2_K_m3_2"] = 0
        with self.assertRaises(ValueError):
            compare_contracts(a, changed)

    def test_comparison_plot_replay_and_rehashed_tampering(self):
        d = self.root / self._testMethodName
        args = [
            "--calculix",
            str(self.archives["CalculiX"]),
            "--elmer",
            str(self.archives["Elmer"]),
            "--output",
            str(d),
        ]
        self.assertEqual(compare_main(args), 0)
        self.assertEqual(compare_main(args + ["--replay"]), 0)
        p = d / "convergence-plot-data.json"
        r = json.loads(p.read_text())
        r["series"]["Elmer"][0]["temperature_L2_K_m3_2"] = 0
        write_json(p, r)
        write_json(d / "manifest.json", manifest(d))
        with self.assertRaises(ValueError):
            compare_main(args + ["--replay"])


if __name__ == "__main__":
    unittest.main()
