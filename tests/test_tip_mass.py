"""Tip-mass physics, coupling and archive tests; all solver data is real."""

import copy
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest

from open_cae_evidence.benchmark import (main, parse_dat, replay, run_mesh, sha256,
                                         verify_modes, write_json)
from open_cae_evidence.frd import verify_frd
from open_cae_evidence.model import area_weights, audit_deck, make_deck, validate
from open_cae_evidence.reference import (BETAS, characteristic, mass_ratio, references,
                                         shape, tip_roots)

CONFIG = Path(__file__).resolve().parents[1] / "benchmarks/cantilever_tip_mass.json"


class TipInputAndReferenceTests(unittest.TestCase):
    def setUp(self):
        self.config = json.loads(CONFIG.read_text())

    def test_mass_units_and_values_rejected(self):
        for key, value in [("unit", "g"), ("unit", "kg/m^3"), ("unit", "N"),
                           ("mass", -1), ("mass", 1170), ("mass", True),
                           ("mass", float("nan")), ("mass", float("inf"))]:
            with self.subTest(key=key, value=value):
                config = copy.deepcopy(self.config)
                config["tip_mass"][key] = value
                with self.assertRaises(ValueError):
                    validate(config)

    def test_mass_placement_and_coupling_rejected(self):
        for key, value in [("position", [0, 0, 0]), ("position", [1, 0.001, 0]),
                           ("position", [1000, 0, 0]), ("position", [1, 0, 0.001]),
                           ("coupling", "single-node"), ("coupling", "rigid"),
                           ("rotary_inertia", "included")]:
            with self.subTest(key=key, value=value):
                config = copy.deepcopy(self.config)
                config["tip_mass"][key] = value
                with self.assertRaises(ValueError):
                    validate(config)

    def test_q8_weights_integrate_translation_rotation_and_quadratics(self):
        for mesh in self.config["meshes"]:
            _, nodes, audit = make_deck(self.config, mesh)
            weights = area_weights(nodes, self.config, mesh)
            self.assertEqual(audit["coupled_face_nodes"], len(weights))
            self.assertAlmostEqual(sum(weights.values()), 1, places=12)
            for axis, side in ((1, 0.015), (2, 0.01)):
                self.assertAlmostEqual(sum(w * nodes[n][axis] for n, w in weights.items()), 0, places=12)
                self.assertAlmostEqual(sum(w * nodes[n][axis]**2 for n, w in weights.items()), side**2/12, places=12)
            self.assertAlmostEqual(sum(w*nodes[n][1]*nodes[n][2] for n, w in weights.items()), 0, places=12)
            self.assertAlmostEqual(audit["total_mass_kg"], 2.34, places=10)

    def test_serialized_mass_and_equations_are_audited(self):
        deck, nodes, _ = make_deck(self.config, [12, 2, 2])
        ref = max(nodes)
        prefix, equations = deck.split("*EQUATION\n", 1)
        block, suffix = equations.split("*BOUNDARY", 1)
        rows = block.splitlines()
        first = rows[1].split(",")
        wrong_weights = rows.copy(); fields = first.copy(); fields[2] = str(float(fields[2])*2)
        wrong_weights[1] = ",".join(fields)
        wrong_axis = rows.copy(); fields = first.copy(); fields[1] = "2"; wrong_axis[1] = ",".join(fields)
        mass_dependent = rows.copy()
        # Swap the first term with the MASS term, retaining identical mathematical coefficients.
        count = int(rows[0]); last_row = (count + 3)//4
        terms = [v.strip() for row in rows[1:last_row+1] for v in row.split(",")]
        terms[:3], terms[-3:] = terms[-3:], terms[:3]
        mass_dependent[1:last_row+1] = [", ".join(terms[i:i+12]) for i in range(0, len(terms), 12)]
        corruptions = [deck.replace(f"49, {ref}", "49, 2"),
                       deck.replace("*MASS, ELSET=TIPMASS\n1.17", "*MASS, ELSET=TIPMASS\n1170"),
                       deck.replace(f"{ref}, 1, 0, 0", f"{ref}, 1, 0.001, 0"),
                       deck.replace("*EQUATION\n" + block, ""),
                       *(prefix + "*EQUATION\n" + "\n".join(r) + "\n*BOUNDARY" + suffix
                         for r in (wrong_weights, wrong_axis, mass_dependent))]
        for corrupt in corruptions:
            with self.subTest(length=len(corrupt)):
                with self.assertRaises(ValueError):
                    audit_deck(corrupt, self.config, [12, 2, 2])

    def test_characteristic_and_dynamic_tip_conditions(self):
        mu = mass_ratio(self.config)
        for beta in tip_roots(mu):
            self.assertAlmostEqual(characteristic(beta, mu), 0, places=11)
        # Independently check the ODE end conditions rather than the determinant only.
        for beta in tip_roots(mu)[:3]:
            gamma = (math.cosh(beta)+math.cos(beta))/(math.sinh(beta)+math.sin(beta))
            moment = math.cosh(beta)+math.cos(beta)-gamma*(math.sinh(beta)+math.sin(beta))
            shear_plus_inertia = (math.sinh(beta)-math.sin(beta)-gamma*(math.cosh(beta)+math.cos(beta))
                                  + mu*beta*shape(beta, 1))
            self.assertAlmostEqual(moment, 0, places=10)
            self.assertAlmostEqual(shear_plus_inertia, 0, places=9)
            self.assertEqual(shape(beta, 0), 0)

    def test_zero_mass_limit_and_monotonic_roots(self):
        self.assertEqual(tip_roots(0), BETAS)
        for a, b in zip(tip_roots(1e-9), BETAS):
            self.assertAlmostEqual(a/b, 1, delta=1e-8)
        for a, b, c in zip(tip_roots(0.5), tip_roots(1), tip_roots(2)):
            self.assertGreater(a, b); self.assertGreater(b, c)
        # For a very heavy mass the first mode approaches sqrt(3EI/(Mt L^3)).
        beta = tip_roots(10)[0]
        self.assertAlmostEqual(beta**4 * 10 / 3, 1, delta=0.03)


class RealTipMassTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.config = json.loads(CONFIG.read_text())
        cls.solver = shutil.which(os.environ.get("CCX", "ccx"))
        if not cls.solver:
            raise RuntimeError("Real CalculiX required, including tip mass and archive tests")
        root = Path(os.environ.get("TEST_EVIDENCE_DIR", "evidence/tests")); root.mkdir(parents=True, exist_ok=True)
        cls.path = Path(tempfile.mkdtemp(prefix="tip-mass-", dir=root)).resolve()
        output = cls.path / "archive"
        if main(["--config", str(CONFIG), "--solver", cls.solver, "--output", str(output)]) != 0:
            raise RuntimeError("Tip-mass benchmark failed during integration setup")
        cls.results = json.loads((output / "results.json").read_text())
        cls.base = cls.results["runs"][0]
        cls.mesh_path = output / "mesh-12x2x2"
        cls.deck, cls.nodes, _ = make_deck(cls.config, [12, 2, 2])
        cls.dat = (cls.mesh_path / "beam.dat").read_text()
        cls.frd = (cls.mesh_path / "beam.frd").read_text()
        cls.f, cls.lam, cls.modes = parse_dat(cls.dat, cls.nodes)

    def scaling(self, label, config, expected):
        run = run_mesh(config, [12, 2, 2], self.path / label, self.solver)
        records = []
        for original, changed, target in zip(self.base["modes"], run["modes"], expected):
            actual = changed["fem_frequency_hz"] / original["fem_frequency_hz"]
            records.append({"mode": changed["mode"], "actual_ratio": actual, "reference_ratio": target,
                            "relative_error": abs(actual/target-1)})
        write_json(self.path / f"{label}-scaling.json", {"config": config, "ratios": records})
        return run, records

    def test_all_six_modes_shapes_meshes_and_coupling(self):
        self.assertEqual(self.results["status"], "passed")
        self.assertEqual(len(self.results["runs"]), 3)
        for run in self.results["runs"]:
            self.assertEqual(len(run["modes"]), 6)
            self.assertTrue(run["frd"]["passed"])
            for mode in run["modes"]:
                self.assertTrue(mode["checks"]["tip_coupling"])
                self.assertTrue(mode["passed"])

    def test_tip_mass_doubling_matches_analytical_ratios(self):
        config = copy.deepcopy(self.config); config["tip_mass"]["mass"] *= 2
        expected = [b["frequency_hz"]/a["frequency_hz"] for a, b in zip(references(self.config), references(config))]
        run, records = self.scaling("double-tip-mass", config, expected)
        self.assertTrue(run["passed"])
        for record in records:
            self.assertLess(record["actual_ratio"], 1)
            self.assertLess(record["relative_error"], 0.01)

    def test_all_mass_scaling_is_inverse_square_root(self):
        config = copy.deepcopy(self.config)
        config["tip_mass"]["mass"] *= 2; config["material"]["density"] *= 2
        _, records = self.scaling("double-all-masses", config, [1/math.sqrt(2)]*6)
        for record in records:
            self.assertAlmostEqual(record["actual_ratio"], record["reference_ratio"], delta=2e-5)

    def test_modulus_scaling_with_tip_mass(self):
        config = copy.deepcopy(self.config); config["material"]["young_modulus"] *= 4
        _, records = self.scaling("quadruple-modulus", config, [2]*6)
        for record in records:
            self.assertAlmostEqual(record["actual_ratio"], 2, delta=2e-5)

    def test_zero_tip_mass_recovers_original_benchmark(self):
        config = copy.deepcopy(self.config); config["tip_mass"]["mass"] = 0
        no_tip = copy.deepcopy(config); del no_tip["tip_mass"]
        self.assertEqual(make_deck(config, [12, 2, 2])[0], make_deck(no_tip, [12, 2, 2])[0])
        run = run_mesh(config, [12, 2, 2], self.path / "zero-tip-mass", self.solver)
        self.assertTrue(run["passed"])
        for a, b in zip(self.base["modes"], run["modes"]):
            self.assertLess(a["fem_frequency_hz"], b["fem_frequency_hz"])

    def test_coupling_residual_detects_corrupted_reference_displacement(self):
        changed = copy.deepcopy(self.modes)
        changed[0][max(self.nodes)] = tuple(v * 1.1 for v in changed[0][max(self.nodes)])
        records = verify_modes(self.config, [12, 2, 2], self.nodes, self.f, self.lam, changed)
        self.assertFalse(records[0]["checks"]["tip_coupling"])

    def test_alternative_face_dependent_dof_is_physically_equivalent(self):
        prefix, equations = self.deck.split("*EQUATION\n", 1)
        block, suffix = equations.split("*BOUNDARY", 1)
        rows = block.splitlines(); cursor = 0
        for _ in range(3):
            count = int(rows[cursor]); cursor += 1
            # Change the algebraic elimination pivot, retaining all coefficients.
            terms = rows[cursor].split(",")
            terms[:3], terms[3:6] = terms[3:6], terms[:3]
            rows[cursor] = ",".join(terms)
            cursor += (count+3)//4
        deck = prefix + "*EQUATION\n" + "\n".join(rows) + "\n*BOUNDARY" + suffix
        audit_deck(deck, self.config, [12, 2, 2])
        path = self.path / "alternative-pivot"; path.mkdir()
        (path / "beam.inp").write_text(deck)
        env = os.environ.copy(); env.update(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1")
        with (path / "solver.log").open("w") as log:
            subprocess.run([self.solver, "-i", "beam"], cwd=path, env=env, stdout=log, stderr=subprocess.STDOUT, check=True, timeout=300)
        self.assertNotIn("*ERROR", (path / "solver.log").read_text())
        f, lam, modes = parse_dat((path / "beam.dat").read_text(), self.nodes)
        self.assertTrue(all(m["passed"] for m in verify_modes(self.config, [12, 2, 2], self.nodes, f, lam, modes)))
        for a, b in zip(f, self.f):
            self.assertAlmostEqual(a/b, 1, delta=2e-5)
        verify_frd((path / "beam.frd").read_text(), deck, self.nodes, f, lam, modes, self.base["solver_version"], {max(self.nodes)})

    def test_frd_semantics_reject_actual_data_corruption(self):
        lines = self.frd.splitlines()
        coord = next(i for i, line in enumerate(lines) if line.strip().startswith("2C")) + 1
        dataset = next(i for i, line in enumerate(lines) if line.strip().startswith("100CL"))
        altered_coords = lines.copy(); altered_coords[coord] = altered_coords[coord][:13] + " 1.00000E+00" + altered_coords[coord][25:]
        altered_frequency = lines.copy(); altered_frequency[dataset] = altered_frequency[dataset].replace(str(self.f[0])[:5], "9.999", 1)
        connectivity = lines.copy()
        conn = next(i for i, line in enumerate(lines) if line[:3] == " -2")
        connectivity[conn] = connectivity[conn].replace("         1", "         2", 1)
        displacement = lines.copy(); pos = dataset + 7
        displacement[pos] = displacement[pos][:13] + " 1.00000E+00" + displacement[pos][25:]
        duplicate = lines.copy(); duplicate.insert(coord, duplicate[coord])
        bad = ["\n".join(s) for s in (altered_coords, altered_frequency, connectivity, displacement, duplicate)]
        bad += [self.frd.replace(" -5  D1", " -5  D2", 1),
                self.frd.replace("1PMODE                         1", "1PMODE                         2", 1),
                self.frd[:-100], self.frd + "unexpected\n"]
        for corrupt in bad:
            with self.subTest(length=len(corrupt)):
                with self.assertRaises(ValueError):
                    verify_frd(corrupt, self.deck, self.nodes, self.f, self.lam, self.modes,
                               self.base["solver_version"], {max(self.nodes)})

    def test_replay_checks_semantics_even_with_updated_manifest(self):
        replay(self.path / "archive")
        corrupt = self.path / "rehashed-corrupt-archive"
        shutil.copytree(self.path / "archive", corrupt)
        frd = corrupt / "mesh-12x2x2/beam.frd"
        frd.write_text(frd.read_text().replace(" -5  D1", " -5  D2", 1))
        manifest = json.loads((corrupt / "manifest.json").read_text())
        manifest["mesh-12x2x2/beam.frd"] = sha256(frd)
        write_json(corrupt / "manifest.json", manifest)
        with self.assertRaisesRegex(ValueError, "component semantics"):
            replay(corrupt)
