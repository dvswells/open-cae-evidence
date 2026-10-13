"""Required integration tests. Missing CalculiX is an error, never a skip/mock."""

import copy
import json
import math
import os
from pathlib import Path
import shutil
import tempfile
import unittest

from open_cae_evidence.benchmark import evaluate, parse_dat, run_mesh
from open_cae_evidence.model import make_deck


class RealSolverTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.config = json.loads((Path(__file__).resolve().parents[1] / "benchmarks/cantilever.json").read_text())
        cls.solver = shutil.which(os.environ.get("CCX", "ccx"))
        if not cls.solver:
            raise RuntimeError("Real CalculiX required for integration tests")
        root = Path(os.environ.get("TEST_EVIDENCE_DIR", "evidence/tests"))
        root.mkdir(parents=True, exist_ok=True)
        cls.path = Path(tempfile.mkdtemp(prefix="real-solver-", dir=root)).resolve()
        cls.base = run_mesh(cls.config, [12, 2, 2], cls.path / "base", cls.solver)
        _, cls.nodes, _ = make_deck(cls.config, [12, 2, 2])
        cls.dat = (cls.path / "base/beam.dat").read_text()

    def test_real_modes_meet_independent_checks(self):
        self.assertTrue(self.base["passed"])
        self.assertEqual(len(self.base["modes"]), 6)

    def test_real_modulus_scaling(self):
        config = copy.deepcopy(self.config)
        config["material"]["young_modulus"] *= 4
        changed = run_mesh(config, [12, 2, 2], self.path / "modulus-times-four", self.solver)
        for a, b in zip(self.base["modes"], changed["modes"]):
            self.assertAlmostEqual(b["fem_frequency_hz"] / a["fem_frequency_hz"], 2, delta=2e-5)

    def test_real_density_scaling(self):
        config = copy.deepcopy(self.config)
        config["material"]["density"] *= 2
        changed = run_mesh(config, [12, 2, 2], self.path / "density-times-two", self.solver)
        for a, b in zip(self.base["modes"], changed["modes"]):
            self.assertAlmostEqual(b["fem_frequency_hz"] / a["fem_frequency_hz"], 1 / math.sqrt(2), delta=2e-5)

    def test_corrupted_real_output_rejected(self):
        # Work from actual .dat output, not a mocked solver response.
        lines = self.dat.splitlines()
        header = next(i for i, s in enumerate(lines) if "E I G E N V A L U E   O U T P U T" in s)
        first = next(i for i in range(header + 1, len(lines)) if len(lines[i].split()) == 5 and lines[i].split()[0] == "1")
        bad_frequency = lines.copy()
        values = bad_frequency[first].split()
        values[3] = "0.0"
        bad_frequency[first] = " ".join(values)
        for corrupt in ["\n".join(bad_frequency), self.dat[:len(self.dat) // 2],
                        self.dat.replace("E I G E N V A L U E    N U M B E R     1", "E I G E N V A L U E    N U M B E R     2", 1)]:
            with self.subTest(corruption_length=len(corrupt)):
                with self.assertRaises(ValueError):
                    parse_dat(corrupt, self.nodes)

    def test_successful_execution_alone_does_not_verify(self):
        frequencies, eigenvalues, modes = parse_dat(self.dat, self.nodes)
        wrong = [f * 1.1 for f in frequencies]
        self.assertFalse(all(r["passed"] for r in evaluate(self.config, self.nodes, wrong, eigenvalues, modes)))
        swapped = copy.deepcopy(modes)
        swapped[0], swapped[2] = swapped[2], swapped[0]
        records = evaluate(self.config, self.nodes, frequencies, eigenvalues, swapped)
        self.assertFalse(records[0]["checks"]["analytical_shape"])
        corrupted_clamp = copy.deepcopy(modes)
        fixed = next(n for n, xyz in self.nodes.items() if xyz[0] == 0)
        corrupted_clamp[0][fixed] = (1, 1, 1)
        records = evaluate(self.config, self.nodes, frequencies, eigenvalues, corrupted_clamp)
        self.assertFalse(records[0]["checks"]["clamped_displacements"])
        zero_modes = [{n: (0.0, 0.0, 0.0) for n in self.nodes} for _ in range(6)]
        with self.assertRaises(ValueError):
            evaluate(self.config, self.nodes, frequencies, eigenvalues, zero_modes)


if __name__ == "__main__":
    unittest.main()
