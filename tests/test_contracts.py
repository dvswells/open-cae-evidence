"""Real-solver archive version compatibility and adapter-independent field tests."""

import copy
import json
import os
from pathlib import Path
import shutil
import tempfile
import unittest

from open_cae_evidence.evidence import manifest, write_json
from open_cae_evidence.results import validate_contract
from open_cae_evidence.runner import replay, run_benchmark
from open_cae_evidence.specification import load_problem
from open_cae_evidence.verification import ModalVerifier


class ContractRegressionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        config_file = (
            Path(__file__).resolve().parents[1] / "benchmarks/cantilever_tip_mass.json"
        )
        cls.config = json.loads(config_file.read_text())
        cls.problem = load_problem(cls.config)
        solver = shutil.which(os.environ.get("CCX", "ccx"))
        if not solver:
            raise RuntimeError("Real solver required for contract regression")
        root = Path(os.environ.get("TEST_EVIDENCE_DIR", "evidence/tests"))
        root.mkdir(parents=True, exist_ok=True)
        cls.path = Path(tempfile.mkdtemp(prefix="contract-", dir=root)).resolve()
        cls.archive = cls.path / "current"
        if run_benchmark(cls.config, cls.archive, solver) != 0:
            raise RuntimeError("Contract regression source run failed")
        cls.contract = json.loads((cls.archive / "normalized-results.json").read_text())

    def test_current_version_and_arbitrary_solver_labels(self):
        validate_contract(self.contract)
        normalized = copy.deepcopy(self.contract["runs"][0])
        # Renumber a solid node above the MASS ID. This transforms only labels
        # on real output, not physics; verification must follow roles, not labels.
        old = normalized["mesh"]["nodes"][0]["id"]
        new = 100000
        for n in normalized["mesh"]["nodes"]:
            if n["id"] == old:
                n["id"] = new
        for e in normalized["mesh"]["elements"]:
            e["connectivity"] = [new if n == old else n for n in e["connectivity"]]
        for q in normalized["quantities"].values():
            for r in q["values"]:
                if r.get("node") == old:
                    r["node"] = new
        normalized["solver"]["name"] = "label-renaming-interop-test"
        result = ModalVerifier().verify(self.problem, normalized)
        self.assertEqual(result["status"], "passed")

    def test_legacy_v1_v2_archives_replay_real_fields(self):
        # Compatibility fixtures retain actual solver fields and original
        # verification structure, removing only newly introduced contract files.
        for version in (1, 2):
            with self.subTest(version=version):
                path = self.path / f"legacy-{version}"
                shutil.copytree(self.archive, path)
                (path / "normalized-results.json").unlink()
                for p in path.glob("mesh-*/solver-results.json"):
                    p.unlink()
                result = json.loads((path / "results.json").read_text())
                result["schema_version"] = version
                write_json(path / "results.json", result)
                write_json(path / "manifest.json", manifest(path))
                replay(path)

    def test_unknown_archive_schema_is_not_treated_as_latest(self):
        path = self.path / "unknown-version"
        shutil.copytree(self.archive, path)
        result = json.loads((path / "results.json").read_text())
        result["schema_version"] = 99
        write_json(path / "results.json", result)
        write_json(path / "manifest.json", manifest(path))
        with self.assertRaisesRegex(ValueError, "Unsupported evidence archive version"):
            replay(path)

    def test_changed_tolerances_do_not_replay_as_passed(self):
        path = self.path / "changed-tolerance"
        shutil.copytree(self.archive, path)
        result = json.loads((path / "results.json").read_text())
        result["tolerances"]["analytical_relative"] = 0.5
        write_json(path / "results.json", result)
        write_json(path / "manifest.json", manifest(path))
        with self.assertRaisesRegex(ValueError, "tolerances"):
            replay(path)

    def test_incomplete_modes_and_inconsistent_eigenvalues_fail(self):
        normalized = copy.deepcopy(self.contract["runs"][0])
        normalized["quantities"]["displacement"]["values"].pop()
        with self.assertRaises(ValueError):
            ModalVerifier().verify(self.problem, normalized)
        normalized = copy.deepcopy(self.contract["runs"][0])
        normalized["quantities"]["eigenvalue"]["values"][0]["value"] *= 2
        with self.assertRaises(ValueError):
            ModalVerifier().verify(self.problem, normalized)
