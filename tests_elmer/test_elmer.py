"""Adversarial tests using real native Elmer solves and archived output.

Missing Elmer is a failure, never a skipped or mocked integration test.
"""

import copy
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest
import xml.etree.ElementTree as ET

from open_cae_evidence.adapters.elmer import ElmerAdapter
from open_cae_evidence.adapters import elmer_input, elmer_output
from open_cae_evidence.comparison import (
    compare_archives,
    compare_contracts,
    main as compare_main,
)
from open_cae_evidence.evidence import manifest, write_json
from open_cae_evidence.results import validate_contract
from open_cae_evidence.runner import replay, run_benchmark
from open_cae_evidence.specification import load_problem
from open_cae_evidence.verification import ThermalVerifier


class ElmerIntegrationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temp = tempfile.TemporaryDirectory()
        cls.root = Path(os.environ.get("CAE_ELMER_TEST_EVIDENCE", cls.temp.name))
        cls.root.mkdir(parents=True, exist_ok=True)
        cls.executable = shutil.which(os.environ.get("ELMER_SOLVER", "ElmerSolver"))
        if cls.executable is None:
            raise RuntimeError("Real pinned ElmerSolver is required")
        cls.config = json.loads(Path("benchmarks/slab_steady_heat.json").read_text())
        cls.problem = load_problem(cls.config)
        cls.adapter = ElmerAdapter()
        cls.archive = cls.root / "valid-elmer"
        if run_benchmark(cls.config, cls.archive, cls.executable, cls.adapter):
            raise RuntimeError("Real Elmer baseline failed")
        cls.contract = json.loads((cls.archive / "normalized-results.json").read_text())
        cls.prepared = cls.adapter.prepare(cls.problem, [2, 2, 2])
        cls.meshdir = cls.archive / "mesh-2x2x2"
        cls.native = next(cls.meshdir.glob("slab*.vtu")).name
        cls.native_text = (cls.meshdir / cls.native).read_text()
        cls.ccx_archive = cls.root / "valid-calculix"
        ccx = shutil.which(os.environ.get("CCX", "ccx"))
        if ccx is None or run_benchmark(cls.config, cls.ccx_archive, ccx):
            raise RuntimeError("Real CalculiX baseline failed")

    @classmethod
    def tearDownClass(cls):
        cls.temp.cleanup()

    def audit(self, **changed):
        files = dict(self.prepared.files)
        files.update(changed)
        return elmer_input.audit(files, self.config, [2, 2, 2])

    def native_change(self, name, modify):
        xml = ET.fromstring(self.native_text)
        field = xml.find(f'.//DataArray[@Name="{name}"]')
        modify(field)
        return ET.tostring(xml, encoding="unicode")

    def parsed_normalized(self, text):
        directory = self.root / self._testMethodName
        shutil.copytree(self.meshdir, directory)
        (directory / self.native).write_text(text)
        return self.adapter.extract(self.prepared, directory, "9.0")

    def copied_archive(self):
        directory = self.root / self._testMethodName
        shutil.copytree(self.archive, directory)
        return directory

    def test_three_meshes_common_contract_and_shared_verifier(self):
        validate_contract(self.contract)
        self.assertEqual(self.contract["schema_version"], "1.0.0")
        self.assertEqual(len(self.contract["runs"]), 3)
        for run in self.contract["runs"]:
            checks = ThermalVerifier().verify(self.problem, run)
            self.assertEqual(checks["status"], "passed")
            self.assertEqual(
                run["quantities"]["heat_flux"]["source"]["kind"],
                "postprocessed_from_solver",
            )
        study = json.loads((self.archive / "results.json").read_text())[
            "mesh_consistency"
        ]
        self.assertTrue(study["passed"])
        self.assertIsNone(study["demonstrated_convergence_order"])

    def test_native_and_normalized_archive_replay(self):
        replay(self.archive, self.adapter)

    def test_modal_is_explicitly_unsupported(self):
        modal = load_problem(json.loads(Path("benchmarks/cantilever.json").read_text()))
        with self.assertRaises(ValueError):
            self.adapter.prepare(modal, [12, 2, 2])

    def test_wrong_conductivity_units(self):
        config = copy.deepcopy(self.config)
        config["units"]["conductivity"] = "W/(mm*K)"
        with self.assertRaises(ValueError):
            load_problem(config)
        with self.assertRaises(ValueError):
            self.audit(
                **{
                    "units.json": self.prepared.files["units.json"].replace(
                        "W/(m*K)", "W/(mm*K)"
                    )
                }
            )

    def test_nonpositive_conductivity(self):
        for k in (0, -10):
            config = copy.deepcopy(self.config)
            config["material"]["conductivity"] = k
            with self.assertRaises(ValueError):
                load_problem(config)

    def test_wrong_serialized_conductivity(self):
        # The independently parsed audit also fails without relying on the
        # canonical text comparison used to reject unsupported keywords.
        with self.assertRaises(ValueError):
            elmer_input.audit_sif(
                self.prepared.files["slab.sif"].replace(
                    "Heat Conductivity = Real 10", "Heat Conductivity = Real 20"
                ),
                self.config,
            )
        with self.assertRaises(ValueError):
            self.audit(
                **{
                    "slab.sif": self.prepared.files["slab.sif"].replace(
                        "Heat Conductivity = Real 10", "Heat Conductivity = Real 20"
                    )
                }
            )

    def test_swapped_or_missing_temperature_boundaries(self):
        sif = self.prepared.files["slab.sif"]
        for changed in (
            sif.replace("Temperature = Real 400", "Temperature = Real 300"),
            sif.replace("  Temperature = Real 400\n", ""),
        ):
            with self.assertRaises(ValueError):
                self.audit(**{"slab.sif": changed})

    def test_inverted_native_cell(self):
        rows = self.prepared.files["mesh.elements"].splitlines()
        fields = rows[0].split()
        fields[3], fields[4] = fields[4], fields[3]
        rows[0] = " ".join(fields)
        with self.assertRaises(ValueError):
            self.audit(**{"mesh.elements": "\n".join(rows) + "\n"})

    def test_duplicate_native_connectivity(self):
        rows = self.prepared.files["mesh.elements"].splitlines()
        fields = rows[0].split()
        fields[-1] = fields[-2]
        rows[0] = " ".join(fields)
        with self.assertRaises(ValueError):
            self.audit(**{"mesh.elements": "\n".join(rows) + "\n"})

    def test_wrong_boundary_assignment(self):
        rows = self.prepared.files["mesh.boundary"].splitlines()
        fields = rows[0].split()
        fields[1] = "2" if fields[1] == "1" else "1"
        rows[0] = " ".join(fields)
        with self.assertRaises(ValueError):
            self.audit(**{"mesh.boundary": "\n".join(rows) + "\n"})

    def test_incorrect_coordinate_mapping(self):
        with self.assertRaises(ValueError):
            self.audit(
                **{
                    "slab.sif": self.prepared.files["slab.sif"].replace(
                        "Coordinate Mapping(3) = 1 2 3", "Coordinate Mapping(3) = 3 2 1"
                    )
                }
            )

    def test_incomplete_real_temperature_output(self):
        text = self.native_change(
            "temperature", lambda a: setattr(a, "text", " ".join(a.text.split()[:-1]))
        )
        with self.assertRaises(ValueError):
            self.parsed_normalized(text)

    def test_reversed_real_native_flux(self):
        text = self.native_change(
            "temperature flux",
            lambda a: setattr(
                a, "text", " ".join(str(-float(x)) for x in a.text.split())
            ),
        )
        run = self.parsed_normalized(text)
        result = ThermalVerifier().verify(self.problem, run)
        write_json(
            self.root / self._testMethodName / "failure-verification.json", result
        )
        self.assertEqual(result["status"], "failed")
        self.assertFalse(result["checks"]["solver_flux_direction_magnitude"])

    def test_wrong_native_temperature(self):
        text = self.native_change(
            "temperature",
            lambda a: setattr(
                a, "text", " ".join(str(float(x) + 5) for x in a.text.split())
            ),
        )
        self.assertEqual(
            ThermalVerifier().verify(self.problem, self.parsed_normalized(text))[
                "status"
            ],
            "failed",
        )

    def test_wrong_native_coordinate_units(self):
        xml = ET.fromstring(self.native_text)
        a = xml.find(".//Points/DataArray")
        a.text = " ".join(str(float(x) * 1000) for x in a.text.split())
        with self.assertRaises(ValueError):
            self.parsed_normalized(ET.tostring(xml, encoding="unicode"))

    def test_corrupt_real_output_connectivity(self):
        def change(a):
            fields = a.text.split()
            fields[0], fields[1] = fields[1], fields[0]
            a.text = " ".join(fields)

        with self.assertRaises(ValueError):
            self.parsed_normalized(self.native_change("connectivity", change))

    def test_missing_native_output(self):
        directory = self.copied_archive()
        (directory / "mesh-2x2x2" / self.native).unlink()
        write_json(directory / "manifest.json", manifest(directory))
        with self.assertRaises(ValueError):
            replay(directory, self.adapter)

    def test_nonconverged_archived_log(self):
        log = (self.meshdir / "solver.log").read_text()
        for changed in (
            log.replace("ELMER SOLVER FINISHED", "INTERRUPTED"),
            log + "\nLinear iteration did not converge to tolerance\n",
            log.replace("HeatSolve: Assembly done", "HeatSolve: assembly missing"),
        ):
            with self.assertRaises(RuntimeError):
                elmer_output.check_log(changed)

    def forced_solve(self, sif):
        directory = self.root / self._testMethodName
        directory.mkdir()
        for name, text in self.prepared.files.items():
            (directory / name).write_text(sif if name == "slab.sif" else text)
        env = os.environ.copy()
        home = Path(self.executable).resolve().parent.parent
        env.update(
            ELMER_HOME=str(home),
            LD_LIBRARY_PATH=str(home / "lib")
            + os.pathsep
            + env.get("LD_LIBRARY_PATH", ""),
            OMP_NUM_THREADS="1",
            OPENBLAS_NUM_THREADS="1",
        )
        with (directory / "solver.log").open("w") as handle:
            run = subprocess.run(
                [self.executable, "slab.sif"],
                cwd=directory,
                env=env,
                stdout=handle,
                stderr=subprocess.STDOUT,
                timeout=120,
                check=False,
            )
        return directory, run.returncode

    def test_real_solver_failure_nonconvergence(self):
        directory, _ = self.forced_solve(
            self.prepared.files["slab.sif"].replace(
                "Linear System Max Iterations = 1000",
                "Linear System Max Iterations = 1",
            )
        )
        with self.assertRaises(RuntimeError):
            elmer_output.check_log((directory / "solver.log").read_text())

    def test_real_swapped_boundaries_converge_but_fail_physics(self):
        sif = self.prepared.files["slab.sif"].replace(
            "Temperature = Real 400", "Temperature = Real TEMP"
        )
        sif = sif.replace("Temperature = Real 300", "Temperature = Real 400").replace(
            "Temperature = Real TEMP", "Temperature = Real 300"
        )
        directory, status = self.forced_solve(sif)
        self.assertEqual(status, 0)
        elmer_output.check_log((directory / "solver.log").read_text())
        # Input audit rejects it. Bypass the adapter only to independently expose
        # the converged actual wrong physical solution to the unchanged verifier.
        with self.assertRaises(ValueError):
            self.audit(**{"slab.sif": sif})
        data = elmer_output.parse_vtu(
            next(directory.glob("slab*.vtu")).read_text(),
            self.prepared.nodes,
            self.prepared.elements,
        )
        del data["nodal_flux"]
        from open_cae_evidence.thermal_verification import verify

        result = verify(self.config, self.prepared.nodes, self.prepared.elements, data)
        write_json(directory / "failure-verification.json", result)
        self.assertFalse(result["passed"])
        self.assertFalse(result["checks"]["nodal_temperature"])
        self.assertFalse(result["checks"]["solver_flux_direction_magnitude"])
        self.assertTrue(result["checks"]["energy_conservation"])

    def test_rehashed_native_flux_tampering(self):
        directory = self.copied_archive()
        text = self.native_change(
            "temperature flux",
            lambda a: setattr(
                a, "text", " ".join(str(-float(x)) for x in a.text.split())
            ),
        )
        (directory / "mesh-2x2x2" / self.native).write_text(text)
        write_json(directory / "manifest.json", manifest(directory))
        with self.assertRaises(ValueError):
            replay(directory, self.adapter)

    def test_rehashed_normalized_contract_tampering(self):
        directory = self.copied_archive()
        value = copy.deepcopy(self.contract)
        value["runs"][0]["quantities"]["temperature"]["values"][0]["value"] += 5
        write_json(directory / "normalized-results.json", value)
        write_json(directory / "manifest.json", manifest(directory))
        with self.assertRaises(ValueError):
            replay(directory, self.adapter)

    def test_incomplete_or_wrong_runtime_provenance(self):
        identity = copy.deepcopy(self.contract["solver"]["provenance"])
        identity["upstream_source_commit"] = "0" * 40
        with self.assertRaises(ValueError):
            self.adapter.validate_provenance(identity)
        identity = copy.deepcopy(self.contract["solver"]["provenance"])
        del identity["runtime_library_sha256"]["share/elmersolver/lib/FluxSolver.so"]
        with self.assertRaises(ValueError):
            self.adapter.validate_provenance(identity)

    def test_incorrect_normalized_physical_units(self):
        for name, unit in (
            ("temperature", "degC"),
            ("heat_flux", "W/mm^2"),
            ("nodal_heat_rate", "kW"),
        ):
            run = copy.deepcopy(self.contract["runs"][0])
            run["quantities"][name]["unit"] = unit
            with self.assertRaises(ValueError):
                ThermalVerifier().verify(self.problem, run)

    def test_coordinate_comparison_and_independent_gates(self):
        ccx = json.loads((self.ccx_archive / "normalized-results.json").read_text())
        changed = copy.deepcopy(self.contract)
        for run in changed["runs"]:
            # Output record order and native IDs already differ between solvers.
            for q in run["quantities"].values():
                q["values"].reverse()
        result = compare_contracts(ccx, changed)
        self.assertEqual(result["status"], "passed")
        changed["runs"][0]["quantities"]["heat_flux"]["values"][0]["value"][0] *= -1
        with self.assertRaises(ValueError):
            compare_contracts(ccx, changed)

    def test_comparison_archive_replay_and_rehashed_tampering(self):
        directory = self.root / self._testMethodName
        args = [
            "--calculix",
            str(self.ccx_archive),
            "--elmer",
            str(self.archive),
            "--output",
            str(directory),
        ]
        self.assertEqual(compare_main(args), 0)
        self.assertEqual(compare_main(args + ["--replay"]), 0)
        value = json.loads((directory / "comparison.json").read_text())
        value["runs"][0]["temperature_max_difference_K"] = 100
        write_json(directory / "comparison.json", value)
        write_json(directory / "manifest.json", manifest(directory))
        with self.assertRaises(ValueError):
            compare_main(args + ["--replay"])


if __name__ == "__main__":
    unittest.main()
