"""Small orchestration layer: problem -> adapter -> verifier -> evidence."""

import argparse
import csv
import json
from pathlib import Path
import shutil
import subprocess
import sys

from .adapters.calculix import CalculiXAdapter
from .evidence import manifest, provenance, sha256, verify_manifest, write_json
from .modal_report import report as modal_report
from .results import (
    ARCHIVE_VERSION,
    build_contract,
    modal_fields,
    thermal_fields,
    validate_contract,
)
from .specification import load_problem
from .thermal_report import report as thermal_report
from .verification import derived_quantities, for_problem


def export_csv(prepared, normalized, directory):
    if prepared.problem.analysis == "modal":
        _, _, modes = modal_fields(normalized)
        for index, mode in enumerate(modes, 1):
            with (directory / f"mode-{index}.csv").open("w", newline="") as handle:
                writer = csv.writer(handle)
                writer.writerow(["node", "x_m", "y_m", "z_m", "ux", "uy", "uz"])
                for n, xyz in sorted(prepared.nodes.items()):
                    writer.writerow([n, *xyz, *mode[n]])
    else:
        data = thermal_fields(normalized)
        with (directory / "temperature.csv").open("w", newline="") as handle:
            w = csv.writer(handle)
            w.writerow(["node", "x_m", "y_m", "z_m", "temperature_K", "heat_input_W"])
            for n, xyz in sorted(prepared.nodes.items()):
                w.writerow([n, *xyz, data["temperature"][n], data["heat_rate"][n]])
        with (directory / "heat-flux.csv").open("w", newline="") as handle:
            w = csv.writer(handle)
            w.writerow(
                [
                    "element",
                    "point",
                    "x_m",
                    "y_m",
                    "z_m",
                    "qx_W_m2",
                    "qy_W_m2",
                    "qz_W_m2",
                ]
            )
            for key, q in sorted(data["heat_flux"].items()):
                w.writerow([*key, *data["integration_coordinates"][key], *q])


def legacy_run(prepared, normalized):
    value = {
        "mesh": prepared.counts,
        "audit": prepared.audit,
        "solver_version": normalized["solver"]["version"],
        "frd": normalized["raw_output_verification"],
        "passed": normalized["verification"]["status"] == "passed",
    }
    if prepared.problem.analysis == "modal":
        value["modes"] = normalized["verification"]["metrics"]
    else:
        value["thermal"] = normalized["verification"]["metrics"]
    return value


def run_mesh(config, mesh, directory, solver, adapter=None):
    problem = load_problem(config)
    adapter = adapter or CalculiXAdapter()
    verifier = for_problem(problem)
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=False)
    prepared = adapter.prepare(problem, mesh)
    if prepared.audit.get("passed") is not True:
        raise ValueError("Serialized solver input audit failed")
    for name, content in prepared.files.items():
        if Path(name).name != name:
            raise ValueError("Prepared filenames must be safe relative leaf names")
        (directory / name).write_text(content)
    write_json(directory / "input-audit.json", prepared.audit)
    version = adapter.execute(prepared, directory, solver)
    normalized = adapter.extract(prepared, directory, version)
    normalized["verification"] = verifier.verify(problem, normalized)
    normalized["verification"]["checks"]["serialized_input_audit"] = True
    derived_quantities(problem, normalized)
    export_csv(prepared, normalized, directory)
    write_json(directory / "solver-results.json", normalized)
    value = legacy_run(prepared, normalized)
    write_json(directory / "verification.json", value)
    return value


def run_benchmark(config, output, executable, adapter=None):
    problem = load_problem(config)
    adapter = adapter or CalculiXAdapter()
    verifier = for_problem(problem)
    output = Path(output).resolve()
    output.mkdir(parents=True, exist_ok=False)
    write_json(output / "config.json", config)
    identity = provenance(executable, adapter)
    write_json(output / "provenance.json", identity)
    result = {
        "schema_version": ARCHIVE_VERSION,
        "status": "failed",
        "config": config,
        "tolerances": verifier.tolerances,
        "runs": [],
    }
    normalized_runs = []
    study = None
    try:
        for counts in config["meshes"]:
            path = output / ("mesh-" + "x".join(map(str, counts)))
            run = run_mesh(config, counts, path, executable, adapter)
            result["runs"].append(run)
            normalized_runs.append(
                json.loads((path / "solver-results.json").read_text())
            )
            if problem.analysis == "modal":
                print(
                    f"Mesh {counts}: {[m['fem_frequency_hz'] for m in run['modes']]}",
                    flush=True,
                )
            else:
                print(
                    f"Mesh {counts}: q={run['thermal']['mean_heat_flux_W_m2']}, Q_out={run['thermal']['surface_outward_heat_flow_W']['x_max']} W",
                    flush=True,
                )
        study = verifier.study(result["runs"])
        result["convergence" if problem.analysis == "modal" else "mesh_consistency"] = (
            study
        )
        study_passed = (
            all(c["passed"] for c in study)
            if isinstance(study, list)
            else study["passed"]
        )
        if all(r["passed"] for r in result["runs"]) and study_passed:
            result["status"] = "passed"
    except (ValueError, RuntimeError, OSError, subprocess.TimeoutExpired) as exc:
        result["error"] = f"{type(exc).__name__}: {exc}"
        print(result["error"], file=sys.stderr)
    contract = build_contract(
        problem, normalized_runs, identity, result["status"], verifier.tolerances, study
    )
    write_json(output / "normalized-results.json", contract)
    write_json(output / "results.json", result)
    (output / "report.md").write_text(
        (modal_report if problem.analysis == "modal" else thermal_report)(result)
    )
    write_json(output / "manifest.json", manifest(output))
    print(f"Benchmark {result['status']}; evidence: {output}")
    return 0 if result["status"] == "passed" else 1


def replay(directory, adapter=None):
    directory = Path(directory).resolve()
    verify_manifest(directory)
    original = json.loads((directory / "results.json").read_text())
    version = original.get("schema_version")
    if type(version) is not int or version not in (1, 2, ARCHIVE_VERSION):
        raise ValueError(f"Unsupported evidence archive version: {version}")
    problem = load_problem(original["config"])
    adapter = adapter or CalculiXAdapter()
    verifier = for_problem(problem)
    if version < 3 and problem.analysis != "modal":
        raise ValueError("Legacy evidence versions support modal analysis only")
    if original["status"] != "passed":
        raise ValueError("Archived scientific checks failed")
    if json.loads((directory / "config.json").read_text()) != problem.parameters:
        raise ValueError("Archived physical specifications disagree")
    if len(original["runs"]) != len(problem.parameters["meshes"]):
        raise ValueError("Incomplete archived mesh runs")
    identity = json.loads((directory / "provenance.json").read_text())
    if identity.get("solver_name", adapter.name) != adapter.name:
        raise ValueError("Archived solver requires an explicit matching adapter")
    normalized_runs = []
    runs = []
    for counts in problem.parameters["meshes"]:
        path = directory / ("mesh-" + "x".join(map(str, counts)))
        prepared = adapter.prepare(problem, counts)
        if prepared.audit.get("passed") is not True:
            raise ValueError("Regenerated serialized input audit failed")
        for name, content in prepared.files.items():
            if (path / name).read_text() != content:
                raise ValueError("Archived deck differs from reproducible input")
        if json.loads((path / "input-audit.json").read_text()) != prepared.audit:
            raise ValueError("Archived input audit differs")
        solver_version = original["runs"][len(runs)]["solver_version"]
        normalized = adapter.extract(prepared, path, solver_version)
        normalized["verification"] = verifier.verify(problem, normalized)
        normalized["verification"]["checks"]["serialized_input_audit"] = True
        derived_quantities(problem, normalized)
        run = legacy_run(prepared, normalized)
        if not run["passed"]:
            raise ValueError("Archived scientific checks failed")
        stored = original["runs"][len(runs)]
        keys = (
            ("modes", "mesh", "audit", "solver_version")
            if problem.analysis == "modal"
            else ("thermal", "mesh", "audit", "solver_version")
        )
        if any(run[k] != stored[k] for k in keys):
            raise ValueError("Stored verification does not match recomputed checks")
        if version >= 2 and run["frd"] != stored["frd"]:
            raise ValueError(
                "Stored FRD verification differs from recomputed semantics"
            )
        if (
            version >= 3
            and json.loads((path / "solver-results.json").read_text()) != normalized
        ):
            raise ValueError("Normalized mesh result differs from raw solver evidence")
        normalized_runs.append(normalized)
        runs.append(run)
    study = verifier.study(runs)
    key = "convergence" if problem.analysis == "modal" else "mesh_consistency"
    if study != original[key] or not (
        all(v["passed"] for v in study) if isinstance(study, list) else study["passed"]
    ):
        raise ValueError("Archived mesh study failed or differs")
    if version >= 3:
        saved = json.loads((directory / "normalized-results.json").read_text())
        validate_contract(saved)
        expected = build_contract(
            problem, normalized_runs, identity, "passed", verifier.tolerances, study
        )
        if saved != expected or original["tolerances"] != verifier.tolerances:
            raise ValueError(
                "Common contract, tolerances or verification differs from raw evidence"
            )
    print(
        "Evidence hashes, input audits, DAT/FRD semantics, normalized fields, independent checks and mesh study verified"
    )


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--config", type=Path, default=Path("benchmarks/cantilever.json")
    )
    parser.add_argument("--output", type=Path, default=Path("evidence/benchmark"))
    parser.add_argument("--solver", default="ccx")
    parser.add_argument("--verify-evidence", type=Path)
    args = parser.parse_args(argv)
    if args.verify_evidence:
        replay(args.verify_evidence)
        return 0
    config = json.loads(args.config.read_text())
    load_problem(config)
    solver = shutil.which(args.solver)
    if not solver:
        raise RuntimeError("Real solver executable required; install calculix-ccx")
    return run_benchmark(config, args.output, solver)
