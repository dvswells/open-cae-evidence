"""Run real CalculiX, verify its output, and retain replayable evidence."""

import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import re
import shutil
import subprocess
import sys

from .model import active_tip, area_weights, audit_deck, make_deck, validate
from .reference import mac, mass_ratio, references, shape
from .frd import verify_frd

# Declared acceptance bounds, not fitted to a stored FEM answer.
TOLERANCES = {"analytical_relative": 0.01, "final_mesh_relative": 0.005,
              "shape_mac_min": 0.99, "direction_energy_min": 0.95,
              "fixed_displacement_relative": 1e-8, "refinement_noise_floor": 5e-5,
              "tip_coupling_relative": 2e-6}


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def parse_dat(text, expected_nodes):
    """Strictly parse eigenvalues and explicitly numbered nodal mode blocks."""
    lines = text.splitlines()
    frequencies, eigenvalues, modes = [], [], []
    in_frequencies = False
    for line in lines:
        if "E I G E N V A L U E   O U T P U T" in line:
            if in_frequencies or frequencies:
                raise ValueError("Duplicate eigenvalue table")
            in_frequencies = True
        if "P A R T I C I P A T I O N" in line:
            in_frequencies = False
        parts = line.split()
        if in_frequencies and len(parts) == 5 and parts[0].isdigit():
            index = int(parts[0])
            lam, omega, frequency, imaginary = map(float, parts[1:])
            if index != len(frequencies) + 1 or not all(math.isfinite(v) for v in (lam, omega, frequency, imaginary)):
                raise ValueError("Malformed eigenvalue table")
            if lam <= 0 or frequency <= 0 or imaginary != 0:
                raise ValueError("Nonpositive/complex mode; check constraints and material")
            if not math.isclose(omega, 2 * math.pi * frequency, rel_tol=2e-6) or not math.isclose(lam, omega**2, rel_tol=2e-6):
                raise ValueError("Eigenvalue/frequency unit consistency failure")
            eigenvalues.append(lam)
            frequencies.append(frequency)
    pending_mode = None
    i = 0
    while i < len(lines):
        line = lines[i]
        if "E I G E N V A L U E    N U M B E R" in line:
            pending_mode = int(line.split()[-1])
        if "displacements (vx,vy,vz) for set ALL" in line:
            if pending_mode != len(modes) + 1:
                raise ValueError("Missing/out-of-order mode shape number")
            displacements = {}
            i += 1
            while i < len(lines) and not lines[i].strip():
                i += 1
            while i < len(lines) and lines[i].strip():
                row = lines[i].split()
                if len(row) != 4 or not row[0].isdigit():
                    raise ValueError("Malformed displacement row")
                n, values = int(row[0]), tuple(float(v) for v in row[1:])
                if n in displacements or not all(math.isfinite(v) for v in values):
                    raise ValueError("Duplicate/nonfinite displacement")
                displacements[n] = values
                i += 1
            if set(displacements) != set(expected_nodes):
                raise ValueError("Incomplete nodal mode shape output")
            modes.append(displacements)
            pending_mode = None
        i += 1
    if len(frequencies) != 6 or len(modes) != 6 or any(b <= a for a, b in zip(frequencies, frequencies[1:])):
        raise ValueError("Exactly six distinct, increasing complete modes required")
    return frequencies, eigenvalues, modes


def evaluate(config, nodes, frequencies, eigenvalues, modes):
    length = config["geometry"]["length"]
    auxiliary = {max(nodes)} if active_tip(config) else set()
    probes = sorted((n for n, (_, y, z) in nodes.items() if n not in auxiliary and abs(y) < 1e-12 and abs(z) < 1e-12), key=lambda n: nodes[n][0])
    if len(probes) < 10 or nodes[probes[0]][0] != 0 or not math.isclose(nodes[probes[-1]][0], length):
        raise ValueError("Insufficient centerline samples")
    records = []
    for index, (reference, frequency, lam, mode) in enumerate(zip(references(config), frequencies, eigenvalues, modes), 1):
        axis = {"y": 1, "z": 2}[reference["direction"]]
        energy = [sum(v[a]**2 for v in mode.values()) for a in range(3)]
        if sum(energy) <= 0 or not math.isfinite(sum(energy)):
            raise ValueError("Mode shape has zero/nonfinite magnitude")
        fraction = energy[axis] / sum(energy)
        peak = max(abs(v) for xyz in mode.values() for v in xyz)
        clamp = max(abs(v) for n, xyz in mode.items() if abs(nodes[n][0]) < 1e-12 for v in xyz) / peak
        samples = [mode[n][axis] for n in probes]
        analytical = [shape(reference["beta"], nodes[n][0] / length) for n in probes]
        assurance = mac(samples, analytical)
        relative = abs(frequency / reference["frequency_hz"] - 1)
        checks = {"analytical_frequency": relative <= TOLERANCES["analytical_relative"],
                  "bending_direction": fraction >= TOLERANCES["direction_energy_min"],
                  "analytical_shape": assurance >= TOLERANCES["shape_mac_min"],
                  "clamped_displacements": clamp <= TOLERANCES["fixed_displacement_relative"]}
        records.append({"mode": index, **reference, "fem_frequency_hz": frequency, "eigenvalue_rad2_s2": lam,
                        "analytical_relative_error": relative, "centerline_mac": assurance,
                        "direction_energy_fraction": fraction, "fixed_relative_displacement": clamp,
                        "checks": checks, "passed": all(checks.values())})
    return records


def verify_modes(config, mesh, nodes, frequencies, eigenvalues, modes):
    records = evaluate(config, nodes, frequencies, eigenvalues, modes)
    if active_tip(config):
        weights = area_weights(nodes, config, mesh)
        for record, mode in zip(records, modes):
            peak = max(abs(v) for xyz in mode.values() for v in xyz)
            residual = max(abs(mode[max(nodes)][a] - sum(w*mode[n][a] for n, w in weights.items())) for a in range(3)) / peak
            record["tip_coupling_relative_residual"] = residual
            # DAT prints seven significant digits. Signed sums need rounding allowance.
            record["checks"]["tip_coupling"] = residual <= TOLERANCES["tip_coupling_relative"]
            record["passed"] = all(record["checks"].values())
    return records


def check_solver_log(log):
    versions = re.findall(r"CalculiX Version ([\d.]+)", log)
    if "*ERROR" in log.upper() or "Job finished" not in log or len(versions) != 1:
        raise RuntimeError("Solver log has errors, is incomplete or lacks a unique version")
    return versions[0]


def run_mesh(config, mesh, directory, solver):
    directory.mkdir(parents=True, exist_ok=False)
    deck, nodes, audit = make_deck(config, mesh)
    (directory / "beam.inp").write_text(deck)
    write_json(directory / "input-audit.json", audit)
    env = os.environ.copy()
    env.update(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", CCX_NPROC_RESULTS="1")
    with (directory / "solver.log").open("w") as log:
        result = subprocess.run([solver, "-i", "beam"], cwd=directory, env=env,
                                stdout=log, stderr=subprocess.STDOUT, timeout=300, check=False)
    log = (directory / "solver.log").read_text(errors="replace")
    if result.returncode != 0 or "*ERROR" in log.upper() or "Job finished" not in log:
        raise RuntimeError(f"CalculiX failed (exit {result.returncode}); see {directory / 'solver.log'}")
    version = check_solver_log(log)
    if not (directory / "beam.frd").is_file() or (directory / "beam.frd").stat().st_size == 0:
        raise ValueError("Missing raw mode-shape FRD output")
    frequencies, eigenvalues, modes = parse_dat((directory / "beam.dat").read_text(), nodes)
    records = verify_modes(config, mesh, nodes, frequencies, eigenvalues, modes)
    frd = verify_frd((directory / "beam.frd").read_text(), deck, nodes, frequencies, eigenvalues, modes,
                     version, {max(nodes)} if active_tip(config) else set())
    for index, mode in enumerate(modes, 1):
        with (directory / f"mode-{index}.csv").open("w", newline="") as handle:
            writer = csv.writer(handle)
            writer.writerow(["node", "x_m", "y_m", "z_m", "ux", "uy", "uz"])
            for n, xyz in sorted(nodes.items()):
                writer.writerow([n, *xyz, *mode[n]])
    value = {"mesh": mesh, "audit": audit, "solver_version": version, "modes": records, "frd": frd,
             "passed": all(r["passed"] for r in records)}
    write_json(directory / "verification.json", value)
    return value


def convergence(runs):
    results = []
    for i in range(6):
        values = [run["modes"][i]["fem_frequency_hz"] for run in runs]
        changes = [abs(b - a) / b for a, b in zip(values, values[1:])]
        passed = (changes[-1] <= TOLERANCES["final_mesh_relative"]
                  and changes[-1] <= max(changes[0], TOLERANCES["refinement_noise_floor"]))
        results.append({"mode": i + 1, "frequencies_hz": values,
                        "relative_changes": changes, "passed": passed})
    return results


def report(result):
    tip = "tip_mass" in result["config"]
    title = "Cantilever with concentrated tip mass" if tip else "Cantilever modal verification"
    equation = ("sech(beta)+cos(beta)+mu*beta*(cos(beta)*tanh(beta)-sin(beta))=0."
                if tip else "cos(beta)*cosh(beta)=-1.")
    lines = [f"# {title}", "", f"Status: **{result['status']}**", "",
             "This is verification of one linear elastic benchmark, not validation against experiments",
             "or evidence of general replacement suitability for a commercial CAE platform.", "",
             "## Model and assumptions", "", "SI: m, kg, s, Pa. Uniform isotropic steel; small undamped free vibrations.",
             "Solid C3D20R quadratic reduced-integration hexahedra; consistent solver mass matrix.",
             "Every x=0 face node is fixed in DOFs 1–3; other faces have no prescribed displacement.",
             "The deck is audited independently after serialization, including material, volume, connectivity and clamp.", "",
             "```json", json.dumps(result["config"], indent=2), "```", "",
             "## Independent checks", "", "Euler–Bernoulli: f = beta²/(2*pi*L²) * sqrt(E*I/(rho*A)); " + equation,
             "Both transverse planes are considered and sorted. Shape checks use the analytical cantilever",
             "eigenfunction at centerline nodes and sign/scale-invariant sampled MAC.",
             "Plane classification uses the fraction of nodal displacement squared in the expected direction;",
             "this is a diagnostic, not a physical modal energy calculation.", "",
             "Acceptance limits are specified before execution: 1% analytical frequency error, MAC >= 0.99,",
             "direction fraction >= 0.95, clamp residual <= 1e-8 of peak displacement, final mesh change <= 0.5%.",
             "The final change must also not exceed the first change (5e-5 relative noise floor).",
             "The 1% bound allows 3D end constraint effects, shear deformation and rotary inertia absent from",
             "slender-beam theory; it is a benchmark acceptance criterion, not a universal error bound.",
             "MAC >= 0.99 permits small cross-section and clamp effects; 0.5% is the selected discretization budget.",
             "Exact fixed displacements should be zero; 1e-8 allows negligible output roundoff.", "",
             "## Calculated frequencies and shapes", "",
             "| Mode | Plane/order | FEM Hz | EB Hz | Error % | MAC | Pass |",
             "|---|---|---:|---:|---:|---:|---|"]
    if result["runs"]:
        for mode in result["runs"][-1]["modes"]:
            lines.append(f"| {mode['mode']} | {mode['direction']}/{mode['order']} | {mode['fem_frequency_hz']:.6f} | {mode['frequency_hz']:.6f} | {100 * mode['analytical_relative_error']:.4f} | {mode['centerline_mac']:.8f} | {mode['passed']} |")
    lines += ["", "## Mesh convergence", "", "| Mode | Frequencies coarse → fine (Hz) | Final change % | Pass |", "|---|---|---:|---|"]
    for item in result.get("convergence", []):
        lines.append(f"| {item['mode']} | {' → '.join(f'{v:.6f}' for v in item['frequencies_hz'])} | {100 * item['relative_changes'][-1]:.5f} | {item['passed']} |")
    lines += ["", "Refinement increases all three element counts. Mode direction and shape checks guard",
              "against matching unrelated eigenmodes by index. No Richardson extrapolation or asymptotic",
              "convergence order is claimed from these nonuniform refinement ratios.", "", "## Evidence and limits", "",
              "Each mesh retains the .inp, .dat, .frd, complete solver log, input audit, verification JSON",
              "and six full nodal displacement CSV files. Eigenvectors have arbitrary sign and amplitude;",
              "CSV displacements are solver normalization values, not physical vibration amplitudes.",
              "manifest.json records SHA-256 digests. provenance.json records solver and Python environment.",
              "Hashes detect alteration, not authorship. Rounded .dat values limit replay precision.",
              "Only six bending modes of this geometry are checked: no torsional/axial reference verification,",
              "experimental validation, damping, prestress, nonlinear response or industrial qualification.",
              "Positive plausible but falsely labeled units cannot be inferred from numbers alone.",
              "See docs/methodology.md for reference derivation, scope and reproducibility instructions."]
    if tip:
        lines += ["", "## Tip inertia and coupling", "",
                  f"Tip/beam mass ratio mu = {mass_ratio(result['config']):.8g}.",
                  "A centered auxiliary MASS node remains independent in three translation DOFs.",
                  "Three MPCs equate its displacement to the exact Q8 area-average of the end face.",
                  "Signed interpolation weights integrate to one and have zero first moments.",
                  "They distribute inertia forces by virtual work without tying every face node rigidly.",
                  "Rotary inertia, eccentricity and attachment stiffness are deliberately neglected.",
                  "Each mode must satisfy the coupling within 2e-6 of peak displacement (DAT rounding).",
                  "CalculiX may warn that its internal MASS material has no elastic constants; this",
                  "inertia-only element is intentionally stiffness-free. All *ERROR messages are fatal.",
                  "See docs/tip-mass.md for the dynamic tip boundary condition and mass-scaling checks."]
        lines += ["", "| Mesh | Total mass kg | Coupled face nodes | Maximum MPC residual |",
                  "|---|---:|---:|---:|"]
        for run in result["runs"]:
            if "tip_mass_kg" in run["audit"]:
                lines.append(f"| {'x'.join(map(str, run['mesh']))} | {run['audit']['total_mass_kg']:.8g} | {run['audit']['coupled_face_nodes']} | {max(m['tip_coupling_relative_residual'] for m in run['modes']):.3e} |")
        lines += ["", "Required CI integration tests separately exercise tip-only and total-mass scaling,",
                  "modulus scaling, the zero-mass limit, and equivalent solid-face elimination pivots.",
                  "Their real inputs, outputs and six-mode ratio JSON files are in evidence/tests/tip-mass-*/.",
                  "This standalone report does not assert that those separate tests have executed."]
    lines += ["", "## Semantic FRD verification", "",
              "For every mesh, ASCII FRD geometry, HEX20 connectivity, solver version, six DISP datasets,",
              "mode IDs, frequencies, eigenvalues and all signed displacement components are checked",
              "against the audited deck and DAT. The component comparison uses 1e-5 relative tolerance",
              "and 1e-12 of peak displacement absolute tolerance for FRD/DAT printed precision.",
              "MASS elements and their auxiliary nodes are omitted by the FRD exporter; their full",
              "mode displacements are preserved in DAT/CSV and independently checked by the MPC residual."]
    if result.get("error"):
        lines += ["", "## Execution failure", "", result["error"]]
    return "\n".join(lines) + "\n"


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def replay(directory):
    """Recompute scientific checks from archived solver data without running AI."""
    manifest = json.loads((directory / "manifest.json").read_text())
    files = {str(p.relative_to(directory)) for p in directory.rglob("*") if p.is_file() and p.name != "manifest.json"}
    if files != set(manifest):
        raise ValueError("Evidence file set differs from manifest")
    for name, digest in manifest.items():
        if sha256(directory / name) != digest:
            raise ValueError(f"Evidence hash mismatch: {name}")
    original = json.loads((directory / "results.json").read_text())
    config = original["config"]
    validate(config)
    runs = []
    for mesh in config["meshes"]:
        path = directory / ("mesh-" + "x".join(map(str, mesh)))
        deck, nodes, _ = make_deck(config, mesh)
        actual = (path / "beam.inp").read_text()
        audit_deck(actual, config, mesh)
        if actual != deck:
            raise ValueError("Archived deck differs from reproducible input")
        f, lam, modes = parse_dat((path / "beam.dat").read_text(), nodes)
        version = check_solver_log((path / "solver.log").read_text())
        frd = verify_frd((path / "beam.frd").read_text(), actual, nodes, f, lam, modes,
                         version, {max(nodes)} if active_tip(config) else set())
        runs.append({"modes": verify_modes(config, mesh, nodes, f, lam, modes), "frd": frd})
    if original["status"] != "passed" or not all(m["passed"] for run in runs for m in run["modes"]) or not all(c["passed"] for c in convergence(runs)):
        raise ValueError("Archived scientific checks failed")
    if [r["modes"] for r in runs] != [r["modes"] for r in original["runs"]] or convergence(runs) != original["convergence"]:
        raise ValueError("Stored verification does not match recomputed checks")
    if original.get("schema_version", 1) >= 2 and [r["frd"] for r in runs] != [r["frd"] for r in original["runs"]]:
        raise ValueError("Stored FRD verification differs from recomputed semantics")
    print("Evidence hashes, decks, DAT/FRD semantics, six modes per mesh, analytical checks, coupling and convergence verified")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=Path("benchmarks/cantilever.json"))
    parser.add_argument("--output", type=Path, default=Path("evidence/benchmark"))
    parser.add_argument("--solver", default="ccx")
    parser.add_argument("--verify-evidence", type=Path)
    args = parser.parse_args(argv)
    if args.verify_evidence:
        replay(args.verify_evidence.resolve())
        return 0
    config = json.loads(args.config.read_text())
    validate(config)
    solver = shutil.which(args.solver)
    if not solver:
        raise RuntimeError("Real CalculiX executable required; install calculix-ccx")
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)  # No stale result reuse or overwriting evidence.
    write_json(output / "config.json", config)
    revision = subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True, text=True)
    dirty = subprocess.run(["git", "status", "--porcelain"], capture_output=True, text=True)
    provenance = {"python": sys.version, "platform": platform.platform(), "solver_path": solver,
                  "solver_binary_sha256": sha256(Path(solver)), "git_revision": revision.stdout.strip() or None,
                  "working_tree_dirty": bool(dirty.stdout.strip()), "github_run_id": os.environ.get("GITHUB_RUN_ID"),
                  "github_sha": os.environ.get("GITHUB_SHA"), "threads": 1}
    write_json(output / "provenance.json", provenance)
    result = {"schema_version": 2, "status": "failed", "config": config, "tolerances": TOLERANCES, "runs": []}
    try:
        for mesh in config["meshes"]:
            run = run_mesh(config, mesh, output / ("mesh-" + "x".join(map(str, mesh))), solver)
            result["runs"].append(run)
            print(f"Mesh {mesh}: {[m['fem_frequency_hz'] for m in run['modes']]}", flush=True)
        result["convergence"] = convergence(result["runs"])
        if all(r["passed"] for r in result["runs"]) and all(c["passed"] for c in result["convergence"]):
            result["status"] = "passed"
    except (ValueError, RuntimeError, OSError, subprocess.TimeoutExpired) as exc:
        result["error"] = f"{type(exc).__name__}: {exc}"
        print(result["error"], file=sys.stderr)
    write_json(output / "results.json", result)
    (output / "report.md").write_text(report(result))
    write_json(output / "manifest.json", {str(p.relative_to(output)): sha256(p) for p in sorted(output.rglob("*")) if p.is_file()})
    print(f"Benchmark {result['status']}; evidence: {output}")
    return 0 if result["status"] == "passed" else 1


if __name__ == "__main__":
    raise SystemExit(main())
