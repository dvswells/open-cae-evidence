"""Corroboration after independent verification and replay of both archives."""

import argparse
import bisect
import json
import math
from pathlib import Path

from .adapters.calculix import CalculiXAdapter
from .adapters.elmer import ElmerAdapter
from .evidence import manifest, sha256, verify_manifest, write_json
from .results import mesh_maps, thermal_fields, validate_contract
from .runner import replay
from .specification import load_problem
from .verification import ThermalVerifier

# Declared before final comparison. Triangle bounds from unchanged individual
# 1e-3 K / 2e-6 relative budgets. Coordinate budget covers CCX DAT 7-digit output.
TOLERANCES = {
    "temperature_absolute_K": 2e-3,
    "flux_relative": 4e-6,
    "heat_flow_relative": 4e-6,
    "energy_balance_difference": 4e-6,
    "coordinate_absolute_m": 1e-7,
}


def matching(left, right, tolerance):
    """Bijective physical-coordinate matching, never native IDs or record order."""
    ordered = sorted(right.items())
    xs = [xyz[0] for xyz, _ in ordered]
    used, pairs = set(), []
    for xyz, value in sorted(left.items()):
        lo, hi = bisect.bisect_left(xs, xyz[0] - tolerance), bisect.bisect_right(
            xs, xyz[0] + tolerance
        )
        candidates = [
            (coord, v)
            for coord, v in ordered[lo:hi]
            if all(abs(a - b) <= tolerance for a, b in zip(xyz, coord))
        ]
        if len(candidates) != 1 or candidates[0][0] in used:
            raise ValueError("Missing, ambiguous or reused physical field coordinate")
        used.add(candidates[0][0])
        pairs.append((value, candidates[0][1]))
    if len(used) != len(right):
        raise ValueError("Cross-solver field coverage differs")
    return pairs


def fields(run):
    nodes, _ = mesh_maps(run)
    data = thermal_fields(run)
    temperatures = {nodes[n]: (t,) for n, t in data["temperature"].items()}
    fluxes = {
        data["integration_coordinates"][key]: q for key, q in data["heat_flux"].items()
    }
    if len(temperatures) != len(nodes) or len(fluxes) != len(data["heat_flux"]):
        raise ValueError("Duplicate physical field locations")
    return temperatures, fluxes


def compare_contracts(left, right):
    for contract in (left, right):
        validate_contract(contract)
        if contract["verification"]["status"] != "passed":
            raise ValueError("Both solvers must independently pass before comparison")
    if (
        left["problem"] != right["problem"]
        or left["problem"]["analysis"] != "steady_heat_conduction"
    ):
        raise ValueError("Comparison requires the same thermal ProblemSpec")
    if {left["solver"]["name"], right["solver"]["name"]} != {"CalculiX", "Elmer"}:
        raise ValueError("Comparison requires two distinct supported solvers")
    problem = load_problem(left["problem"]["parameters"])
    verifier = ThermalVerifier()
    indexed = [
        {tuple(r["mesh"]["counts"]): r for r in c["runs"]} for c in (left, right)
    ]
    runs = []
    for counts in problem.parameters["meshes"]:
        a, b = (v[tuple(counts)] for v in indexed)
        checks = [verifier.verify(problem, r) for r in (a, b)]
        if any(c["status"] != "passed" for c in checks):
            raise ValueError("Actual fields fail independent physical verification")
        af, bf = fields(a), fields(b)
        pairs_t = matching(af[0], bf[0], 1e-12)
        pairs_q = matching(af[1], bf[1], TOLERANCES["coordinate_absolute_m"])
        t_error = max(abs(x[0] - y[0]) for x, y in pairs_t)
        scale_q = max(math.sqrt(sum(v * v for v in x)) for x, _ in pairs_q)
        q_difference = max(
            math.sqrt(sum((v - w) ** 2 for v, w in zip(x, y))) for x, y in pairs_q
        )
        metrics = [c["metrics"] for c in checks]
        scale_Q = abs(metrics[0]["surface_outward_heat_flow_W"]["x_max"])
        flows = [m["surface_outward_heat_flow_W"] for m in metrics]
        flow_difference = {face: flows[1][face] - flows[0][face] for face in flows[0]}
        rates = [m["boundary_heat_input_W"] for m in metrics]
        rate_difference = {face: rates[1][face] - rates[0][face] for face in rates[0]}
        energy = [m["normalized_energy_imbalance"] for m in metrics]
        passed = (
            t_error <= TOLERANCES["temperature_absolute_K"]
            and q_difference / scale_q <= TOLERANCES["flux_relative"]
            and max(
                abs(v) for v in (*flow_difference.values(), *rate_difference.values())
            )
            / scale_Q
            <= TOLERANCES["heat_flow_relative"]
            and abs(energy[0] - energy[1]) <= TOLERANCES["energy_balance_difference"]
        )
        runs.append(
            {
                "counts": counts,
                "temperature_pairs": len(pairs_t),
                "flux_pairs": len(pairs_q),
                "temperature_max_difference_K": t_error,
                "flux_max_vector_difference_W_m2": q_difference,
                "flux_max_relative_difference": q_difference / scale_q,
                "outward_heat_flow_difference_W": flow_difference,
                "boundary_heat_input_difference_W": rate_difference,
                "normalized_energy_imbalances": energy,
                "independent_verification": {
                    c["solver"]["name"]: m for c, m in zip((left, right), metrics)
                },
                "passed": passed,
            }
        )
    return {
        "schema_name": "open-cae-evidence.cross-solver-comparison",
        "schema_version": "1.0.0",
        "status": "passed" if all(r["passed"] for r in runs) else "failed",
        "solvers": [left["solver"]["name"], right["solver"]["name"]],
        "problem": left["problem"],
        "tolerances": TOLERANCES,
        "runs": runs,
        "field_semantics": "Coordinates: bijective matching; CalculiX native IP flux versus Elmer native L2 nodal flux interpolated to the same physical Gauss locations. No temperature interpolation.",
        "limitations": [
            "Cross-solver agreement corroborates independent verification; it is not proof of correctness.",
            "Affine Q1 linear patch has no demonstrated discretization convergence order.",
            "Native integration-point versus recovered nodal flux semantics differ for nonconstant fields.",
            "Printed native output limits distinguishable numerical errors.",
        ],
    }


def compare_archives(calculix, elmer):
    for path, adapter in ((calculix, CalculiXAdapter()), (elmer, ElmerAdapter())):
        replay(path, adapter)
    paths = [calculix / "normalized-results.json", elmer / "normalized-results.json"]
    result = compare_contracts(*(json.loads(p.read_text()) for p in paths))
    result["source_contract_sha256"] = [sha256(p) for p in paths]
    result["source_manifest_sha256"] = [
        sha256(p / "manifest.json") for p in (calculix, elmer)
    ]
    result["native_evidence_replayed"] = True
    return result


def report(result):
    lines = [
        "# CalculiX / Elmer thermal corroboration",
        "",
        f"Status: **{result['status']}**",
        "",
        result["field_semantics"],
        "",
        "Each solver passed the unchanged ThermalVerifier and native archive replay before comparison.",
        "",
        "| Mesh | Max T difference K | Max flux vector difference W/m² | Cold flow difference W | Energy imbalances |",
        "|---|---:|---:|---:|---|",
    ]
    for r in result["runs"]:
        lines.append(
            f"| {'x'.join(map(str,r['counts']))} | {r['temperature_max_difference_K']:.6g} | {r['flux_max_vector_difference_W_m2']:.6g} | {r['outward_heat_flow_difference_W']['x_max']:.6g} | {r['normalized_energy_imbalances']} |"
        )
    lines += [
        "",
        "Acceptance budgets (declared before final comparison):",
        json.dumps(result["tolerances"], indent=2),
        "",
        "Individual analytical errors are retained in comparison.json. Modeling and affine Q1 discretization are shared;",
        "Elmer uses nodal L2 flux recovery and CalculiX writes native Gauss flux. Recovery changes general fields,",
        "but preserves a constant flux in this case up to numerical error. Native ASCII output rounds values to",
        "about seven significant digits; zero printed temperature error does not establish exact arithmetic.",
        "Mesh consistency across three levels is a linear patch check, not a convergence-order measurement.",
        "",
        *result["limitations"],
        "",
    ]
    return "\n".join(lines)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--calculix", type=Path, required=True)
    parser.add_argument("--elmer", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--replay", action="store_true")
    args = parser.parse_args(argv)
    result = compare_archives(args.calculix, args.elmer)
    if args.replay:
        verify_manifest(args.output)
        if json.loads((args.output / "comparison.json").read_text()) != result or (
            args.output / "report.md"
        ).read_text() != report(result):
            raise ValueError("Archived comparison differs from native replay")
    else:
        args.output.mkdir(parents=True, exist_ok=False)
        write_json(args.output / "comparison.json", result)
        (args.output / "report.md").write_text(report(result))
        write_json(args.output / "manifest.json", manifest(args.output))
    print(f"Cross-solver comparison {result['status']}; native evidence replayed")
    return 0 if result["status"] == "passed" else 1


if __name__ == "__main__":
    raise SystemExit(main())
