"""Compare common Q1-based definitions after both native archives pass replay."""

import argparse
import copy
import json
import math
from pathlib import Path

from .adapters.calculix import CalculiXAdapter
from .adapters.elmer import ElmerAdapter
from .comparison import matching
from .convergence_plot import write as plot_write, replay as plot_replay
from .evidence import manifest, write_json, verify_manifest, sha256
from .results import validate_contract, mesh_maps, thermal_fields
from .runner import replay
from .source_quantities import add_quantities
from .specification import load_problem
from .thermal_problem import sourced
from .verification import for_problem

TOLERANCES = {
    "temperature_absolute_K": 2e-3,
    "reaction_relative": 4e-6,
    "coordinate_absolute_m": 1e-12,
}
NORMS = (
    "temperature_L2_K_m3_2",
    "gradient_H1_K_m1_2",
    "raw_flux_L2_W_m_minus1_2",
    "native_flux_L2_W_m_minus1_2",
    "recovered_flux_L2_W_m_minus1_2",
)


def compare_contracts(left, right):
    for c in (left, right):
        validate_contract(c)
        if c["verification"]["status"] != "passed":
            raise ValueError("Independent source verifications must pass first")
    if left["problem"] != right["problem"] or not sourced(
        left["problem"]["parameters"]
    ):
        raise ValueError("Identical sourced physical specifications required")
    if [left["solver"]["name"], right["solver"]["name"]] != ["CalculiX", "Elmer"]:
        raise ValueError("Two independent adapters required")
    problem = load_problem(left["problem"]["parameters"])
    verifier = for_problem(problem)
    rows = []
    series = {"CalculiX": [], "Elmer": []}
    studies = {}
    for ra, rb in zip(left["runs"], right["runs"]):
        counts = ra["mesh"]["counts"]
        if counts != rb["mesh"]["counts"]:
            raise ValueError("Physical mesh resolutions differ")
        checked = []
        for c, run in zip((left, right), (ra, rb)):
            v = verifier.verify(problem, run)
            v["checks"]["serialized_input_audit"] = True
            if v["status"] != "passed" or v != run["verification"]:
                raise ValueError(
                    "Physical or stored convergence metrics fail independent recomputation"
                )
            # Derived common fields also must be reconstructed from actual DOFs.
            expected = copy.deepcopy(run)
            add_quantities(problem, expected)
            for name in (
                "raw_temperature_gradient",
                "raw_fourier_heat_flux",
                "heat_flow",
                "boundary_reaction_heat_flow",
            ):
                if run["quantities"].get(name) != expected["quantities"][name]:
                    raise ValueError(
                        "Derived physical field contradicts numerical DOFs"
                    )
            checked.append(v["metrics"])
            series[c["solver"]["name"]].append(v["metrics"])
        maps = []
        for run in (ra, rb):
            nodes, _ = mesh_maps(run)
            data = thermal_fields(run)
            T = {nodes[n]: (v,) for n, v in data["temperature"].items()}
            q = {
                tuple(r["coordinates"]): tuple(r["value"])
                for r in run["quantities"]["raw_fourier_heat_flux"]["values"]
            }
            maps.append((T, q))
        pairs = matching(maps[0][0], maps[1][0], TOLERANCES["coordinate_absolute_m"])
        diffT = max(abs(a[0] - b[0]) for a, b in pairs)
        pairsq = matching(maps[0][1], maps[1][1], TOLERANCES["coordinate_absolute_m"])
        diffq = max(
            math.sqrt(sum((x - y) ** 2 for x, y in zip(a, b))) for a, b in pairsq
        )
        g = problem.parameters["geometry"]
        V = g["length"] * g["width"] * g["height"]
        k = problem.parameters["material"]["conductivity"]
        scales = [g[name] / n for name, n in zip(("length", "width", "height"), counts)]
        boundT = diffT * math.sqrt(V)
        boundG = 2 * diffT * math.sqrt(sum(1 / d**2 for d in scales)) * math.sqrt(V)
        normdiff = {name: checked[1][name] - checked[0][name] for name in NORMS}
        reactiondiff = {
            f: checked[1]["boundary_reaction_outward_heat_flow_W"][f] - v
            for f, v in checked[0]["boundary_reaction_outward_heat_flow_W"].items()
        }
        checks = {
            "nodal_coordinate_temperature_agreement": diffT
            <= TOLERANCES["temperature_absolute_K"],
            "common_temperature_norm_triangle_bound": abs(normdiff[NORMS[0]])
            <= boundT + 1e-10,
            "common_gradient_norm_triangle_bound": abs(normdiff[NORMS[1]])
            <= boundG + 1e-10,
            "common_flux_norm_triangle_bound": abs(normdiff[NORMS[2]])
            <= k * boundG + 1e-9,
            "common_raw_flux_field_consistency": diffq
            <= 2 * k * diffT * math.sqrt(sum(1 / d**2 for d in scales)) + 1e-7,
            "boundary_reaction_agreement": max(abs(v) for v in reactiondiff.values())
            / checked[0]["actual_source_integral_W"]
            <= TOLERANCES["reaction_relative"],
        }
        rows.append(
            {
                "counts": counts,
                "h_x_m": checked[0]["h_x_m"],
                "independent_metrics": {"CalculiX": checked[0], "Elmer": checked[1]},
                "temperature_max_coordinate_difference_K": diffT,
                "temperature_pairs": len(pairs),
                "raw_flux_pairs": len(pairsq),
                "raw_flux_max_coordinate_vector_difference_W_m2": diffq,
                "norm_differences_Elmer_minus_CalculiX": normdiff,
                "temperature_L2_triangle_bound_K_m3_2": boundT,
                "gradient_H1_triangle_bound_K_m1_2": boundG,
                "reaction_flow_differences_W": reactiondiff,
                "checks": checks,
                "passed": all(checks.values()),
            }
        )
    for c in (left, right):
        name = c["solver"]["name"]
        study = verifier.study([{"passed": True, "thermal": m} for m in series[name]])
        if not study["passed"] or study != c["verification"]["mesh_study"]:
            raise ValueError("Stored or measured convergence study fails")
        studies[name] = study
    return {
        "schema_name": "open-cae-evidence.source-convergence-comparison",
        "schema_version": "1.0.0",
        "status": "passed" if all(r["passed"] for r in rows) else "failed",
        "problem": left["problem"],
        "tolerances": TOLERANCES,
        "runs": rows,
        "convergence_studies": studies,
        "semantics": "Common raw flux is -k grad(T_h) from genuine DOFs. Native IP flux and differently recovered flux remain separate diagnostics. Coordinates, not native IDs, define matching.",
    }


def compare_archives(calculix, elmer):
    replay(calculix, CalculiXAdapter())
    replay(elmer, ElmerAdapter())
    files = [p / "normalized-results.json" for p in (calculix, elmer)]
    r = compare_contracts(*(json.loads(p.read_text()) for p in files))
    r["source_contract_sha256"] = [sha256(p) for p in files]
    r["source_manifest_sha256"] = [
        sha256(p / "manifest.json") for p in (calculix, elmer)
    ]
    r["native_archives_replayed"] = True
    return r


def report(r):
    lines = [
        "# Source-driven CalculiX / Elmer convergence comparison",
        "",
        f"Status: **{r['status']}**",
        "",
        r["semantics"],
        "",
        "Both individual physical gates, source-aware conservation and native replay passed before comparison.",
        "![Measured convergence](convergence.svg)",
        "",
        "| Mesh | CalculiX T L2 | Elmer T L2 | CalculiX H1 | Elmer H1 | CalculiX recovered flux L2 | Elmer recovered flux L2 | Max DOF temperature difference K |",
        "|---|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for row in r["runs"]:
        a, b = (row["independent_metrics"][name] for name in ("CalculiX", "Elmer"))
        lines.append(
            f"| {'x'.join(map(str,row['counts']))} | {a[NORMS[0]]:.9g} | {b[NORMS[0]]:.9g} | {a[NORMS[1]]:.9g} | {b[NORMS[1]]:.9g} | {a[NORMS[-1]]:.9g} | {b[NORMS[-1]]:.9g} | {row['temperature_max_coordinate_difference_K']:.7g} |"
        )
    lines += [
        "",
        "Units: T L2 K m^(3/2), H1 seminorm K m^(1/2), flux L2 W m^(-1/2). All machine-readable norms retain full computed precision.",
        "",
        "## Actual observed orders",
        "",
        "```json",
        json.dumps(r["convergence_studies"], indent=2),
        "```",
        "",
        "Nodal coincidence is not convergence. Volume quadrature measures the nonquadratic Q1 field, and off-center",
        "points prevent the center-only cancellation of this quadratic-reference gradient error. Native/recovered",
        "flux differences are not averaged away or required to mimic the raw gradient order. Reaction balances",
        "are independent of one-sided boundary gradient reconstruction, whose coarse errors remain visible.",
        "CalculiX prints E13.6 temperature; its rounding is differentiated in reconstructed gradients.",
        "Elmer ASCII temperature has eight significant digits. Quantization bounds and null precision-limited orders",
        "are explicit. Agreement corroborates independent verification, not industrial qualification or experimental validation.",
        "",
    ]
    return "\n".join(lines)


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    for name in ("calculix", "elmer", "output"):
        p.add_argument("--" + name, type=Path, required=True)
    p.add_argument("--replay", action="store_true")
    args = p.parse_args(argv)
    result = compare_archives(args.calculix, args.elmer)
    series = {
        name: [r["independent_metrics"][name] for r in result["runs"]]
        for name in ("CalculiX", "Elmer")
    }
    if args.replay:
        verify_manifest(args.output)
        if json.loads((args.output / "comparison.json").read_text()) != result or (
            args.output / "report.md"
        ).read_text() != report(result):
            raise ValueError("Sourced comparison/metrics/report differs from replay")
        plot_replay(args.output, series)
    else:
        args.output.mkdir(parents=True, exist_ok=False)
        write_json(args.output / "comparison.json", result)
        (args.output / "report.md").write_text(report(result))
        plot_write(args.output, series)
        write_json(args.output / "manifest.json", manifest(args.output))
    print("Source convergence comparison " + result["status"])
    return 0 if result["status"] == "passed" else 1


if __name__ == "__main__":
    raise SystemExit(main())
