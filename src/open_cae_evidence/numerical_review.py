"""Reproducible adversarial audit of apparently excellent sourced-slab evidence."""

import argparse
import json
import math
from pathlib import Path
from . import q1, source_reference
from .evidence import write_json, manifest, verify_manifest, sha256
from .results import mesh_maps, thermal_fields
from .source_comparison import compare_archives
from .adapters.elmer import ElmerAdapter
from .adapters.elmer_output import check_log
from .specification import load_problem
from .verification import for_problem


def review(calculix, elmer, failure_root):
    comparison = compare_archives(calculix, elmer)
    entries = []
    for name, path in (("CalculiX", calculix), ("Elmer", elmer)):
        c = json.loads((path / "normalized-results.json").read_text())
        config = c["problem"]["parameters"]
        run = c["runs"][0]
        nodes, elements = mesh_maps(run)
        data = thermal_fields(run)
        wrong_center_integral = 0.0
        for conn in elements.values():
            start, span = q1.cell(nodes, conn)
            _, grad = q1.evaluate(
                tuple(data["temperature"][n] for n in conn), span, (0.0, 0.0, 0.0)
            )
            exact = source_reference.analytical(config, start[0] + span[0] / 2)[
                "gradient_K_m"
            ]
            wrong_center_integral += sum(
                (a - b) ** 2 for a, b in zip(grad, exact)
            ) * math.prod(span)
        metrics = [r["verification"]["metrics"] for r in c["runs"]]
        checks = {
            "volume_error_not_confused_with_nodal_error": metrics[0][
                "temperature_L2_K_m3_2"
            ]
            > 100 * metrics[0]["precision_floors"]["temperature_L2_K_m3_2"],
            "center_only_false_zero_detected": math.sqrt(wrong_center_integral) < 1e-8
            and metrics[0]["gradient_H1_K_m1_2"] > 1,
            "full_27_point_volume_quadrature": all(
                m["volume_quadrature_samples"] == 27 * len(r["mesh"]["elements"])
                for m, r in zip(metrics, c["runs"])
            ),
            "raw_boundary_flux_is_not_mislabelled_reaction": abs(
                metrics[0]["raw_net_outward_W"] - metrics[0]["reaction_net_outward_W"]
            )
            > 100,
            "measured_orders_resolved_and_supported": c["verification"]["mesh_study"][
                "passed"
            ],
            "native_and_raw_fields_separate": all(
                r["quantities"]["raw_fourier_heat_flux"]["source"]["kind"]
                == "postprocessed_from_solver"
                for r in c["runs"]
            ),
        }
        entries.append(
            {
                "solver": name,
                "checks": checks,
                "coarse_nodal_error_K": metrics[0]["temperature_max_error_K"],
                "coarse_volume_temperature_L2_K_m3_2": metrics[0][
                    "temperature_L2_K_m3_2"
                ],
                "incorrect_center_only_gradient_norm": math.sqrt(wrong_center_integral),
                "correct_volume_gradient_H1": metrics[0]["gradient_H1_K_m1_2"],
                "coarse_reaction_net_W": metrics[0]["reaction_net_outward_W"],
                "coarse_raw_flux_net_W": metrics[0]["raw_net_outward_W"],
                "observed_rates": c["verification"]["mesh_study"]["observed_rates"],
            }
        )
    failed_directory = (
        failure_root / "test_real_wrong_source_converges_conserves_but_fails_reference"
    )
    actual = load_problem(
        json.loads((failed_directory / "actual-physical-spec.json").read_text())
    )
    expected = load_problem(
        json.loads((failed_directory / "expected-physical-spec.json").read_text())
    )
    adapter = ElmerAdapter()
    prepared = adapter.prepare(actual, actual.parameters["meshes"][0])
    version = check_log((failed_directory / "solver.log").read_text())
    numerical = adapter.extract(prepared, failed_directory, version)
    bad = for_problem(expected).verify(expected, numerical)
    if bad != json.loads((failed_directory / "failure-verification.json").read_text()):
        raise ValueError(
            "Wrong-source failure evidence differs from native numerical replay"
        )
    wrong_source = (
        bad["status"] == "failed"
        and bad["checks"]["reaction_source_energy_conservation"]
        and not bad["checks"]["audited_source_matches_specification"]
    )
    return {
        "review_version": "1.0.0",
        "status": (
            "passed"
            if wrong_source
            and comparison["status"] == "passed"
            and all(all(e["checks"].values()) for e in entries)
            else "failed"
        ),
        "review_scope": "Adversarial numerical evidence audit, not an external human endorsement or industrial qualification",
        "entries": entries,
        "converged_balanced_wrong_source_rejected": wrong_source,
        "wrong_source_native_output_reextracted": True,
        "failure_native_files_sha256": {
            p.name: sha256(p)
            for p in sorted(failed_directory.iterdir())
            if p.is_file()
        },
        "conclusion": "Excellent nodal values result from this uniform Q1 nodal patch property. Genuine volume T, gradient and raw-flux errors are nonzero and converge; recovered flux and reaction estimates retain distinct semantics.",
    }


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    for name in ("calculix", "elmer", "failure-root", "output"):
        p.add_argument("--" + name, type=Path, required=True)
    p.add_argument("--replay", action="store_true")
    args = p.parse_args(argv)
    r = review(args.calculix, args.elmer, args.failure_root)
    text = (
        "# Adversarial numerical evidence review\n\n"
        + r["review_scope"]
        + "\n\nStatus: **"
        + r["status"]
        + "**\n\n"
        + r["conclusion"]
        + "\n\n```json\n"
        + json.dumps(r, indent=2)
        + "\n```\n"
    )
    if args.replay:
        verify_manifest(args.output)
        if (
            json.loads((args.output / "review.json").read_text()) != r
            or (args.output / "report.md").read_text() != text
        ):
            raise ValueError("Adversarial review differs from native evidence")
    else:
        args.output.mkdir(parents=True, exist_ok=False)
        write_json(args.output / "review.json", r)
        (args.output / "report.md").write_text(text)
        write_json(args.output / "manifest.json", manifest(args.output))
    print("Adversarial numerical review " + r["status"])
    return 0 if r["status"] == "passed" else 1


if __name__ == "__main__":
    raise SystemExit(main())
