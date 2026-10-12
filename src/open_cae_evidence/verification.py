"""Numerical verification receives normalized physical fields, never solver files."""

import math

from . import modal_verification, thermal_verification
from .results import (
    mesh_maps,
    modal_fields,
    thermal_fields,
    quantity,
    require_quantities,
)


class ModalVerifier:
    tolerances = modal_verification.TOLERANCES

    def verify(self, problem, normalized):
        require_quantities(problem.analysis, normalized)
        nodes, _ = mesh_maps(normalized)
        f, lam, modes = modal_fields(normalized)
        auxiliary = [
            n["id"] for n in normalized["mesh"]["nodes"] if n["role"] == "point_mass"
        ]
        active = problem.parameters.get("tip_mass", {}).get("mass", 0) > 0
        if len(auxiliary) != (1 if active else 0):
            raise ValueError("Modal mass-node roles differ from physical problem")
        # Existing verified Q8 transfer routines use a highest-ID reference.
        # Adapt arbitrary external labels locally; preserve public result IDs.
        if active and auxiliary[0] != max(nodes):
            old = auxiliary[0]
            new = max(nodes) + 1
            nodes[new] = nodes.pop(old)
            modes = [
                {new if n == old else n: v for n, v in mode.items()} for mode in modes
            ]
        if any(not math.isfinite(v) or v <= 0 for v in (*f, *lam)) or any(
            b <= a for a, b in zip(f, f[1:])
        ):
            raise ValueError(
                "Normalized eigenvalues/frequencies must be finite, positive and ordered"
            )
        if any(
            not math.isclose(l, (2 * math.pi * v) ** 2, rel_tol=2e-6)
            for l, v in zip(lam, f)
        ):
            raise ValueError("Normalized eigenvalue/frequency units are inconsistent")
        records = modal_verification.verify_modes(
            problem.parameters, normalized["mesh"]["counts"], nodes, f, lam, modes
        )
        checks = {
            f"mode_{r['mode']}_{name}": passed
            for r in records
            for name, passed in r["checks"].items()
        }
        return {
            "status": "passed" if all(checks.values()) else "failed",
            "checks": checks,
            "metrics": records,
        }

    def study(self, runs):
        return modal_verification.convergence(runs)


class ThermalVerifier:
    tolerances = thermal_verification.TOLERANCES

    def verify(self, problem, normalized):
        require_quantities(problem.analysis, normalized)
        nodes, elements = mesh_maps(normalized)
        data = thermal_fields(normalized)
        result = thermal_verification.verify(problem.parameters, nodes, elements, data)
        return {
            "status": "passed" if result["passed"] else "failed",
            "checks": result["checks"],
            "metrics": result,
        }

    def study(self, runs):
        return thermal_verification.mesh_consistency(runs)


def for_problem(problem):
    return (
        ThermalVerifier()
        if problem.analysis == "steady_heat_conduction"
        else ModalVerifier()
    )


def derived_quantities(problem, normalized):
    if problem.analysis != "steady_heat_conduction":
        return
    metrics = normalized["verification"]["metrics"]
    normalized["quantities"]["heat_flow"] = quantity(
        "W",
        "boundary",
        ["outward_heat_flow"],
        [
            {"face": face, "value": v}
            for face, v in metrics["surface_outward_heat_flow_W"].items()
        ],
        "postprocessed_from_solver",
        normalized["quantities"]["heat_flux"]["source"]["files"],
        "Tensor Gauss extrapolation/interpolation and 2x2 face integration of actual solver IP flux on affine HEX8 cells; positive is outward.",
    )
