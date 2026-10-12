"""Independent slab verification on physical mesh/fields, without CalculiX imports."""

import itertools
import math

from .thermal_problem import analytical

TOLERANCES = {
    "temperature_absolute_K": 1e-3,
    "gradient_relative": 1e-5,
    "flux_relative": 2e-6,
    "heat_flow_relative": 2e-6,
    "conservation_relative": 2e-6,
    "mesh_consistency_relative": 2e-6,
}
SIGNS = (
    (-1, -1, -1),
    (1, -1, -1),
    (1, 1, -1),
    (-1, 1, -1),
    (-1, -1, 1),
    (1, -1, 1),
    (1, 1, 1),
    (-1, 1, 1),
)


def verify(config, nodes, elements, data):
    ref = analytical(config, 0)
    g = config["geometry"]
    lengths = [g["length"], g["width"], g["height"]]
    qscale = abs(ref["heat_flux_W_m2"][0])
    gradscale = abs(ref["gradient_K_m"][0])
    flow = ref["heat_flow_W"]
    expected = {(e, p) for e in elements for p in range(1, 9)}
    if (
        set(data)
        != {"temperature", "heat_rate", "heat_flux", "integration_coordinates"}
        or set(data["temperature"]) != set(nodes)
        or set(data["heat_rate"]) != set(nodes)
        or set(data["heat_flux"]) != expected
        or set(data["integration_coordinates"]) != expected
    ):
        raise ValueError("Incomplete normalized thermal fields")
    all_values = [
        *data["temperature"].values(),
        *data["heat_rate"].values(),
        *(v for q in data["heat_flux"].values() for v in q),
        *(v for xyz in data["integration_coordinates"].values() for v in xyz),
    ]
    if not all(math.isfinite(v) for v in all_values):
        raise ValueError("Nonfinite normalized thermal result")
    temp_error = max(
        abs(t - analytical(config, nodes[n][0])["temperature_K"])
        for n, t in data["temperature"].items()
    )
    flux_error = max(
        abs(q[a] - ref["heat_flux_W_m2"][a]) / qscale
        for q in data["heat_flux"].values()
        for a in range(3)
    )
    gradients = []
    flux_closure = 0
    gradient_error = 0
    surfaces = {
        f"{axis}_{side}": 0.0 for axis in ("x", "y", "z") for side in ("min", "max")
    }
    root = 1 / math.sqrt(3)
    for e, conn in elements.items():
        start = nodes[conn[0]]
        end = nodes[conn[6]]
        span = [b - a for a, b in zip(start, end)]
        if any(v <= 0 for v in span):
            raise ValueError("Invalid normalized thermal element")
        samples = {}
        for p in range(1, 9):
            key = e, p
            xyz = data["integration_coordinates"][key]
            local = [2 * (x - a) / d - 1 for x, a, d in zip(xyz, start, span)]
            signs = tuple(1 if v > 0 else -1 for v in local)
            if signs in samples or any(
                not math.isclose(abs(v), root, rel_tol=3e-5, abs_tol=1e-6)
                for v in local
            ):
                raise ValueError(
                    "Integration coordinates do not cover the eight Gauss points"
                )
            samples[signs] = data["heat_flux"][key]
            # Q1 FE derivative computed from actual nodal output, independently
            # of both solver HFL and the closed-form reference.
            grad = []
            exact_local = [s * root for s in signs]
            for a in range(3):
                derivative = 0.0
                for n, s in zip(conn, SIGNS):
                    basis = (
                        s[a]
                        / 8
                        * math.prod(
                            1 + s[b] * exact_local[b] for b in range(3) if b != a
                        )
                        * 2
                        / span[a]
                    )
                    derivative += basis * data["temperature"][n]
                grad.append(derivative)
            gradients.append(grad)
            gradient_error = max(
                gradient_error,
                max(
                    abs(grad[a] - ref["gradient_K_m"][a]) / gradscale for a in range(3)
                ),
            )
            flux_closure = max(
                flux_closure,
                max(
                    abs(
                        samples[signs][a] + config["material"]["conductivity"] * grad[a]
                    )
                    / qscale
                    for a in range(3)
                ),
            )
        # Integrate actual solver flux on every exterior face. Tensor-product
        # interpolation from the 8 Gauss points to 2x2 face quadrature points
        # is exact for Q1 gradients on these affine rectangular elements.
        for a, axis in enumerate(("x", "y", "z")):
            area = math.prod(span[b] for b in range(3) if b != a)
            for side, sign, pos in (("min", -1, start[a]), ("max", 1, end[a])):
                bound = 0 if sign < 0 else lengths[a]
                if not math.isclose(pos, bound, abs_tol=1e-12):
                    continue
                q_integral = 0.0
                other = [b for b in range(3) if b != a]
                for pair in itertools.product((-root, root), repeat=2):
                    point = [0.0] * 3
                    point[a] = sign
                    for b, v in zip(other, pair):
                        point[b] = v
                    value = sum(
                        q[a]
                        * math.prod((1 + s[b] * point[b] / root) / 2 for b in range(3))
                        for s, q in samples.items()
                    )
                    q_integral += sign * value * area / 4
                surfaces[f"{axis}_{side}"] += q_integral
    rates = data["heat_rate"]
    hot = sum(
        v for n, v in rates.items() if math.isclose(nodes[n][0], 0, abs_tol=1e-12)
    )
    cold = sum(
        v
        for n, v in rates.items()
        if math.isclose(nodes[n][0], g["length"], abs_tol=1e-12)
    )
    free = max(
        (
            abs(v)
            for n, v in rates.items()
            if not math.isclose(nodes[n][0], 0, abs_tol=1e-12)
            and not math.isclose(nodes[n][0], g["length"], abs_tol=1e-12)
        ),
        default=0,
    )
    conservation = max(abs(sum(rates.values())), abs(sum(surfaces.values()))) / flow
    checks = {
        "nodal_temperature": temp_error <= TOLERANCES["temperature_absolute_K"],
        "gradient_direction_magnitude": gradient_error
        <= TOLERANCES["gradient_relative"],
        "solver_flux_direction_magnitude": flux_error <= TOLERANCES["flux_relative"],
        "fourier_closure_from_nodal_output": flux_closure
        <= TOLERANCES["gradient_relative"],
        "integrated_heat_flow": max(
            abs(surfaces["x_min"] + flow), abs(surfaces["x_max"] - flow)
        )
        / flow
        <= TOLERANCES["heat_flow_relative"],
        "adiabatic_sides": max(
            abs(v) for k, v in surfaces.items() if not k.startswith("x")
        )
        / flow
        <= TOLERANCES["heat_flow_relative"],
        "boundary_heat_rates": max(abs(hot - flow), abs(cold + flow), free) / flow
        <= TOLERANCES["heat_flow_relative"],
        "flux_reaction_agreement": max(
            abs(hot + surfaces["x_min"]), abs(cold + surfaces["x_max"])
        )
        / flow
        <= TOLERANCES["heat_flow_relative"],
        "energy_conservation": conservation <= TOLERANCES["conservation_relative"],
    }
    return {
        "passed": all(checks.values()),
        "checks": checks,
        "temperature_max_error_K": temp_error,
        "gradient_max_relative_error": gradient_error,
        "flux_max_relative_error": flux_error,
        "fourier_closure_relative_error": flux_closure,
        "mean_gradient_K_m": [
            sum(v[a] for v in gradients) / len(gradients) for a in range(3)
        ],
        "mean_heat_flux_W_m2": [
            sum(v[a] for v in data["heat_flux"].values()) / len(data["heat_flux"])
            for a in range(3)
        ],
        "surface_outward_heat_flow_W": surfaces,
        "boundary_heat_input_W": {"x_min": hot, "x_max": cold},
        "normalized_energy_imbalance": conservation,
        "reference": ref,
    }


def mesh_consistency(runs):
    values = [r["thermal"]["surface_outward_heat_flow_W"]["x_max"] for r in runs]
    spread = (max(values) - min(values)) / max(abs(v) for v in values)
    return {
        "kind": "mesh_consistency",
        "demonstrated_convergence_order": None,
        "relative_heat_flow_spread": spread,
        "passed": all(r["passed"] for r in runs)
        and spread <= TOLERANCES["mesh_consistency_relative"],
        "reason": "The analytical linear field is in the affine DC3D8/Q1 space; no monotonic error reduction is required.",
    }
