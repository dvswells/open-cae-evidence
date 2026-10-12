"""Source-aware physical checks and measured volume errors; solver-independent.

Budgets fixed before final convergence solves. Source-free gates are untouched.
"""

import math
from . import q1, source_reference
from .results import mesh_maps, thermal_fields

TOLERANCES = {
    "nodal_temperature_absolute_K": 1e-3,
    "reaction_relative": 2e-6,
    "conservation_relative": 2e-6,
    "native_flux_consistency_relative": 2e-6,
    "source_integral_relative": 1e-12,
    "temperature_order_deviation": 0.1,
    "gradient_order_deviation": 0.1,
    "precision_floor_multiplier": 10.0,
}


def nodal_vectors(quantity, nodes):
    if (
        quantity["unit"] != "W/m^2"
        or quantity["location"] != "node"
        or quantity["components"] != ["x", "y", "z"]
    ):
        raise ValueError("Recovered flux units/location/components differ")
    rows = {r["node"]: tuple(r["value"]) for r in quantity["values"]}
    if (
        len(rows) != len(quantity["values"])
        or set(rows) != set(nodes)
        or any(len(v) != 3 or not all(map(math.isfinite, v)) for v in rows.values())
    ):
        raise ValueError("Incomplete/nonfinite recovered native flux")
    return rows


def verify(config, normalized):
    nodes, elements = mesh_maps(normalized)
    counts = normalized["mesh"]["counts"]
    lengths = [config["geometry"][name] for name in ("length", "width", "height")]
    if (
        len(nodes) != math.prod(n + 1 for n in counts)
        or len(elements) != math.prod(counts)
        or any(n["role"] != "solid" for n in normalized["mesh"]["nodes"])
        or any(
            e["topology"] != "hexahedron8" or e["integration"] != "full"
            for e in normalized["mesh"]["elements"]
        )
    ):
        raise ValueError("Sourced slab requires the complete physical HEX8 mesh")
    lattice = set()
    for xyz in nodes.values():
        key = tuple(round(x / l * c) for x, l, c in zip(xyz, lengths, counts))
        if key in lattice or any(
            i < 0 or i > c or not math.isclose(x, l * i / c, rel_tol=0, abs_tol=1e-12)
            for i, x, l, c in zip(key, xyz, lengths, counts)
        ):
            raise ValueError(
                "Normalized thermal coordinates, units or node coverage differ"
            )
        lattice.add(key)
    cells = set()
    for conn in elements.values():
        start, span = q1.cell(nodes, conn)
        key = tuple(round(x / l * c) for x, l, c in zip(start, lengths, counts))
        if (
            key in cells
            or any(i < 0 or i >= c for i, c in zip(key, counts))
            or any(
                not math.isclose(d, l / c, rel_tol=1e-12)
                for d, l, c in zip(span, lengths, counts)
            )
        ):
            raise ValueError(
                "Overlapping, incomplete or nonuniform normalized source cells"
            )
        cells.add(key)
    data = thermal_fields(normalized)
    if (
        set(data["temperature"]) != set(nodes)
        or set(data["heat_rate"]) != set(nodes)
        or set(data["heat_flux"]) != {(e, p) for e in elements for p in range(1, 9)}
    ):
        raise ValueError("Incomplete sourced numerical fields")
    if not all(
        math.isfinite(v)
        for key in ("temperature", "heat_rate")
        for v in data[key].values()
    ) or not all(math.isfinite(v) for q in data["heat_flux"].values() for v in q):
        raise ValueError("Nonfinite sourced numerical fields")
    model = normalized["raw_output_verification"].get("volumetric_source")
    if not model or model.get("unit") != "W/m^3" or model.get("region") != "all_volume":
        raise ValueError("Missing/ambiguous audited volumetric source and units")
    source = {r["element"]: r["value"] for r in model["element_values"]}
    if (
        len(source) != len(model["element_values"])
        or set(source) != set(elements)
        or not all(math.isfinite(v) for v in source.values())
    ):
        raise ValueError("Incomplete/invalid audited source region coverage")
    representation = normalized["raw_output_verification"].get(
        "thermal_flux_representation"
    )
    if representation not in ("native_integration_point", "global_l2_nodal_projection"):
        raise ValueError("Unsupported native flux representation")
    kind = normalized["quantities"]["heat_flux"]["source"]["kind"]
    if kind != (
        "solver"
        if representation == "native_integration_point"
        else "postprocessed_from_solver"
    ):
        raise ValueError("Mislabelled native/recovered flux provenance")
    recovered = nodal_vectors(
        normalized["quantities"]["recovered_nodal_heat_flux"], nodes
    )
    g = config["geometry"]
    L = g["length"]
    k = config["material"]["conductivity"]
    balance = source_reference.balance(config)
    generation = balance["generation_W"]
    qscale = max(
        abs(source_reference.analytical(config, x)["heat_flux_W_m2"][0]) for x in (0, L)
    )
    temp_error = max(
        abs(t - source_reference.analytical(config, nodes[n][0])["temperature_K"])
        for n, t in data["temperature"].items()
    )
    sums = {
        key: 0.0
        for key in (
            "temperature",
            "gradient",
            "raw_flux",
            "native_flux",
            "recovered_flux",
        )
    }
    maxima = {
        key: 0.0 for key in ("temperature", "raw_flux", "native_flux", "recovered_flux")
    }
    surfaces = {f"{a}_{s}": 0.0 for a in ("x", "y", "z") for s in ("min", "max")}
    native_surfaces = dict(surfaces)
    projection = {n: [0.0, 0.0, 0.0] for n in nodes}
    test_volume = {n: 0.0 for n in nodes}
    native_closure = 0.0
    native_excess = 0.0
    sampling_closure = 0.0
    volume = 0.0
    actual_generation = 0.0
    raw_mean = [0.0, 0.0, 0.0]
    resolution = (
        normalized["raw_output_verification"]["temperature_print_resolution_K"] / 2
    )
    for e, conn in elements.items():
        start, span = q1.cell(nodes, conn)
        V = math.prod(span)
        volume += V
        actual_generation += source[e] * V
        values = tuple(data["temperature"][n] for n in conn)
        rec = tuple(recovered[n] for n in conn)
        samples = q1.gp_samples(data, nodes, conn, e)
        for local, weight in q1.G3:
            xyz = tuple(a + d * (x + 1) / 2 for a, d, x in zip(start, span, local))
            T, grad = q1.evaluate(values, span, local)
            raw = tuple(-k * v for v in grad)
            exact = source_reference.analytical(config, xyz[0])
            native = q1.gp_vector(samples, local)
            recq = q1.vector(rec, local)
            for a in range(3):
                raw_mean[a] += raw[a] * weight * V / 8
            errT = T - exact["temperature_K"]
            errG = sum((v - w) ** 2 for v, w in zip(grad, exact["gradient_K_m"]))
            errors = {
                "temperature": errT * errT,
                "gradient": errG,
                "raw_flux": sum(
                    (v - w) ** 2 for v, w in zip(raw, exact["heat_flux_W_m2"])
                ),
                "native_flux": sum(
                    (v - w) ** 2 for v, w in zip(native, exact["heat_flux_W_m2"])
                ),
                "recovered_flux": sum(
                    (v - w) ** 2 for v, w in zip(recq, exact["heat_flux_W_m2"])
                ),
            }
            for key, value in errors.items():
                sums[key] += value * weight * V / 8
                if key in maxima:
                    maxima[key] = max(maxima[key], math.sqrt(value))
        # 2x2x2 native sample consistency and the global Galerkin projection residual.
        for local in q1.G2:
            _, grad = q1.evaluate(values, span, local)
            raw = tuple(-k * v for v in grad)
            native = samples[local]
            native_closure = max(
                native_closure, max(abs(v - w) for v, w in zip(native, raw)) / qscale
            )
            native_excess = max(
                native_excess,
                max(
                    max(0.0, abs(native[a] - raw[a]) - 2 * k * resolution / span[a])
                    / qscale
                    for a in range(3)
                ),
            )
            if representation == "global_l2_nodal_projection":
                recq = q1.vector(rec, local)
                sampling_closure = max(
                    sampling_closure,
                    max(abs(v - w) for v, w in zip(native, recq)) / qscale,
                )
            weights, _ = q1.basis(local)
            for n, w in zip(conn, weights):
                test_volume[n] += w * V / 8
                for a in range(3):
                    projection[n][a] += (native[a] - raw[a]) * w * V / 8
        # Two-point face quadrature is exact for these Q1 derivatives/recovered fields.
        for a, axis in enumerate(("x", "y", "z")):
            area = V / span[a]
            for sign, side, pos in (
                (-1, "min", start[a]),
                (1, "max", start[a] + span[a]),
            ):
                if not math.isclose(
                    pos,
                    0 if sign < 0 else g[("length", "width", "height")[a]],
                    rel_tol=0,
                    abs_tol=1e-12,
                ):
                    continue
                others = [b for b in range(3) if b != a]
                for pair in (
                    (-1 / math.sqrt(3), -1 / math.sqrt(3)),
                    (-1 / math.sqrt(3), 1 / math.sqrt(3)),
                    (1 / math.sqrt(3), -1 / math.sqrt(3)),
                    (1 / math.sqrt(3), 1 / math.sqrt(3)),
                ):
                    local = [0.0, 0.0, 0.0]
                    local[a] = sign
                    for b, x in zip(others, pair):
                        local[b] = x
                    _, grad = q1.evaluate(values, span, local)
                    surfaces[f"{axis}_{side}"] += sign * (-k * grad[a]) * area / 4
                    native_surfaces[f"{axis}_{side}"] += (
                        sign * q1.gp_vector(samples, local)[a] * area / 4
                    )
    if not math.isclose(volume, L * g["width"] * g["height"], rel_tol=1e-12):
        raise ValueError("Incomplete/overlapping physical source volume")
    rates = data["heat_rate"]
    hot = sum(v for n, v in rates.items() if abs(nodes[n][0]) < 1e-12)
    cold = sum(v for n, v in rates.items() if abs(nodes[n][0] - L) < 1e-12)
    free = max(
        (
            abs(v)
            for n, v in rates.items()
            if abs(nodes[n][0]) >= 1e-12 and abs(nodes[n][0] - L) >= 1e-12
        ),
        default=0.0,
    )
    reaction = {"x_min": -hot, "x_max": -cold}
    conservation = (
        max(
            abs(-hot - cold - actual_generation),
            abs(-sum(rates.values()) - actual_generation),
        )
        / generation
    )
    projection_residual = max(
        abs(v) / test_volume[n] / qscale for n, vec in projection.items() for v in vec
    )
    h = L / normalized["mesh"]["counts"][0]
    theory = source_reference.interpolation_norms(config, h)
    norms = {key: math.sqrt(value) for key, value in sums.items()}
    # Triangle inequality bounds from the ACTUAL nodal perturbation, not fitted errors.
    perturbation_T = temp_error * math.sqrt(volume)
    counts = normalized["mesh"]["counts"]
    spans = [g[name] / c for name, c in zip(("length", "width", "height"), counts)]
    perturbation_G = (
        2 * temp_error * math.sqrt(sum(1 / d**2 for d in spans)) * math.sqrt(volume)
    )
    quadrature_consistency = abs(
        norms["temperature"] - theory["temperature_L2_K_m3_2"]
    ) <= perturbation_T + 1e-10 * max(theory["temperature_L2_K_m3_2"], 1) and abs(
        norms["gradient"] - theory["gradient_H1_K_m1_2"]
    ) <= perturbation_G + 1e-10 * max(
        theory["gradient_H1_K_m1_2"], 1
    )
    boundary_error = {
        face: surfaces[face] - target
        for face, target in balance["outward_heat_flow_W"].items()
    }
    expected_boundary_error = (
        config["internal_heat_generation"] * h / 2 * g["width"] * g["height"]
    )
    checks = {
        "audited_source_matches_specification": all(
            v == config["internal_heat_generation"] for v in source.values()
        ),
        "source_volume_integral": abs(actual_generation - generation) / generation
        <= TOLERANCES["source_integral_relative"],
        "nodal_temperature": temp_error <= TOLERANCES["nodal_temperature_absolute_K"],
        "nonnodal_quadrature_matches_independent_interpolation_identity": quadrature_consistency,
        "native_flux_consistency": (
            native_excess
            if representation == "native_integration_point"
            else max(projection_residual, sampling_closure)
        )
        <= TOLERANCES["native_flux_consistency_relative"],
        "boundary_reaction_reference": max(
            abs(reaction[f] - v) for f, v in balance["outward_heat_flow_W"].items()
        )
        / generation
        <= TOLERANCES["reaction_relative"],
        "free_node_equilibrium": free / generation <= TOLERANCES["reaction_relative"],
        "reaction_source_energy_conservation": conservation
        <= TOLERANCES["conservation_relative"],
        "adiabatic_sides": max(
            abs(v) for f, v in surfaces.items() if not f.startswith("x")
        )
        / generation
        <= TOLERANCES["reaction_relative"],
        "raw_boundary_flux_discretization_bound": max(
            abs(v) for v in boundary_error.values()
        )
        <= expected_boundary_error
        + 2 * k * temp_error / h * g["width"] * g["height"]
        + 1e-8,
    }
    return {
        "passed": all(checks.values()),
        "checks": checks,
        "h_x_m": h,
        "volume_m3": volume,
        "temperature_max_error_K": temp_error,
        "temperature_L2_K_m3_2": norms["temperature"],
        "temperature_RMS_K": norms["temperature"] / math.sqrt(volume),
        "gradient_H1_K_m1_2": norms["gradient"],
        "raw_flux_L2_W_m_minus1_2": norms["raw_flux"],
        "native_flux_L2_W_m_minus1_2": norms["native_flux"],
        "recovered_flux_L2_W_m_minus1_2": norms["recovered_flux"],
        "maximum_quadrature_errors": maxima,
        "sampling": "3x3x3 Gauss per affine cell (degree 5), includes off-center points; maxima over these samples, not global suprema",
        "volume_quadrature_samples": 27 * len(elements),
        "native_flux_representation": representation,
        "native_pointwise_fourier_difference_relative": native_closure,
        "native_fourier_difference_excess_after_print_bound_relative": native_excess,
        "native_global_projection_residual_relative": projection_residual,
        "surface_outward_heat_flow_W": surfaces,
        "native_surface_outward_heat_flow_W": native_surfaces,
        "boundary_heat_input_W": {"x_min": hot, "x_max": cold},
        "boundary_reaction_outward_heat_flow_W": reaction,
        "actual_source_integral_W": actual_generation,
        "reaction_net_outward_W": sum(reaction.values()),
        "raw_boundary_flux_error_W": boundary_error,
        "raw_boundary_heat_flow_error_W": max(abs(v) for v in boundary_error.values()),
        "raw_net_outward_W": sum(surfaces.values()),
        "normalized_energy_imbalance": conservation,
        "independent_interpolation_norms": theory,
        "mean_heat_flux_W_m2": [v / volume for v in raw_mean],
        "reference": source_reference.analytical(config, 0),
        "precision_floors": {
            "raw_boundary_heat_flow_error_W": 2 * k * resolution / h
            * g["width"] * g["height"],
            "temperature_L2_K_m3_2": resolution * math.sqrt(volume),
            "gradient_H1_K_m1_2": 2
            * resolution
            * math.sqrt(sum(1 / d**2 for d in spans))
            * math.sqrt(volume),
            "raw_flux_L2_W_m_minus1_2": 2
            * k
            * resolution
            * math.sqrt(sum(1 / d**2 for d in spans))
            * math.sqrt(volume),
            "recovered_flux_L2_W_m_minus1_2": normalized["raw_output_verification"][
                "recovered_flux_print_resolution_W_m2"
            ]
            / 2
            * math.sqrt(volume),
        },
    }


def study(runs):
    metrics = [r["thermal"] for r in runs]
    keys = (
        "temperature_L2_K_m3_2",
        "gradient_H1_K_m1_2",
        "raw_flux_L2_W_m_minus1_2",
        "recovered_flux_L2_W_m_minus1_2",
        "raw_boundary_heat_flow_error_W",
    )
    rates = []
    for a, b in zip(metrics, metrics[1:]):
        row = {
            "h_coarse_m": a["h_x_m"],
            "h_fine_m": b["h_x_m"],
            "orders": {},
            "precision_limited": [],
        }
        for key in keys:
            if min(a[key], b[key]) <= TOLERANCES["precision_floor_multiplier"] * max(
                a["precision_floors"][key], b["precision_floors"][key]
            ):
                row["orders"][key] = None
                row["precision_limited"].append(key)
            else:
                row["orders"][key] = math.log(a[key] / b[key]) / math.log(
                    a["h_x_m"] / b["h_x_m"]
                )
        rates.append(row)
    resolved = all(
        row["orders"]["temperature_L2_K_m3_2"] is not None
        and row["orders"]["gradient_H1_K_m1_2"] is not None
        for row in rates
    )
    supported = resolved and all(
        abs(row["orders"]["temperature_L2_K_m3_2"] - 2)
        <= TOLERANCES["temperature_order_deviation"]
        and abs(row["orders"]["gradient_H1_K_m1_2"] - 1)
        <= TOLERANCES["gradient_order_deviation"]
        and abs(row["orders"]["raw_flux_L2_W_m_minus1_2"] - 1)
        <= TOLERANCES["gradient_order_deviation"]
        for row in rates
    )
    decreasing = all(
        b[key] < a[key] for a, b in zip(metrics, metrics[1:]) for key in keys[:3]
    )
    return {
        "kind": "measured_discretization_convergence",
        "passed": len(runs) >= 3
        and all(r["passed"] for r in runs)
        and supported
        and decreasing,
        "observed_rates": rates,
        "resolved_primary_orders": resolved,
        "theoretical_orders_supported": supported,
        "reason": "Measured nonnodal errors; nodal coincidence is not convergence. Recovered flux orders are diagnostic, not constrained to raw FE orders.",
    }
