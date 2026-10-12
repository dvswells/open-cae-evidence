"""Explicit derived physical quantities from actual normalized nodal output."""

from . import q1
from .results import mesh_maps, thermal_fields, quantity


def add_quantities(problem, normalized):
    nodes, elements = mesh_maps(normalized)
    data = thermal_fields(normalized)
    gradients = []
    fluxes = []
    k = problem.parameters["material"]["conductivity"]
    for e, conn in sorted(elements.items()):
        start, span = q1.cell(nodes, conn)
        values = tuple(data["temperature"][n] for n in conn)
        for p, local in enumerate(q1.G2, 1):
            _, grad = q1.evaluate(values, span, local)
            xyz = [a + d * (x + 1) / 2 for a, d, x in zip(start, span, local)]
            gradients.append(
                {"element": e, "point": p, "coordinates": xyz, "value": list(grad)}
            )
            fluxes.append(
                {
                    "element": e,
                    "point": p,
                    "coordinates": xyz,
                    "value": [-k * v for v in grad],
                }
            )
    files = normalized["quantities"]["temperature"]["source"]["files"]
    normalized["quantities"]["raw_temperature_gradient"] = quantity(
        "K/m",
        "integration_point",
        ["x", "y", "z"],
        gradients,
        "postprocessed_from_solver",
        files,
        "Unrecovered Q1 derivative of actual numerical temperature DOFs at canonical physical Gauss locations; no native gradient claim.",
    )
    normalized["quantities"]["raw_fourier_heat_flux"] = quantity(
        "W/m^2",
        "integration_point",
        ["x", "y", "z"],
        fluxes,
        "postprocessed_from_solver",
        files,
        "-k grad(T_h) from actual numerical nodal output and audited conductivity. Common cross-solver physical definition, distinct from native/projected heat_flux.",
    )
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
        files,
        "Exterior one-sided raw Q1-gradient face quadrature. Discretization estimate, not a Dirichlet reaction or an energy-conserving flux reconstruction.",
    )
    normalized["quantities"]["boundary_reaction_heat_flow"] = quantity(
        "W",
        "boundary",
        ["outward_heat_flow"],
        [
            {"face": face, "value": v}
            for face, v in metrics["boundary_reaction_outward_heat_flow_W"].items()
        ],
        "postprocessed_from_solver",
        normalized["quantities"]["nodal_heat_rate"]["source"]["files"],
        "Negative sum of normalized actual boundary heat-input reactions, with audited consistent body-source removal from CalculiX RFL where required.",
    )
