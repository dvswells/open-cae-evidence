"""Physical thermal specification and independent Fourier reference; no solver I/O."""

import math

UNITS = {
    "length": "m",
    "temperature": "K",
    "conductivity": "W/(m*K)",
    "heat_flux": "W/m^2",
    "heat_flow": "W",
}


def sourced(config):
    return config.get("problem_id") == "homogeneous_slab_volumetric_source"


def validate(config):
    keys = {
        "spec_version",
        "problem_id",
        "analysis",
        "units",
        "geometry",
        "material",
        "boundary",
        "internal_heat_generation",
        "meshes",
    }
    if sourced(config):
        keys.add("source_region")
    if (
        set(config) != keys
        or config["spec_version"] != ("1.1" if sourced(config) else "1.0")
        or config["analysis"] != "steady_heat_conduction"
    ):
        raise ValueError("Unknown or incomplete thermal problem specification")
    if config["problem_id"] not in (
        "homogeneous_slab_steady_heat",
        "homogeneous_slab_volumetric_source",
    ) or config["units"] != (
        {**UNITS, "volumetric_heat_generation": "W/m^3"} if sourced(config) else UNITS
    ):
        raise ValueError("This slab benchmark requires explicit SI thermal units")
    g, m, b = config["geometry"], config["material"], config["boundary"]
    if (
        set(g) != {"length", "width", "height"}
        or set(m) != {"conductivity"}
        or set(b) != {"x_min_temperature", "x_max_temperature", "other_faces"}
    ):
        raise ValueError("Invalid thermal material, geometry or boundary keys")
    values = [
        *g.values(),
        m["conductivity"],
        b["x_min_temperature"],
        b["x_max_temperature"],
        config["internal_heat_generation"],
    ]
    if any(type(v) not in (int, float) or not math.isfinite(v) for v in values):
        raise ValueError("Thermal physical inputs must be finite numbers")
    if (
        not all(0.001 <= v <= 10 for v in g.values())
        or not 0 < m["conductivity"] <= 10000
    ):
        raise ValueError("Positive SI dimensions and conductivity are required")
    if not 0 < b["x_max_temperature"] < b["x_min_temperature"] <= 1000:
        raise ValueError(
            "Assign the hotter temperature to x_min, colder to x_max, in K"
        )
    if b["other_faces"] != "adiabatic":
        raise ValueError("Adiabatic sides are required")
    if sourced(config):
        if (
            config["source_region"] != "all_volume"
            or not 0 < config["internal_heat_generation"] <= 1e7
        ):
            raise ValueError("Source requires positive finite W/m^3 over all volume")
    elif config["internal_heat_generation"] != 0:
        raise ValueError("Reference requires adiabatic sides and zero heat generation")
    if len(config["meshes"]) < 3:
        raise ValueError("At least three meshes are required for consistency")
    previous = None
    for mesh in config["meshes"]:
        if (
            not isinstance(mesh, list)
            or len(mesh) != 3
            or any(type(n) is not int or not 2 <= n <= 32 for n in mesh)
            or math.prod(mesh) > 4096
        ):
            raise ValueError("Thermal meshes require three bounded integer counts")
        if previous and any(a <= b for a, b in zip(mesh, previous)):
            raise ValueError("Refine every direction for the consistency study")
        previous = mesh


def analytical(config, x):
    """Closed-form Fourier solution, not FE interpolation or solver-derived data."""
    if sourced(config):
        from .source_reference import analytical as source_analytical

        return source_analytical(config, x)
    g, b = config["geometry"], config["boundary"]
    gradient = (b["x_max_temperature"] - b["x_min_temperature"]) / g["length"]
    flux = -config["material"]["conductivity"] * gradient
    return {
        "temperature_K": b["x_min_temperature"] + gradient * x,
        "gradient_K_m": [gradient, 0.0, 0.0],
        "heat_flux_W_m2": [flux, 0.0, 0.0],
        "heat_flow_W": flux * g["width"] * g["height"],
    }
