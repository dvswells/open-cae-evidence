"""Physical modal benchmark validation independent of solver input syntax."""

import math

UNITS = {"length": "m", "mass": "kg", "time": "s", "stress": "Pa", "density": "kg/m^3"}


def validate(config):
    """Reject invalid/mislabeled inputs, not silently convert them.

    Bounds define this metallic slender-beam benchmark's domain, not universal
    material laws. A plausible number with a false unit label cannot be detected.
    """
    base_keys = {"units", "geometry", "material", "boundary", "meshes"}
    if set(config) not in (base_keys, base_keys | {"tip_mass"}):
        raise ValueError("Unexpected or missing configuration keys")
    if config["units"] != UNITS:
        raise ValueError("Explicit SI units are required")
    g, m = config["geometry"], config["material"]
    if set(g) != {"length", "width", "height"} or set(m) != {
        "young_modulus",
        "poisson_ratio",
        "density",
    }:
        raise ValueError("Unexpected geometry/material keys")
    for value in (*g.values(), *m.values()):
        if (
            isinstance(value, bool)
            or not isinstance(value, (int, float))
            or not math.isfinite(value)
        ):
            raise ValueError("Physical inputs must be finite numbers")
    if not all(0.001 <= v <= 10 for v in g.values()):
        raise ValueError("Dimensions outside benchmark SI domain [0.001, 10] m")
    if g["length"] / max(g["width"], g["height"]) < 50:
        raise ValueError("Beam must be slender: L/max(b,h) >= 50")
    if not 1e9 <= m["young_modulus"] <= 1e12 or not 1000 <= m["density"] <= 30000:
        raise ValueError("Metallic benchmark expects E in Pa and density in kg/m^3")
    if not -1 < m["poisson_ratio"] < 0.5:
        raise ValueError("Isotropic elastic stability requires -1 < nu < 0.5")
    if config["boundary"] != {"plane": "x=0", "dofs": [1, 2, 3]} or any(
        type(d) is not int for d in config["boundary"]["dofs"]
    ):
        raise ValueError("Clamp x=0 in all three translational DOFs")
    if "tip_mass" in config:
        tip = config["tip_mass"]
        if not isinstance(tip, dict) or set(tip) != {
            "mass",
            "unit",
            "position",
            "coupling",
            "rotary_inertia",
        }:
            raise ValueError("Complete tip mass specification required")
        mass = tip["mass"]
        beam_mass = math.prod(g.values()) * m["density"]
        if (
            type(mass) not in (int, float)
            or not math.isfinite(mass)
            or not 0 <= mass <= 10 * beam_mass
        ):
            raise ValueError("Tip mass must be finite and in [0,10*beam mass] kg")
        if tip["unit"] != "kg":
            raise ValueError("Tip mass is a mass in kg, not density or force")
        if tip["position"] != [g["length"], 0, 0] or any(
            type(v) not in (int, float) for v in tip["position"]
        ):
            raise ValueError("Tip mass must be centered at [L,0,0]")
        if tip["coupling"] != "area-average-q8" or tip["rotary_inertia"] != "neglected":
            raise ValueError(
                "Use area-average-q8 translation coupling with neglected rotary inertia"
            )
    meshes = config["meshes"]
    if not isinstance(meshes, list) or len(meshes) < 3:
        raise ValueError("At least three refinement levels are required")
    previous = None
    for mesh in meshes:
        if (
            not isinstance(mesh, list)
            or len(mesh) != 3
            or any(type(n) is not int or n < 2 or n > 128 for n in mesh)
            or any(n % 2 for n in mesh)
        ):
            raise ValueError("Meshes require three even integer counts in [2,128]")
        if math.prod(mesh) > 10000:
            raise ValueError("Mesh exceeds small benchmark resource budget")
        if previous and any(a <= b for a, b in zip(mesh, previous)):
            raise ValueError("Refine all three directions monotonically")
        previous = mesh
