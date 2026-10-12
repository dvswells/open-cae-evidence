"""Select and validate physical problems without selecting a solver."""

from .interfaces import ProblemSpec
from .modal_problem import validate as validate_modal
from .thermal_problem import validate as validate_thermal
from .thermal_problem import sourced


def load_problem(config):
    if config.get("analysis") == "steady_heat_conduction":
        validate_thermal(config)
        b = config["boundary"]
        conditions = [
            {
                "type": "prescribed_temperature",
                "face": face,
                "unit": "K",
                "value": b[key],
            }
            for face, key in (
                ("x_min", "x_min_temperature"),
                ("x_max", "x_max_temperature"),
            )
        ]
        conditions += [
            {"type": "adiabatic", "face": face, "unit": "W/m^2", "value": 0}
            for face in ("y_min", "y_max", "z_min", "z_max")
        ]
        return ProblemSpec(
            config["problem_id"],
            config["analysis"],
            config,
            config["units"],
            conditions,
            ("temperature", "heat_flux", "nodal_heat_rate"),
            (
                "transient_heat",
                "temperature_dependent_conductivity",
                *(
                    ("nonuniform_heat_generation",)
                    if sourced(config)
                    else ("internal_heat_generation",)
                ),
                "experimental_validation",
            ),
        )
    validate_modal(config)
    tip = "tip_mass" in config
    return ProblemSpec(
        "cantilever_tip_mass_modal" if tip else "cantilever_modal",
        "modal",
        config,
        config["units"],
        [
            {
                "type": "prescribed_displacement",
                "face": "x_min",
                "components": ["x", "y", "z"],
                "value": [0, 0, 0],
                "unit": "m",
            }
        ],
        ("frequency", "eigenvalue", "displacement"),
        ("torsional_axial_reference", "damping", "prestress", "experimental_validation")
        + (("tip_rotary_inertia", "attachment_compliance") if tip else ()),
    )
