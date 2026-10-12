"""Versioned solver-independent result contract and physical field helpers."""

from functools import lru_cache
from importlib.resources import files
import json
import math

from jsonschema import Draft202012Validator

SCHEMA_VERSION = "1.0.0"
ARCHIVE_VERSION = 3
REQUIRED = {
    "modal": {
        "frequency": ("Hz", "mode", ["frequency"]),
        "eigenvalue": ("rad^2/s^2", "mode", ["eigenvalue"]),
        "displacement": ("arbitrary", "node_mode", ["x", "y", "z"]),
    },
    "steady_heat_conduction": {
        "temperature": ("K", "node", ["temperature"]),
        "heat_flux": ("W/m^2", "integration_point", ["x", "y", "z"]),
        "nodal_heat_rate": ("W", "node", ["heat_input"]),
    },
}


def require_quantities(analysis, normalized):
    for name, (unit, location, components) in REQUIRED[analysis].items():
        q = normalized["quantities"].get(name)
        if q is None or (q["unit"], q["location"], q["components"]) != (
            unit,
            location,
            components,
        ):
            raise ValueError(
                f"Required {analysis} quantity {name} with physical units/components unavailable"
            )


def quantity(unit, location, components, values, source, files_used, note=""):
    return {
        "unit": unit,
        "location": location,
        "components": components,
        "values": values,
        "source": {"kind": source, "files": files_used, "note": note},
    }


def physical_mesh(prepared):
    mesh = {
        "coordinate_unit": "m",
        "counts": prepared.counts,
        "nodes": [],
        "elements": [],
    }
    for n, xyz in sorted(prepared.nodes.items()):
        mesh["nodes"].append(
            {
                "id": n,
                "coordinates": list(xyz),
                "role": "point_mass" if n in prepared.auxiliary_nodes else "solid",
            }
        )
    for e, conn in sorted(prepared.elements.items()):
        mesh["elements"].append(
            {
                "id": e,
                "topology": prepared.topology,
                "connectivity": list(conn),
                "integration": prepared.integration,
            }
        )
    if prepared.auxiliary_nodes:
        mesh["elements"].append(
            {
                "id": max(prepared.elements) + 1,
                "topology": "point_mass",
                "connectivity": sorted(prepared.auxiliary_nodes),
                "integration": "not_applicable",
            }
        )
    return mesh


def mesh_maps(normalized):
    nodes = {r["id"]: tuple(r["coordinates"]) for r in normalized["mesh"]["nodes"]}
    elements = {
        r["id"]: tuple(r["connectivity"])
        for r in normalized["mesh"]["elements"]
        if r["topology"] != "point_mass"
    }
    if len(nodes) != len(normalized["mesh"]["nodes"]) or len(elements) != sum(
        r["topology"] != "point_mass" for r in normalized["mesh"]["elements"]
    ):
        raise ValueError("Duplicate normalized mesh IDs")
    return nodes, elements


def modal_fields(normalized):
    q = normalized["quantities"]
    freq = q["frequency"]["values"]
    lam = q["eigenvalue"]["values"]
    if [v["mode"] for v in freq] != list(range(1, 7)) or [
        v["mode"] for v in lam
    ] != list(range(1, 7)):
        raise ValueError("Six numbered normalized eigenpairs required")
    nodes, _ = mesh_maps(normalized)
    modes = [{} for _ in range(6)]
    for row in q["displacement"]["values"]:
        n, m = row["node"], row["mode"]
        if not 1 <= m <= 6 or n in modes[m - 1] or len(row["value"]) != 3:
            raise ValueError("Duplicate/invalid normalized mode shape")
        modes[m - 1][n] = tuple(row["value"])
    if any(set(m) != set(nodes) for m in modes):
        raise ValueError("Incomplete normalized mode shapes")
    return [r["value"] for r in freq], [r["value"] for r in lam], modes


def thermal_fields(normalized):
    q = normalized["quantities"]
    data = {
        "temperature": {},
        "heat_rate": {},
        "heat_flux": {},
        "integration_coordinates": {},
    }
    for name, key in (("temperature", "temperature"), ("nodal_heat_rate", "heat_rate")):
        for r in q[name]["values"]:
            if r["node"] in data[key]:
                raise ValueError("Duplicate normalized thermal node")
            data[key][r["node"]] = r["value"]
    for r in q["heat_flux"]["values"]:
        key = r["element"], r["point"]
        if key in data["heat_flux"] or len(r["value"]) != 3:
            raise ValueError("Duplicate/invalid normalized thermal flux")
        data["heat_flux"][key] = tuple(r["value"])
        data["integration_coordinates"][key] = tuple(r["coordinates"])
    return data


@lru_cache(maxsize=1)
def validator():
    schema = json.loads(
        files("open_cae_evidence").joinpath("schemas/results-v1.json").read_text()
    )
    Draft202012Validator.check_schema(schema)
    return Draft202012Validator(schema)


def validate_contract(contract):
    if contract.get("schema_version") != SCHEMA_VERSION:
        raise ValueError(
            f"Unsupported results schema version: {contract.get('schema_version')}"
        )

    def finite(value):
        if isinstance(value, float) and not math.isfinite(value):
            raise ValueError("Nonfinite results contract value")
        if isinstance(value, dict):
            for v in value.values():
                finite(v)
        elif isinstance(value, list):
            for v in value:
                finite(v)

    finite(contract)
    errors = validator().iter_errors(contract)
    first = next(errors, None)
    if first is not None:
        raise ValueError(
            f"Invalid results contract at {list(first.path)}: {first.message}"
        )
    from .specification import load_problem

    problem = load_problem(contract["problem"]["parameters"])
    for key, value in (
        ("id", problem.id),
        ("analysis", problem.analysis),
        ("units", problem.units),
        ("boundary_conditions", problem.boundary_conditions),
    ):
        if contract["problem"][key] != value:
            raise ValueError(
                "Normalized physical specification disagrees with its parameters"
            )
    analysis = contract["problem"]["analysis"]
    for run in contract["runs"]:
        nodes, elements = mesh_maps(run)
        require_quantities(analysis, run)
        if run["verification"]["status"] == "passed" and (
            not all(run["verification"]["checks"].values())
            or not run["raw_output_verification"]["passed"]
        ):
            raise ValueError("Passed normalized run contains failed checks")
        if analysis == "modal":
            modal_fields(run)
        else:
            data = thermal_fields(run)
            expected = {(e, p) for e in elements for p in range(1, 9)}
            if (
                set(data["temperature"]) != set(nodes)
                or set(data["heat_rate"]) != set(nodes)
                or set(data["heat_flux"]) != expected
            ):
                raise ValueError("Incomplete normalized thermal result coverage")
    if contract["verification"]["status"] == "passed" and (
        len(contract["runs"]) < 3
        or any(r["verification"]["status"] != "passed" for r in contract["runs"])
    ):
        raise ValueError("Passed contract requires all mesh runs verified")
    if (
        contract["verification"]["status"] == "passed"
        and [r["mesh"]["counts"] for r in contract["runs"]]
        != problem.parameters["meshes"]
    ):
        raise ValueError(
            "Passed contract must contain exactly the prescribed mesh study"
        )
    return contract


def build_contract(problem, runs, provenance, status, tolerances, study):
    # Tolerance numbers retain their original values. Annotate physical dimensions.
    limits = {
        k: {"value": v, "unit": "K" if k.endswith("_K") else "1"}
        for k, v in tolerances.items()
    }
    contract = {
        "schema_name": "open-cae-evidence.results",
        "schema_version": SCHEMA_VERSION,
        "problem": {
            "id": problem.id,
            "analysis": problem.analysis,
            "units": problem.units,
            "parameters": problem.parameters,
            "boundary_conditions": problem.boundary_conditions,
        },
        "solver": {
            "name": provenance["solver_name"],
            "adapter_version": provenance["adapter_version"],
            "provenance": provenance,
        },
        "runs": runs,
        "verification": {"status": status, "tolerances": limits, "mesh_study": study},
        "unsupported_capabilities": list(problem.unsupported_capabilities),
    }
    return validate_contract(contract)
