"""Audited native Elmer 808 brick mesh and the supported SI SIF subset.

No analytical temperature or heat-flow values are used by serialization.
"""

import itertools
import json
import math

from ..thermal_problem import validate, sourced

CORNERS = (
    (0, 0, 0),
    (1, 0, 0),
    (1, 1, 0),
    (0, 1, 0),
    (0, 0, 1),
    (1, 0, 1),
    (1, 1, 1),
    (0, 1, 1),
)
# Outward oriented quads in Elmer 404 ordering, from audited 808 corners.
FACES = (
    (0, 4, 7, 3),
    (1, 2, 6, 5),
    (0, 1, 5, 4),
    (3, 7, 6, 2),
    (0, 3, 2, 1),
    (4, 5, 6, 7),
)


def sif(config):
    b = config["boundary"]
    k = config["material"]["conductivity"]
    text = f"""! SI: coordinates m; temperature K; conductivity W/(m*K); flux W/m^2; loads W
Header
  CHECK KEYWORDS Warn
  Mesh DB "." "."
End
Simulation
  Max Output Level = 5
  Coordinate System = Cartesian 3D
  Coordinate Mapping(3) = 1 2 3
  Simulation Type = Steady State
  Steady State Max Iterations = 1
  Output Intervals = 1
End
Body 1
  Target Bodies(1) = 1
  Equation = 1
  Material = 1
End
Equation 1
  Active Solvers(3) = 1 2 3
End
Material 1
  Heat Conductivity = Real {k:.16g}
  Density = Real 1
End
Solver 1
  Equation = Heat Equation
  Procedure = "HeatSolve" "HeatSolver"
  Variable = Temperature
  Variable Dofs = 1
  Linear System Solver = Iterative
  Linear System Iterative Method = BiCGStab
  Linear System Preconditioning = ILU0
  Linear System Max Iterations = 1000
  Linear System Convergence Tolerance = 1.0e-12
  Linear System Abort Not Converged = True
  Linear System Residual Output = 1
  Steady State Convergence Tolerance = 1.0e-10
  Nonlinear System Max Iterations = 1
  Nonlinear System Convergence Tolerance = 1.0e-10
  Calculate Loads = True
  Stabilize = False
End
Solver 2
  Equation = ComputeFlux
  Procedure = "FluxSolver" "FluxSolver"
  Calculate Flux = True
  Target Variable = String Temperature
  Flux Coefficient = String "Heat Conductivity"
  Linear System Solver = Direct
  Linear System Direct Method = Umfpack
End
Solver 3
  Equation = ResultOutput
  Procedure = "ResultOutputSolve" "ResultOutputSolver"
  Exec Solver = After Saving
  Output File Name = slab
  Vtu Format = True
  Binary Output = False
  Single Precision = False
  Save Bulk Only = True
  Save Geometry Ids = True
  Scalar Field 1 = String Temperature
  Scalar Field 2 = String "Temperature Loads"
  Vector Field 1 = String "Temperature Flux"
End
Boundary Condition 1
  Target Boundaries(1) = 1
  Temperature = Real {b['x_min_temperature']:.16g}
End
Boundary Condition 2
  Target Boundaries(1) = 2
  Temperature = Real {b['x_max_temperature']:.16g}
End
Boundary Condition 3
  Target Boundaries(4) = 3 4 5 6
  Heat Flux = Real 0
End
"""
    if sourced(config):
        text = text.replace("  Material = 1\n", "  Material = 1\n  Body Force = 1\n", 1)
        text = text.replace(
            "Solver 1\n",
            f'! Uniform volumetric source W/m^3, not a density-scaled specific source\nBody Force 1\n  Volumetric Heat Source = Real {config["internal_heat_generation"]:.16g}\nEnd\nSolver 1\n',
            1,
        )
    return text


def generate(config, counts):
    validate(config)
    if counts not in config["meshes"]:
        raise ValueError("Elmer mesh must be a declared physical mesh resolution")
    lengths = [config["geometry"][x] for x in ("length", "width", "height")]
    grid, nodes, elements = {}, {}, {}
    # Different numbering from CalculiX, deliberately. The order is native 808.
    keys = list(itertools.product(*(range(c + 1) for c in counts)))
    for n, key in enumerate(reversed(keys), 1):
        grid[key] = n
        nodes[n] = tuple(
            v * length / count for v, length, count in zip(key, lengths, counts)
        )
    boundary = []
    for start in itertools.product(*(range(c) for c in counts)):
        e = len(elements) + 1
        conn = tuple(
            grid[tuple(a + b for a, b in zip(start, corner))] for corner in CORNERS
        )
        elements[e] = conn
        for face, local in enumerate(FACES):
            axis, upper = divmod(face, 2)
            if start[axis] == (counts[axis] - 1 if upper else 0):
                boundary.append(
                    (len(boundary) + 1, face + 1, e, 0, 404, *(conn[i] for i in local))
                )
    files = {
        "slab.sif": sif(config),
        "mesh.header": f"{len(nodes)} {len(elements)} {len(boundary)}\n2\n808 {len(elements)}\n404 {len(boundary)}\n",
        "mesh.nodes": "".join(
            f"{n} -1 " + " ".join(f"{x:.16g}" for x in xyz) + "\n"
            for n, xyz in nodes.items()
        ),
        "mesh.elements": "".join(
            " ".join(map(str, (e, 1, 808, *conn))) + "\n"
            for e, conn in elements.items()
        ),
        "mesh.boundary": "".join(" ".join(map(str, row)) + "\n" for row in boundary),
        "units.json": json.dumps(config["units"], sort_keys=True, indent=2) + "\n",
    }
    return files, nodes, elements


def audit_sif(text, config):
    """Parse physical assignments independently of the SIF text generator."""
    sections, active = {}, None
    for raw in text.splitlines():
        line = raw.split("!", 1)[0].strip().lower()
        if not line:
            continue
        if line == "end":
            if active is None:
                raise ValueError("Unexpected SIF End")
            active = None
        elif active is None:
            if line in sections:
                raise ValueError("Duplicate SIF section")
            active = line
            sections[active] = {}
        else:
            if active == "header" and line == "check keywords warn":
                sections[active]["check keywords"] = "warn"
                continue
            if active == "header" and line.startswith("mesh db "):
                if "mesh db" in sections[active]:
                    raise ValueError("Duplicate native mesh location")
                sections[active]["mesh db"] = line[len("mesh db ") :]
                continue
            if "=" not in line:
                raise ValueError("Malformed SIF assignment")
            key, value = (v.strip() for v in line.split("=", 1))
            if key in sections[active]:
                raise ValueError("Duplicate SIF keyword")
            sections[active][key] = value
    if active is not None:
        raise ValueError("Unclosed SIF section")
    expected = {
        "header": {"mesh db": '"." "."'},
        "simulation": {
            "coordinate system": "cartesian 3d",
            "coordinate mapping(3)": "1 2 3",
            "simulation type": "steady state",
            "steady state max iterations": "1",
        },
        "body 1": {"target bodies(1)": "1", "material": "1", "equation": "1"},
        "equation 1": {"active solvers(3)": "1 2 3"},
        "boundary condition 1": {"target boundaries(1)": "1"},
        "boundary condition 2": {"target boundaries(1)": "2"},
        "boundary condition 3": {"target boundaries(4)": "3 4 5 6"},
        "solver 1": {
            "procedure": '"heatsolve" "heatsolver"',
            "variable": "temperature",
            "variable dofs": "1",
            "calculate loads": "true",
            "stabilize": "false",
        },
        "solver 2": {
            "procedure": '"fluxsolver" "fluxsolver"',
            "target variable": "string temperature",
            "flux coefficient": 'string "heat conductivity"',
            "calculate flux": "true",
        },
        "solver 3": {"procedure": '"resultoutputsolve" "resultoutputsolver"'},
    }
    for section, assignments in expected.items():
        if any(
            sections.get(section, {}).get(key) != value
            for key, value in assignments.items()
        ):
            raise ValueError("Parsed Elmer physical assignments/procedure differ")
    numerical = [
        ("material 1", "heat conductivity", config["material"]["conductivity"]),
        ("material 1", "density", 1),
        (
            "boundary condition 1",
            "temperature",
            config["boundary"]["x_min_temperature"],
        ),
        (
            "boundary condition 2",
            "temperature",
            config["boundary"]["x_max_temperature"],
        ),
        ("boundary condition 3", "heat flux", 0),
    ]
    if sourced(config):
        if sections.get("body 1", {}).get("body force") != "1":
            raise ValueError("Source is not applied over body 1/all cells")
        numerical.append(
            (
                "body force 1",
                "volumetric heat source",
                config["internal_heat_generation"],
            )
        )
    for section, key, expected_value in numerical:
        raw = sections.get(section, {}).get(key, "")
        try:
            value = float(raw.removeprefix("real "))
        except ValueError as exc:
            raise ValueError("Missing/malformed native physical property") from exc
        if not math.isfinite(value) or not math.isclose(
            value, expected_value, rel_tol=1e-14, abs_tol=1e-14
        ):
            raise ValueError(
                "Parsed Elmer conductivity/temperature/flux differs from ProblemSpec"
            )


def audit(files, config, counts):
    """Independently parse geometry, cells, oriented faces and physical SIF."""
    validate(config)
    audit_sif(files["slab.sif"], config)
    if (
        files["units.json"]
        != json.dumps(config["units"], sort_keys=True, indent=2) + "\n"
    ):
        raise ValueError("Elmer SI unit mapping differs")
    # Canonical supported keyword subset forbids hidden sources, scaling or BC overrides.
    if files["slab.sif"] != sif(config):
        raise ValueError(
            "Elmer material, boundary, coordinate or procedure audit failed"
        )
    lengths = [config["geometry"][x] for x in ("length", "width", "height")]
    nodes, grid = {}, {}
    for row in files["mesh.nodes"].splitlines():
        n, partition, *xyz = row.split()
        n, xyz = int(n), tuple(map(float, xyz))
        if (
            partition != "-1"
            or len(xyz) != 3
            or n in nodes
            or not all(map(math.isfinite, xyz))
        ):
            raise ValueError("Invalid native node")
        key = tuple(round(x / length * c) for x, length, c in zip(xyz, lengths, counts))
        if key in grid or any(
            i < 0 or i > c or not math.isclose(x, i * length / c, abs_tol=1e-12)
            for x, i, length, c in zip(xyz, key, lengths, counts)
        ):
            raise ValueError("Native mesh coordinates or length units differ")
        nodes[n], grid[key] = xyz, n
    if len(nodes) != math.prod(c + 1 for c in counts):
        raise ValueError("Incomplete native nodes")
    elements, cells, expected_faces = {}, set(), set()
    for row in files["mesh.elements"].splitlines():
        e, body, kind, *conn = map(int, row.split())
        if (
            body != 1
            or kind != 808
            or e in elements
            or len(conn) != 8
            or len(set(conn)) != 8
            or not set(conn) <= nodes.keys()
        ):
            raise ValueError("Invalid native connectivity/type/body")
        start = tuple(
            round(x / length * c)
            for x, length, c in zip(nodes[conn[0]], lengths, counts)
        )
        if start in cells or any(i >= c for i, c in zip(start, counts)):
            raise ValueError("Repeated/outside cell")
        for n, corner in zip(conn, CORNERS):
            if n != grid.get(tuple(a + b for a, b in zip(start, corner))):
                raise ValueError("Inverted or non-affine Elmer 808 cell")
        elements[e], _ = tuple(conn), cells.add(start)
        for face, local in enumerate(FACES):
            axis, upper = divmod(face, 2)
            if start[axis] == (counts[axis] - 1 if upper else 0):
                expected_faces.add((face + 1, e, 0, 404, *(conn[i] for i in local)))
    if len(cells) != math.prod(counts):
        raise ValueError("Incomplete native cells")
    faces, ids = set(), set()
    for row in files["mesh.boundary"].splitlines():
        n, *entry = map(int, row.split())
        if n in ids or tuple(entry) in faces:
            raise ValueError("Duplicate boundary face")
        ids.add(n)
        faces.add(tuple(entry))
    if faces != expected_faces:
        raise ValueError("Boundary assignment/orientation/parent differs")
    header = f"{len(nodes)} {len(elements)} {len(faces)}\n2\n808 {len(elements)}\n404 {len(faces)}\n"
    if files["mesh.header"] != header:
        raise ValueError("Native mesh header differs")
    result = {
        "passed": True,
        "nodes": len(nodes),
        "elements": len(elements),
        "boundary_faces": len(faces),
        "volume_m3": math.prod(lengths),
        "positive_affine_jacobians": True,
        "native_element": 808,
        "native_boundary_element": 404,
        "units": config["units"],
    }
    if sourced(config):
        row = next(
            line
            for line in files["slab.sif"].splitlines()
            if line.strip().startswith("Volumetric Heat Source =")
        )
        value = float(row.split("Real ", 1)[1])
        result["volumetric_source"] = {
            "unit": "W/m^3",
            "region": "all_volume",
            "source_files": ["slab.sif", "mesh.elements", "units.json"],
            "element_values": [
                {"element": e, "value": value} for e in sorted(elements)
            ],
        }
    return result
