"""CalculiX DC3D8 input generation and serialized input audit."""

import math

from .thermal_problem import validate

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


def make_deck(config, mesh):
    validate(config)
    nx, ny, nz = mesh
    g = config["geometry"]
    nodes, grid, elements = {}, {}, {}
    for i in range(nx + 1):
        for j in range(ny + 1):
            for k in range(nz + 1):
                n = len(nodes) + 1
                grid[i, j, k] = n
                nodes[n] = (
                    i * g["length"] / nx,
                    j * g["width"] / ny,
                    k * g["height"] / nz,
                )
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                elements[len(elements) + 1] = tuple(
                    grid[i + a, j + b, k + c] for a, b, c in CORNERS
                )
    lines = [
        "*HEADING",
        "Open CAE Evidence steady heat slab; SI m W K",
        "*NODE, NSET=ALL",
    ]
    lines += [
        f"{n}, " + ", ".join(f"{v:.16g}" for v in xyz) for n, xyz in nodes.items()
    ]
    lines += ["*ELEMENT, TYPE=DC3D8, ELSET=SLAB"]
    lines += [", ".join(map(str, [e, *conn])) for e, conn in elements.items()]
    for name, value in (("HOT", 0), ("COLD", g["length"])):
        ids = [n for n, xyz in nodes.items() if xyz[0] == value]
        lines += [f"*NSET, NSET={name}"]
        lines += [", ".join(map(str, ids[i : i + 16])) for i in range(0, len(ids), 16)]
    b = config["boundary"]
    lines += [
        "*MATERIAL, NAME=SOLID",
        "*CONDUCTIVITY",
        f'{config["material"]["conductivity"]:.16g}',
        "*SOLID SECTION, ELSET=SLAB, MATERIAL=SOLID",
        "*INITIAL CONDITIONS, TYPE=TEMPERATURE",
        f'ALL, {b["x_max_temperature"]:.16g}',
        "*STEP",
        "*HEAT TRANSFER, STEADY STATE",
        "1, 1",
        "*BOUNDARY",
        f'HOT, 11, 11, {b["x_min_temperature"]:.16g}',
        f'COLD, 11, 11, {b["x_max_temperature"]:.16g}',
        "*NODE PRINT, NSET=ALL",
        "NT, RFL",
        "*EL PRINT, ELSET=SLAB, GLOBAL=YES",
        "HFL, COORD",
        "*NODE FILE",
        "NT, RFL",
        "*EL FILE",
        "HFL",
        "*END STEP",
    ]
    deck = "\n".join(lines) + "\n"
    return deck, nodes, elements, audit_deck(deck, config, mesh)


def audit_deck(deck, config, mesh):
    validate(config)
    sections = []
    for line in deck.splitlines():
        if not line.strip() or line.startswith("**"):
            continue
        if line.startswith("*"):
            sections.append((line, []))
        elif not sections:
            raise ValueError("Thermal input data before keyword")
        else:
            sections[-1][1].append(line)
    expected = [
        "*HEADING",
        "*NODE, NSET=ALL",
        "*ELEMENT, TYPE=DC3D8, ELSET=SLAB",
        "*NSET, NSET=HOT",
        "*NSET, NSET=COLD",
        "*MATERIAL, NAME=SOLID",
        "*CONDUCTIVITY",
        "*SOLID SECTION, ELSET=SLAB, MATERIAL=SOLID",
        "*INITIAL CONDITIONS, TYPE=TEMPERATURE",
        "*STEP",
        "*HEAT TRANSFER, STEADY STATE",
        "*BOUNDARY",
        "*NODE PRINT, NSET=ALL",
        "*EL PRINT, ELSET=SLAB, GLOBAL=YES",
        "*NODE FILE",
        "*EL FILE",
        "*END STEP",
    ]
    if [s[0] for s in sections] != expected:
        raise ValueError(
            "Thermal keyword contract violated; no loads/sources permitted"
        )
    data = dict(sections)
    nodes, elements = {}, {}
    for row in data[expected[1]]:
        fields = row.split(",")
        n = int(fields[0])
        if n in nodes or len(fields) != 4:
            raise ValueError("Duplicate/malformed thermal node")
        nodes[n] = tuple(float(v) for v in fields[1:])
    g, b = config["geometry"], config["boundary"]
    lengths = [g["length"], g["width"], g["height"]]
    if len(nodes) != math.prod(n + 1 for n in mesh) or any(
        not math.isfinite(v) for xyz in nodes.values() for v in xyz
    ):
        raise ValueError("Invalid thermal mesh node coverage")
    grid = {}
    for n, xyz in nodes.items():
        key = tuple(round(x / l * c) for x, l, c in zip(xyz, lengths, mesh))
        if key in grid or any(
            k < 0 or k > c or not math.isclose(x, l * k / c, abs_tol=1e-12)
            for x, l, k, c in zip(xyz, lengths, key, mesh)
        ):
            raise ValueError("Thermal coordinates differ from structured slab")
        grid[key] = n
    cells, used = set(), set()
    for row in data[expected[2]]:
        e, *conn = [int(v) for v in row.split(",")]
        if (
            e in elements
            or len(conn) != 8
            or len(set(conn)) != 8
            or not set(conn) <= set(nodes)
        ):
            raise ValueError("Bad thermal connectivity")
        start = tuple(
            round(x / l * c) for x, l, c in zip(nodes[conn[0]], lengths, mesh)
        )
        if start in cells or any(k >= c for k, c in zip(start, mesh)):
            raise ValueError("Duplicate/outside thermal cell")
        for n, offset in zip(conn, CORNERS):
            if n != grid.get(tuple(a + b for a, b in zip(start, offset))):
                raise ValueError("Inverted/incorrect DC3D8 element")
        elements[e] = tuple(conn)
        cells.add(start)
        used.update(conn)
    if len(cells) != math.prod(mesh) or used != set(nodes):
        raise ValueError("Incomplete slab cells")
    for name, x in (("HOT", 0), ("COLD", g["length"])):
        ids = [int(v) for row in data[f"*NSET, NSET={name}"] for v in row.split(",")]
        if len(ids) != len(set(ids)) or set(ids) != {
            n for n, xyz in nodes.items() if math.isclose(xyz[0], x, abs_tol=1e-12)
        }:
            raise ValueError("Thermal boundary face membership incorrect")
    exact = {
        "*CONDUCTIVITY": [f'{config["material"]["conductivity"]:.16g}'],
        "*BOUNDARY": [
            f'HOT, 11, 11, {b["x_min_temperature"]:.16g}',
            f'COLD, 11, 11, {b["x_max_temperature"]:.16g}',
        ],
        "*INITIAL CONDITIONS, TYPE=TEMPERATURE": [
            f'ALL, {b["x_max_temperature"]:.16g}'
        ],
        "*HEAT TRANSFER, STEADY STATE": ["1, 1"],
        "*NODE PRINT, NSET=ALL": ["NT, RFL"],
        "*EL PRINT, ELSET=SLAB, GLOBAL=YES": ["HFL, COORD"],
        "*NODE FILE": ["NT, RFL"],
        "*EL FILE": ["HFL"],
    }
    if any(data[key] != value for key, value in exact.items()):
        raise ValueError(
            "Serialized conductivity, temperatures, procedure or output differs from specification"
        )
    return {
        "passed": True,
        "nodes": len(nodes),
        "elements": len(elements),
        "volume_m3": math.prod(lengths),
        "hot_nodes": (mesh[1] + 1) * (mesh[2] + 1),
        "cold_nodes": (mesh[1] + 1) * (mesh[2] + 1),
    }
