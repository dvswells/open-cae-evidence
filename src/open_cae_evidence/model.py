"""SI contract, structured C3D20R mesh, and independently audited solver deck."""

import math

UNITS = {"length": "m", "mass": "kg", "time": "s", "stress": "Pa", "density": "kg/m^3"}
OFFSETS = (
    (0, 0, 0), (2, 0, 0), (2, 2, 0), (0, 2, 0),
    (0, 0, 2), (2, 0, 2), (2, 2, 2), (0, 2, 2),
    (1, 0, 0), (2, 1, 0), (1, 2, 0), (0, 1, 0),
    (1, 0, 2), (2, 1, 2), (1, 2, 2), (0, 1, 2),
    (0, 0, 1), (2, 0, 1), (2, 2, 1), (0, 2, 1),
)


def validate(config):
    """Reject invalid/mislabeled inputs, not silently convert them.

    Bounds define this metallic slender-beam benchmark's domain, not universal
    material laws. A plausible number with a false unit label cannot be detected.
    """
    if set(config) != {"units", "geometry", "material", "boundary", "meshes"}:
        raise ValueError("Unexpected or missing configuration keys")
    if config["units"] != UNITS:
        raise ValueError("Explicit SI units are required")
    g, m = config["geometry"], config["material"]
    if set(g) != {"length", "width", "height"} or set(m) != {"young_modulus", "poisson_ratio", "density"}:
        raise ValueError("Unexpected geometry/material keys")
    for value in (*g.values(), *m.values()):
        if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
            raise ValueError("Physical inputs must be finite numbers")
    if not all(0.001 <= v <= 10 for v in g.values()):
        raise ValueError("Dimensions outside benchmark SI domain [0.001, 10] m")
    if g["length"] / max(g["width"], g["height"]) < 50:
        raise ValueError("Beam must be slender: L/max(b,h) >= 50")
    if not 1e9 <= m["young_modulus"] <= 1e12 or not 1000 <= m["density"] <= 30000:
        raise ValueError("Metallic benchmark expects E in Pa and density in kg/m^3")
    if not -1 < m["poisson_ratio"] < 0.5:
        raise ValueError("Isotropic elastic stability requires -1 < nu < 0.5")
    if (config["boundary"] != {"plane": "x=0", "dofs": [1, 2, 3]}
            or any(type(d) is not int for d in config["boundary"]["dofs"])):
        raise ValueError("Clamp x=0 in all three translational DOFs")
    meshes = config["meshes"]
    if not isinstance(meshes, list) or len(meshes) < 3:
        raise ValueError("At least three refinement levels are required")
    previous = None
    for mesh in meshes:
        if (not isinstance(mesh, list) or len(mesh) != 3
                or any(type(n) is not int or n < 2 or n > 128 for n in mesh)
                or any(n % 2 for n in mesh)):
            raise ValueError("Meshes require three even integer counts in [2,128]")
        if math.prod(mesh) > 10000:
            raise ValueError("Mesh exceeds small benchmark resource budget")
        if previous and any(a <= b for a, b in zip(mesh, previous)):
            raise ValueError("Refine all three directions monotonically")
        previous = mesh


def make_deck(config, mesh):
    validate(config)
    nx, ny, nz = mesh
    g, m = config["geometry"], config["material"]
    nodes, ids, elements = {}, {}, []
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                connectivity = []
                for a, b, c in OFFSETS:
                    key = (2 * i + a, 2 * j + b, 2 * k + c)
                    if key not in ids:
                        node = len(ids) + 1
                        ids[key] = node
                        nodes[node] = (key[0] * g["length"] / (2 * nx),
                                       key[1] * g["width"] / (2 * ny) - g["width"] / 2,
                                       key[2] * g["height"] / (2 * nz) - g["height"] / 2)
                    connectivity.append(ids[key])
                elements.append(connectivity)
    fixed = [n for n, (x, _, _) in nodes.items() if x == 0]
    lines = ["*HEADING", "Open CAE Evidence; SI m kg s Pa; C3D20R cantilever", "*NODE, NSET=ALL"]
    lines.extend(f"{n}, {x:.16g}, {y:.16g}, {z:.16g}" for n, (x, y, z) in nodes.items())
    lines.append("*ELEMENT, TYPE=C3D20R, ELSET=BEAM")
    for e, conn in enumerate(elements, 1):
        lines.append(", ".join(map(str, [e, *conn[:15]])) + ",")
        lines.append(", ".join(map(str, conn[15:])))
    lines.append("*NSET, NSET=FIXED")
    for i in range(0, len(fixed), 16):
        lines.append(", ".join(map(str, fixed[i:i + 16])))
    lines += ["*MATERIAL, NAME=STEEL", "*ELASTIC",
              f'{m["young_modulus"]:.16g}, {m["poisson_ratio"]:.16g}',
              "*DENSITY", f'{m["density"]:.16g}',
              "*SOLID SECTION, ELSET=BEAM, MATERIAL=STEEL", "*BOUNDARY", "FIXED, 1, 3",
              "*STEP", "*FREQUENCY", "6", "*NODE PRINT, NSET=ALL", "U", "*NODE FILE", "U", "*END STEP"]
    deck = "\n".join(lines) + "\n"
    audit = audit_deck(deck, config, mesh)
    return deck, nodes, audit


def audit_deck(deck, config, mesh):
    """Parse the serialized deck rather than trusting the mesh generator's state.

    Strictly accepts this benchmark's keyword contract. Checks all C3D20R
    connectivity/edge midpoints, positive volume, dimensions, material and clamp.
    """
    sections = []
    for raw in deck.splitlines():
        line = raw.strip()
        if not line or line.startswith("**"):
            continue
        if line.startswith("*"):
            sections.append((line.upper(), []))
        else:
            if not sections:
                raise ValueError("Data before keyword")
            sections[-1][1].append(line)
    expected = ["*HEADING", "*NODE, NSET=ALL", "*ELEMENT, TYPE=C3D20R, ELSET=BEAM",
                "*NSET, NSET=FIXED", "*MATERIAL, NAME=STEEL", "*ELASTIC", "*DENSITY",
                "*SOLID SECTION, ELSET=BEAM, MATERIAL=STEEL", "*BOUNDARY", "*STEP",
                "*FREQUENCY", "*NODE PRINT, NSET=ALL", "*NODE FILE", "*END STEP"]
    if [s[0] for s in sections] != expected:
        raise ValueError("Solver deck keyword contract violated")
    data = {key: rows for key, rows in sections}
    nodes = {}
    for row in data["*NODE, NSET=ALL"]:
        n, *xyz = row.split(",")
        if len(xyz) != 3 or int(n) in nodes:
            raise ValueError("Invalid/duplicate node")
        nodes[int(n)] = tuple(float(v) for v in xyz)
    if not nodes or not all(math.isfinite(v) for xyz in nodes.values() for v in xyz):
        raise ValueError("Invalid coordinates")
    g, m = config["geometry"], config["material"]
    limits = [(0, g["length"]), (-g["width"] / 2, g["width"] / 2), (-g["height"] / 2, g["height"] / 2)]
    for axis, (low, high) in enumerate(limits):
        coords = [xyz[axis] for xyz in nodes.values()]
        if not math.isclose(min(coords), low, abs_tol=1e-12) or not math.isclose(max(coords), high, abs_tol=1e-12):
            raise ValueError("Mesh dimensions differ from physical inputs")
    fixed = [int(n) for row in data["*NSET, NSET=FIXED"] for n in row.split(",")]
    if len(fixed) != len(set(fixed)) or set(fixed) != {n for n, xyz in nodes.items() if abs(xyz[0]) < 1e-12}:
        raise ValueError("Clamp must contain exactly all x=0 face nodes")
    if data["*BOUNDARY"] != ["FIXED, 1, 3"]:
        raise ValueError("Incorrect constrained DOFs")
    if data["*FREQUENCY"] != ["6"] or data["*NODE PRINT, NSET=ALL"] != ["U"] or data["*NODE FILE"] != ["U"]:
        raise ValueError("Six modes and complete displacement output required")
    elastic = [float(v) for v in data["*ELASTIC"][0].split(",")]
    if elastic != [m["young_modulus"], m["poisson_ratio"]] or data["*DENSITY"] != [f'{m["density"]:.16g}']:
        raise ValueError("Serialized material differs from inputs")
    rows = data["*ELEMENT, TYPE=C3D20R, ELSET=BEAM"]
    if len(rows) != 2 * math.prod(mesh):
        raise ValueError("Incorrect element count")
    used, volume, element_ids = set(), 0.0, set()
    cells = set()
    for i in range(0, len(rows), 2):
        values = [int(v.strip()) for v in (rows[i] + rows[i + 1]).split(",")]
        e, *conn = values
        if e in element_ids or len(conn) != 20 or len(set(conn)) != 20:
            raise ValueError("Invalid/duplicate quadratic element")
        element_ids.add(e)
        try:
            points = [nodes[n] for n in conn]
        except KeyError as exc:
            raise ValueError("Missing node") from exc
        start = points[0]
        end = points[6]
        lengths = [b - a for a, b in zip(start, end)]
        if any(not math.isclose(d, (high - low) / count, rel_tol=1e-10) for d, (low, high), count in zip(lengths, limits, mesh)):
            raise ValueError("Element inverted or not a structured rectangular cell")
        cell = tuple(round((v - low) / d) for v, (low, _), d in zip(start, limits, lengths))
        if cell in cells or any(c < 0 or c >= count for c, count in zip(cell, mesh)):
            raise ValueError("Overlapping/out-of-domain element")
        cells.add(cell)
        for point, offset in zip(points, OFFSETS):
            if any(not math.isclose(v, a + d * q / 2, abs_tol=1e-12) for v, a, d, q in zip(point, start, lengths, offset)):
                raise ValueError("Incorrect C3D20R connectivity/midside coordinates")
        used.update(conn)
        volume += math.prod(lengths)
    if used != set(nodes) or not math.isclose(volume, math.prod(g.values()), rel_tol=1e-10):
        raise ValueError("Mesh volume or node coverage mismatch")
    return {"nodes": len(nodes), "elements": len(element_ids), "fixed_nodes": len(fixed),
            "volume_m3": volume, "mass_kg": volume * m["density"], "passed": True}
