"""Strict ASCII VTU semantics; native L2 nodal flux -> Q1 Gauss samples.

The transformation uses real FluxSolver output, never the analytical reference.
"""

import itertools
import math
import re
import xml.etree.ElementTree as ET

from .elmer_input import CORNERS


def check_log(log):
    versions = re.findall(r"ELMER SOLVER \(v ([^)]+)\) STARTED", log)
    bad = (
        "ERROR",
        "FATAL",
        "did not converge",
        "Too many iterations",
        "diverged",
        "Iteration halted",
        "WARNING",
    )
    if (
        len(versions) != 1
        or log.count("ELMER SOLVER FINISHED") != 1
        or any(x.lower() in log.lower() for x in bad)
    ):
        raise RuntimeError(
            "Elmer failed, non-converged, warning-bearing or incomplete log"
        )
    if not all(x in log for x in ("HeatSolve", "FluxSolver", "ResultOutput")):
        raise RuntimeError("Missing thermal, flux or output solver execution")
    if "MAIN: Version: 9.0 (Rev: 9352634," not in log:
        raise RuntimeError(
            "Native Elmer version/source revision differs from pinned runtime"
        )
    if log.count("HeatSolve: Assembly done") != 1:
        raise RuntimeError("Missing/ambiguous heat linear solve")
    block = log.split("HeatSolve: Assembly done", 1)[-1].split("ComputeChange:", 1)[0]
    residuals = re.findall(r"^\s+\d+\s+([\d.]+E[+-]\d+)\s*$", block, re.M | re.I)
    if not residuals or float(residuals[-1]) > 1e-12:
        raise RuntimeError("Heat linear system lacks a converged residual")
    return versions[0].strip()


def array(element, components, integer=False):
    if (
        element is None
        or (element.get("format") or element.get("FORMAT")) != "ascii"
        or int(element.get("NumberOfComponents", "1")) != components
    ):
        raise ValueError("Missing/unsupported native VTU field encoding/components")
    dtype = element.get("type")
    if dtype not in (("Int32", "UInt8") if integer else ("Float64",)):
        raise ValueError("Native field precision/type differs")
    try:
        data = [int(x) if integer else float(x) for x in (element.text or "").split()]
    except ValueError as exc:
        raise ValueError("Malformed native numerical field") from exc
    if len(data) % components or not all(map(math.isfinite, data)):
        raise ValueError("Incomplete/nonfinite native field")
    return [tuple(data[i : i + components]) for i in range(0, len(data), components)]


def parse_vtu(text, nodes, elements):
    try:
        root = ET.fromstring(text)
    except ET.ParseError as exc:
        raise ValueError("Malformed native VTU XML") from exc
    pieces = root.findall("./UnstructuredGrid/Piece")
    if root.get("type") != "UnstructuredGrid" or len(pieces) != 1:
        raise ValueError("Unexpected native mesh type/partitioning")
    piece = pieces[0]
    nn, ne = len(nodes), len(elements)
    if (
        int(piece.get("NumberOfPoints", "-1")) != nn
        or int(piece.get("NumberOfCells", "-1")) != ne
    ):
        raise ValueError("Native result mesh coverage differs")
    coordinates = array(piece.find("./Points/DataArray"), 3)
    if len(coordinates) != nn:
        raise ValueError("Incomplete native coordinates")
    # Coordinate matching, deliberately independent of input or output node IDs.
    expected = {tuple(round(x, 12) for x in xyz): n for n, xyz in nodes.items()}
    mapping = []
    for xyz in coordinates:
        key = tuple(round(x, 12) for x in xyz)
        if key not in expected or any(
            abs(a - b) > 1e-12 for a, b in zip(xyz, nodes[expected[key]])
        ):
            raise ValueError("Native coordinate system/units differ")
        mapping.append(expected[key])
    if len(set(mapping)) != nn:
        raise ValueError("Duplicate/incomplete native physical coordinates")
    cells = {a.get("Name"): a for a in piece.findall("./Cells/DataArray")}
    if len(cells) != len(piece.findall("./Cells/DataArray")):
        raise ValueError("Duplicate native topology arrays")
    offsets = [x[0] for x in array(cells.get("offsets"), 1, True)]
    kinds = [x[0] for x in array(cells.get("types"), 1, True)]
    conn = [x[0] for x in array(cells.get("connectivity"), 1, True)]
    if (
        offsets != list(range(8, 8 * ne + 1, 8))
        or kinds != [12] * ne
        or len(conn) != 8 * ne
        or any(x < 0 or x >= nn for x in conn)
    ):
        raise ValueError("Native output is not complete VTK HEX8 topology")
    actual = [
        tuple(mapping[x] for x in conn[i : i + 8]) for i in range(0, len(conn), 8)
    ]
    if len(set(actual)) != ne or set(actual) != set(elements.values()):
        raise ValueError(
            "Native result connectivity/orientation differs from audited mesh"
        )
    bodies = piece.findall('./CellData/DataArray[@Name="GeometryIds"]')
    if len(bodies) != 1 or array(bodies[0], 1, True) != [(1,)] * ne:
        raise ValueError("Native result body assignments differ")
    fields = {a.get("Name"): a for a in piece.findall("./PointData/DataArray")}
    if len(fields) != len(piece.findall("./PointData/DataArray")):
        raise ValueError("Duplicate native fields")
    temperature = array(fields.get("temperature"), 1)
    loads = array(fields.get("temperature loads"), 1)
    flux = array(fields.get("temperature flux"), 3)
    if any(len(x) != nn for x in (temperature, loads, flux)):
        raise ValueError("Incomplete native temperature/loads/flux")
    data = {
        "temperature": dict(zip(mapping, (v[0] for v in temperature))),
        "heat_rate": dict(zip(mapping, (v[0] for v in loads))),
        "nodal_flux": dict(zip(mapping, flux)),
        "heat_flux": {},
        "integration_coordinates": {},
    }
    # Native FluxSolver recovers q by global L2 projection into continuous Q1.
    # Sample that native recovered field at the common 2x2x2 physical Gauss points.
    # These samples are postprocessed_from_solver, NOT native IP quantities.
    gp = tuple(itertools.product((-1 / math.sqrt(3), 1 / math.sqrt(3)), repeat=3))
    signs = [tuple(2 * x - 1 for x in c) for c in CORNERS]
    for e, connectivity in elements.items():
        for p, local in enumerate(gp, 1):
            weights = [
                math.prod(1 + s * x for s, x in zip(sign, local)) / 8 for sign in signs
            ]
            data["heat_flux"][e, p] = tuple(
                sum(w * data["nodal_flux"][n][a] for w, n in zip(weights, connectivity))
                for a in range(3)
            )
            data["integration_coordinates"][e, p] = tuple(
                sum(w * nodes[n][a] for w, n in zip(weights, connectivity))
                for a in range(3)
            )
    return data
