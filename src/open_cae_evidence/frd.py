"""Semantic verification of CalculiX's ASCII FRD, independent of file hashes.

This deliberately supports only the C3D20R/DISP modal export used here. MASS
auxiliary nodes/elements are absent from CalculiX FRD; their output is checked
in DAT. Unsupported binary or other datasets fail rather than silently pass.
"""

import math
import re

FRD_RELATIVE = 1e-5  # FRD prints six significant digits; DAT prints seven.


def vector_row(line):
    if line[:3] != " -1" or len(line.rstrip()) != 49:
        raise ValueError("Malformed fixed-width FRD vector record")
    node = int(line[3:13])
    values = tuple(float(line[i:i+12]) for i in (13, 25, 37))
    if not all(math.isfinite(v) for v in values):
        raise ValueError("Nonfinite FRD data")
    return node, values


def solid_connectivity(deck):
    rows = deck.split("*ELEMENT, TYPE=C3D20R, ELSET=BEAM\n", 1)[1].split("*", 1)[0].strip().splitlines()
    result = {}
    # FRD HEX20 stores vertical edge midsides before top-face edge midsides.
    permutation = (*range(12), *range(16, 20), *range(12, 16))
    for i in range(0, len(rows), 2):
        e, *conn = [int(v) for v in (rows[i] + rows[i+1]).split(",")]
        result[e] = tuple(conn[k] for k in permutation)
    return result


def verify_frd(text, deck, nodes, frequencies, eigenvalues, dat_modes, solver_version, auxiliary=None):
    lines = text.splitlines()
    expected_nodes = set(nodes) - set(auxiliary or ())
    expected_elements = solid_connectivity(deck)
    coordinates, elements, modes = {}, {}, []
    seen_coordinates = seen_elements = False
    mode_number = stiffness = frequency = None
    version = None
    i = 0
    ended = False
    while i < len(lines):
        line = lines[i]
        tag = line.strip().split()[0] if line.strip() else ""
        if line.startswith("    1UVERSION"):
            if version is not None:
                raise ValueError("Duplicate FRD solver version")
            version = line.split("Version", 1)[1].strip()
        elif tag == "2C":
            if seen_coordinates or line.split()[-1] != "1":
                raise ValueError("Duplicate or unsupported FRD coordinate block")
            count = int(line.split()[1]); seen_coordinates = True
            i += 1
            while i < len(lines) and lines[i][:3] == " -1":
                n, xyz = vector_row(lines[i])
                if n in coordinates:
                    raise ValueError("Duplicate FRD node")
                coordinates[n] = xyz; i += 1
            if i >= len(lines) or lines[i].strip() != "-3" or len(coordinates) != count:
                raise ValueError("Incomplete FRD coordinates")
        elif tag == "3C":
            if seen_elements or line.split()[-1] != "1":
                raise ValueError("Duplicate or unsupported FRD topology block")
            count = int(line.split()[1]); seen_elements = True
            i += 1
            while i < len(lines) and lines[i][:3] == " -1":
                parts = lines[i].split()
                if len(parts) != 5 or parts[2:] != ["4", "0", "1"]:
                    raise ValueError("Expected FRD HEX20 element type")
                e = int(parts[1]); conn = []
                i += 1
                while i < len(lines) and lines[i][:3] == " -2":
                    conn.extend(int(v) for v in lines[i][3:].split()); i += 1
                if e in elements or len(conn) != 20:
                    raise ValueError("Duplicate/incomplete FRD element")
                elements[e] = tuple(conn)
            if i >= len(lines) or lines[i].strip() != "-3" or len(elements) != count:
                raise ValueError("Incomplete FRD topology")
        elif tag == "1PMODE":
            if mode_number is not None:
                raise ValueError("Duplicate FRD mode metadata")
            mode_number = int(line.split()[1])
        elif tag == "1PGK":
            if stiffness is not None:
                raise ValueError("Duplicate FRD eigenvalue metadata")
            stiffness = float(line.split()[1])
        elif tag.startswith("100C"):
            fields = line.split()
            number = len(modes) + 1
            if (len(fields) != 7 or fields[0] != "100CL" or fields[4] != "2"
                    or fields[5] != f"{number}MODAL" or fields[6] != "1"
                    or int(fields[1]) != 100 + number or mode_number != number
                    or int(fields[3]) != len(expected_nodes) or number > 6):
                raise ValueError("FRD dataset type/count/mode identity mismatch")
            frequency = float(fields[2])
            if (stiffness is None or not math.isfinite(frequency) or not math.isfinite(stiffness)
                    or not math.isclose(frequency, frequencies[number-1], rel_tol=2e-6)
                    or not math.isclose(stiffness, eigenvalues[number-1], rel_tol=2e-6)):
                raise ValueError("FRD eigenvalue/frequency differs from DAT")
            i += 1
            metadata = [" -4  DISP        4    1", " -5  D1          1    2    1    0",
                        " -5  D2          1    2    2    0", " -5  D3          1    2    3    0",
                        " -5  ALL         1    2    0    0    1ALL"]
            if lines[i:i+5] != metadata:
                raise ValueError("FRD displacement component semantics differ")
            i += 5
            values = {}
            while i < len(lines) and lines[i][:3] == " -1":
                n, xyz = vector_row(lines[i])
                if n in values:
                    raise ValueError("Duplicate FRD displacement node")
                values[n] = xyz; i += 1
            if i >= len(lines) or lines[i].strip() != "-3" or set(values) != expected_nodes:
                raise ValueError("Incomplete FRD mode shape node coverage")
            peak = max(abs(v) for xyz in dat_modes[number-1].values() for v in xyz)
            for n, xyz in values.items():
                if any(not math.isclose(a, b, rel_tol=FRD_RELATIVE, abs_tol=1e-12*peak)
                       for a, b in zip(xyz, dat_modes[number-1][n])):
                    raise ValueError("FRD displacement differs from DAT (sign, scale, component or mode)")
            modes.append(values)
            mode_number = stiffness = frequency = None
        elif tag == "9999":
            if i != len(lines)-1:
                raise ValueError("Data after FRD end marker")
            ended = True
        elif (not line.strip() or tag == "1C" or tag.startswith("1U")
              or tag in {"1PSTEP", "1PGM", "1PHID", "1PSUBC"}):
            pass
        else:
            raise ValueError(f"Unexpected FRD record: {tag}")
        i += 1
    if (not ended or version != solver_version or set(coordinates) != expected_nodes
            or elements != expected_elements or len(modes) != 6 or mode_number is not None or stiffness is not None):
        raise ValueError("FRD version, geometry, topology, modal coverage or termination mismatch")
    for n, xyz in coordinates.items():
        if any(not math.isclose(a, b, rel_tol=FRD_RELATIVE, abs_tol=1e-12) for a, b in zip(xyz, nodes[n])):
            raise ValueError("FRD coordinates differ from audited deck")
    return {"passed": True, "format": "CalculiX ASCII FRD", "solver_version": version,
            "solid_nodes": len(coordinates), "solid_elements": len(elements), "modal_datasets": len(modes),
            "auxiliary_nodes_checked_in_dat": sorted(auxiliary or ()), "dat_frd_relative_tolerance": FRD_RELATIVE}
