"""CalculiX thermal DAT extraction and ASCII FRD cross-file semantics."""

import math

from .frd import vector_row


def parse_dat(text, nodes, elements):
    names = {
        "temperatures for set ALL": ("temperature", 2),
        "heat generation for set ALL": ("heat_rate", 2),
        "heat flux (elem, integ.pnt.,qx,qy,qz)": ("heat_flux", 5),
        "global coordinates (elem, integ.pnt.,x,y,z)": ("integration_coordinates", 5),
    }
    data = {}
    lines = text.splitlines()
    i = 0
    while i < len(lines):
        line = lines[i].strip()
        matched = [info for prefix, info in names.items() if line.startswith(prefix)]
        if not line:
            i += 1
            continue
        if len(matched) != 1 or "and time" not in line or float(line.split()[-1]) != 1:
            raise ValueError("Unexpected thermal DAT field, time or trailing data")
        name, width = matched[0]
        if name in data:
            raise ValueError(
                "Repeated thermal output; only the single final steady state is supported"
            )
        records = {}
        i += 1
        while i < len(lines) and not lines[i].strip():
            i += 1
        while i < len(lines) and lines[i].strip():
            row = lines[i].split()
            if len(row) != width:
                raise ValueError("Malformed thermal DAT row")
            key = int(row[0]) if width == 2 else (int(row[0]), int(row[1]))
            values = [float(v) for v in row[1 if width == 2 else 2 :]]
            if key in records or not all(math.isfinite(v) for v in values):
                raise ValueError("Duplicate/nonfinite thermal DAT values")
            records[key] = values[0] if width == 2 else tuple(values)
            i += 1
        data[name] = records
    if set(data) != {
        "temperature",
        "heat_rate",
        "heat_flux",
        "integration_coordinates",
    }:
        raise ValueError("Incomplete thermal DAT quantities")
    if set(data["temperature"]) != set(nodes) or set(data["heat_rate"]) != set(nodes):
        raise ValueError("Incomplete thermal nodal output")
    expected = {(e, p) for e in elements for p in range(1, 9)}
    if (
        set(data["heat_flux"]) != expected
        or set(data["integration_coordinates"]) != expected
    ):
        raise ValueError("Incomplete eight-point DC3D8 output")
    return data


def verify_frd(text, nodes, elements, data, version):
    lines = text.splitlines()
    i = 0
    xyz = {}
    topology = {}
    fields = {}
    ended = False
    versions = []
    while i < len(lines):
        line = lines[i]
        tag = line.strip().split()[0] if line.strip() else ""
        if line.startswith("    1UVERSION"):
            versions.append(line.split("Version", 1)[1].strip())
        elif tag in ("2C", "3C"):
            count = int(line.split()[1])
            i += 1
            target = xyz if tag == "2C" else topology
            if target or line.split()[-1] != "1":
                raise ValueError("Duplicate/unsupported thermal FRD mesh block")
            while i < len(lines) and lines[i][:3] == " -1":
                if tag == "2C":
                    n, values = vector_row(lines[i])
                    i += 1
                else:
                    row = lines[i].split()
                    if len(row) != 5 or row[2:] != ["1", "0", "1"]:
                        raise ValueError("Thermal FRD requires HEX8")
                    n = int(row[1])
                    i += 1
                    if i >= len(lines) or lines[i][:3] != " -2":
                        raise ValueError("Missing thermal FRD connectivity")
                    values = tuple(int(v) for v in lines[i][3:].split())
                    i += 1
                if n in target:
                    raise ValueError("Duplicate thermal FRD mesh ID")
                target[n] = values
            if len(target) != count or i >= len(lines) or lines[i].strip() != "-3":
                raise ValueError("Incomplete thermal FRD mesh")
        elif tag == "100CL":
            row = line.split()
            if (
                len(row) != 7
                or row[1] != "101"
                or float(row[2]) != 1
                or int(row[3]) != len(nodes)
                or row[4:] != ["0", "1", "1"]
            ):
                raise ValueError("Thermal FRD step/time/count mismatch")
            i += 1
            header = lines[i].split()
            i += 1
            if header not in (
                ["-4", "NDTEMP", "1", "1"],
                ["-4", "RFL", "1", "1"],
                ["-4", "FLUX", "4", "1"],
            ):
                raise ValueError("Unknown thermal FRD dataset")
            name = header[1]
            labels = {
                "NDTEMP": [["-5", "T", "1", "1", "0", "0"]],
                "RFL": [["-5", "RFL", "1", "1", "0", "0"]],
                "FLUX": [
                    ["-5", "F1", "1", "2", "1", "0"],
                    ["-5", "F2", "1", "2", "2", "0"],
                    ["-5", "F3", "1", "2", "3", "0"],
                    ["-5", "ALL", "1", "2", "0", "0", "1ALL"],
                ],
            }
            if [s.split() for s in lines[i : i + len(labels[name])]] != labels[name]:
                raise ValueError("Thermal FRD component semantics differ")
            i += len(labels[name])
            values = {}
            while i < len(lines) and lines[i][:3] == " -1":
                if name == "FLUX":
                    n, v = vector_row(lines[i])
                else:
                    if len(lines[i].rstrip()) != 25:
                        raise ValueError("Bad thermal FRD scalar width")
                    n = int(lines[i][3:13])
                    v = float(lines[i][13:25])
                    if not math.isfinite(v):
                        raise ValueError("Nonfinite thermal FRD scalar")
                if n in values:
                    raise ValueError("Duplicate thermal FRD value")
                values[n] = v
                i += 1
            if (
                name in fields
                or set(values) != set(nodes)
                or i >= len(lines)
                or lines[i].strip() != "-3"
            ):
                raise ValueError("Incomplete/duplicate thermal FRD dataset")
            fields[name] = values
        elif tag == "9999":
            if i != len(lines) - 1:
                raise ValueError("Trailing thermal FRD data")
            ended = True
        elif not line.strip() or tag == "1C" or tag.startswith("1U") or tag == "1PSTEP":
            pass
        else:
            raise ValueError(f"Unsupported thermal FRD record {tag}")
        i += 1
    if (
        not ended
        or versions != [version]
        or set(xyz) != set(nodes)
        or topology != elements
        or set(fields) != {"NDTEMP", "FLUX", "RFL"}
    ):
        raise ValueError("Thermal FRD mesh, version or quantity coverage mismatch")
    for n, coord in xyz.items():
        if any(
            not math.isclose(a, b, rel_tol=1e-5, abs_tol=1e-12)
            for a, b in zip(coord, nodes[n])
        ):
            raise ValueError("Thermal FRD coordinates differ")
    for field, key in (("NDTEMP", "temperature"), ("RFL", "heat_rate")):
        scale = max(abs(v) for v in data[key].values())
        if any(
            not math.isclose(
                fields[field][n], value, rel_tol=1e-5, abs_tol=1e-12 * scale
            )
            for n, value in data[key].items()
        ):
            raise ValueError(f"Thermal FRD {field} differs from DAT")
    # FLUX is extrapolated/averaged, not the raw integration-point quantity.
    # CalculiX HEX8 extrapolation gives about 4e-5 relative bias for constant flux.
    average = [
        sum(q[a] for q in data["heat_flux"].values()) / len(data["heat_flux"])
        for a in range(3)
    ]
    scale = max(abs(v) for v in average)
    if any(
        not math.isclose(q[a], average[a], rel_tol=1e-4, abs_tol=1e-10 * scale)
        for q in fields["FLUX"].values()
        for a in range(3)
    ):
        raise ValueError("Thermal FRD extrapolated FLUX differs from uniform DAT field")
    return {
        "passed": True,
        "format": "CalculiX ASCII FRD",
        "solver_version": version,
        "nodes": len(nodes),
        "elements": len(elements),
        "datasets": ["NDTEMP", "FLUX", "RFL"],
        "nodal_flux_is_extrapolated": True,
        "nodal_flux_relative_tolerance": 1e-4,
    }
