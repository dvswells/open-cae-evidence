"""Strict CalculiX modal DAT extraction."""

import math


def parse_dat(text, expected_nodes):
    """Strictly parse eigenvalues and explicitly numbered nodal mode blocks."""
    lines = text.splitlines()
    frequencies, eigenvalues, modes = [], [], []
    in_frequencies = False
    for line in lines:
        if "E I G E N V A L U E   O U T P U T" in line:
            if in_frequencies or frequencies:
                raise ValueError("Duplicate eigenvalue table")
            in_frequencies = True
        if "P A R T I C I P A T I O N" in line:
            in_frequencies = False
        parts = line.split()
        if in_frequencies and len(parts) == 5 and parts[0].isdigit():
            index = int(parts[0])
            lam, omega, frequency, imaginary = map(float, parts[1:])
            if index != len(frequencies) + 1 or not all(
                math.isfinite(v) for v in (lam, omega, frequency, imaginary)
            ):
                raise ValueError("Malformed eigenvalue table")
            if lam <= 0 or frequency <= 0 or imaginary != 0:
                raise ValueError(
                    "Nonpositive/complex mode; check constraints and material"
                )
            if not math.isclose(
                omega, 2 * math.pi * frequency, rel_tol=2e-6
            ) or not math.isclose(lam, omega**2, rel_tol=2e-6):
                raise ValueError("Eigenvalue/frequency unit consistency failure")
            eigenvalues.append(lam)
            frequencies.append(frequency)
    pending_mode = None
    i = 0
    while i < len(lines):
        line = lines[i]
        if "E I G E N V A L U E    N U M B E R" in line:
            pending_mode = int(line.split()[-1])
        if "displacements (vx,vy,vz) for set ALL" in line:
            if pending_mode != len(modes) + 1:
                raise ValueError("Missing/out-of-order mode shape number")
            displacements = {}
            i += 1
            while i < len(lines) and not lines[i].strip():
                i += 1
            while i < len(lines) and lines[i].strip():
                row = lines[i].split()
                if len(row) != 4 or not row[0].isdigit():
                    raise ValueError("Malformed displacement row")
                n, values = int(row[0]), tuple(float(v) for v in row[1:])
                if n in displacements or not all(math.isfinite(v) for v in values):
                    raise ValueError("Duplicate/nonfinite displacement")
                displacements[n] = values
                i += 1
            if set(displacements) != set(expected_nodes):
                raise ValueError("Incomplete nodal mode shape output")
            modes.append(displacements)
            pending_mode = None
        i += 1
    if (
        len(frequencies) != 6
        or len(modes) != 6
        or any(b <= a for a, b in zip(frequencies, frequencies[1:]))
    ):
        raise ValueError("Exactly six distinct, increasing complete modes required")
    return frequencies, eigenvalues, modes
