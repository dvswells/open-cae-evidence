"""Reproducible pointwise evidence at every nonnodal volume quadrature location."""

import csv
import hashlib
import io
from . import q1, source_reference
from .results import mesh_maps, thermal_fields


def rows(config, normalized):
    nodes, elements = mesh_maps(normalized)
    data = thermal_fields(normalized)
    k = config["material"]["conductivity"]
    yield [
        "element",
        "sample",
        "x_m",
        "y_m",
        "z_m",
        "volume_weight_m3",
        "T_FE_K",
        "T_exact_K",
        "T_error_K",
        "grad_x_FE_K_m",
        "grad_y_FE_K_m",
        "grad_z_FE_K_m",
        "grad_x_exact_K_m",
        "qx_raw_W_m2",
        "qy_raw_W_m2",
        "qz_raw_W_m2",
        "qx_exact_W_m2",
    ]
    for e, conn in sorted(elements.items()):
        start, span = q1.cell(nodes, conn)
        Tvalues = tuple(data["temperature"][n] for n in conn)
        for p, (local, weight) in enumerate(q1.G3, 1):
            xyz = tuple(a + d * (v + 1) / 2 for a, d, v in zip(start, span, local))
            T, grad = q1.evaluate(Tvalues, span, local)
            exact = source_reference.analytical(config, xyz[0])
            import math

            yield [
                e,
                p,
                *xyz,
                weight * math.prod(span) / 8,
                T,
                exact["temperature_K"],
                T - exact["temperature_K"],
                *grad,
                exact["gradient_K_m"][0],
                *(-k * v for v in grad),
                exact["heat_flux_W_m2"][0],
            ]


def write(path, config, normalized):
    with path.open("w", newline="") as handle:
        csv.writer(handle).writerows(rows(config, normalized))


def replay(path, config, normalized):
    # Hash the deterministic reconstructed CSV stream, without keeping a large
    # volume dataset in memory. This checks semantics even after manifest rehash.
    digest = hashlib.sha256()
    buffer = io.StringIO(newline="")
    writer = csv.writer(buffer)
    for row in rows(config, normalized):
        writer.writerow(row)
        digest.update(buffer.getvalue().encode())
        buffer.seek(0)
        buffer.truncate(0)
    if hashlib.sha256(path.read_bytes()).digest() != digest.digest():
        raise ValueError("Nonnodal sample evidence differs from numerical replay")
