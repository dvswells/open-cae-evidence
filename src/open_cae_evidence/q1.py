"""Affine rectangular HEX8 evaluation from genuine degrees of freedom.

No analytical reference is evaluated here. Coordinates/order are physical.
"""

import itertools
import math
from functools import lru_cache

SIGNS = (
    (-1, -1, -1),
    (1, -1, -1),
    (1, 1, -1),
    (-1, 1, -1),
    (-1, -1, 1),
    (1, -1, 1),
    (1, 1, 1),
    (-1, 1, 1),
)
G2 = tuple(itertools.product((-1 / math.sqrt(3), 1 / math.sqrt(3)), repeat=3))
G3 = tuple(
    (p, math.prod(w))
    for p, w in (
        (tuple(v[0] for v in row), tuple(v[1] for v in row))
        for row in itertools.product(
            ((-math.sqrt(3 / 5), 5 / 9), (0.0, 8 / 9), (math.sqrt(3 / 5), 5 / 9)),
            repeat=3,
        )
    )
)


@lru_cache(maxsize=128)
def basis(local):
    weights = tuple(
        math.prod(1 + s * x for s, x in zip(sign, local)) / 8 for sign in SIGNS
    )
    derivative = tuple(
        tuple(
            sign[a] * math.prod(1 + sign[b] * local[b] for b in range(3) if b != a) / 8
            for sign in SIGNS
        )
        for a in range(3)
    )
    return weights, derivative


def cell(nodes, conn):
    start = nodes[conn[0]]
    span = tuple(b - a for a, b in zip(start, nodes[conn[6]]))
    if len(conn) != 8 or len(set(conn)) != 8 or any(d <= 0 for d in span):
        raise ValueError("Invalid affine HEX8 orientation/topology")
    for n, sign in zip(conn, SIGNS):
        if any(
            not math.isclose(
                nodes[n][a],
                start[a] + span[a] * (sign[a] + 1) / 2,
                rel_tol=0,
                abs_tol=1e-12,
            )
            for a in range(3)
        ):
            raise ValueError("Q1 evaluator requires canonical affine bricks")
    return start, span


def evaluate(temperatures, span, local):
    weights, derivatives = basis(tuple(local))
    value = sum(w * t for w, t in zip(weights, temperatures))
    gradient = tuple(
        2 / span[a] * sum(d * t for d, t in zip(derivatives[a], temperatures))
        for a in range(3)
    )
    return value, gradient


def vector(values, local):
    weights, _ = basis(tuple(local))
    return tuple(sum(w * v[a] for w, v in zip(weights, values)) for a in range(3))


def gp_samples(data, nodes, conn, e):
    start, span = cell(nodes, conn)
    samples = {}
    for p in range(1, 9):
        xyz = data["integration_coordinates"][e, p]
        local = tuple(2 * (x - a) / d - 1 for x, a, d in zip(xyz, start, span))
        point = tuple((1 if x > 0 else -1) / math.sqrt(3) for x in local)
        if point in samples or any(
            not math.isclose(x, y, rel_tol=3e-5, abs_tol=1e-6)
            for x, y in zip(local, point)
        ):
            raise ValueError("Incomplete/wrong physical Gauss mapping")
        samples[point] = data["heat_flux"][e, p]
    return samples


@lru_cache(maxsize=128)
def gp_weights(local):
    root = 1 / math.sqrt(3)
    return tuple(
        math.prod((1 + p[a] * local[a] / (root * root)) / 2 for a in range(3))
        for p in G2
    )


def gp_vector(samples, local):
    weights = gp_weights(tuple(local))
    return tuple(sum(w * samples[p][a] for w, p in zip(weights, G2)) for a in range(3))


def consistent_source_load(nodes, elements, source):
    """Integrate N_i*qdot on actual affine cells. No boundary/reference flow."""
    rows = {row["element"]: row["value"] for row in source["element_values"]}
    if (
        len(rows) != len(source["element_values"])
        or set(rows) != set(elements)
        or source["unit"] != "W/m^3"
    ):
        raise ValueError("Incomplete/ambiguous source coverage/units")
    loads = {n: 0.0 for n in nodes}
    for e, conn in elements.items():
        _, span = cell(nodes, conn)
        load = rows[e] * math.prod(span) / 8
        for n in conn:
            loads[n] += load
    return loads
