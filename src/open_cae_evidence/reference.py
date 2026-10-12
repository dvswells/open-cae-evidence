"""Independent Euler-Bernoulli equations; no FEM data used to fit tolerances."""

import math

# Positive roots of cos(beta)*cosh(beta) = -1 (dimensionless).
BETAS = (1.875104068711961, 4.694091132974174, 7.854757438237613,
         10.995540734875467, 14.13716839104647, 17.278759532088237)


def references(config):
    g, m = config["geometry"], config["material"]
    length, b, h = g["length"], g["width"], g["height"]
    area = b * h
    modes = []
    for direction, inertia in (("z", b * h**3 / 12), ("y", h * b**3 / 12)):
        for order, beta in enumerate(BETAS, 1):
            frequency = beta**2 / (2 * math.pi * length**2) * math.sqrt(m["young_modulus"] * inertia / (m["density"] * area))
            modes.append({"direction": direction, "order": order, "beta": beta, "frequency_hz": frequency})
    return sorted(modes, key=lambda mode: mode["frequency_hz"])[:6]


def shape(beta, xi):
    ratio = (math.cosh(beta) + math.cos(beta)) / (math.sinh(beta) + math.sin(beta))
    return math.cosh(beta * xi) - math.cos(beta * xi) - ratio * (math.sinh(beta * xi) - math.sin(beta * xi))


def mac(a, b):
    """Unweighted sampled modal assurance criterion, invariant to sign/scale."""
    denominator = sum(v * v for v in a) * sum(v * v for v in b)
    if len(a) != len(b) or denominator <= 0:
        raise ValueError("MAC needs equal, nonzero vectors")
    return sum(x * y for x, y in zip(a, b)) ** 2 / denominator
