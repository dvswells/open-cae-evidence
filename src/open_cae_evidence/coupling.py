"""Physical Q8 area-average transfer; no solver keyword syntax."""

import math


def active_tip(config):
    return config.get("tip_mass", {}).get("mass", 0) > 0


def area_weights(nodes, config, mesh):
    """Exact area integral of end-face Q8 interpolation divided by total area.

    Each face corner integrates to -1/12 and each midside node to 1/3 of
    its cell area. Signed weights are consistent interpolation integrals;
    the added mass matrix m*w*w^T remains positive semidefinite.
    """
    _, ny, nz = mesh
    g = config["geometry"]
    if not active_tip(config):
        return {}
    face = {
        n: xyz
        for n, xyz in nodes.items()
        if math.isclose(xyz[0], g["length"], abs_tol=1e-12) and n != max(nodes)
    }
    # The auxiliary mass node has the largest ID; mesh face nodes precede it.
    grid = {
        (
            round((y / g["width"] + 0.5) * 2 * ny),
            round((z / g["height"] + 0.5) * 2 * nz),
        ): n
        for n, (_, y, z) in face.items()
    }
    weights = {n: 0.0 for n in face}
    for j in range(ny):
        for k in range(nz):
            for a, b, factor in (
                (0, 0, -1 / 12),
                (2, 0, -1 / 12),
                (2, 2, -1 / 12),
                (0, 2, -1 / 12),
                (1, 0, 1 / 3),
                (2, 1, 1 / 3),
                (1, 2, 1 / 3),
                (0, 1, 1 / 3),
            ):
                try:
                    n = grid[(2 * j + a, 2 * k + b)]
                except KeyError as exc:
                    raise ValueError("Incomplete Q8 end face") from exc
                weights[n] += factor / (ny * nz)
    return weights
