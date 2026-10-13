"""Independent integration of -k T'' = qdot; no FE fields or solver I/O."""


def analytical(config, x):
    g, b = config["geometry"], config["boundary"]
    L, k, s = (
        g["length"],
        config["material"]["conductivity"],
        config["internal_heat_generation"],
    )
    hot, cold = b["x_min_temperature"], b["x_max_temperature"]
    return {
        "temperature_K": hot + (cold - hot) * x / L + s * x * (L - x) / (2 * k),
        "gradient_K_m": [(cold - hot) / L + s * (L - 2 * x) / (2 * k), 0.0, 0.0],
        "heat_flux_W_m2": [k * (hot - cold) / L + s * (x - L / 2), 0.0, 0.0],
        "heat_flow_W": (k * (hot - cold) / L + s * (x - L / 2))
        * g["width"]
        * g["height"],
    }


def balance(config):
    g, b = config["geometry"], config["boundary"]
    L, A = g["length"], g["width"] * g["height"]
    k, s = config["material"]["conductivity"], config["internal_heat_generation"]
    base = k * (b["x_min_temperature"] - b["x_max_temperature"]) / L
    return {
        "outward_heat_flow_W": {
            "x_min": -A * (base - s * L / 2),
            "x_max": A * (base + s * L / 2),
        },
        "generation_W": s * L * A,
        "net_outward_W": s * L * A,
    }


def interpolation_norms(config, h):
    """Closed-form Q1 interpolant error, an independent quadrature cross-check.

    Not a solver result or a substitute for measured nonnodal norms.
    Integral_0^h [s t(h-t)/(2k)]² dt = (s/2k)² h^5/30.
    Integral_0^h [s(t-h/2)/k]² dt = (s/k)² h^3/12.
    """
    import math

    g = config["geometry"]
    V = g["length"] * g["width"] * g["height"]
    s, k = config["internal_heat_generation"], config["material"]["conductivity"]
    return {
        "temperature_L2_K_m3_2": s / (2 * k) * h * h * math.sqrt(V / 30),
        "gradient_H1_K_m1_2": s / k * h * math.sqrt(V / 12),
        "raw_flux_L2_W_m_minus1_2": s * h * math.sqrt(V / 12),
    }


def definition(config):
    L = config["geometry"]["length"]
    return {
        "reference_version": "1.0.0",
        "kind": "independent_analytical_reference",
        "equation": "-k Laplacian(T) = qdot",
        "temperature": "Th + (Tc-Th)*x/L + qdot*x*(L-x)/(2*k)",
        "gradient_x": "(Tc-Th)/L + qdot*(L-2*x)/(2*k)",
        "flux_x": "k*(Th-Tc)/L + qdot*(x-L/2)",
        "parameters": config,
        "midpoint": analytical(config, L / 2),
        "balance": balance(config),
        "assumptions": [
            "homogeneous positive isotropic k",
            "uniform volumetric source",
            "steady conduction",
            "opposing fixed temperatures",
            "adiabatic transverse faces",
            "affine rectangular slab",
        ],
        "interpolant_norm_identities": {
            "L2_T": "(qdot/(2*k))*h_x^2*sqrt(V/30)",
            "H1_T": "(qdot/k)*h_x*sqrt(V/12)",
            "L2_raw_flux": "qdot*h_x*sqrt(V/12)",
        },
    }
