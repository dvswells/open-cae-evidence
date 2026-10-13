"""Technical report for measured sourced-slab errors; no computed solver values."""

import json


def report(result):
    lines = [
        "# Source-driven 3D thermal slab",
        "",
        f"Status: **{result['status']}**; solver: {result.get('solver_name')}",
        "",
        "-k Laplacian(T)=qdot; T=Th+(Tc-Th)x/L+qdot*x*(L-x)/(2k).",
        "Q1 nodal temperatures can be exact while the piecewise field is not quadratic. All primary error norms",
        "are measured from actual DOFs at 3x3x3 volume Gauss points, including off-center points.",
        "L2(T): K m^(3/2); gradient H1 seminorm: K m^(1/2); flux L2: W m^(-1/2).",
        "Native flux and recovered flux are separate diagnostics, never the common unrecovered Fourier field.",
        "",
        "| Mesh | h_x m | Nodal T error K | Volume T L2 | T RMS K | Gradient H1 | Raw flux L2 | Recovered flux L2 |",
        "|---|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for r in result["runs"]:
        m = r["thermal"]
        lines.append(
            f"| {'x'.join(map(str,r['mesh']))} | {m['h_x_m']:.6g} | {m['temperature_max_error_K']:.6g} | {m['temperature_L2_K_m3_2']:.9g} | {m['temperature_RMS_K']:.9g} | {m['gradient_H1_K_m1_2']:.9g} | {m['raw_flux_L2_W_m_minus1_2']:.9g} | {m['recovered_flux_L2_W_m_minus1_2']:.9g} |"
        )
    lines += [
        "",
        "## Boundary/source balance",
        "",
        "Positive reaction input supplies heat to the body; outward boundary flow has the opposite sign.",
        "Source integral uses the audited serialized source and actual cell volumes. CalculiX native RFL includes",
        "body loading: subtract integral N_i*qdot dV to obtain the boundary residual. Elmer Temperature Loads",
        "already include that subtraction. No analytical boundary flow is substituted for numerical output.",
        "Raw one-sided Q1 boundary gradients are discretization estimates, not boundary reactions. Their",
        "net outward integral need not equal the source on a coarse mesh; its error must decrease.",
        "",
        "| Mesh | Source W | Reaction hot outward W | Reaction cold outward W | Reaction net W | Raw hot outward W | Raw cold outward W | Relative balance residual |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for r in result["runs"]:
        m = r["thermal"]
        a = m["boundary_reaction_outward_heat_flow_W"]
        b = m["surface_outward_heat_flow_W"]
        lines.append(
            f"| {'x'.join(map(str,r['mesh']))} | {m['actual_source_integral_W']:.9g} | {a['x_min']:.9g} | {a['x_max']:.9g} | {m['reaction_net_outward_W']:.9g} | {b['x_min']:.9g} | {b['x_max']:.9g} | {m['normalized_energy_imbalance']:.5g} |"
        )
    lines += [
        "",
        "## Measured convergence and declared budgets",
        "",
        "![Measured convergence](convergence.svg)",
        "",
        "```json",
        json.dumps(result.get("mesh_consistency", {}), indent=2),
        "```",
        "",
        "```json",
        json.dumps(result["tolerances"], indent=2),
        "```",
        "",
        "Rates below the recorded output-resolution floor are null, not invented. Recovered-flux orders are",
        "reported separately and are not forced to equal the raw-gradient order. Source-free gates are unchanged.",
        "Independent closed-form interpolant norm identities check the quadrature; they do not replace measured errors.",
        "Native files, full DOFs, source audits, normalized raw/native/recovered fields, per-mesh checks and manifests",
        "are archived and re-extracted on replay without solver execution. See docs/source-thermal.md.",
        "This ideal uniform-source affine slab is numerical verification, not experimental validation or industrial qualification.",
        "",
    ]
    if result.get("error"):
        lines += ["## Failure", "", result["error"]]
    return "\n".join(lines)
