"""Thermal technical report built from numerical checks, not solver completion."""

import json


def report(result):
    elmer = result.get("solver_name") == "Elmer"
    lines = [
        "# Homogeneous slab: steady 3D heat conduction",
        "",
        f"Status: **{result['status']}**",
        "",
        "This is numerical verification of an ideal conduction problem, not experimental validation or industrial qualification.",
        "",
        "## Physical model",
        "",
        (
            "Elmer 808 (eight-node trilinear thermal brick), constant isotropic conductivity,"
            if elmer
            else "DC3D8 (eight-node, fully integrated thermal brick), constant isotropic conductivity,"
        ),
        "prescribed opposing-face temperatures, adiabatic remaining faces, zero internal generation.",
        "The steady procedure has no physical transient time.",
        "",
        "```json",
        json.dumps(result["config"], indent=2),
        "```",
        "",
        "## Independent reference and numerical checks",
        "",
        "From div(k grad T)=0: T(x)=T_hot+(T_cold-T_hot)*x/L; grad T=[(T_cold-T_hot)/L,0,0].",
        "Fourier: q=-k grad T; Q=q_x*width*height. The reference uses no FE interpolation or solver values.",
        (
            "Numerical Q1 derivatives use actual nodal temperatures. Native L2 nodal flux, sampled at Gauss points, is"
            if elmer
            else "Numerical Q1 derivatives use actual nodal temperatures. Actual solver integration-point HFL is"
        ),
        "compared to both these derivatives and the closed-form reference. Flux is integrated on all six",
        "exterior faces by Gauss-point interpolation/extrapolation to 2x2 face quadrature, exact for the",
        (
            "affine rectangular Q1 recovered flux field. Rates are independently summed from native Temperature Loads."
            if elmer
            else "affine rectangular Q1 gradient field. Rates are also independently summed from actual solver RFL."
        ),
        "Positive nodal load means heat supplied to the body; positive surface flow is outward. They have opposite signs.",
        "",
        "Acceptance: nodal temperature error <=1e-3 K; normalized gradient/Fourier closure <=1e-5;",
        "flux magnitude/direction, integrated flow, side insulation, reactions, conservation and mesh spread <=2e-6.",
        "These limits allow printed DAT rounding and roundoff; they are case-specific verification budgets.",
        "",
        "| Mesh | Max T error K | Mean q_x W/m² | Cold outward flow W | Hot input W | Cold input W | Energy imbalance / Q | Pass |",
        "|---|---:|---:|---:|---:|---:|---:|---|",
    ]
    for run in result["runs"]:
        t = run["thermal"]
        lines.append(
            f"| {'x'.join(map(str,run['mesh']))} | {t['temperature_max_error_K']:.5g} | {t['mean_heat_flux_W_m2'][0]:.9g} | {t['surface_outward_heat_flow_W']['x_max']:.9g} | {t['boundary_heat_input_W']['x_min']:.9g} | {t['boundary_heat_input_W']['x_max']:.9g} | {t['normalized_energy_imbalance']:.4g} | {t['passed']} |"
        )
    lines += [
        "",
        "## Mesh study",
        "",
        json.dumps(result.get("mesh_consistency", {}), indent=2),
        "",
        "This is mesh consistency, not demonstrated discretization convergence. The exact linear solution is",
        "represented in the affine Q1 space; no monotonic improvement or convergence order is required/claimed.",
        "",
        "## Evidence, extraction limits and scope",
        "",
    ]
    if elmer:
        lines += [
            "Native SIF, mesh.header/nodes/elements/boundary, ASCII VTU, solver log, SI mapping, normalized fields,",
            "audits, numerical verification, provenance and hashes are retained. Elmer FluxSolver projects",
            "-k grad(T) into a continuous nodal Q1 field through global L2 projection. Interpolation to physical",
            "Gauss points is explicitly postprocessed_from_solver; it is not native integration-point output.",
            "Temperature Loads are independently computed by Elmer from A_bulk*T-b_bulk. VTU units are implicit:",
            "the audited SIF and units.json enforce m/K/W throughout, identity Cartesian coordinate mapping.",
            "Density=1 is a required inert Elmer material field, unused by this steady source-free conduction problem.",
            "Only the heat linear residual is an iterative convergence criterion; flux projection uses direct UMFPACK.",
            "See docs/elmer.md for headless source build, recovery semantics, output precision and limitations.",
        ]
    else:
        lines += [
            "Raw .inp, .dat, .frd, logs, serialized material/BC audit, normalized fields, CSV, JSON verification,",
            "solver provenance and hashes are retained. normalized-results.json uses the common versioned contract.",
            "FRD temperature/RFL and mesh are checked against DAT/audited input. FRD nodal FLUX is extrapolated",
            "and averaged: the observed constant-field extrapolation bias is about 4e-5. Its diagnostic consistency",
            "allowance is 1e-4; all physical flux/energy gates use the raw DAT integration-point values at 2e-6.",
            "There is no native surface-flux output in this adapter. Integrated flow is explicitly labeled",
            "postprocessed_from_solver, not fabricated or presented as a native solver quantity.",
        ]
    lines += [
        "Only homogeneous isotropic steady conduction on affine rectangular HEX8 meshes is covered.",
        "No radiation, convection, internal sources, transient capacity, nonlinear properties or coupled physics",
        "is verified. JSON/hashes provide reproducibility and integrity, not authorship or malicious forgery protection.",
        "See docs/thermal.md and docs/architecture.md. CI runs the negative tests and both preserved modal benchmarks.",
    ]
    if result.get("error"):
        lines += ["", "## Failure", "", result["error"]]
    return "\n".join(lines) + "\n"
