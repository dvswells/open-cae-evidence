"""Human-readable modal evidence report (legacy field compatibility)."""

import json
from .reference import mass_ratio


def report(result):
    tip = "tip_mass" in result["config"]
    title = (
        "Cantilever with concentrated tip mass"
        if tip
        else "Cantilever modal verification"
    )
    equation = (
        "sech(beta)+cos(beta)+mu*beta*(cos(beta)*tanh(beta)-sin(beta))=0."
        if tip
        else "cos(beta)*cosh(beta)=-1."
    )
    lines = [
        f"# {title}",
        "",
        f"Status: **{result['status']}**",
        "",
        "This is verification of one linear elastic benchmark, not validation against experiments",
        "or evidence of general replacement suitability for a commercial CAE platform.",
        "",
        "## Model and assumptions",
        "",
        "SI: m, kg, s, Pa. Uniform isotropic steel; small undamped free vibrations.",
        "Solid C3D20R quadratic reduced-integration hexahedra; consistent solver mass matrix.",
        "Every x=0 face node is fixed in DOFs 1–3; other faces have no prescribed displacement.",
        "The deck is audited independently after serialization, including material, volume, connectivity and clamp.",
        "",
        "```json",
        json.dumps(result["config"], indent=2),
        "```",
        "",
        "## Independent checks",
        "",
        "Euler–Bernoulli: f = beta²/(2*pi*L²) * sqrt(E*I/(rho*A)); " + equation,
        "Both transverse planes are considered and sorted. Shape checks use the analytical cantilever",
        "eigenfunction at centerline nodes and sign/scale-invariant sampled MAC.",
        "Plane classification uses the fraction of nodal displacement squared in the expected direction;",
        "this is a diagnostic, not a physical modal energy calculation.",
        "",
        "Acceptance limits are specified before execution: 1% analytical frequency error, MAC >= 0.99,",
        "direction fraction >= 0.95, clamp residual <= 1e-8 of peak displacement, final mesh change <= 0.5%.",
        "The final change must also not exceed the first change (5e-5 relative noise floor).",
        "The 1% bound allows 3D end constraint effects, shear deformation and rotary inertia absent from",
        "slender-beam theory; it is a benchmark acceptance criterion, not a universal error bound.",
        "MAC >= 0.99 permits small cross-section and clamp effects; 0.5% is the selected discretization budget.",
        "Exact fixed displacements should be zero; 1e-8 allows negligible output roundoff.",
        "",
        "## Calculated frequencies and shapes",
        "",
        "| Mode | Plane/order | FEM Hz | EB Hz | Error % | MAC | Pass |",
        "|---|---|---:|---:|---:|---:|---|",
    ]
    if result["runs"]:
        for mode in result["runs"][-1]["modes"]:
            lines.append(
                f"| {mode['mode']} | {mode['direction']}/{mode['order']} | {mode['fem_frequency_hz']:.6f} | {mode['frequency_hz']:.6f} | {100 * mode['analytical_relative_error']:.4f} | {mode['centerline_mac']:.8f} | {mode['passed']} |"
            )
    lines += [
        "",
        "## Mesh convergence",
        "",
        "| Mode | Frequencies coarse → fine (Hz) | Final change % | Pass |",
        "|---|---|---:|---|",
    ]
    for item in result.get("convergence", []):
        lines.append(
            f"| {item['mode']} | {' → '.join(f'{v:.6f}' for v in item['frequencies_hz'])} | {100 * item['relative_changes'][-1]:.5f} | {item['passed']} |"
        )
    lines += [
        "",
        "Refinement increases all three element counts. Mode direction and shape checks guard",
        "against matching unrelated eigenmodes by index. No Richardson extrapolation or asymptotic",
        "convergence order is claimed from these nonuniform refinement ratios.",
        "",
        "## Evidence and limits",
        "",
        "Each mesh retains the .inp, .dat, .frd, complete solver log, input audit, verification JSON",
        "and six full nodal displacement CSV files. Eigenvectors have arbitrary sign and amplitude;",
        "CSV displacements are solver normalization values, not physical vibration amplitudes.",
        "manifest.json records SHA-256 digests. provenance.json records solver and Python environment.",
        "Hashes detect alteration, not authorship. Rounded .dat values limit replay precision.",
        "Only six bending modes of this geometry are checked: no torsional/axial reference verification,",
        "experimental validation, damping, prestress, nonlinear response or industrial qualification.",
        "Positive plausible but falsely labeled units cannot be inferred from numbers alone.",
        "See docs/methodology.md for reference derivation, scope and reproducibility instructions.",
    ]
    if tip:
        lines += [
            "",
            "## Tip inertia and coupling",
            "",
            f"Tip/beam mass ratio mu = {mass_ratio(result['config']):.8g}.",
            "A centered auxiliary MASS node remains independent in three translation DOFs.",
            "Three MPCs equate its displacement to the exact Q8 area-average of the end face.",
            "Signed interpolation weights integrate to one and have zero first moments.",
            "They distribute inertia forces by virtual work without tying every face node rigidly.",
            "Rotary inertia, eccentricity and attachment stiffness are deliberately neglected.",
            "Each mode must satisfy the coupling within 2e-6 of peak displacement (DAT rounding).",
            "CalculiX may warn that its internal MASS material has no elastic constants; this",
            "inertia-only element is intentionally stiffness-free. All *ERROR messages are fatal.",
            "See docs/tip-mass.md for the dynamic tip boundary condition and mass-scaling checks.",
        ]
        lines += [
            "",
            "| Mesh | Total mass kg | Coupled face nodes | Maximum MPC residual |",
            "|---|---:|---:|---:|",
        ]
        for run in result["runs"]:
            if "tip_mass_kg" in run["audit"]:
                lines.append(
                    f"| {'x'.join(map(str, run['mesh']))} | {run['audit']['total_mass_kg']:.8g} | {run['audit']['coupled_face_nodes']} | {max(m['tip_coupling_relative_residual'] for m in run['modes']):.3e} |"
                )
        lines += [
            "",
            "Required CI integration tests separately exercise tip-only and total-mass scaling,",
            "modulus scaling, the zero-mass limit, and equivalent solid-face elimination pivots.",
            "Their real inputs, outputs and six-mode ratio JSON files are in evidence/tests/tip-mass-*/.",
            "This standalone report does not assert that those separate tests have executed.",
        ]
    lines += [
        "",
        "## Semantic FRD verification",
        "",
        "For every mesh, ASCII FRD geometry, HEX20 connectivity, solver version, six DISP datasets,",
        "mode IDs, frequencies, eigenvalues and all signed displacement components are checked",
        "against the audited deck and DAT. The component comparison uses 1e-5 relative tolerance",
        "and 1e-12 of peak displacement absolute tolerance for FRD/DAT printed precision.",
        "MASS elements and their auxiliary nodes are omitted by the FRD exporter; their full",
        "mode displacements are preserved in DAT/CSV and independently checked by the MPC residual.",
    ]
    if result.get("error"):
        lines += ["", "## Execution failure", "", result["error"]]
    return "\n".join(lines) + "\n"
