# Source-driven thermal verification

This benchmark extends the source-free patch test without changing it. Both
adapters use `benchmarks/slab_volumetric_source.json`, the common `ProblemSpec`
and `ThermalVerifier`. No analytical temperatures or fluxes are imposed inside
the volume. The initial development results are not a substitute for the
commit-specific Actions artifacts, which contain the actual measured values.

## Physical reference and sign conventions

The homogeneous slab occupies `[0,0.1] × [0,0.2] × [0,0.3]` m. Conductivity is
10 W/(m K), uniform volumetric generation is 100000 W/m³, opposing x faces are
400 and 300 K, and the other four faces are adiabatic. Input specification 1.1
requires explicit `volumetric_heat_generation: W/m^3` and `source_region:
all_volume`. This deliberately restricted initial case requires finite positive
generation (at most 10⁷ W/m³); zero belongs to the unchanged specification 1.0
source-free case. Unknown units, versions and region assignments fail.

Integrating `-k T'' = s` gives `T' = C - s x/k`. Applying the two temperatures
gives

```
T(x) = Th + (Tc-Th)x/L + s x(L-x)/(2k)
grad_x(T) = (Tc-Th)/L + s(L-2x)/(2k)
q_x(x) = -k grad_x(T) = k(Th-Tc)/L + s(x-L/2).
```

The transverse derivatives vanish. Thus the governing equation, both Dirichlet
conditions and adiabatic sides hold. Integration of `div(q)=s` yields
`∮ q·n dA = ∫ s dV`. Here T(L/2)=362.5 K, q_x(0)=5000 and q_x(L)=15000 W/m².
Outward flow at the hot face is **−300 W**, at the cold face **+900 W**, and
generation/net outward flow is **+600 W**. Positive nodal heat input is opposite
to outward flow. These values are reference predictions, never replacement
solver output. Independent finite-difference tests check the PDE and Fourier
identity separately from the finite-element implementation.

## Native loading, input audits and reactions

CalculiX uses full-integration DC3D8 and `*HEAT TRANSFER, STEADY STATE` with
`*DFLUX` on the full SLAB element set, label **BF**. BF is power per unit volume,
as specified in the official [CalculiX 2.21 manual](https://www.dhondt.de/ccx_2.21.pdf),
section 7.39. The independent deck audit checks actual value, sign, unit comment,
set coverage, connectivity, material, temperatures and unsupported load cards.
It rejects replacing the source with prescribed internal temperatures.

Elmer uses affine native 808 HEX8 elements and its HeatSolver with body force
`Volumetric Heat Source = Real 100000`. This is the explicit volumetric branch
of the pinned upstream `fem/src/HeatSolve.F90`, rather than its density-scaled
`Heat Source` fallback. The independently parsed SIF audit checks the actual
Body/Body Force assignment and scalar, unit metadata, native mesh orientation,
boundary IDs, equations and material. Native mechanisms assemble the load in
both solvers; neither adapter evaluates the reference to form numerical DOFs.

For this CalculiX source problem the requested native RFL is **K T**, containing
the body-load contribution as well as boundary residuals. It is not the separated
Dirichlet reaction. For each affine constant-source cell the consistent nodal
load is independently integrated as `f_i = ∫ N_i s dV = s Vcell/8`. Subtracting
these audited actual loads from native RFL gives `r = K T - f`. Free-node
residuals must vanish. The original combined nodal rates are retained as
`native_total_nodal_heat_rate`; corrected `nodal_heat_rate` is explicitly marked
`postprocessed_from_solver`. Elmer's native **Temperature Loads** already provide
this residual, and remain marked as solver output. Neither uses predicted
analytical boundary reactions to fill missing data.

Source integration uses actual cell volumes and independently audited serialized
source coverage. Reaction conservation is checked against that actual integral,
then the source and reactions are checked separately against the requested
physical specification. Consequently a wrong-source solve can conserve energy
and still fail the physical gates.

## Nonnodal errors and quadrature

On this uniform mesh the exact quadratic solution can be reproduced at every
node by Q1 FE. Zero nodal error is therefore insufficient. `q1.py` evaluates the
trilinear field and its physical gradient from **actual exported numerical DOFs**.
Whole-volume error norms are

```
E_T = sqrt(∫ (T_h-T_exact)^2 dV)                 [K m^(3/2)]
E_g = sqrt(∫ |grad(T_h)-grad(T_exact)|² dV)      [K m^(1/2)]
E_q = sqrt(∫ |-k grad(T_h)-q_exact|² dV)        [W m^(-1/2)].
```

Each affine cell uses 3×3×3 Gauss quadrature, exact through degree five in each
coordinate. It integrates the squared quadratic temperature error and gradient
error, and includes off-center points. Center-only gradient sampling would give
a false zero for the nodally exact quadratic case. An adversarial review explicitly
demonstrates that trap with actual solver DOFs. Per-cell maxima are maxima over
these 27 physical samples, not certified continuum suprema. The CSV retains every
sample's coordinates, volume weight, FE/exact temperatures, gradients and fluxes.
Nodal error, volume RMS error and full dimensioned norms are all retained.

As an independent integration check, exact quadratic nodal interpolation has
`E_T = (s/2k) h_x² sqrt(V/30)`, `E_g = (s/k) h_x sqrt(V/12)` and
`E_q = s h_x sqrt(V/12)`. These identities follow from direct polynomial
integration on an interval; they do not compute or replace reported errors.
Measured DOF deviations bound the difference from those identities by triangle
inequalities. Four production meshes have counts 2³, 4³, 8³ and 16³; h_x halves
from 0.05 to 0.00625 m while all physical dimensions stay fixed.

Measured successive order is `log(E_coarse/E_fine)/log(h_coarse/h_fine)`.
Temperature, gradient and raw Fourier flux must decrease and support the stated
Q1 expectations 2, 1 and 1 within 0.1. This is a case-specific scientific gate,
not a general prediction for arbitrary geometry or material. An error below ten
times the propagated native output-resolution floor has a null order. Recovery
and boundary-reconstruction orders are diagnostics and are not fitted to the
primary expectations. No meaningful order is claimed for the unchanged exact
linear patch test.

## Distinct flux semantics and boundary estimates

The normalized archive retains five distinct representations:

| Quantity | Source/meaning |
|---|---|
| `raw_temperature_gradient` | Q1 derivative of actual exported T, sampled at 2³ Gauss locations |
| `raw_fourier_heat_flux` | −k times that gradient, common cross-solver field |
| CalculiX `heat_flux` | Native DAT integration-point HFL |
| Elmer `heat_flux` | Q1 sample of native L2-recovered nodal flux; explicitly postprocessed |
| `recovered_nodal_heat_flux` | CalculiX FRD extrapolated/averaged flux or Elmer native L2 nodal recovery |
| `boundary_reaction_heat_flow` | Outward sum of actual/normalized boundary residuals |
| `heat_flow` | Independently reconstructed outward raw Q1 boundary flux estimate |

Native CalculiX HFL is checked against the Q1 Fourier field with propagated
printing uncertainty. FRD recovered values are checked against the native DAT
extrapolation/averaging. Elmer recovery is checked with the actual global L2
projection residual; it is not forced to equal the discontinuous raw gradient
pointwise. Volume norms of both native-represented and recovered fields are
reported separately. Recovery rates around 1.5 in this case differ between
solvers. An O(h) boundary error over an O(h) layer can explain such an L2 rate;
this case does not establish a universal recovery order or solver ranking.

One-sided Q1 boundary gradients are discretization estimates. Their individual
heat-flow errors are approximately `s h_x A/2` (150, 75, 37.5, 18.75 W), whereas
boundary reactions satisfy the balance at native output accuracy. The raw net
outward estimate approaches 600 W from about 300, 450, 525 and 562.5 W. Demanding
roundoff agreement between these representations would be scientifically wrong.
Both estimates and their measured errors/orders remain visible in evidence.

## Declared accuracy budgets and output precision

Source-case budgets were fixed before final comparative runs: nodal T error
≤0.001 K; reaction/reference and reaction/source balance relative error ≤2×10⁻⁶;
native flux consistency excess after print propagation ≤2×10⁻⁶; integrated source
relative error ≤10⁻¹². Reaction residuals scale by the expected 600 W generation;
native flux consistency scales by the reference peak flux. Recovery is checked
against its own reconstruction/projection equations, not analytical equality.
These budgets allow native decimal output and assembly/linear-solve roundoff,
while remaining far smaller than measured discretization errors. The original
source-free tolerances and all modal physics/gates remain unchanged.

CalculiX DAT NT prints E13.6 (seven significant digits); near 300–400 K its step
is 10⁻⁴ K. Elmer ASCII VTU prints eight significant digits (10⁻⁵ K there).
CalculiX resolution is parsed from the actual native tokens; Elmer precision
metadata records the pinned native VTU format at this benchmark's field scale.
Different field scales/formats would need a renewed precision audit. Differentiation amplifies
temperature half-ULP δ: componentwise gradient uncertainty ≤2δ/h_component,
Fourier uncertainty ≤2kδ/h_component. Early CalculiX 8³ pointwise comparisons
exceeded a naive fixed consistency threshold because DAT rounded DOFs while HFL
was computed from internal values. The final check removes only this explicitly
propagated printing bound before applying the original 2×10⁻⁶ consistency budget.
It does not weaken source-free gates or conceal discretization error. Small
differences between cloud CalculiX 2.20 and CI 2.21, and between Elmer builds, may
reflect assembly, recovery and printing; use artifact values for each exact run.

Cross-solver comparison first replays both native archives and repeats independent
physical gates. It matches coordinates, not node IDs, compares common Q1 fields
and norm triangle bounds, and checks reactions with a 4×10⁻⁶ combined relative
budget and T with 0.002 K. Differently recovered flux norms are reported without
averaging or claiming identical numerical semantics. Agreement is corroboration.

## Reproduction, adversarial gates and evidence

The Ubuntu Actions workflow builds pinned official Elmer headlessly, installs
signed CalculiX 2.21, and executes 48 original +26 Elmer +24 source tests. The
source tests retain real native corruption and wrong-source solve evidence.
They cover missing/wrong/unit/sign/region sources, numerical T corruption,
reversed flux, coordinate transforms, incomplete native output, fabricated yet
balanced reactions, altered norms and rehashed archive tampering including CSV
samples and plot data. Missing solvers fail tests rather than skipping them.

The owner need only use GitHub on Android. Contributors can independently use
the workflow or the following commands after its documented runtime setup:

```sh
python -m pip install '.[plots]'
python -m unittest discover -s tests_source -v
python -m open_cae_evidence.benchmark --config benchmarks/slab_volumetric_source.json --output evidence/source-calculix
python -m open_cae_evidence.benchmark --adapter elmer --config benchmarks/slab_volumetric_source.json --output evidence/source-elmer
python -m open_cae_evidence.source_comparison --calculix evidence/source-calculix --elmer evidence/source-elmer --output evidence/source-comparison
python -m open_cae_evidence.source_comparison --calculix evidence/source-calculix --elmer evidence/source-elmer --output evidence/source-comparison --replay
```

Native replay needs Python and pinned plot dependencies, but **no solver execution
or AI**. Use `benchmark --verify-evidence` (with `--adapter elmer` for Elmer).
`numerical_review` recomputes the nodal/center trap and inspects the converged
wrong-source rejection before PR creation; CI repeats it and its replay.
All output directories must be new. Reports include machine-precision measured
tables, semantic provenance, tolerances and deterministic log-log SVG plots.
`analytical-reference.json` and per-mesh `nonnodal-samples.csv` supplement
the existing native files, contracts, audits, verification JSON and manifests.
Replay recomputes their semantics, not just hashes. SHA-256 integrity alone is
not independent scientific correctness. Artifacts expire after 30 days.

Results contract **1.0.0 remains unchanged**; existing arbitrary named optional
quantities carry units/locations/provenance. The input source is audited metadata,
not a fabricated solver field. Archive index version 3 and legacy source-free
and modal replay are retained. Downstream readers must inspect quantity source
and recovery semantics rather than interpreting an integration-point sample
location as native provenance. No MCP service is required or implemented.

## Limits

This is an ideal homogeneous, uniform-source, affine structured 3D Q1 problem
whose exact physics is one-dimensional. It tests true volume discretization error
but not arbitrary distorted meshes, spatially varying conductivity/source,
nonlinear/transient/coupled physics, radiation, contact or real material data.
Exported decimal DOFs limit very fine-mesh orders. CalculiX's separated reaction
is a defensible residual reconstruction, not a claimed native separate reaction.
Recovered flux is solver-specific. Energy conservation alone cannot establish
correct physical inputs. No experimental validation, certification, industrial
qualification or general ANSYS equivalence follows from this benchmark.
