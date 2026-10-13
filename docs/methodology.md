# Cantilever modal benchmark specification

## Physical contract

The beam axis is x ∈ [0,1] m, width y ∈ [-0.0075,0.0075] m and
height z ∈ [-0.005,0.005] m. E = 210 GPa, nu = 0.3,
rho = 7800 kg/m³. These are benchmark definitions, not a certified steel grade.
Volume = 1.5e-4 m³; total geometric mass = 1.17 kg. No loads or prestress.
The complete x=0 face (including midside nodes) is fixed in translations 1,2,3.
Solids have no independent rotational degrees of freedom. All remaining faces
are free. This ideal clamp is not a model of a real fixture.

CalculiX is unitless: its input numbers must form a consistent m–kg–s system.
The input JSON explicitly declares length, mass, time, stress and density units.
No implicit conversion is performed. Nonfinite inputs, nonpositive dimensions,
unstable isotropic Poisson ratios, missing DOFs and wrong unit labels are rejected.
Benchmark-specific bounds (dimensions 0.001–10 m, E 1e9–1e12 Pa,
rho 1000–30000 kg/m³, L/max(b,h) >= 50) detect common SI mistakes and define
the applicable domain. They do not validate every real material or catch a false
unit label attached to plausible numbers. Input provenance remains essential.

## FEM formulation and mesh

Structured C3D20R elements have 20 nodes and reduced integration (eight integration
points). The standard CalculiX linear elastic solid formulation and solver mass
matrix are used. The generated meshes contain 12×2×2 (48), 24×4×4 (384) and
48×6×6 (1728) elements. No mesher, proprietary reference or fitted FEM answer is
required. A separate serialized-deck reader verifies each node, element corner
and midside position, nonoverlapping cells, geometric volume, material, keyword
contract, all fixed nodes and all displacement outputs before launching CalculiX.
Positive frequencies plus analytical shapes provide checks against spurious low
modes; the workflow does not claim a universal proof against hourglassing.

Six eigenpairs are requested with `*FREQUENCY`. `.dat` eigenvalues must be positive
and real, ordered with explicit mode numbers, and consistent with omega=2*pi*f and
lambda=omega² within 2e-6 relative tolerance (seven significant digits printed by
the solver). Every mode must include all nodes, and fixed displacements must be
zero within 1e-8 of peak displacement. `.frd` and the full `.dat` are retained.
CalculiX's eigenvector normalization and sign are arbitrary: CSV displacement
values are not physical vibration amplitudes.

## Independent analytical reference

For uniform Euler–Bernoulli bending, E I phi'''' = rho A omega² phi. The clamp
requires phi(0)=phi'(0)=0; the free tip requires phi''(L)=phi'''(L)=0.
Solving the four-coefficient characteristic system gives cos(beta) cosh(beta)=-1.
For each positive root beta:

```
f = beta²/(2*pi*L²) * sqrt(E*I/(rho*A))
phi(xi) = cosh(beta*xi) - cos(beta*xi)
          - gamma * (sinh(beta*xi) - sin(beta*xi))
gamma = (cosh(beta) + cos(beta)) / (sinh(beta) + sin(beta))
xi = x/L
```

A = b h. Displacement in z bends about y (I_y=b h³/12); displacement in y
bends about z (I_z=h b³/12). Both families are sorted together. For this geometry
the first six are z1, y1, z2, y2, z3, y3. Higher torsional and axial modes are not
analytically verified. If a changed configuration brings such a mode into the
first six, the bending direction/shape checks should fail; it needs a new benchmark
specification rather than relaxed tolerances.

The reference roots are recorded to double precision in `reference.py` and tested
against cos(beta)+sech(beta)=0. Analytical free-tip conditions and the known 1.5
plane frequency ratio are tested separately. The derivation above is sufficient
to reproduce the reference without third-party data.

The shape comparison uses unweighted centerline node samples:
MAC = (u·v)² / ((u·u)(v·v)). It is invariant to sign/scale and is not a mass-weighted
orthogonality or full-field residual check. Nodal squared displacement fractions
identify the expected bending direction; they are not physical modal energies.

## Acceptance criteria and limits

| Check | Required limit | Basis |
|---|---|---|
| Frequency vs EB, every mesh/mode | ≤1% relative | Slender beam (L/h=100, L/b≈66.7); room for shear, rotary inertia and local clamp effects omitted by EB |
| Analytical centerline MAC | ≥0.99 | Detect wrong shapes/orders while permitting small 3D clamp and cross-section effects |
| Expected direction fraction | ≥0.95 | Well-separated two-plane bending; significant mixing indicates wrong classification/model |
| Clamp nodal displacement / peak | ≤1e-8 | Constraints should be exact; allow negligible numerical/output roundoff |
| Final two-mesh frequency change | ≤0.5% | Chosen discretization budget below the analytical comparison budget |
| Change decreases with refinement | final ≤max(first,5e-5) | Reject divergent trends while allowing a 0.005% numerical noise floor |
| Real E×4 / rho×2 scaling | ratios 2 / 1/sqrt(2), absolute error ≤2e-5 | Exact eigenvalue scaling for unchanged geometry and nu; accommodates printed precision |

These are acceptance criteria for this case, not certified uncertainty bounds.
Three anisotropic refinement levels do not establish an asymptotic order or a
Richardson-extrapolated exact solution. No experiment, damping, geometric
nonlinearity, prestress, contact, real fixture compliance or manufacturing variation is included.

## Negative checks and evidence replay

Tests reject negative/zero E, negative density, suspect unit-scaled material
values, wrong units, invalid nu, NaN/Inf, non-slender geometry, bad refinement,
corrupted connectivity and incomplete/misplaced clamps. Additional tests alter
**actual solver output** to demonstrate detection of zero frequencies, truncated
displacement output and wrong mode numbers. Verification rejects a 10% frequency
change, switched first/second bending shapes and a nonzero fixed-node displacement.
No mocked output replaces the actual three-mesh benchmark.

Two further real solves check f proportional to sqrt(E) and to 1/sqrt(rho).
All integration tests fail when the solver is missing. CI captures the number of
executed tests and preserves the process exit status when piping logs.

To replay a downloaded `benchmark` directory without rerunning the solver:

```sh
python -m open_cae_evidence.benchmark --verify-evidence path/to/benchmark
```

This verifies hashes, regenerates inputs, audits the archived decks and recomputes
analytical, shape, boundary and convergence checks from `.dat`. SHA-256 manifests
detect accidental alteration, not identity or malicious replacement of the entire
manifest. Rerunning CI is the independent solver execution. The run's source SHA,
solver binary hash, solver version, Python/platform and library versions identify
its provenance. Ubuntu's solver is pinned at 2.21-1; transitive libraries and
GitHub's ubuntu-24.04 image are recorded rather than claimed immutable. Minor
platform-dependent roundoff is expected, not byte-for-byte deterministic outputs.

## Open documentation

- [CalculiX project](https://www.calculix.de/) (GPL solver).
- [CalculiX 2.21 manual](https://www.dhondt.de/ccx_2.21.pdf): C3D20R, `*FREQUENCY`,
  `*BOUNDARY`, `*NODE PRINT`, `*NODE FILE` and eigenvalue output.
- The equations, root definitions, configuration and tests in this repository
  are original open reference definitions under MIT; no restricted dataset is used.
