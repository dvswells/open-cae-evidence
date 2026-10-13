# Cantilever with concentrated tip mass

The unloaded cantilever geometry, steel properties, fixed face and three C3D20R
meshes are unchanged. A **1.17 kg** ideal point mass is located at (L,0,0), giving tip/beam
mass ratio mu=1 and total model mass 2.34 kg. This is an explicit mathematical
idealization: no size, eccentricity, rotational inertia or attachment compliance.
It is not evidence that a finite industrial payload can neglect those effects.

## Translational mass and section transfer

CalculiX's native one-node `MASS` element carries three independent translation
DOFs at an auxiliary node, distinct from every solid mesh node, at the tip center.
Three `*EQUATION` constraints impose:

```
u_mass,a = sum_i w_i * u_i,a = (1/A) integral_tip u_a dA,   a = 1,2,3
```

Here the integral is the **exact area average of the Q8 end-face interpolation**,
assembled over every solid element face on x=L. Each face corner contributes
-A_cell/(12*A), each midside node +A_cell/(3*A); contributions at shared nodes
are accumulated. All face nodes participate. Signed weights are interpolation
integrals, not negative physical point masses. Tests independently check their
constant, linear, quadratic and mixed moments: sum(w)=1, sum(w*y)=sum(w*z)=0,
sum(w*y²)=b²/12, sum(w*z²)=h²/12 and sum(w*y*z)=0.

The added kinetic energy is T_tip = Mt/2 * |u_dot_mass|². Substitution gives the
positive semidefinite translational mass contribution Mt * w*w^T in each axis.
Virtual work transfers the inertia force as f_i = w_i*f_tip. The vanishing first
moments preserve a centered, rotation-free resultant. Only three average
translations are constrained: face nodes are not rigidly tied together and
section rotations and warping remain possible. This is an ideal massless
distributing attachment, not a modeled rigid payload body. In the beam limit,
the mean transverse tip displacement equals the beam centerline displacement.

The first term of each CalculiX equation is a solid face DOF (the largest absolute
face coefficient is chosen); the MASS node remains independent. Moving the mass
DOF to the dependent side caused `*ERROR in add_sm_ei` during an exploratory
CalculiX 2.20 run despite exit code zero. The input auditor rejects that ordering,
and every solver error is fatal. A real-solver test changes the dependent solid
node while retaining the identical constraint, then checks all six frequencies
within 2e-5 in ratio, shapes and coupling. The choice is algebraic elimination,
not an arbitrary physical attachment to that solid node.

CalculiX may emit a warning that its internally created MASS material has no
elastic constants. This inertia-only element intentionally contributes no
stiffness; the beam retains its audited elastic material. The warning is retained
in the log. No errors are permitted.

## Independent eigenvalue and shape reference

For the same uniform Euler–Bernoulli beam, E I phi'''' = rho A omega² phi.
At the clamp phi=phi'=0. At the mass end the bending moment is zero and the
shear balances the point-mass inertia:

```
phi''(L) = 0
EI*phi'''(L) + Mt*omega²*phi(L) = 0
mu = Mt/(rho*A*L)
beta^4 = rho*A*omega²*L^4/(EI)
```

The general clamped-end solution is C*(cosh(beta*xi)-cos(beta*xi)) +
D*(sinh(beta*xi)-sin(beta*xi)), xi=x/L. Zero moment implies D/C=-gamma,
gamma=(cosh(beta)+cos(beta))/(sinh(beta)+sin(beta)). Substituting the dynamic
shear condition gives:

```
1 + cos(beta)*cosh(beta)
  + mu*beta*(cos(beta)*sinh(beta) - sin(beta)*cosh(beta)) = 0
```

The implemented, algebraically equivalent equation divides by cosh(beta):

```
sech(beta) + cos(beta)
  + mu*beta*(cos(beta)*tanh(beta) - sin(beta)) = 0
```

Its first six positive roots are independently bracketed on a pi/64 grid and
bisected 60 times. The supported ratio is 0≤mu≤10; the scan must find exactly
six roots. No FEM values are used to solve this scalar equation. For mu=0 the
unloaded cantilever roots and input deck are recovered exactly. Root residuals and the
unscaled moment/shear boundary conditions are tested independently.

For each plane f=beta²/(2*pi*L²)*sqrt(E*I/(rho*A)); the two families are sorted.
For mu=1 the first six are z1, y1, z2, y2, z3, y3. Their shape is:

```
phi(xi) = cosh(beta*xi) - cos(beta*xi)
          - gamma*(sinh(beta*xi) - sin(beta*xi))
```

The beta values differ from the unloaded-tip case; using unloaded cantilever
frequencies or multiplying all modes by a common tip-mass correction factor is incorrect.
The shape MAC excludes the auxiliary node to avoid a duplicated centerline tip
sample. DAT/CSV still retain it for every mode and axis.

## Preserved quality gates and scaling

Every mesh/mode uses the unchanged 1% analytical frequency, MAC≥0.99, direction
fraction≥0.95, clamp≤1e-8 and final-mesh-change≤0.5% gates from the unloaded
cantilever. The refinement trend check remains active. All six modes must also satisfy the
three MPCs within **2e-6 of peak nodal displacement**, allowing signed sums of
seven-significant-digit DAT values. This is not an attachment displacement budget.

Real solver tests additionally require:

- Doubling **only the tip mass** lowers all six tracked frequencies; each ratio
  agrees with the independently recomputed analytical ratio within 1%. This
  is a verification tolerance, not a claim that every frequency scales as 1/sqrt(Mt).
- Doubling **both beam density and tip mass** gives 1/sqrt(2) in every mode,
  within absolute ratio tolerance 2e-5 for printed precision.
- Multiplying E by four gives frequency ratio two, with the same 2e-5 allowance.
- Zero tip mass produces the original deck exactly and recovers the original
  real-solver benchmark within its existing analytical gates.

Each scaling case retains real inputs, outputs and a `*-scaling.json` file
containing all six actual ratios and independent reference ratios under
`evidence/tests/tip-mass-*/`. These are required integration tests in CI; the
standalone three-mesh CLI report does not by itself assert that the separate
scaling test suite has executed. CI retains the test log and its actual count.

## Negative inputs and semantic FRD verification

Negative tests cover g/N/density units for a mass, negative/nonfinite/implausible
mass, wrong axial location, eccentric placement, missing equations, wrong weights,
wrong DOFs, direct attachment to a solid node and elimination of the mass DOF.
One test corrupts the auxiliary displacement in actual solver data to prove the
coupling residual can fail despite successful solver execution.

FRD is checked semantically on **both** benchmarks during execution and replay:
coordinates and node IDs, all solid HEX20 connectivity (with the documented FRD
midside permutation), six DISP datasets, component labels, mode IDs, frequencies,
eigenvalues, solver version and all signed nodal components against DAT. The
FRD component tolerance is 1e-5 relative plus 1e-12 of peak displacement absolute,
because FRD prints six significant digits versus seven in DAT. Frequency and
eigenvalue agreement use the existing 2e-6 precision allowance. Only this ASCII
modal export is supported; unknown datasets, binary output, duplicated/missing
records and missing termination fail explicitly.

CalculiX omits one-node MASS elements and their auxiliary nodes from FRD. Requiring
them there would incorrectly reject a valid export. The audited deck and DAT
must contain them, DAT/CSV preserve their six full displacement vectors, and the
coupling residual covers their physics. FRD must cover **every solid node**.

Tests corrupt actual FRD geometry, topology, frequency, mode numbering, component
metadata, displacement values and termination. A replay test deliberately updates
the hash manifest after corrupting FRD components and still must fail the semantic
check. Hashes alone would miss this discrepancy. These checks establish cross-file
consistency, not authenticity against replacement of an entire coordinated dataset.
Unloaded cantilever archives remain replayable with the additional semantic checks.

## Reproduction and scope

The GitHub Actions workflow runs both cases, publishes both reports and retains
complete evidence for 30 days.
Contributors can use the same workflow or:

```sh
python -m unittest discover -s tests -v
python -m open_cae_evidence.benchmark --config benchmarks/cantilever_tip_mass.json --output evidence/tip-mass
python -m open_cae_evidence.benchmark --verify-evidence evidence/tip-mass
```

Use a new output directory. The existing Ubuntu/CalculiX installation and pins
in README apply. All code, configuration, reference equations and tests are MIT
licensed; CalculiX remains GPL licensed. No new dependency or solver is introduced.

This is numerical verification of a specified linear, undamped beam plus ideal
translational tip inertia. No experimental validation, finite payload rotary
inertia, attachment dynamics, torsional/axial analytical verification, convergence
order proof or industrial qualification is claimed.

Reference documentation: [CalculiX 2.21 manual](https://www.dhondt.de/ccx_2.21.pdf),
sections 6.2.44 (`MASS`), 7.55 (`*EQUATION`) and 7.81 (`*MASS`). The equations and
derivation above are open independent references, not proprietary solver results.
