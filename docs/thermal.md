# Independent 3D steady heat-conduction benchmark

## Definition and analytical derivation

A homogeneous slab occupies x∈[0,0.1] m, y∈[0,0.2] m, z∈[0,0.3] m.
Its isotropic conductivity is **10 W/(m K)**. The entire x=0 face is held at
400 K and the entire x=0.1 face at 300 K. Other faces are adiabatic. There is
no internal generation, thermal capacity, convection, radiation or mechanical
coupling in this reference. Conductivity and temperatures are definitions, not
measurements of a qualified material/system.

For constant positive k and no heat sources, div(q)=0 and Fourier q=-k grad(T)
give Laplace(T)=0. The opposing constant temperatures and insulated sides admit
the unique linear solution:

```
T(x) = 400 K - (1000 K/m)*x
grad(T) = [-1000, 0, 0] K/m
q = [10000, 0, 0] W/m²
A = 0.2*0.3 = 0.06 m²
Q = q_x*A = 600 W
```

The outward flow is -600 W at x_min, +600 W at x_max and zero on each side.
Heat supplied to the body is +600 W at the hot face and -600 W at the cold face.
At steady state the total outward flow and total supplied heat are separately zero.
The reference expressions in `thermal_problem.analytical` use only the physical
configuration. They never read numerical data or the finite-element interpolation.

## Real solver and independently audited input

CalculiX uses **DC3D8** eight-node, fully integrated thermal bricks (the documented
thermal alias of C3D8) under `*HEAT TRANSFER, STEADY STATE`. Temperature DOF 11 is
prescribed on every node of the two faces. Omitting other heat loads/film/radiation
imposes natural zero-normal-flux conditions on the sides. Steady pseudo-time is 1;
it does not represent transient heating time. Initial temperatures initialize the
solver but are not a heat-capacity model.

The three meshes are 2³, 4³ and 8³ cells: 8/64/512 elements and 27/125/729 nodes.
Serialization is audited for complete structured connectivity, orientation,
coordinates, volume, exact face membership, conductivity, prescribed temperatures,
procedure and output quantities. Additional loads/sources or missing/changed BCs
are rejected. SI unit labels are explicit; unknown/mislabelled units, nonfinite,
zero/negative conductivity, inverted or equal hot/cold assignments and unsupported
internal generation are rejected. Plausible numbers attached to false unit labels
cannot be discovered from numbers alone.

## Measured numerical quantities and conservation

The adapter extracts actual solver **NT**, **HFL**, **COORD** and **RFL** from DAT.
Temperature/rate output must include every node and HFL/COORD every one of the
eight integration points per element. Global integration coordinates are checked
against the complete Gauss-point set in each affine rectangular cell.

The independent verifier differentiates Q1 shape functions applied to the actual
nodal temperature output. These numerical gradients are compared both to the
closed-form gradient and to the solver HFL through q+k*grad(T)=0. This is a
numerical consistency check separate from the independent analytical reference;
its derivative implementation is not used to generate solver flux or reference T.

There is no native surface-integrated flux output in this adapter. Instead it
interpolates/extrapolates the actual eight-point flux field onto 2×2 quadrature
points on each exterior face, using tensor-product linear basis functions on the
Gauss lattice, then integrates q·n with actual cell-face areas. This is exact for
the Q1 gradient field on the affine rectangular cells used here. It is explicitly
postprocessed solver data, not an invented native output or an analytical value.
Curved cells, nonuniform k and other element types require different verified
integration rules and are unsupported.

All six surfaces are evaluated separately. RFL is independently summed on the hot
and cold faces and checked to be negligible on free-temperature nodes. Because no
imposed heat sources are present, the DAT label "heat generation" represents the
nodal heat input/reaction here, not volumetric generation. The verifier compares
each face flow to both its Fourier reference and the oppositely signed supplied
heat, checks side insulation, and checks both global balances. Swapped BCs can
conserve energy and still be physically wrong; direction/reference checks are
mandatory alongside conservation.

## Tolerances and mesh study

| Check | Limit | Rationale |
|---|---|---|
| Nodal T vs linear reference | 1e-3 K absolute | Conservative DAT printed-temperature allowance over the supported K domain |
| Numerical gradient and q+k grad(T) closure | 1e-5 relative to axial reference scale | Differentiating rounded nodal output amplifies its last digits |
| Actual HFL vector vs reference | 2e-6 of reference q_x | Seven-significant-digit DAT and solver roundoff; includes transverse components and sign |
| Each integrated end flow / reactions / side insulation | 2e-6 of reference Q | Actual numerical area integration and rounded nodal-rate summation |
| Global balances | 2e-6 of reference Q | Detect inconsistent energy, not merely balanced analytical constants |
| Across-mesh heat-flow spread | 2e-6 relative | Consistency of three independently solved meshes |

The exact linear solution lies in the affine Q1 space. All meshes can therefore
be essentially exact, with roundoff/printed-data error fluctuating. No monotonic
convergence requirement, inferred convergence order or discretization-error
extrapolation is appropriate. This study demonstrates **mesh consistency**;
it is not a meaningful general convergence-rate test for DC3D8.

FRD mesh, temperature and reaction fields are checked against the input and DAT.
Its nodal FLUX is an extrapolated/averaged diagnostic, not the raw solver HFL.
The constant-flux vertical slice showed approximately 4e-5 relative extrapolation
bias; a separate 1e-4 diagnostic agreement allowance covers this effect. This is
not a relaxation of any modal gate or the 2e-6 thermal physical flux gates. All
thermal physics/energy checks use the actual DAT integration-point field. A
future nonuniform-flux benchmark needs a different FRD comparison policy.

## Negative evidence, reproduction and limitations

Tests reject conductivity units/values, wrong/missing temperatures, source loads,
procedure/connectivity corruption and incomplete fields. They corrupt real DAT,
FRD, normalized temperature, heat-flux sign and reaction output. One deliberately
swapped-boundary **real solve** finishes and conserves energy but fails the
independent temperature/flux checks. Double-conductivity execution doubles the
actual flow without changing the temperature solution. Rehashed inconsistent
contracts still fail replay. All negative and positive test evidence is retained
with its expected outcome, and all previous 29 structural tests remain unchanged.

The existing pinned Ubuntu Actions workflow runs all three benchmarks, validates
the common results schema, replays all evidence and uploads raw and normalized
files, reports, provenance and test logs. The owner needs no terminal or server.
Contributors can independently reproduce without an AI:

```sh
python -m pip install .
python -m unittest discover -s tests -v
python -m open_cae_evidence.benchmark --config benchmarks/slab_steady_heat.json --output evidence/thermal
python -m open_cae_evidence.benchmark --verify-evidence evidence/thermal
```

The CalculiX/Ubuntu installation in README applies; output directories must be
new. Project inputs, references, code and generated benchmark data are MIT
licensed; the separately distributed solver retains its GPL license. Pinned
dependencies, raw files and transitive runtime versions are retained, but hosted
runner images are not immutable and hashes do not establish authenticity.

No nonlinear, anisotropic, transient, conjugate or coupled thermal problem is
verified. There is no experiment, hardware comparison, industrial qualification
or ANSYS-equivalence claim. See [architecture](architecture.md) and
[results contract](results-schema.md) for extension and interoperability limits.

Documentation: [CalculiX 2.21 manual](https://www.dhondt.de/ccx_2.21.pdf), `*ELEMENT`
thermal aliases, `*CONDUCTIVITY`, `*HEAT TRANSFER`, `*BOUNDARY`, `*NODE PRINT`
(NT/RFL), `*EL PRINT` (HFL/COORD). The Fourier derivation above is independent,
open reference material and does not use proprietary result data.
