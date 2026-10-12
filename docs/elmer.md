# Independent Elmer thermal verification

## Runtime and trust

The Ubuntu 24.04 signed package indices queried in Actions contain no package
matching `^elmer`. We therefore build the official
[ElmerCSC/elmerfem](https://github.com/ElmerCSC/elmerfem) **release-9.0** source at
commit **9352634829dfde24439000571d64438d8c683264** (the annotated release tag resolves
to this commit, not to the tag object's hash). `scripts/build-elmer.sh` checks the
fetched commit and absence of tracked modifications. Git HTTPS verification and
signed Ubuntu APT validation remain enabled. No executable is obtained from a
third party, no shell installer is fetched, and no proprietary library is used.

The official README documents the Ubuntu CMake/C/C++/Fortran/BLAS/LAPACK build.
We explicitly disable GUI, ElmerPost, ElmerIce, contributed modules, MPI, OpenMP,
MUMPS, Hypre and Lua, and select generic BLAS. Two build workers limit memory use.
Only Solver_TGT, ElmerGrid, HeatSolve, FluxSolver and ResultOutputSolve and their
upstream dependencies are compiled. An initial all-module build failed in an
unrelated DCRComplexSolve electromagnetic module on modern gfortran; compiling
only the required targets avoids that module without editing upstream source.
The global upstream installer expects all optional modules, so our script copies
only the compiled thermal runtime and its required shared libraries/data files.
Build RPATHs are removed; ELMER_HOME and LD_LIBRARY_PATH select that runtime.
ElmerGrid is exercised on the official HeatUniso example, although the audited
slab generator writes native mesh files directly and needs no mesh conversion.

The workflow pins GCC/G++/gfortran 13.3.0 Ubuntu packages, CMake 3.28.3 and BLAS/
LAPACK development packages. It retains package versions, compiler output,
CMakeCache, source commit, configure/build/install logs, ldd dependency lists,
binary hashes, upstream license texts and a compiled-runtime archive. A clean
build is performed on each run, without cache-dependent claims. The build/job
budget is 35 minutes; each benchmark solver process is limited to 120 seconds.
The common Python version, package dependencies and Action revisions are pinned.
Hosted runner images, transitive system packages and timestamps can still change:
this establishes repeatable numerical execution, **not bitwise-reproducible binaries**.
ElmerSolver's library is LGPL; other Elmer components are GPL under upstream's
license policy. Project code/generated benchmark data remain MIT. Licenses are
copied with the compiled artifacts; upstream source is publicly available at the
recorded commit. We do not relicense or vendor upstream code in the project.

## Same physical problem, different implementation

Both adapters consume `load_problem(benchmarks/slab_steady_heat.json)` unchanged:
0.1 x 0.2 x 0.3 m; k=10 W/(m K); x_min=400 K; x_max=300 K; zero heat flux on the
other faces; no sources; meshes 2³, 4³ and 8³. Elmer's scalar HeatSolver solves
the Cartesian steady diffusion equation with standard 808 trilinear bricks.
The audited native 808 corner order is read from the pinned `elements.def`:
(-,-,-),(+,-,-),(+,+,-),(-,+,-),(-,-,+),(+,-,+),(+,+,+),(-,+,+).
It has the same affine Q1 approximation space and eight-point volume quadrature
as CalculiX DC3D8. Native node IDs are deliberately reversed relative to
CalculiX; sharing an FE approximation space is not sharing the solver executable
or its assembly/linear solution implementation.

Mesh files contain native 808 bodies and outward-oriented 404 boundary quads.
Audits independently parse coordinates, lattice membership, complete cells,
connectivity, positive affine Jacobians, face orientation/parent membership and
six physical face assignments. The SIF audit accepts only our canonical supported
keyword subset, rejecting property changes, missing/swapped boundary temperatures,
coordinate permutations, sources and unsupported extra keywords. SI units are
enforced by the physical specification, SIF and archived units.json. Elmer VTU
does not embed physical units; we explicitly declare that limitation and apply
no undisclosed conversion or coordinate transformation.

Elmer requests a Density field even for this steady problem. Density=1 is an
inert implementation requirement: no transient capacity, convection or source
term uses it here. It does not define another physical slab or enter the
analytical reference. The actual native temperature distribution, flux and
loads still must pass every original physical gate. The heat linear system uses
BiCGStab/ILU0 with 1e-12 relative residual and abort on nonconvergence. One steady
iteration solves this linear problem; relative change from a zero initial
field is not a nonlinear/transient convergence criterion. Flux projection uses
direct UMFPACK, avoiding ill-defined relative iterative residuals for nearly
zero transverse components. No solver warnings are accepted by the adapter.

## Native numerical quantities and normalization

ASCII Float64 VTU preserves temperature, Temperature Loads and Temperature Flux,
the full physical coordinates, VTK HEX8 connectivity and body GeometryIds.
The output is checked for one complete mesh, unique fields/coordinates, finite
values, exact field component counts, positive ordering and body assignment.
Output nodes are matched to the audited mesh by **physical coordinates** with
1e-12 m tolerance, never by native ordering. The `FORMAT` attribute spelling on
Elmer's native VTK cell-type array is explicitly supported. Unknown encodings,
precision, field mappings, missing fields and incomplete output fail extraction.
Native boundary surfaces are audited from serialized mesh.boundary and verified
through actual boundary temperatures/heat loads; VTU bulk output alone does not
claim to measure or store boundary-condition metadata.

- Temperature is a native nodal scalar, interpreted as K by the audited SI model.
- Temperature Loads are native `A_bulk T - b_bulk` before essential conditions;
  positive means heat supplied to the body. They are not calculated from the
  normalized recovered flux or the analytical reference.
- FluxSolver computes -k grad(T) and globally projects it by L2 into a continuous
  nodal Q1 field. Native Cartesian components are retained as optional
  `recovered_nodal_heat_flux` (source=solver).
- The required common `heat_flux` quantity samples that **actual native recovered
  field** by Q1 interpolation at the eight physical Gauss locations. Its source
  is **postprocessed_from_solver**, never native integration-point output. It is
  not a gradient substituted from analytical expectations.
- Common boundary heat-flow integration is unchanged: evaluate/extrapolate the
  recorded eight-point field at 2x2 face quadrature. For the constant slab flux
  it is exact up to numerical/recovery/output errors. Native loads provide a
  separate reaction-based conservation check.

No schema or ThermalVerifier change is needed. Contract 1.0.0 already supports
postprocessed field sources and optional quantities; every source/location note
must be respected. The existing analytical Fourier reference, all nine physical
checks, mesh study and tolerances are unchanged. Elmer's seven-significant-digit
ASCII writer limits distinguishable errors even though its arrays are Float64.
Zero printed temperature error therefore does not establish exact arithmetic.

## Cross-solver comparison and adversarial gates

`comparison.py` first replays **both native archives** and independently invokes
the same ThermalVerifier on both normalized contracts. Only then does it compare
complete temperature and flux fields, six outward face integrals, input reactions
and energy balances on all three mesh resolutions. Temperatures match bijectively
at coordinates without interpolation (1e-12 m). Gauss samples match bijectively
within 1e-7 m, covering seven-digit CalculiX DAT coordinate rounding in this slab.
No node/element/Gauss label equivalence is presumed. Elmer flux is interpolated
from native nodal recovery; CalculiX flux is native Gauss output. They represent
the same constant Cartesian vector at comparable physical points for this case;
their semantics must not be treated as identical on a general nonconstant field.

Comparison budgets were declared before the final comparison: 2e-3 K for T,
4e-6 relative for flux/heat flow/energy-balance difference, from the triangle
bounds of the unchanged individual 1e-3 K/2e-6 gates. All individual gradient,
Fourier closure, flux, reaction, conservation and mesh gates remain necessary.
The comparison separates physical-model sameness and affine Q1 discretization
from Elmer recovery and native output rounding. Agreement corroborates the
independent checks; agreement alone cannot establish correctness.

The separate `tests_elmer` suite runs real solvers and intentionally corrupts
SIF conductivity/units/temperatures, face assignments, topology/orientation,
coordinates, native VTU temperatures/flux/connectivity and completion/convergence
logs. Actual one-iteration failure and a real **converged swapped-temperature
solve** are retained. The latter conserves energy but fails the analytical
temperature/gradient/flux/flow gates, demonstrating why conservation and successful
execution alone are insufficient. Rehashed native/normalized/comparison archive
tampering is rejected by semantic replay, not only by checksum disagreement.
CI retains these failing inputs/outputs/logs and the failure verdicts alongside
the passing evidence. Hash manifests are integrity aids, not authentication of
a dataset against coordinated wholesale malicious replacement.

## Reproduction and current scope

From GitHub/ChatGPT Android, inspect **Independent Elmer thermal verification**
in the automatic push/PR checks. workflow_dispatch becomes available after the
workflow reaches the default branch with owner approval. No owner-operated machine is needed.
The Actions artifact contains the build/runtime proof, 48 original tests, Elmer
tests/failure cases, both modal archives, both thermal archives and comparison.
For independent contributors on an already provisioned open-source runtime:

```sh
python -m unittest discover -s tests -v
python -m unittest discover -s tests_elmer -v
python -m open_cae_evidence.benchmark --adapter elmer --config benchmarks/slab_steady_heat.json --output evidence/elmer-thermal
python -m open_cae_evidence.benchmark --adapter elmer --verify-evidence evidence/elmer-thermal
python -m open_cae_evidence.comparison --calculix evidence/thermal --elmer evidence/elmer-thermal --output evidence/comparison
python -m open_cae_evidence.comparison --calculix evidence/thermal --elmer evidence/elmer-thermal --output evidence/comparison --replay
```

Replay needs Python and archived files, not the solver/compiler or AI inference.
The adapter is intentionally restricted to the current affine homogeneous steady
slab. Modal Elmer, tetrahedra, curved meshes, sources, transients and nonlinear
conductivity are unsupported. No general recovery-error bound, experimental
validation, industrial qualification or ANSYS equivalence follows from this patch
test. A nonconstant Poisson/source reference would materially strengthen testing
and measure actual discretization/recovery error, but requires new physical
source/energy-balance verification. It is deferred to Mission 5 to preserve the
unchanged verifier and quality gates of this mission. No convergence order is
claimed for the exactly representable linear slab.
