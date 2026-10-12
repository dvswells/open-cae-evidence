# Open CAE Evidence

Open-source numerical engineering benchmarks with independently replayable
evidence. The aim is to support industrial evaluation of open CAE tools through
reproducible verification. A passing benchmark is not a claim of equivalence to
ANSYS or proof of suitability for an industrial design.

## Mission 1: a verified modal cantilever

The first vertical slice runs **real CalculiX** on a 1 m × 15 mm × 10 mm
steel cantilever, computes its first six frequencies and complete nodal mode
shapes, and checks three successively refined C3D20R solid meshes. Independent
Euler–Bernoulli frequency and shape comparisons, serialized input audits,
physical input rejection and real-solver material scaling tests are required.
See [methodology](docs/methodology.md) for equations, assumptions and tolerances.

## Mission 2: concentrated tip mass

A second benchmark adds a centered 1.17 kg ideal point mass, transferred through
the exact area-average motion of the tip face with three multipoint equations.
It retains all Mission 1 gates, adds an independently solved tip-mass characteristic
equation, verifies mass-scaling and coupling, and checks archived FRD semantics
against the audited mesh and DAT. See [tip-mass formulation](docs/tip-mass.md).
The original benchmark remains unchanged and both cases run in CI.

## Mission 3: solver-independent evidence and steady thermal conduction

The two modal cases now run through separate physical specifications, a CalculiX
adapter, independent numerical verifiers and common evidence generation/replay.
All three benchmarks emit the versioned solver-independent
`normalized-results.json` contract, including units, mesh, BCs, provenance,
quantities, tolerances, physical outcomes and unsupported capabilities.
See [architecture and adapter interfaces](docs/architecture.md) and
[results schema](docs/results-schema.md).

The third benchmark solves a real 3D homogeneous slab with opposing prescribed
temperatures and insulated sides. Fourier references, nodal temperatures,
gradients, actual solver flux, integrated flow and independent reaction-based
energy conservation are checked on three meshes. This is a **mesh-consistency**
test: the exact linear solution belongs to the affine Q1 element space, so no
monotonic error reduction is required. See [thermal derivation](docs/thermal.md).

## Mission 4: independent Elmer interoperability

The same slab specification now runs with **CalculiX and Elmer 9.0**, using the
unchanged `ThermalVerifier` and results contract 1.0.0. Elmer is compiled headlessly
from a pinned official source commit. Audited native 808 meshes/SIF, complete VTU
temperatures, native heat loads and L2-recovered nodal flux are retained. Gauss
sampling of the recovered field is explicitly labeled as postprocessing.
Physical-coordinate comparisons require both independent verification and native
archive replay to pass first. See [Elmer runtime, semantics and limits](docs/elmer.md).

[The Elmer Actions workflow](.github/workflows/elmer.yml) runs all **48 preserved
tests plus 26 Elmer/interoperability tests**, all four solver/benchmark executions,
comparison and replay. It retains build provenance, complete native evidence and
actual converged-but-wrong boundary-condition failure evidence. Agreement is
corroboration, not proof of correctness or industrial qualification.

## Mission 5: volumetric source and measured convergence

Both solvers also solve the same slab with native 100000 W/m³ heat generation
on four meshes. Actual nodal DOFs define nonnodal volume temperature, gradient
and raw Fourier-flux errors; quadrature prevents exact nodal interpolation from
masquerading as zero discretization error. Reaction/source balance, separate
native/recovered flux diagnostics, coordinate-based comparison, convergence
plots and semantic archive replay are required. CI retains 74 existing tests
and adds 24 source/adversarial tests. See [source derivation, accuracy budgets
and reproduction](docs/source-thermal.md).

## Cloud and Android workflow

Development happens through Codex Cloud and GitHub. No local installation or
terminal is required from the owner. Every push and pull request automatically
runs [the Ubuntu workflow](.github/workflows/modal.yml). From GitHub on Android,
open **Actions → Verified cantilever modal benchmark → run** to read the technical
report in the run summary and download `cantilever-evidence-<run>-<attempt>`.
The workflow also supports **Run workflow** once it exists on the default branch.
PRs remain unmerged until owner approval.

Also inspect **Actions → Independent Elmer thermal verification** for the
`elmer-evidence-<run>-<attempt>` artifact and the common/cross-solver reports.
Both workflows run automatically on pushes and PRs; dispatch becomes available
from the default branch once the workflow is merged with owner approval.

Artifacts retain each solver input, raw `.dat` and `.frd`, solver logs, nodal CSV
for all six modes, audits, JSON verification, convergence results, SHA-256 manifest
and Markdown report. They are uploaded even after a failed check, when available,
and expire after 30 days. Download evidence for longer-term review.

## Independent reproduction

Scientists and contributors can reproduce this without an AI agent using the
same workflow, or run these commands on Ubuntu 24.04 with Python 3.12.14:

```sh
sudo apt-get update
sudo apt-get install --no-install-recommends -y calculix-ccx=2.21-1
python -m pip install .
python -m unittest discover -s tests -v
python -m open_cae_evidence.benchmark --output evidence/benchmark
python -m open_cae_evidence.benchmark --verify-evidence evidence/benchmark
python -m open_cae_evidence.benchmark --config benchmarks/cantilever_tip_mass.json --output evidence/tip-mass
python -m open_cae_evidence.benchmark --verify-evidence evidence/tip-mass
python -m open_cae_evidence.benchmark --config benchmarks/slab_steady_heat.json --output evidence/thermal
python -m open_cae_evidence.benchmark --verify-evidence evidence/thermal
```

Output directories must be new, so stale results cannot silently validate a run.
Simulation and verification need no AI service. Source convergence plots/replay
use the pinned open-source `.[plots]` extra. JSON Schema validation uses pinned
open-source `jsonschema` and pinned transitive dependencies; the numerical code
uses Python's standard library. The build backend and Actions
are pinned; Ubuntu package signatures and TLS verification stay enabled. Transitive
system libraries and hosted runner images can change; their versions are recorded.

## Scope and license

Only the two ideal linear elastic cantilevers and homogeneous isotropic steady
conduction on affine rectangular elements are covered. No
FreeCAD, graphical interface or optimizer is included. This is
numerical verification, not experimental validation or industrial qualification.
See [contributor guidelines](CONTRIBUTING.md) and [agent instructions](AGENTS.md).

Project code, documentation, input and analytical reference definitions are
MIT licensed. CalculiX is an independently distributed GPL-licensed dependency;
we do not relicense or vendor it. Elmer retains upstream GPL/LGPL terms;
compiled artifacts include the upstream license texts. Generated benchmark data is offered under MIT.

No MCP server or coupled-physics workflow is implemented. Elmer supports only the
verified steady slab here. The architecture documents extension points and limits.
