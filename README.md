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

## Cloud and Android workflow

Development happens through Codex Cloud and GitHub. No local installation or
terminal is required from the owner. Every push and pull request automatically
runs [the Ubuntu workflow](.github/workflows/modal.yml). From GitHub on Android,
open **Actions → Verified cantilever modal benchmark → run** to read the technical
report in the run summary and download `cantilever-evidence-<run>-<attempt>`.
The workflow also supports **Run workflow** once it exists on the default branch.
PRs remain unmerged until owner approval.

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
```

Output directories must be new, so stale results cannot silently validate a run.
Runtime Python code uses only the standard library. The build backend and Actions
are pinned; Ubuntu package signatures and TLS verification stay enabled. Transitive
system libraries and hosted runner images can change; their versions are recorded.

## Scope and license

Only a uniform, slender, isotropic, undamped, unprestressed, linear elastic
cantilever, with or without ideal translational tip inertia, is covered. No
FreeCAD, extra solver or optimizer is included. This is
numerical verification, not experimental validation or industrial qualification.
See [contributor guidelines](CONTRIBUTING.md) and [agent instructions](AGENTS.md).

Project code, documentation, input and analytical reference definitions are
MIT licensed. CalculiX is an independently distributed GPL-licensed dependency;
we do not relicense or vendor it. Generated benchmark data is offered under MIT.
