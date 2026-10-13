# Architecture and public interfaces

The architecture separates a physical problem/verification boundary and a
solver/serialization boundary. CalculiX and Elmer are separate
adapters; Elmer currently implements only the steady slab. Both use the same
ProblemSpec and ThermalVerifier, with no duplicated analytical reference.
The uniform-source slab uses the same interfaces; ThermalVerifier
dispatches by physical problem, not adapter. Source-free verification is unchanged.

```mermaid
flowchart LR
    P[Physical ProblemSpec] --> I[Adapter prepare: audited input]
    I --> X[Adapter execute: real solver]
    X --> N[Adapter extract: normalized fields]
    P --> V[Independent physical verifier]
    N --> V
    V --> E[Common JSON contract and evidence]
    R[Archived raw files] --> N
    E --> R
```

## Responsibility map

| Component | Responsibility | Must not do |
|---|---|---|
| `specification.py`, `modal_problem.py`, `thermal_problem.py` | Physical configuration, units, BCs, supported domain | Invoke a solver |
| `interfaces.py` | Small dataclasses and structural Python protocols | Own physics or a plugin framework |
| `adapters/calculix.py` | Generate input, execute, extract, normalize, validate raw output semantics | Fit or evaluate analytical references |
| `adapters/elmer.py`, `elmer_input.py`, `elmer_output.py` | Audit native SIF/808 mesh, execute pinned Elmer, extract VTU and explicitly normalize recovered flux | Change the physical verifier or label recovered flux as native integration-point data |
| `comparison.py` | Replay both archives, recheck independent gates and compare physical coordinates | Treat solver agreement as correctness |
| `model.py`, `thermal_input.py` | CalculiX input serializers/auditors | Claim successful output verification |
| `adapters/modal_dat.py`, `frd.py`, `thermal_output.py` | CalculiX DAT/FRD codecs and cross-file checks | Substitute analytical numbers for missing solver results |
| `reference.py`, `thermal_problem.analytical` | Independent analytical equations | Read FEM output |
| `source_reference.py`, `source_verification.py` | Source PDE/reference, volume error norms, source-aware balance and observed orders | Depend on solver file formats |
| `q1.py`, `source_quantities.py` | Physical Q1 interpolation/gradient, consistent source loads, common derived fields | Substitute reference fields for DOFs |
| `source_samples.py`, `source_comparison.py`, `convergence_plot.py`, `numerical_review.py` | Replayable point samples, common-field comparison, measured plots and adversarial numerical review | Treat nodal coincidence or conservation alone as correctness |
| `modal_verification.py`, `thermal_verification.py`, `verification.py` | Check normalized physical quantities, constraints, conservation, mesh studies | Parse solver files or select a solver |
| `results.py`, packaged JSON Schema | Versioned, unit-bearing interoperability contract | Treat schema validity as physical verification |
| `evidence.py`, report modules | Persistence, hashes, provenance, reports | Start user-operated services |
| `runner.py` | Connect the components, CLI, compatibility replay | Embed solver-specific reference physics |
| `benchmark.py` | Preserve legacy imports and CLI | Duplicate the implementation |

The public imports from `benchmark`, `model` and `reference` remain available.
Q8 transfer math is isolated
in `coupling.py`; the serialized equations remain the CalculiX adapter's responsibility.

## Minimal interfaces

`ProblemSpec` contains `id`, `analysis`, validated `parameters`, explicit `units`,
physical `boundary_conditions`, `required_quantities` and `unsupported_capabilities`.
`load_problem(config)` preserves unversioned modal physical input files;
the source-free thermal input explicitly declares `spec_version: "1.0"`.
The separate sourced input uses 1.1 with explicit W/m³ source and region.

`SolverAdapter` is a Python structural protocol:

1. `prepare(problem, counts) -> PreparedInput`: return relative text input files,
   an audited serialized deck, physical node coordinates/connectivity, canonical
   topology/integration information and explicit auxiliary-node roles. No run yet.
2. `execute(prepared, directory, executable) -> solver_version`: run the real solver
   noninteractively, preserve logs, inspect success/errors and fail on missing output
   or execution errors. A zero process exit is insufficient.
3. `extract(prepared, directory, solver_version) -> normalized_run`: read actual
   files, validate coverage/semantics, and return unit-bearing quantities in physical
   coordinates. Verify relevant logs/version on replay as well. Do not run a solver
   or read analytical references inside extraction.

`PreparedInput` is a frozen dataclass; its mapping payloads are ordinary Python data,
not deeply immutable objects. Input file names must be safe relative leaf names.
The adapter audits its serialized material/BCs against the physical specification.
The runner retains raw files and the audit and compares regenerated input on replay.

`NumericalVerifier.verify(problem, normalized_run)` returns `status`, Boolean
`checks` and numerical `metrics`. `study(runs)` verifies the mesh study. The modal
study retains its frequency-convergence criteria; the slab study is mesh consistency.
No verifier imports the CalculiX adapter or DAT/FRD codecs. Modal compatibility
routines adapt arbitrary external node labels using explicit `point_mass` roles;
the public normalized IDs are not changed.

`run_benchmark(config, output, executable, adapter=None)` accepts an explicit
adapter instance (default CalculiX). `replay(directory, adapter=None)` regenerates
input, re-extracts raw fields, repeats physical checks, validates the schema and
compares the contract with the archived one. No solver runs during replay.
An archive for another adapter requires its matching explicit adapter instance.
CLI selection is `--adapter calculix` (default) or `--adapter elmer`, with optional
`--solver` executable override. Elmer archives require `--adapter elmer` on replay.
Optional adapter hooks `runtime_provenance(executable)` and
`validate_provenance(identity)` record/check source and library identity; they
do not enter physical verification and are not required by the three-method
protocol. [Elmer](elmer.md) is the concrete second-adapter example.

## Results contract and evidence

See [results schema](results-schema.md) and the packaged
[`results-v1.json`](../src/open_cae_evidence/schemas/results-v1.json).
`normalized-results.json` is the public, solver-independent contract (1.0.0).
`results.json` is the internal evidence index (archive version 3) with retained
legacy modal fields for compatibility. Each mesh has `solver-results.json`, holding
its normalized fields and verification. Keeping these layers explicit avoids
silently assigning new meanings to legacy fields.

JSON Schema is validated using pinned, open-source `jsonschema`; physical coverage,
units and verdict consistency are checked additionally. This library is not a
simulation or orchestration framework. All schema resources ship in the Python wheel.

The common archive layout for all six solver/benchmark executions is:

```
config.json
provenance.json
results.json                 # evidence index v3; legacy fields retained
normalized-results.json      # public contract v1.0.0
report.md
manifest.json                # SHA-256 of every other file
mesh-<nx>x<ny>x<nz>/
  <beam|slab>.inp/.dat/.frd
  # Elmer instead: slab.sif, mesh.header/nodes/elements/boundary,
  # units.json, slab_t0001.vtu (native fields and physical topology)
  solver.log                 # plus other native solver files
  input-audit.json
  solver-results.json
  verification.json
  <mode-N|temperature|heat-flux>.csv
```

Sourced archives also contain `analytical-reference.json`, `convergence.svg`,
`convergence-plot-data.json`, and per-mesh `nonnodal-samples.csv`.
Their semantic replay recomputes all samples, norms, plots and source-aware
checks from native DOFs. Plotting is an optional pinned dependency extra needed
for sourced benchmark generation/replay. See [source case](source-thermal.md).

Only new output directories are accepted. Evidence with symlinks, unsafe manifest
paths, changed/incomplete file sets, unsupported versions or disagreeing metadata
fails replay. Hashes do not authenticate source or prevent coordinated replacement
of an entire dataset. Raw numerical re-extraction and independent physical checks
are necessary even if every hash matches. Raw output/saved fields are compared,
not rerun-solver eigenvectors, because signs and platform roundoff can differ.

Legacy archive versions **1 and 2** replay their real modal files through the same
adapter and physical checks. They lack the common contract and are not rewritten
in place. Automated tests construct legacy metadata layouts from real current
solver fields. Unknown archive or contract versions fail explicitly;
there is no `version >= latest` fallback.

## Implementing another adapter

Implement the three protocol methods and normalize to the existing physical contract.
For the current problems a second solver must preserve SI geometry/material/BCs,
canonical element ordering and physical field meanings, including input audits.
Map native eigenpairs to Hz, rad²/s² and sign/scale-arbitrary displacement vectors.
Map thermal quantities to K, W/m² and W, with global integration-point coordinates
and positive-input nodal heat-rate convention. Different native result quantities
must not be silently aliased. Omit unsupported quantities and declare the limitation;
the relevant benchmark cannot pass if its required quantity is unavailable.

Call `run_benchmark(..., adapter=YourAdapter())`. The verifier consumes the same
fields, so no physical verification logic needs to change. Add real execution,
negative parser/physical tests, source redistribution/licensing information and
CI evidence. The existing node-label interoperability test transforms identifiers
in real CalculiX data; it is not execution or validation of a second solver.

The present verifier/topology contract is intentionally narrow: structured affine
HEX8 slab fields and the current HEX20 modal meshes plus an optional point mass.
An adapter can normalize another solver on these meshes now. Supporting tetrahedra,
curved cells or different analyses requires an explicit schema/verifier extension;
changing the adapter alone must not make unsupported physics appear verified.

Elmer demonstrates the second adapter without modifying the physical verifier
or the common schema. Native nodal L2 flux is retained separately and sampled
at common physical Gauss coordinates as `postprocessed_from_solver`. It must
not be described as a native integration-point field. See [Elmer](elmer.md).

Future, unimplemented integrations:

- **Code_Aster** could supply a solver adapter with the same physical
  problem specifications, normalized fields and raw evidence obligations.
- **Gmsh** is a mesh provider, not a solver; it could supply mesh generation with
  independent geometry/topology audits before a solver adapter serializes the input.
- **FreeCAD** is a geometry/preprocessing source, not an analytical reference;
  its future headless use would enter before audited physical input generation.

No such additional tool, coupled-physics analysis or optimizer is implemented here.

## Downstream interoperability

External integrations should consume the static versioned JSON/artifacts,
honor units, location/component semantics, optional quantities, provenance and
unsupported-capability declarations. Schema validity is not a scientific verdict.
Integrations must inspect physical outcomes and raw-evidence replay, reject unsupported
schema versions and enforce size/path/access limits when handling untrusted artifacts.
Missing quantities must not be converted to zero, and benchmark outcomes must
not be presented as certification.
