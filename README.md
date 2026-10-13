# Open CAE Evidence

[![Modal verification](https://github.com/dvswells/open-cae-evidence/actions/workflows/modal.yml/badge.svg?branch=feat%2Fsource-driven-thermal-convergence)](https://github.com/dvswells/open-cae-evidence/actions/workflows/modal.yml?query=branch%3Afeat%2Fsource-driven-thermal-convergence)
[![Thermal verification](https://github.com/dvswells/open-cae-evidence/actions/workflows/elmer.yml/badge.svg?branch=feat%2Fsource-driven-thermal-convergence)](https://github.com/dvswells/open-cae-evidence/actions/workflows/elmer.yml?query=branch%3Afeat%2Fsource-driven-thermal-convergence)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

Reproducible numerical engineering benchmarks with independently replayable
evidence. Open CAE Evidence combines real CalculiX and Elmer simulations with
analytical references, physical checks and archived results to support evaluation
of open-source CAE tools.

## Benchmarks

| Physical problem | Solver | Independent checks | Methodology |
|---|---|---|---|
| Steel cantilever, 1 m × 15 mm × 10 mm | CalculiX | First six frequencies and mode shapes, Euler–Bernoulli references, three-mesh refinement and material scaling | [Modal analysis](docs/methodology.md) |
| Cantilever with a centered 1.17 kg tip mass | CalculiX | Full-face mass transfer, dynamic-tip characteristic equation, mode shapes and mass scaling | [Tip-mass formulation](docs/tip-mass.md) |
| Steady homogeneous thermal slab | CalculiX, Elmer | Fourier solution, temperature/flux fields, boundary flow, energy conservation and three-mesh consistency | [Source-free conduction](docs/thermal.md) |
| Thermal slab with uniform 100000 W/m³ generation | CalculiX, Elmer | Nonnodal volume errors, four-mesh convergence, source/reaction balance and separate native/recovered flux diagnostics | [Source-driven conduction](docs/source-thermal.md) |

The source-free linear slab is a patch-exact mesh-consistency test. The sourced
slab measures genuine discretization error: exact nodal values alone cannot
establish convergence. Cross-solver comparisons use matching physical coordinates
and explicitly documented field semantics; agreement is corroborating evidence.

## Verification and evidence

Physical specifications, solver adapters, independent verification and evidence
generation are separate components. All six solver/benchmark executions emit
the versioned `normalized-results.json` contract with units, mesh information,
boundary conditions, provenance, quantities, tolerances and verification outcomes.
See [architecture and interfaces](docs/architecture.md),
[results contract](docs/results-schema.md) and [Elmer runtime and semantics](docs/elmer.md).

GitHub Actions runs 98 tests, all benchmarks, cross-solver comparisons and archived
evidence replay. Negative tests reject corrupted inputs, outputs and rehashed
evidence, including numerically converged cases with incorrect physical inputs.

Open a workflow run to read its technical reports and test logs. Full
artifacts are not uploaded during normal CI runs. Repository maintainers can
select **Run workflow** in GitHub Actions with `upload_evidence` enabled to
retain the native solver inputs/outputs, logs, complete numerical fields, input
audits, normalized results, analytical definitions, convergence tables/plots,
failure evidence, build provenance and SHA-256 manifests for 30 days.
Independent users can reproduce these results from the source and pinned
dependencies. Replay re-extracts numerical fields and repeats physical checks;
hash integrity alone does not establish scientific correctness.

## Reproduction

Both workflows run automatically on pull requests and pushes to `main`. For independent
execution, use Ubuntu 24.04 and Python 3.12.14. Install CalculiX and follow the
[pinned headless Elmer setup](docs/elmer.md) before running all test suites:

```sh
sudo apt-get update
sudo apt-get install --no-install-recommends -y calculix-ccx=2.21-1
python -m pip install '.[plots]'
python -m unittest discover -s tests -v
python -m unittest discover -s tests_elmer -v
python -m unittest discover -s tests_source -v
python -m open_cae_evidence.benchmark --output evidence/benchmark
python -m open_cae_evidence.benchmark --verify-evidence evidence/benchmark
```

Use new output directories. Benchmark-specific commands and acceptance criteria
are in the methodology links above. Archived replay requires Python and the
recorded dependencies, but does not execute a solver. Dependencies and Actions
are pinned where practical; transitive system libraries and hosted runner images
can change, so their versions are recorded with the evidence.

## Scope and contribution

The current scope is ideal linear elastic cantilevers and homogeneous isotropic
steady conduction on affine rectangular elements. Assumptions and reference
limitations are explicit in each methodology. These benchmarks provide numerical
verification, not experimental validation, certification or industrial qualification.

See [contributor guidelines](CONTRIBUTING.md) for scientific quality gates and
pull request requirements.

## License

Project code, documentation, inputs, analytical reference definitions and generated
benchmark data are [MIT licensed](LICENSE). CalculiX is a separately distributed
GPL-licensed dependency. Elmer retains its upstream GPL/LGPL terms; compiled
artifacts include the upstream license texts.
