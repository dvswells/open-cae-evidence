# Contributing

Use a branch and pull request. The owner reviews from ChatGPT Android; CI reports
and artifacts must be sufficient for review without local tools. A merge requires
owner approval. Explain the physical problem, changed behavior, checks actually
executed and remaining limitations in the PR.

Use the automated Ubuntu workflow or the independent commands in README. Changes
to physics, reference equations, meshing or parsers need relevant positive and
negative tests and real-solver evidence. Preserve the failing case during diagnosis.
Do not replace a scientific test with a successful process exit or zero-test run.

Keep physical specifications, solver adapters, independent verification and
evidence responsibilities separate. Read [architecture](docs/architecture.md)
and [the public results contract](docs/results-schema.md) before extending an
adapter or analysis. Required quantities must have explicit units, locations,
components and actual solver/postprocessing provenance. Do not manufacture
missing values. Preserve legacy archive replay and unknown-version rejection.

CI must run both structural benchmarks, the thermal slab with CalculiX and Elmer,
all 48 preserved tests and the 26 tests in `tests_elmer`, comparison and replay.
See [Elmer reproduction and native semantics](docs/elmer.md). Do not require monotonic
convergence for a patch-exact linear thermal solution. The source distribution
includes the JSON Schema; packaging changes must preserve it.

Mission 5 additionally requires all 24 `tests_source` tests, both four-mesh
native sourced benchmarks, source comparison, numerical adversarial review and
their replay. Install the pinned `.[plots]` extra. See [source thermal](docs/source-thermal.md)
for nonnodal quadrature, printing bounds, reaction/source balance and distinct
raw/native/recovered flux semantics. Preserve all 74 previous tests and their
acceptance criteria. Conservation alone cannot verify the requested source.

Reference definitions must have equations, units, assumptions and an open source
or an independently reproducible derivation. Specify tolerances before comparing
new solver output; changes to acceptance criteria need physical justification.
Keep raw evidence and recorded versions. Avoid committing large solver binaries
or regenerable result files; Actions artifacts carry run-specific evidence.

License contributions under this repository's MIT license. Do not include
proprietary solver results or restricted reference datasets. CalculiX retains its
own GPL terms; Elmer retains its GPL/LGPL terms. Never include credentials in reports, commits or artifacts.
