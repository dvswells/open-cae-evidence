# Contributing

Use a branch and pull request. The owner reviews from ChatGPT Android; CI reports
and artifacts must be sufficient for review without local tools. A merge requires
owner approval. Explain the physical problem, changed behavior, checks actually
executed and remaining limitations in the PR.

Use the automated Ubuntu workflow or the independent commands in README. Changes
to physics, reference equations, meshing or parsers need relevant positive and
negative tests and real-solver evidence. Preserve the failing case during diagnosis.
Do not replace a scientific test with a successful process exit or zero-test run.

Reference definitions must have equations, units, assumptions and an open source
or an independently reproducible derivation. Specify tolerances before comparing
new solver output; changes to acceptance criteria need physical justification.
Keep raw evidence and recorded versions. Avoid committing large solver binaries
or regenerable result files; Actions artifacts carry run-specific evidence.

License contributions under this repository's MIT license. Do not include
proprietary solver results or restricted reference datasets. CalculiX retains its
own GPL terms. Never include credentials in reports, commits or artifacts.
