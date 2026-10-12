# Instructions for Codex and contributors

- Use the existing checkout in the isolated Codex Cloud environment. Do not make
  a worktree unless the owner explicitly asks.
- Develop through Codex Cloud and GitHub; the owner uses ChatGPT Android. Never
  require the owner to run terminal commands, install software locally or use a desktop.
- Open a PR for changes. Never merge without owner approval.
- Execute the real solver and scientific checks; no solver mocks as FEM evidence.
  An absent solver must fail integration tests, not skip them.
- Run `python -m unittest discover -s tests -v`, the complete three-mesh benchmark
  for both `benchmarks/cantilever.json` and `benchmarks/cantilever_tip_mass.json`
  in separate new output directories, and `--verify-evidence` on each before
  declaring readiness. Both archived DAT and FRD semantics must pass.
  Test runner exit status must survive log capture (`set -o pipefail`).
- Use the pinned Ubuntu Actions workflow for independent reproduction. Keep solver
  inputs, outputs, logs, provenance, hashes, reports and failed-run evidence.
- Keep numerical assumptions and justified acceptance criteria explicit. Do not
  widen tolerances or drop tests to make a failed run pass. Diagnose failures.
- All delivered source, tests and reference definitions must be open source.
  No proprietary solver, paid API or paid infrastructure. Prefer existing open
  tools. Do not add FreeCAD, MCP, another solver or optimization until this
  benchmark is verified and the owner authorizes expanding scope.
- Treat modal eigenvectors as sign/scale arbitrary. Match physical directions and
  shapes, not only ascending frequency indices. Check units and the serialized clamp.
- Keep the tip MASS node independent; transfer its translation through the full
  Q8 end-face average, not an arbitrary single-node attachment. Audit mass units,
  location, weights and equations, and verify their residuals for every mode.
  Tip-only mass changes do not obey a universal inverse-square-root law; use
  the independent dynamic-tip characteristic equation. See docs/tip-mass.md.
- Do not claim commercial CAE equivalence or experimental validation from this test.
- Inspect status before editing and preserve unrelated owner changes. Never expose
  credentials, disable signatures/TLS, or use remote shell installers.
