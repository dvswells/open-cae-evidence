# Public results contract 1.0.0

The normative structural schema is
[`src/open_cae_evidence/schemas/results-v1.json`](../src/open_cae_evidence/schemas/results-v1.json),
JSON Schema draft 2020-12. Python `results.validate_contract` additionally checks
finite numbers, physical specification consistency, required quantity units,
complete field coverage and verdict consistency. Independent verification is a
separate mandatory operation; schema validation cannot establish physics accuracy.

Root fields:

| Field | Meaning |
|---|---|
| `schema_name`, `schema_version` | `open-cae-evidence.results`, exact supported version `1.0.0` |
| `problem` | ID, analysis, physical units/parameters, explicit boundary conditions |
| `solver` | Adapter identity/version and executable hash, Python/platform, source SHA/dirty state, CI provenance |
| `runs` | Actual normalized fields, physical mesh, solver version, verification and raw-output checks per mesh |
| `verification` | Passed/failed verdict, dimensioned tolerances, mesh convergence or consistency study |
| `unsupported_capabilities` | Physical/reference scope not covered by this benchmark |

Mesh coordinates carry `coordinate_unit: "m"`. Nodes have stable positive integer
IDs, coordinates and `solid`/`point_mass` roles. Elements identify topology,
connectivity and full/reduced/not-applicable integration. HEX8 corner order is
bottom (000,100,110,010), top (001,101,111,011). HEX20 adds bottom/top/vertical
edge midsides in the standard input ordering; the adapter handles FRD's different
ordering. IDs are not required to match another solver's native labels.

Quantities are named maps; unsupported quantities are absent, not padded with zeros.
Each quantity has `unit`, `location`, ordered `components`, `values` and `source`.
Source identifies `solver` versus `postprocessed_from_solver`, file names within the
mesh directory and an explanatory note. Coordinates on integration-point records
use the mesh coordinate unit. Each record contains a scalar/vector `value` and the
IDs appropriate to its location: node, mode, node+mode, element+point+coordinates,
or boundary face.

| Analysis | Required quantity | Unit | Location | Meaning |
|---|---|---|---|---|
| Modal | `frequency` | Hz | mode | Positive cyclic natural frequency |
| Modal | `eigenvalue` | rad²/s² (`rad^2/s^2`) | mode | Squared angular frequency |
| Modal | `displacement` | `arbitrary` | node_mode | Solver-normalized eigenvector; no physical vibration amplitude is claimed |
| Steady heat | `temperature` | K | node | Absolute nodal temperature |
| Steady heat | `heat_flux` | W/m² (`W/m^2`) | integration_point | Global Fourier flux, with actual point coordinates |
| Steady heat | `nodal_heat_rate` | W | node | Positive heat input into the body; boundary reactions in the no-source case |
| Steady heat | `heat_flow` (derived) | W | boundary | Integrated outward flux from actual solver IP data |

Different analyses do not expose identical quantities. The CalculiX adapter does
not pretend that its extrapolated FRD nodal FLUX is the raw integration-point HFL,
or that surface integration is a native solver output. Per-run unsupported
capabilities make those limitations explicit. Required but unavailable quantities
cause failure; another benchmark may define a different defensible evidence policy.

For example, an actual thermal integration-point record has this structure
(the value depends on the solver run; this shape example is not evidence):

```json
{
  "unit": "W/m^2",
  "location": "integration_point",
  "components": ["x", "y", "z"],
  "values": [{"element": 1, "point": 1, "coordinates": [0.01, 0.02, 0.03], "value": [10000, 0, 0]}],
  "source": {"kind": "solver", "files": ["slab.dat"], "note": "Global heat flux"}
}
```

Dimensionless tolerances carry unit `1`; temperature absolute errors carry K.
Every mesh outcome has explicit Boolean checks and measured metrics, not only a
solver exit code. The archive contract includes a failed status for incomplete
execution when extractable evidence is available; a passed contract requires
exactly the prescribed mesh study and all required checks.

Versioning is explicit. The **results contract** version 1.0.0 differs from the
**archive index** versions 1/2/3 and the thermal **input specification** version 1.0.
Unknown versions are rejected pending an implemented compatibility/migration policy.
New topologies or analyses need a documented compatible version and tests. Existing
Mission 1/2 raw files and index fields remain available; legacy replay adds checks
without silently mutating or manufacturing missing contract files.
