"""Elmer 9.0 thermal adapter. Shared ProblemSpec and ThermalVerifier unchanged."""

import os
import subprocess

from ..interfaces import PreparedInput
from ..results import physical_mesh, quantity
from . import elmer_input, elmer_output
from ..thermal_problem import sourced


class ElmerAdapter:
    name = "Elmer"
    version = "1.0.0"

    def validate_provenance(self, identity):
        import re

        required = {
            "lib/libelmersolver.so",
            "lib/libmpi_stubs.so",
            "lib/libmatc.so",
            "lib/libfhuti.so",
            "lib/libarpack.so",
            "share/elmersolver/lib/HeatSolve.so",
            "share/elmersolver/lib/FluxSolver.so",
            "share/elmersolver/lib/ResultOutputSolve.so",
            "share/elmersolver/lib/elements.def",
            "share/elmersolver/lib/SOLVER.KEYWORDS",
        }
        hashes = identity.get("runtime_library_sha256", {})
        if (
            identity.get("upstream_source_commit")
            != "9352634829dfde24439000571d64438d8c683264"
            or identity.get("upstream_repository")
            != "https://github.com/ElmerCSC/elmerfem"
            or identity.get("adapter_version") != self.version
            or set(hashes) != required
            or any(
                not isinstance(h, str) or not re.fullmatch("[0-9a-f]{64}", h)
                for h in hashes.values()
            )
        ):
            raise ValueError("Incomplete/inconsistent Elmer runtime provenance")

    def runtime_provenance(self, executable):
        from ..evidence import sha256

        home = executable.resolve().parent.parent
        commit = (home / "source-commit.txt").read_text().strip()
        if commit != "9352634829dfde24439000571d64438d8c683264":
            raise ValueError("Elmer runtime source commit differs from audited build")
        identity = {
            "upstream_repository": "https://github.com/ElmerCSC/elmerfem",
            "upstream_source_commit": commit,
            "runtime_library_sha256": {
                str(p.relative_to(home)): sha256(p)
                for sub in ("lib", "share/elmersolver/lib")
                for p in sorted((home / sub).glob("*"))
                if p.is_file()
            },
            "runtime_configuration": "serial, no GUI/OpenMP/MPI; HeatSolve/FluxSolver/ResultOutputSolve",
        }
        self.validate_provenance({**identity, "adapter_version": self.version})
        return identity

    def prepare(self, problem, counts):
        if problem.analysis != "steady_heat_conduction":
            raise ValueError(
                "Elmer adapter supports only the verified steady slab; modal unsupported"
            )
        files, nodes, elements = elmer_input.generate(problem.parameters, counts)
        audit = elmer_input.audit(files, problem.parameters, counts)
        return PreparedInput(
            problem,
            counts,
            files,
            "slab",
            nodes,
            elements,
            audit,
            "hexahedron8",
            "full",
        )

    def audit_files(self, prepared, directory):
        files = {name: (directory / name).read_text() for name in prepared.files}
        return elmer_input.audit(files, prepared.problem.parameters, prepared.counts)

    def execute(self, prepared, directory, executable):
        self.audit_files(prepared, directory)
        env = os.environ.copy()
        env.update(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1")
        from pathlib import Path

        home = Path(executable).resolve().parent.parent
        env.update(
            ELMER_HOME=str(home),
            LD_LIBRARY_PATH=str(home / "lib")
            + os.pathsep
            + env.get("LD_LIBRARY_PATH", ""),
        )
        with (directory / "solver.log").open("w") as log:
            result = subprocess.run(
                [executable, "slab.sif"],
                cwd=directory,
                env=env,
                stdout=log,
                stderr=subprocess.STDOUT,
                timeout=120,
                check=False,
            )
        if result.returncode:
            raise RuntimeError(f"Elmer exited {result.returncode}; see solver.log")
        return elmer_output.check_log((directory / "solver.log").read_text())

    def extract(self, prepared, directory, solver_version):
        self.audit_files(prepared, directory)
        if (
            elmer_output.check_log((directory / "solver.log").read_text())
            != solver_version
        ):
            raise ValueError("Native solver version differs from archive")
        files = list(directory.glob("slab*.vtu"))
        if len(files) != 1:
            raise ValueError("Missing/ambiguous native Elmer output")
        native = files[0].name
        data = elmer_output.parse_vtu(
            files[0].read_text(), prepared.nodes, prepared.elements
        )
        quantities = {
            "temperature": quantity(
                "K",
                "node",
                ["temperature"],
                [
                    {"node": n, "value": v}
                    for n, v in sorted(data["temperature"].items())
                ],
                "solver",
                [native],
                "Native unitless scalar interpreted as K through audited SI SIF and units.json; no temperature conversion.",
            ),
            "nodal_heat_rate": quantity(
                "W",
                "node",
                ["heat_input"],
                [{"node": n, "value": v} for n, v in sorted(data["heat_rate"].items())],
                "solver",
                [native],
                "Temperature Loads = A_bulk*T-b_bulk, positive supplied heat; not derived from recovered flux.",
            ),
            "recovered_nodal_heat_flux": quantity(
                "W/m^2",
                "node",
                ["x", "y", "z"],
                [
                    {"node": n, "value": list(v)}
                    for n, v in sorted(data["nodal_flux"].items())
                ],
                "solver",
                [native],
                "Native FluxSolver global L2 projection of -k*grad(T) to continuous Q1 nodal field, Cartesian components.",
            ),
            "heat_flux": quantity(
                "W/m^2",
                "integration_point",
                ["x", "y", "z"],
                [
                    {
                        "element": e,
                        "point": p,
                        "coordinates": list(data["integration_coordinates"][e, p]),
                        "value": list(q),
                    }
                    for (e, p), q in sorted(data["heat_flux"].items())
                ],
                "postprocessed_from_solver",
                [native, "mesh.nodes", "mesh.elements"],
                "Q1 interpolation of native L2-recovered nodal Fourier flux at physical 2x2x2 Gauss points; NOT native integration-point flux. Coordinates m.",
            ),
        }
        normalized = {
            "mesh": physical_mesh(prepared),
            "solver": {
                "name": self.name,
                "version": solver_version,
                "adapter_version": self.version,
            },
            "quantities": quantities,
            "raw_output_verification": {
                "passed": True,
                "native_vtu": native,
                "nodes": len(prepared.nodes),
                "elements": len(prepared.elements),
                "coordinate_matching": True,
                "complete_temperature_loads_recovered_flux": True,
                "audited_boundary_assignments": True,
                "flux_semantics": "native L2 nodal projection; Q1 Gauss sampling",
                "completed_without_numerical_warnings": True,
            },
            "unsupported_capabilities": [
                "modal",
                "native_integration_point_flux",
                "native_boundary_surface_flux",
                "non_affine_meshes",
                "embedded_native_physical_units",
            ],
        }
        if sourced(prepared.problem.parameters):
            audit = self.audit_files(prepared, directory)
            normalized["raw_output_verification"]["volumetric_source"] = audit[
                "volumetric_source"
            ]
            normalized["raw_output_verification"][
                "thermal_flux_representation"
            ] = "global_l2_nodal_projection"
            normalized["raw_output_verification"][
                "temperature_print_resolution_K"
            ] = 1e-5
            normalized["raw_output_verification"][
                "recovered_flux_print_resolution_W_m2"
            ] = 0.001
        return normalized
