"""CalculiX adapter: serialization, execution, extraction and raw-file semantics.

No analytical frequency, temperature or conservation reference is evaluated here.
"""

import os
import re
import subprocess

from .. import frd, model, thermal_input, thermal_output
from ..coupling import active_tip
from ..interfaces import PreparedInput
from ..results import physical_mesh, quantity
from .modal_dat import parse_dat


def check_solver_log(log):
    versions = re.findall(r"CalculiX Version ([\d.]+)", log)
    if "*ERROR" in log.upper() or "Job finished" not in log or len(versions) != 1:
        raise RuntimeError(
            "Solver log has errors, is incomplete or lacks a unique version"
        )
    return versions[0]


class CalculiXAdapter:
    name = "CalculiX"
    version = "1.0.0"

    def prepare(self, problem, counts):
        config = problem.parameters
        if problem.analysis == "steady_heat_conduction":
            deck, nodes, elements, audit = thermal_input.make_deck(config, counts)
            return PreparedInput(
                problem,
                counts,
                {"slab.inp": deck},
                "slab",
                nodes,
                elements,
                audit,
                "hexahedron8",
                "full",
            )
        deck, nodes, audit = model.make_deck(config, counts)
        block = (
            deck.split("*ELEMENT, TYPE=C3D20R, ELSET=BEAM\n", 1)[1]
            .split("*", 1)[0]
            .strip()
            .splitlines()
        )
        elements = {}
        for i in range(0, len(block), 2):
            e, *conn = [int(v) for v in (block[i] + block[i + 1]).split(",")]
            elements[e] = tuple(conn)
        auxiliary = frozenset({max(nodes)}) if active_tip(config) else frozenset()
        return PreparedInput(
            problem,
            counts,
            {"beam.inp": deck},
            "beam",
            nodes,
            elements,
            audit,
            "hexahedron20",
            "reduced",
            auxiliary,
        )

    def execute(self, prepared, directory, executable):
        env = os.environ.copy()
        env.update(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", CCX_NPROC_RESULTS="1")
        with (directory / "solver.log").open("w") as log:
            result = subprocess.run(
                [executable, "-i", prepared.stem],
                cwd=directory,
                env=env,
                stdout=log,
                stderr=subprocess.STDOUT,
                timeout=300,
                check=False,
            )
        if result.returncode != 0:
            raise RuntimeError(
                f"CalculiX failed (exit {result.returncode}); see {directory/'solver.log'}"
            )
        return check_solver_log((directory / "solver.log").read_text(errors="replace"))

    def extract(self, prepared, directory, solver_version):
        if check_solver_log((directory / "solver.log").read_text()) != solver_version:
            raise ValueError("Archived solver version differs from its log")
        stem = prepared.stem
        raw = [stem + ".dat"]
        if prepared.problem.analysis == "modal":
            f, lam, modes = parse_dat(
                (directory / (stem + ".dat")).read_text(), prepared.nodes
            )
            semantics = frd.verify_frd(
                (directory / (stem + ".frd")).read_text(),
                prepared.files[stem + ".inp"],
                prepared.nodes,
                f,
                lam,
                modes,
                solver_version,
                prepared.auxiliary_nodes,
            )
            quantities = {
                "frequency": quantity(
                    "Hz",
                    "mode",
                    ["frequency"],
                    [{"mode": i, "value": v} for i, v in enumerate(f, 1)],
                    "solver",
                    raw,
                ),
                "eigenvalue": quantity(
                    "rad^2/s^2",
                    "mode",
                    ["eigenvalue"],
                    [{"mode": i, "value": v} for i, v in enumerate(lam, 1)],
                    "solver",
                    raw,
                ),
                "displacement": quantity(
                    "arbitrary",
                    "node_mode",
                    ["x", "y", "z"],
                    [
                        {"node": n, "mode": i, "value": list(xyz)}
                        for i, mode in enumerate(modes, 1)
                        for n, xyz in sorted(mode.items())
                    ],
                    "solver",
                    raw,
                    "Solver-normalized eigenvectors have arbitrary sign and amplitude; these are not physical vibration amplitudes.",
                ),
            }
            unsupported = (
                ["physical_modal_amplitudes", "mass_nodes_in_frd"]
                if prepared.auxiliary_nodes
                else ["physical_modal_amplitudes"]
            )
        else:
            data = thermal_output.parse_dat(
                (directory / (stem + ".dat")).read_text(),
                prepared.nodes,
                prepared.elements,
            )
            semantics = thermal_output.verify_frd(
                (directory / (stem + ".frd")).read_text(),
                prepared.nodes,
                prepared.elements,
                data,
                solver_version,
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
                    raw,
                ),
                "nodal_heat_rate": quantity(
                    "W",
                    "node",
                    ["heat_input"],
                    [
                        {"node": n, "value": v}
                        for n, v in sorted(data["heat_rate"].items())
                    ],
                    "solver",
                    raw,
                    "Positive means heat supplied to the body; RFL is boundary reaction input here because all imposed sources are absent.",
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
                    "solver",
                    raw,
                    "Global Fourier heat-flux vector; integration-point coordinates use mesh coordinate_unit m.",
                ),
            }
            unsupported = [
                "native_boundary_surface_flux",
                "exact_nodal_flux_extrapolation",
            ]
        return {
            "mesh": physical_mesh(prepared),
            "solver": {
                "name": self.name,
                "version": solver_version,
                "adapter_version": self.version,
            },
            "quantities": quantities,
            "raw_output_verification": semantics,
            "unsupported_capabilities": unsupported,
        }
