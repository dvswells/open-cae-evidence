"""Independent modal checks on normalized fields; no solver execution or parsing."""

from .coupling import active_tip, area_weights
from .reference import mac, references, shape
import math

TOLERANCES = {
    "analytical_relative": 0.01,
    "final_mesh_relative": 0.005,
    "shape_mac_min": 0.99,
    "direction_energy_min": 0.95,
    "fixed_displacement_relative": 1e-8,
    "refinement_noise_floor": 5e-5,
    "tip_coupling_relative": 2e-6,
}


def evaluate(config, nodes, frequencies, eigenvalues, modes):
    length = config["geometry"]["length"]
    auxiliary = {max(nodes)} if active_tip(config) else set()
    probes = sorted(
        (
            n
            for n, (_, y, z) in nodes.items()
            if n not in auxiliary and abs(y) < 1e-12 and abs(z) < 1e-12
        ),
        key=lambda n: nodes[n][0],
    )
    if (
        len(probes) < 10
        or nodes[probes[0]][0] != 0
        or not math.isclose(nodes[probes[-1]][0], length)
    ):
        raise ValueError("Insufficient centerline samples")
    records = []
    for index, (reference, frequency, lam, mode) in enumerate(
        zip(references(config), frequencies, eigenvalues, modes), 1
    ):
        axis = {"y": 1, "z": 2}[reference["direction"]]
        energy = [sum(v[a] ** 2 for v in mode.values()) for a in range(3)]
        if sum(energy) <= 0 or not math.isfinite(sum(energy)):
            raise ValueError("Mode shape has zero/nonfinite magnitude")
        fraction = energy[axis] / sum(energy)
        peak = max(abs(v) for xyz in mode.values() for v in xyz)
        clamp = (
            max(
                abs(v)
                for n, xyz in mode.items()
                if abs(nodes[n][0]) < 1e-12
                for v in xyz
            )
            / peak
        )
        samples = [mode[n][axis] for n in probes]
        analytical = [shape(reference["beta"], nodes[n][0] / length) for n in probes]
        assurance = mac(samples, analytical)
        relative = abs(frequency / reference["frequency_hz"] - 1)
        checks = {
            "analytical_frequency": relative <= TOLERANCES["analytical_relative"],
            "bending_direction": fraction >= TOLERANCES["direction_energy_min"],
            "analytical_shape": assurance >= TOLERANCES["shape_mac_min"],
            "clamped_displacements": clamp <= TOLERANCES["fixed_displacement_relative"],
        }
        records.append(
            {
                "mode": index,
                **reference,
                "fem_frequency_hz": frequency,
                "eigenvalue_rad2_s2": lam,
                "analytical_relative_error": relative,
                "centerline_mac": assurance,
                "direction_energy_fraction": fraction,
                "fixed_relative_displacement": clamp,
                "checks": checks,
                "passed": all(checks.values()),
            }
        )
    return records


def verify_modes(config, mesh, nodes, frequencies, eigenvalues, modes):
    records = evaluate(config, nodes, frequencies, eigenvalues, modes)
    if active_tip(config):
        weights = area_weights(nodes, config, mesh)
        for record, mode in zip(records, modes):
            peak = max(abs(v) for xyz in mode.values() for v in xyz)
            residual = (
                max(
                    abs(
                        mode[max(nodes)][a]
                        - sum(w * mode[n][a] for n, w in weights.items())
                    )
                    for a in range(3)
                )
                / peak
            )
            record["tip_coupling_relative_residual"] = residual
            # DAT prints seven significant digits. Signed sums need rounding allowance.
            record["checks"]["tip_coupling"] = (
                residual <= TOLERANCES["tip_coupling_relative"]
            )
            record["passed"] = all(record["checks"].values())
    return records


def convergence(runs):
    results = []
    for i in range(6):
        values = [run["modes"][i]["fem_frequency_hz"] for run in runs]
        changes = [abs(b - a) / b for a, b in zip(values, values[1:])]
        passed = changes[-1] <= TOLERANCES["final_mesh_relative"] and changes[
            -1
        ] <= max(changes[0], TOLERANCES["refinement_noise_floor"])
        results.append(
            {
                "mode": i + 1,
                "frequencies_hz": values,
                "relative_changes": changes,
                "passed": passed,
            }
        )
    return results
