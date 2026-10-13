"""Standalone scientific log-log SVG from measured error data, headless Agg."""

import io
import json
from .evidence import write_json

FIELDS = (
    ("temperature_L2_K_m3_2", "Temperature L2 [K m$^{3/2}$]"),
    ("gradient_H1_K_m1_2", "Gradient H1 seminorm [K m$^{1/2}$]"),
    ("raw_flux_L2_W_m_minus1_2", "Raw Fourier flux L2 [W m$^{-1/2}$]"),
    ("recovered_flux_L2_W_m_minus1_2", "Native recovered flux L2 [W m$^{-1/2}$]"),
)


def plot_data(series):
    return {
        "schema_version": "1.0.0",
        "x_unit": "m",
        "x_definition": "axial element length h_x",
        "series": {
            name: [
                {"h_x_m": m["h_x_m"], **{key: m[key] for key, _ in FIELDS}}
                for m in rows
            ]
            for name, rows in series.items()
        },
        "semantics": "Measured volume norms from actual numerical DOFs; recovered flux is separate. No nodal or center-only error surrogate.",
    }


def svg(data):
    import os, tempfile
    from pathlib import Path

    os.environ.setdefault(
        "MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "open-cae-evidence-mpl")
    )
    os.environ.setdefault(
        "XDG_CACHE_HOME", str(Path(tempfile.gettempdir()) / "open-cae-evidence-cache")
    )
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    with plt.rc_context(
        {
            "svg.hashsalt": "open-cae-evidence-convergence-v1",
            "font.family": "DejaVu Sans",
        }
    ):
        fig, axes = plt.subplots(2, 2, figsize=(10, 7), layout="constrained")
        for ax, (key, label) in zip(axes.flat, FIELDS):
            for name, rows in data["series"].items():
                ax.loglog(
                    [r["h_x_m"] for r in rows],
                    [r[key] for r in rows],
                    marker="o",
                    label=name,
                )
            ax.set_xlabel("Axial mesh length h_x [m]")
            ax.set_ylabel(label)
            ax.grid(True, which="both", alpha=0.3)
            ax.legend()
        fig.suptitle("Source-driven slab: measured nonnodal convergence")
        handle = io.StringIO()
        fig.savefig(handle, format="svg", metadata={"Date": None})
        plt.close(fig)
        return handle.getvalue()


def write(directory, series):
    data = plot_data(series)
    write_json(directory / "convergence-plot-data.json", data)
    (directory / "convergence.svg").write_text(svg(data))


def replay(directory, series):
    data = plot_data(series)
    if json.loads((directory / "convergence-plot-data.json").read_text()) != data or (
        directory / "convergence.svg"
    ).read_text() != svg(data):
        raise ValueError(
            "Convergence plot or plotted numerical data differs from replay"
        )
