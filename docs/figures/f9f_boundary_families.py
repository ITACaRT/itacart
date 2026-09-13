"""Generate the F9f figures: absorbing cells and the polar cap.

Run from the repository root::

    python docs/figures/f9f_boundary_families.py

Writes into ``docs/_static/f9/``, where the ``phase_figures`` extension
discovers it. Every number plotted is computed from ``itacart`` itself at run
time; nothing here is transcribed, so a figure that stops matching the package
stops matching on the next run rather than drifting quietly.

Requires ``matplotlib``, which is in the ``notebooks`` extra rather than in
``docs`` — the documentation build reads the PNGs and never runs this file.
"""

from __future__ import annotations

import pathlib

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

import itacart  # noqa: E402
from itacart import boundary, geodesy  # noqa: E402

OUT = pathlib.Path(__file__).resolve().parents[1] / "_static" / "f9"
SIDE_M = 10_000.0
INK, ACCENT, MUTED = "#222222", "#b03a2e", "#8c8c8c"


def _style(ax, title: str, xlabel: str, ylabel: str) -> None:
    ax.set_title(title, fontsize=11, color=INK)
    ax.set_xlabel(xlabel, fontsize=9, color=INK)
    ax.set_ylabel(ylabel, fontsize=9, color=INK)
    ax.tick_params(labelsize=8, colors=INK)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    ax.grid(True, linewidth=0.4, alpha=0.3)


def absorbing_cell_against_its_neighbour() -> None:
    """An absorbing trapezoid beside the interior cell of the same row.

    Drawn on the sinusoidal plane, in metres, where the lattice is regular and
    the shapes are what the construction says they are. In geodetic degrees the
    two cells sit three degrees apart and each renders as a sliver, which shows
    nothing.
    """
    pairs = (
        ("SE(1900/0196)", MUTED, "interior"),
        ("SE(1930/0196)", ACCENT, "absorbing"),
    )

    fig, axes = plt.subplots(
        1, 2, figsize=(7.6, 3.8), dpi=160, sharex=True, sharey=True
    )
    for ax, (index, colour, label) in zip(axes, pairs):
        ring = itacart.cell_to_boundary(index, close=True)
        plane = [geodesy.geodetic_to_sinusoidal(lon, lat) for lon, lat in ring]
        xs = [x - plane[0][0] for x, _ in plane]
        ys = [y - plane[0][1] for _, y in plane]
        ax.plot(xs, ys, color=colour, linewidth=1.6)
        ax.fill(xs, ys, color=colour, alpha=0.14)
        ax.set_aspect("equal")
        share = itacart.effective_cell_area(index) / itacart.nominal_cell_area(1)
        _style(
            ax,
            f"{label} — {boundary.cell_shape(index)}",
            "metres east of the anchor",
            "metres north of the anchor",
        )
        ax.text(
            0.03,
            0.94,
            f"{index}\n{share:.3f} of the nominal area",
            transform=ax.transAxes,
            fontsize=8,
            color=colour,
            va="top",
            linespacing=1.5,
        )

    # Shared limits, or the trapezoid looks the same size as the parallelogram
    # and the figure argues the opposite of what it is for.
    axes[0].set_xlim(-11_000, 11_000)
    axes[0].set_ylim(-11_000, 1_500)
    axes[0].set_xticks([-10_000, -5_000, 0, 5_000, 10_000])
    axes[0].set_yticks([-10_000, -5_000, 0])

    fig.suptitle(
        "Absorption, on the plane where the shapes are regular", fontsize=11, color=INK
    )
    fig.tight_layout()
    fig.savefig(OUT / "f9_22_absorbing_cell_against_interior.png")
    plt.close(fig)


def lattice_narrowing_to_the_pole() -> None:
    """Columns per row, from the equator to the last row that exists.

    The tail carries the point and is invisible at full scale, so it gets an
    inset rather than an annotation that collides with itself.
    """
    rows = list(range(0, 1001, 5))
    columns = [boundary.last_lattice_column("SE", row, SIDE_M) for row in rows]

    fig, ax = plt.subplots(figsize=(7.0, 4.2), dpi=160)
    ax.plot(rows, columns, color=INK, linewidth=1.4)
    ax.fill_between(rows, columns, color=INK, alpha=0.08)
    _style(
        ax,
        "The lattice narrows towards the pole",
        "row index (0 at the equator)",
        "last addressable column",
    )

    tail_rows = list(range(985, 1001))
    tail = [boundary.last_lattice_column("SE", row, SIDE_M) for row in tail_rows]
    inset = ax.inset_axes((0.52, 0.42, 0.44, 0.46))
    inset.step(tail_rows, tail, where="post", color=ACCENT, linewidth=1.3)
    inset.set_title("the last fifteen rows", fontsize=8, color=ACCENT)
    inset.tick_params(labelsize=7, colors=INK)
    inset.grid(True, linewidth=0.4, alpha=0.3)
    for spine in ("top", "right"):
        inset.spines[spine].set_visible(False)
    inset.set_ylim(-4, max(tail) * 1.35)
    inset.annotate(
        "999: six cells in the world\n1000: two, both absorbing",
        xy=(999.2, 2),
        xytext=(985.4, max(tail) * 1.02),
        fontsize=7,
        color=ACCENT,
        linespacing=1.4,
        arrowprops=dict(arrowstyle="-", color=ACCENT, linewidth=0.6),
    )

    fig.tight_layout()
    fig.savefig(OUT / "f9_23_lattice_narrowing_to_the_pole.png")
    plt.close(fig)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    absorbing_cell_against_its_neighbour()
    lattice_narrowing_to_the_pole()
    print(f"wrote 2 figures into {OUT}")


if __name__ == "__main__":
    main()
