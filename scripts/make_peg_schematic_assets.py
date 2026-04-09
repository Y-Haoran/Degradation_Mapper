#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from ase import Atoms
from ase.data import atomic_numbers, covalent_radii
from ase.io import read
from matplotlib.patches import Circle, FancyArrowPatch


ROOT = Path("/nfs/users/nfs_d/ds39/git_repos/Degradation_Mapper")
PEG_DIR = Path(
    "/lustre/scratch127/mave/sge_analysis/team229/ds39/RD/LLaMA-Mesh/formulations_AI/OPoly26/"
    "evaluations/formulation_descriptor_family_expansion_v1/peg"
)

COLORS = {
    "carbon": "#2d2d2d",
    "oxygen": "#d62728",
    "bond": "#777777",
    "hotspot": "#d62728",
    "energy": "#1f77b4",
    "force": "#d62728",
    "panel": "#f6f6f6",
}


def draw_chain(ax, pts, oxy_idx=None, hotspot=None, linewidth=2.6):
    oxy_idx = set(oxy_idx or [])
    hotspot = set(hotspot or [])
    for i in range(len(pts) - 1):
        color = COLORS["hotspot"] if {i, i + 1} == hotspot else COLORS["bond"]
        lw = 4.2 if {i, i + 1} == hotspot else linewidth
        ax.plot(
            [pts[i][0], pts[i + 1][0]],
            [pts[i][1], pts[i + 1][1]],
            color=color,
            lw=lw,
            solid_capstyle="round",
            zorder=1,
        )
    for i, (x, y) in enumerate(pts):
        color = COLORS["oxygen"] if i in oxy_idx else COLORS["carbon"]
        radius = 0.012 if i in oxy_idx else 0.010
        ax.add_patch(Circle((x, y), radius, facecolor=color, edgecolor="white", lw=0.8, zorder=2))


def polymer_panel(ax):
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    top = np.array(
        [
            [0.08, 0.75],
            [0.17, 0.82],
            [0.26, 0.75],
            [0.35, 0.82],
            [0.44, 0.75],
            [0.53, 0.82],
            [0.62, 0.75],
            [0.71, 0.82],
            [0.80, 0.75],
            [0.89, 0.82],
        ]
    )
    draw_chain(ax, top, oxy_idx={2, 5, 8}, hotspot={5, 6})

    # degraded fragments
    left = np.array(
        [
            [0.10, 0.30],
            [0.19, 0.37],
            [0.28, 0.30],
            [0.37, 0.37],
            [0.46, 0.30],
        ]
    )
    right = np.array(
        [
            [0.56, 0.30],
            [0.65, 0.37],
            [0.74, 0.30],
            [0.83, 0.37],
            [0.92, 0.30],
        ]
    )
    draw_chain(ax, left, oxy_idx={2})
    draw_chain(ax, right, oxy_idx={2})

    # carbonyl hints
    ax.plot([0.46, 0.52], [0.30, 0.34], color=COLORS["bond"], lw=2.4)
    ax.plot([0.46, 0.52], [0.27, 0.31], color=COLORS["bond"], lw=1.2)
    ax.add_patch(Circle((0.54, 0.36), 0.013, facecolor=COLORS["oxygen"], edgecolor="white", lw=0.8))
    ax.plot([0.50, 0.56], [0.34, 0.30], color=COLORS["bond"], lw=2.4)
    ax.plot([0.50, 0.56], [0.31, 0.27], color=COLORS["bond"], lw=1.2)
    ax.add_patch(Circle((0.58, 0.28), 0.013, facecolor=COLORS["oxygen"], edgecolor="white", lw=0.8))


def bond_pairs(atoms: Atoms) -> list[tuple[int, int]]:
    pos = atoms.get_positions()
    syms = atoms.get_chemical_symbols()
    pairs = []
    for i in range(len(atoms)):
        ri = covalent_radii[atomic_numbers[syms[i]]]
        for j in range(i + 1, len(atoms)):
            rj = covalent_radii[atomic_numbers[syms[j]]]
            if np.linalg.norm(pos[i] - pos[j]) <= 1.18 * (ri + rj):
                pairs.append((i, j))
    return pairs


def fragment_panel(ax, atoms: Atoms, hotspot=(4, 5)):
    pos = atoms.get_positions().copy()
    pos -= pos.mean(axis=0, keepdims=True)
    rot = np.array([[0.86, -0.42, 0.30], [0.30, 0.89, 0.35], [-0.41, -0.17, 0.90]])
    pos = pos @ rot.T
    # weak perspective
    x = pos[:, 0] + 0.18 * pos[:, 2]
    y = pos[:, 1] + 0.08 * pos[:, 2]
    pts = np.c_[x, y]

    for i, j in bond_pairs(atoms):
        color = COLORS["hotspot"] if {i, j} == set(hotspot) else "#9a9a9a"
        lw = 4.4 if {i, j} == set(hotspot) else 1.8
        ax.plot([pts[i, 0], pts[j, 0]], [pts[i, 1], pts[j, 1]], color=color, lw=lw, solid_capstyle="round")

    sizes = {"H": 24, "C": 62, "O": 78}
    for idx, (xy, sym) in enumerate(zip(pts, atoms.get_chemical_symbols())):
        color = COLORS["oxygen"] if sym == "O" else ("#d9d9d9" if sym == "H" else COLORS["carbon"])
        edge = "#111111" if idx in hotspot else "white"
        lw = 1.4 if idx in hotspot else 0.6
        ax.scatter([xy[0]], [xy[1]], s=sizes.get(sym, 50), c=color, edgecolors=edge, linewidths=lw, zorder=3)

    ax.set_aspect("equal")
    ax.axis("off")


def scan_panel(ax, df: pd.DataFrame):
    x = df["bond_length_A"].to_numpy()
    e = df["relative_energy_eV"].to_numpy()
    f = df["mean_target_force_norm_eV_per_A"].to_numpy()
    e = e / e.max()
    f = f / f.max()
    ax.plot(x, e, color=COLORS["energy"], lw=3.2)
    ax.plot(x, f, color=COLORS["force"], lw=3.2)
    imax_e = int(np.argmax(e))
    imax_f = int(np.argmax(f))
    ax.scatter([x[imax_e]], [e[imax_e]], s=55, color=COLORS["energy"], zorder=4)
    ax.scatter([x[imax_f]], [f[imax_f]], s=55, color=COLORS["force"], zorder=4)
    ax.axvline(x[imax_f], color="#555555", lw=1.5, linestyle="--", alpha=0.7)
    ax.set_facecolor("white")
    for spine in ax.spines.values():
        spine.set_linewidth(1.2)
        spine.set_color("#444444")
    ax.set_xticks([])
    ax.set_yticks([])


def main():
    atoms = read(PEG_DIR / "optimized_fragment.xyz")
    metrics = pd.read_csv(PEG_DIR / "bond_scan_metrics.csv")

    fig = plt.figure(figsize=(13.5, 4.8))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.15, 1.1, 0.95], wspace=0.16)
    ax_left = fig.add_subplot(gs[0, 0])
    ax_mid = fig.add_subplot(gs[0, 1])
    ax_right = fig.add_subplot(gs[0, 2])

    polymer_panel(ax_left)
    fragment_panel(ax_mid, atoms)
    scan_panel(ax_right, metrics)

    fig.patches.extend(
        [
            FancyArrowPatch((0.33, 0.52), (0.39, 0.52), arrowstyle="-|>", mutation_scale=18, lw=2.0, color="#666666", transform=fig.transFigure),
            FancyArrowPatch((0.64, 0.52), (0.71, 0.52), arrowstyle="-|>", mutation_scale=18, lw=2.0, color="#666666", transform=fig.transFigure),
        ]
    )

    out_png = ROOT / "figures" / "peg_descriptor_schematic.png"
    out_pdf = ROOT / "figures" / "peg_descriptor_schematic.pdf"
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    fig.savefig(out_pdf, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
