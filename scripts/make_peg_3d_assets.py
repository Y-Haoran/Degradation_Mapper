#!/usr/bin/env python3
from __future__ import annotations

import shutil
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from ase import Atoms
from ase.data import atomic_numbers, covalent_radii
from ase.io import read


ROOT = Path("/nfs/users/nfs_d/ds39/git_repos/Degradation_Mapper")
SRC_XYZ = Path(
    "/lustre/scratch127/mave/sge_analysis/team229/ds39/RD/LLaMA-Mesh/formulations_AI/OPoly26/"
    "evaluations/formulation_descriptor_family_expansion_v1/peg/optimized_fragment.xyz"
)

HOTSPOT = {4, 5}
ELEMENT_COLORS = {
    "H": "#d9d9d9",
    "C": "#2b2b2b",
    "O": "#d62728",
    "N": "#1f77b4",
}


def bonds(atoms: Atoms) -> list[tuple[int, int]]:
    pos = atoms.get_positions()
    syms = atoms.get_chemical_symbols()
    pairs = []
    for i in range(len(atoms)):
        ri = covalent_radii[atomic_numbers[syms[i]]]
        for j in range(i + 1, len(atoms)):
            rj = covalent_radii[atomic_numbers[syms[j]]]
            cutoff = 1.18 * (ri + rj)
            if np.linalg.norm(pos[i] - pos[j]) <= cutoff:
                pairs.append((i, j))
    return pairs


def draw_panel(ax, atoms: Atoms, elev: float, azim: float, title: str) -> None:
    pos = atoms.get_positions().copy()
    pos -= pos.mean(axis=0, keepdims=True)

    for i, j in bonds(atoms):
        color = "#b8b8b8"
        lw = 1.6
        if {i, j} == HOTSPOT:
            color = "#d62728"
            lw = 4.2
        ax.plot(
            [pos[i, 0], pos[j, 0]],
            [pos[i, 1], pos[j, 1]],
            [pos[i, 2], pos[j, 2]],
            color=color,
            lw=lw,
            solid_capstyle="round",
            zorder=1,
        )

    for idx, (xyz, sym) in enumerate(zip(pos, atoms.get_chemical_symbols())):
        size = 160 if sym != "H" else 55
        edge = "#111111" if idx in HOTSPOT else "white"
        lw = 1.6 if idx in HOTSPOT else 0.6
        ax.scatter(
            [xyz[0]],
            [xyz[1]],
            [xyz[2]],
            s=size,
            c=ELEMENT_COLORS.get(sym, "#7f7f7f"),
            edgecolors=edge,
            linewidths=lw,
            depthshade=False,
            zorder=3,
        )
        if idx in HOTSPOT:
            ax.text(xyz[0] + 0.10, xyz[1] + 0.10, xyz[2] + 0.10, str(idx), color="#b2182b", fontsize=10, weight="bold")

    ax.view_init(elev=elev, azim=azim)
    lim = np.abs(pos).max() * 1.18
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_zlim(-lim, lim)
    ax.set_title(title, fontsize=12, fontweight="bold")
    ax.set_axis_off()


def main() -> None:
    atoms = read(SRC_XYZ)

    fig = plt.figure(figsize=(12.5, 6.4))
    ax1 = fig.add_subplot(1, 2, 1, projection="3d")
    ax2 = fig.add_subplot(1, 2, 2, projection="3d")
    draw_panel(ax1, atoms, elev=18, azim=-58, title="PEG fragment 3D view")
    draw_panel(ax2, atoms, elev=10, azim=32, title="Alternate view")

    fig.suptitle(
        "PEG/PEO representative fragment in 3D\nHotspot probe bond highlighted in red (atoms 4-5, O-C)",
        fontsize=15,
        fontweight="bold",
        y=0.98,
    )
    fig.text(
        0.5,
        0.03,
        "This is the actual optimized fragment used for the mapper PEG descriptor example.",
        ha="center",
        fontsize=11,
        color="#444444",
    )

    out_fig = ROOT / "figures" / "peg_fragment_3d.png"
    fig.savefig(out_fig, dpi=220, bbox_inches="tight")
    plt.close(fig)

    shutil.copy2(SRC_XYZ, ROOT / "results" / "peg_fragment_optimized.xyz")


if __name__ == "__main__":
    main()
