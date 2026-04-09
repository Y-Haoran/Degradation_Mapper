#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from ase import Atoms
from ase.io import read
from ase.data import atomic_numbers, covalent_radii


ROOT = Path("/nfs/users/nfs_d/ds39/git_repos/Degradation_Mapper")
SRC = Path(
    "/lustre/scratch127/mave/sge_analysis/team229/ds39/RD/LLaMA-Mesh/formulations_AI/OPoly26/"
    "evaluations/formulation_descriptor_family_expansion_v1/peg"
)


ELEMENT_COLORS = {
    "H": "#d9d9d9",
    "C": "#222222",
    "O": "#d62728",
    "N": "#1f77b4",
}


def atom_color(symbol: str) -> str:
    return ELEMENT_COLORS.get(symbol, "#7f7f7f")


def bond_pairs(atoms: Atoms) -> list[tuple[int, int]]:
    positions = atoms.get_positions()
    symbols = atoms.get_chemical_symbols()
    pairs: list[tuple[int, int]] = []
    for i in range(len(atoms)):
        zi = atomic_numbers[symbols[i]]
        ri = covalent_radii[zi]
        for j in range(i + 1, len(atoms)):
            zj = atomic_numbers[symbols[j]]
            rj = covalent_radii[zj]
            cutoff = 1.18 * (ri + rj)
            dist = np.linalg.norm(positions[i] - positions[j])
            if dist <= cutoff:
                pairs.append((i, j))
    return pairs


def draw_fragment(ax, atoms: Atoms, atom_i: int, atom_j: int) -> None:
    pos = atoms.get_positions()[:, :2]
    pos = pos - pos.mean(axis=0, keepdims=True)
    pos[:, 0] *= 1.05
    pairs = bond_pairs(atoms)
    for i, j in pairs:
        color = "#c8c8c8"
        lw = 2.0
        zorder = 1
        if {i, j} == {atom_i, atom_j}:
            color = "#d62728"
            lw = 5.0
            zorder = 3
        ax.plot(
            [pos[i, 0], pos[j, 0]],
            [pos[i, 1], pos[j, 1]],
            color=color,
            lw=lw,
            solid_capstyle="round",
            zorder=zorder,
        )
    for idx, (xy, symbol) in enumerate(zip(pos, atoms.get_chemical_symbols())):
        size = 140 if symbol != "H" else 60
        edge = "#000000" if idx in {atom_i, atom_j} else "white"
        linewidth = 2.2 if idx in {atom_i, atom_j} else 0.8
        ax.scatter(
            [xy[0]],
            [xy[1]],
            s=size,
            color=atom_color(symbol),
            edgecolors=edge,
            linewidths=linewidth,
            zorder=4,
        )
        if idx in {atom_i, atom_j}:
            ax.text(
                xy[0] + 0.08,
                xy[1] + 0.08,
                f"{idx}",
                fontsize=10,
                fontweight="bold",
                color="#b2182b",
            )
    ax.set_title("1. PEG fragment and hotspot bond", loc="left", fontsize=12, fontweight="bold")
    ax.text(
        0.01,
        0.97,
        "Representative dimethyl-capped PEG/PEO trimer\nHotspot probe bond: atoms 4-5 (O-C)",
        transform=ax.transAxes,
        va="top",
        fontsize=10,
        bbox=dict(boxstyle="round,pad=0.3", facecolor="#fff7f3", edgecolor="#d62728"),
    )
    ax.set_aspect("equal")
    ax.axis("off")


def draw_scan(ax, df: pd.DataFrame, stretch_tol: float, barrier: float, peak_force: float) -> None:
    x = df["bond_length_A"].to_numpy()
    rel_e = df["relative_energy_eV"].to_numpy()
    force = df["mean_target_force_norm_eV_per_A"].to_numpy()

    ax.plot(x, rel_e, color="#1f77b4", lw=2.6, label="Relative energy")
    ax.set_xlabel("Bond length (A)")
    ax.set_ylabel("Relative energy (eV)", color="#1f77b4")
    ax.tick_params(axis="y", labelcolor="#1f77b4")
    ax.grid(alpha=0.25, linestyle="--")
    ax.axvline(stretch_tol, color="#d62728", lw=1.8, linestyle="--")
    ax.axhline(barrier, color="#1f77b4", lw=1.2, linestyle=":")
    ax2 = ax.twinx()
    ax2.plot(x, force, color="#d62728", lw=2.2, label="Mean target force")
    ax2.set_ylabel("Mean target force (eV/A)", color="#d62728")
    ax2.tick_params(axis="y", labelcolor="#d62728")
    ax.set_title("2. Mapper-based bond scan", loc="left", fontsize=12, fontweight="bold")
    ax.text(
        0.01,
        0.98,
        f"Barrier proxy = {barrier:.3f} eV\nStretch tolerance = {stretch_tol:.3f} A\nPeak force = {peak_force:.3f} eV/A",
        transform=ax.transAxes,
        va="top",
        fontsize=10,
        bbox=dict(boxstyle="round,pad=0.3", facecolor="#f7fbff", edgecolor="#1f77b4"),
    )


def draw_workflow(ax, descriptors: dict, spec: dict, protocol: dict) -> None:
    ax.axis("off")
    ax.set_title("3. Descriptor outputs and ML join", loc="left", fontsize=12, fontweight="bold")
    text = (
        "Simple workflow\n"
        "  a. choose representative fragment\n"
        "  b. choose likely hotspot bond\n"
        "  c. run mapper bond scan\n"
        "  d. extract descriptor values\n\n"
        "PEG descriptor values\n"
        f"  hotspot bond type: {descriptors['hotspot_bond_type']}\n"
        f"  barrier proxy: {descriptors['barrier_proxy_eV']:.3f} eV\n"
        f"  stretch tolerance: {descriptors['stretch_tolerance_A']:.3f} A\n"
        f"  peak force: {descriptors['peak_force_eV_per_A']:.3f} eV/A\n\n"
        "How this enters ML\n"
        "  All PEG/PEO formulation rows\n"
        "  receive the same family-level PEG descriptor\n"
        "  in the current pilot.\n\n"
        "Example formulation row\n"
        "  polymer_family = PEG/PEO\n"
        f"  barrier_proxy_eV = {descriptors['barrier_proxy_eV']:.3f}\n"
        f"  hotspot_bond_type = {descriptors['hotspot_bond_type']}\n"
        f"  stretch_tolerance_A = {descriptors['stretch_tolerance_A']:.3f}\n"
        f"  peak_force_eV_per_A = {descriptors['peak_force_eV_per_A']:.3f}\n"
    )
    ax.text(
        0.02,
        0.98,
        text,
        va="top",
        fontsize=10.5,
        family="monospace",
        bbox=dict(boxstyle="round,pad=0.4", facecolor="#fcfcfc", edgecolor="#444444"),
    )


def main() -> None:
    with open(SRC / "descriptors.json") as f:
        payload = json.load(f)
    metrics = pd.read_csv(SRC / "bond_scan_metrics.csv")
    atoms = read(SRC / "optimized_fragment.xyz")

    descriptors = payload["descriptors"]
    protocol = payload["protocol"]
    spec = payload["spec"]
    atom_i = int(descriptors["target_atom_i"])
    atom_j = int(descriptors["target_atom_j"])

    fig = plt.figure(figsize=(15, 8), constrained_layout=True)
    gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 0.12], width_ratios=[1.0, 1.2])
    ax_frag = fig.add_subplot(gs[:, 0])
    ax_scan = fig.add_subplot(gs[0, 1])
    ax_text = fig.add_subplot(gs[1, 1])

    draw_fragment(ax_frag, atoms, atom_i, atom_j)
    draw_scan(
        ax_scan,
        metrics,
        descriptors["stretch_tolerance_A"],
        descriptors["barrier_proxy_eV"],
        descriptors["peak_force_eV_per_A"],
    )
    draw_workflow(ax_text, descriptors, spec, protocol)

    fig.suptitle(
        "PEG Example: from hotspot bond probe to mapper-derived formulation descriptor",
        fontsize=15,
        fontweight="bold",
        y=1.02,
    )

    out_fig = ROOT / "figures" / "peg_descriptor_example.png"
    fig.savefig(out_fig, dpi=220, bbox_inches="tight")
    plt.close(fig)

    result_df = pd.DataFrame(
        [
            {
                "polymer_key": payload["polymer_key"],
                "polymer_name": payload["polymer_name"],
                "fragment_id": spec["fragment_id"],
                "repeat_units": spec["repeat_units"],
                "end_capping_strategy": spec["end_capping_strategy"],
                "probe_atom_i": atom_i,
                "probe_atom_j": atom_j,
                "hotspot_bond_type": descriptors["hotspot_bond_type"],
                "barrier_proxy_eV": descriptors["barrier_proxy_eV"],
                "stretch_tolerance_A": descriptors["stretch_tolerance_A"],
                "peak_force_eV_per_A": descriptors["peak_force_eV_per_A"],
            }
        ]
    )
    result_df.to_csv(ROOT / "results" / "peg_descriptor_example.csv", index=False)
    metrics.to_csv(ROOT / "results" / "peg_bond_scan_metrics.csv", index=False)


if __name__ == "__main__":
    main()
