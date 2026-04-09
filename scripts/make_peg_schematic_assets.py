#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Circle


ROOT = Path("/nfs/users/nfs_d/ds39/git_repos/Degradation_Mapper")


def add_box(ax, x, y, w, h, title, lines, face="#ffffff", edge="#444444"):
    box = FancyBboxPatch(
        (x, y),
        w,
        h,
        boxstyle="round,pad=0.02,rounding_size=0.02",
        linewidth=1.4,
        facecolor=face,
        edgecolor=edge,
    )
    ax.add_patch(box)
    ax.text(x + w / 2, y + h - 0.05, title, ha="center", va="top", fontsize=12, fontweight="bold")
    top = y + h - 0.12
    step = 0.085
    for i, line in enumerate(lines):
        ax.text(x + 0.04, top - i * step, line, ha="left", va="top", fontsize=11, family="monospace")


def draw_network(ax, center=(0.50, 0.53), scale=0.10):
    xs_left = [center[0] - scale * 1.2] * 3
    ys_left = [center[1] + scale, center[1], center[1] - scale]
    xs_mid = [center[0]] * 4
    ys_mid = [center[1] + scale * 1.45, center[1] + scale * 0.45, center[1] - scale * 0.45, center[1] - scale * 1.45]
    xs_right = [center[0] + scale * 1.2] * 2
    ys_right = [center[1] + scale * 0.5, center[1] - scale * 0.5]

    for xl, yl in zip(xs_left, ys_left):
        for xm, ym in zip(xs_mid, ys_mid):
            ax.plot([xl, xm], [yl, ym], color="#444444", lw=1.1, zorder=1)
    for xm, ym in zip(xs_mid, ys_mid):
        for xr, yr in zip(xs_right, ys_right):
            ax.plot([xm, xr], [ym, yr], color="#444444", lw=1.1, zorder=1)

    for x, y in zip(xs_left, ys_left):
        ax.add_patch(Circle((x, y), 0.012, facecolor="#6a5acd", edgecolor="#333333", lw=0.8, zorder=2))
    for x, y in zip(xs_mid, ys_mid):
        ax.add_patch(Circle((x, y), 0.012, facecolor="#ffcc4d", edgecolor="#333333", lw=0.8, zorder=2))
    for x, y in zip(xs_right, ys_right):
        ax.add_patch(Circle((x, y), 0.012, facecolor="#d9534f", edgecolor="#333333", lw=0.8, zorder=2))


def main():
    fig, ax = plt.subplots(figsize=(14, 7))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    ax.text(
        0.03,
        0.95,
        "PEG Example: chemistry schematic + mapper-derived descriptor outputs",
        fontsize=18,
        fontweight="bold",
        ha="left",
        va="top",
    )

    add_box(
        ax,
        0.03,
        0.54,
        0.28,
        0.28,
        "PEG",
        [
            "A  HO-(CH2-CH2-O)n-(CH2-CH2-O)m-OH",
            "",
            "Likely hotspot family:",
            "central ether backbone O-C bond",
        ],
        face="#fffdf7",
        edge="#444444",
    )

    add_box(
        ax,
        0.03,
        0.16,
        0.28,
        0.28,
        "Degradation Products (schematic)",
        [
            "B  HO-(CH2-CH2-O)n-CHO",
            "   OHC-(CH2-CH2-O)m-OH",
            "",
            "   HO-(CH2-CH2-O)n-COOH",
            "   HOOC-(CH2-CH2-O)m-OH",
        ],
        face="#fffaf9",
        edge="#444444",
    )

    arrow1 = FancyArrowPatch((0.33, 0.52), (0.43, 0.52), arrowstyle="-|>", mutation_scale=18, lw=2.2, color="#555555")
    ax.add_patch(arrow1)
    ax.text(0.38, 0.57, "choose representative\nPEG fragment", ha="center", va="bottom", fontsize=10)

    draw_network(ax, center=(0.50, 0.53), scale=0.085)
    ax.text(0.50, 0.35, "mapper-based O-C bond probe", ha="center", fontsize=12, fontweight="bold")
    ax.text(0.50, 0.30, "rigid bond-stretch scan", ha="center", fontsize=10)

    arrow2 = FancyArrowPatch((0.59, 0.52), (0.68, 0.52), arrowstyle="-|>", mutation_scale=18, lw=2.2, color="#555555")
    ax.add_patch(arrow2)
    ax.text(0.635, 0.57, "extract compact\nbond-fragility descriptors", ha="center", va="bottom", fontsize=10)

    descriptor_box = FancyBboxPatch(
        (0.70, 0.29),
        0.26,
        0.36,
        boxstyle="round,pad=0.03,rounding_size=0.02",
        linewidth=1.6,
        facecolor="#f7fbff",
        edgecolor="#3b6ea8",
    )
    ax.add_patch(descriptor_box)
    ax.text(0.72, 0.61, "PEG descriptor values", ha="left", va="top", fontsize=13, fontweight="bold")
    bullets = [
        "hotspot bond type: O-C",
        "barrier proxy: 4.330 eV",
        "stretch tolerance: 1.728 A",
        "peak force: 3.194 eV/A",
    ]
    for i, line in enumerate(bullets):
        ax.text(0.73, 0.56 - i * 0.07, f"• {line}", ha="left", va="top", fontsize=12.5)

    ax.text(
        0.70,
        0.20,
        "How it enters ML:\nall PEG/PEO formulation rows receive\nthis PEG family-level descriptor in the\ncurrent pilot.",
        ha="left",
        va="top",
        fontsize=11,
        bbox=dict(boxstyle="round,pad=0.3", facecolor="#fcfcfc", edgecolor="#777777"),
    )

    ax.text(
        0.03,
        0.04,
        "Note: the left panel is a schematic chemistry motivation figure. The mapper scan is a local O-C bond probe, not a full claim that one scan uniquely determines the exact final degradation products.",
        fontsize=10,
        ha="left",
        va="bottom",
        color="#444444",
    )

    out = ROOT / "figures" / "peg_descriptor_schematic.png"
    fig.savefig(out, dpi=220, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
