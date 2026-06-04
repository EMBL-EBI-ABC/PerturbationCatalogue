import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np

fig, ax = plt.subplots(figsize=(16, 7), dpi=300)
ax.set_xlim(0, 16)
ax.set_ylim(0, 7)
ax.axis("off")

# Colors
C_TEXT_MAIN = "#0f172a"
C_TEXT_SUB = "#64748b"
C_TEXT_TINY = "#94a3b8"
C_DIVIDER = "#e2e8f0"
C_CELL = "#f8fafc"
C_NUC = "#f1f5f9"
C_RED = "#ef4444"
C_BLUE = "#3b82f6"
C_GREEN = "#10b981"
C_BC_1 = "#8b5cf6"
C_BC_2 = "#f59e0b"
C_HANDLE = "#e2e8f0"

# Fonts
FONT_TITLE = {
    "fontsize": 18,
    "fontweight": "bold",
    "color": C_TEXT_MAIN,
    "ha": "center",
    "family": "sans-serif",
}
FONT_SUB = {"fontsize": 12, "color": C_TEXT_SUB, "ha": "center", "family": "sans-serif"}
FONT_TINY = {
    "fontsize": 9,
    "color": C_TEXT_TINY,
    "family": "sans-serif",
    "weight": "bold",
}
FONT_MONO = {
    "fontsize": 10,
    "family": "monospace",
    "color": "white",
    "ha": "center",
    "va": "center",
    "weight": "bold",
}
FONT_LABEL = {
    "fontsize": 11,
    "color": C_TEXT_MAIN,
    "family": "sans-serif",
    "ha": "center",
    "weight": "bold",
}

# 3 Columns: Centers at 2.66, 8.0, 13.33
ax.plot([5.33, 5.33], [0.5, 6.5], color=C_DIVIDER, linestyle="--", linewidth=1.5)
ax.plot([10.66, 10.66], [0.5, 6.5], color=C_DIVIDER, linestyle="--", linewidth=1.5)

# --- Phase 1: Transduction ---
ax.text(2.66, 6.2, "1. Transduction", **FONT_TITLE)
ax.text(2.66, 5.8, "CRISPR Perturbation Library", **FONT_SUB)
ax.text(2.66, 0.8, "POOL OF SINGLY PERTURBED CELLS", ha="center", **FONT_TINY)


def draw_cell(ax, x, y, radius, color):
    cell = patches.Circle(
        (x, y), radius, facecolor=C_CELL, edgecolor="#cbd5e1", linewidth=2.5, zorder=2
    )
    nuc = patches.Circle(
        (x - radius * 0.1, y + radius * 0.1),
        radius * 0.4,
        facecolor=C_NUC,
        edgecolor="#cbd5e1",
        linewidth=1.5,
        zorder=3,
    )
    ax.add_patch(cell)
    ax.add_patch(nuc)
    cas9 = patches.Circle(
        (x + radius * 0.4, y - radius * 0.3),
        radius * 0.2,
        facecolor="#e2e8f0",
        edgecolor="#94a3b8",
        linewidth=1.5,
        zorder=3,
    )
    ax.add_patch(cas9)
    ax.plot(
        [x + radius * 0.5, x + radius * 0.8],
        [y - radius * 0.4, y - radius * 0.1],
        color=color,
        lw=3,
        solid_capstyle="round",
        zorder=4,
    )
    ax.plot(
        [x + radius * 0.8, x + radius * 0.6],
        [y - radius * 0.1, y + radius * 0.2],
        color=color,
        lw=3,
        solid_capstyle="round",
        zorder=4,
    )


def draw_virus(ax, x, y, color):
    poly = patches.RegularPolygon(
        (x, y), 6, radius=0.25, facecolor="white", edgecolor="#64748b", lw=2, zorder=3
    )
    ax.add_patch(poly)
    ax.plot(
        [x - 0.12, x + 0.12],
        [y - 0.1, y + 0.1],
        color=color,
        lw=2.5,
        solid_capstyle="round",
        zorder=4,
    )


draw_cell(ax, 2.66, 4.2, 0.8, C_RED)
draw_virus(ax, 0.6, 5.0, C_RED)
ax.annotate(
    "",
    xy=(1.7, 4.5),
    xytext=(0.85, 4.9),
    arrowprops=dict(arrowstyle="->", color="#94a3b8", lw=2),
)

draw_cell(ax, 2.66, 2.0, 0.8, C_BLUE)
draw_virus(ax, 4.7, 2.8, C_BLUE)
ax.annotate(
    "",
    xy=(3.6, 2.35),
    xytext=(4.45, 2.7),
    arrowprops=dict(arrowstyle="->", color="#94a3b8", lw=2),
)

# --- Phase 2: Barcoding ---
ax.text(8.0, 6.2, "2. Molecular Barcoding", **FONT_TITLE)
ax.text(8.0, 5.8, "In-Droplet Reverse Transcription", **FONT_SUB)
ax.text(
    8.0, 0.8, "SHARED BARCODE LINKS PHENOTYPE AND GENOTYPE", ha="center", **FONT_TINY
)

# Bead Arc
arc = patches.Arc(
    (5.8, 3.1),
    width=2.5,
    height=4.5,
    angle=0,
    theta1=-90,
    theta2=90,
    color="#7dd3fc",
    lw=3,
)
ax.add_patch(arc)
ax.text(
    6.4, 4.8, "CELL LYSIS", fontsize=10, color="#94a3b8", style="italic", weight="bold"
)
ax.annotate(
    "",
    xy=(6.5, 4.5),
    xytext=(5.9, 4.7),
    arrowprops=dict(
        arrowstyle="-", color="#cbd5e1", linestyle="--", connectionstyle="arc3,rad=-0.2"
    ),
)

# Bracket
ax.plot([6.0, 5.8, 5.8, 6.0], [4.0, 4.0, 2.2, 2.2], color="#cbd5e1", lw=2)
ax.text(
    5.6, 3.1, "SAME CELL BARCODE", rotation=90, va="center", ha="center", **FONT_TINY
)


def draw_oligo(y_pos, capture_name, transcript_name, transcript_color, is_sgrna=False):
    # Backbone (zorder=1 to sit behind rects)
    ax.plot([6.0, 9.8], [y_pos, y_pos], color="#64748b", lw=2, zorder=1)

    # CBC
    ax.add_patch(
        patches.Rectangle((6.2, y_pos - 0.2), 0.8, 0.4, facecolor=C_BC_1, zorder=2)
    )
    ax.text(6.6, y_pos, "CBC", **FONT_MONO, zorder=3)

    # UMI
    ax.add_patch(
        patches.Rectangle((7.1, y_pos - 0.2), 0.6, 0.4, facecolor=C_BC_2, zorder=2)
    )
    ax.text(7.4, y_pos, "UMI", **FONT_MONO, zorder=3)

    # Capture sequence
    cap_width = 1.6 if is_sgrna else 1.3
    ax.add_patch(
        patches.Rectangle(
            (7.8, y_pos - 0.2), cap_width, 0.4, facecolor=C_HANDLE, zorder=2
        )
    )
    ax.text(
        7.8 + cap_width / 2,
        y_pos,
        capture_name,
        fontdict={
            "fontsize": 10,
            "family": "monospace",
            "color": C_TEXT_MAIN,
            "ha": "center",
            "va": "center",
            "weight": "bold",
        },
        zorder=3,
    )

    # Transcript pairing
    start_x = 7.8 + cap_width + 0.1
    if is_sgrna:
        x = np.linspace(start_x, 10.3, 100)
        y = y_pos + 0.2 + 0.15 * np.sin((x - start_x) * 12)
        ax.plot(x, y, color=transcript_color, lw=3.5, zorder=2)
        for bx in np.arange(start_x + 0.05, 9.7, 0.15):
            ax.plot(
                [bx, bx], [y_pos + 0.05, y_pos + 0.2], color="#94a3b8", lw=1.5, zorder=1
            )
    else:
        ax.plot(
            [start_x, 10.3],
            [y_pos + 0.2, y_pos + 0.2],
            color=transcript_color,
            lw=3.5,
            zorder=2,
        )
        for bx in np.arange(start_x + 0.05, 9.7, 0.15):
            ax.plot(
                [bx, bx], [y_pos + 0.05, y_pos + 0.2], color="#94a3b8", lw=1.5, zorder=1
            )

    ax.text(
        10.4,
        y_pos + 0.1,
        transcript_name,
        color=transcript_color if transcript_color != "#cbd5e1" else C_TEXT_MAIN,
        fontsize=12,
        va="center",
        ha="left",
        weight="bold",
    )


draw_oligo(3.9, "Poly(dT)", "Endogenous mRNA", "#94a3b8", False)
draw_oligo(2.5, "Capture Seq", "sgRNA Transcript", C_BLUE, True)

# --- Phase 3: Matrix ---
ax.text(13.33, 6.2, "3. Readout Matrix", **FONT_TITLE)
ax.text(13.33, 5.8, "Genotype-Phenotype Map", **FONT_SUB)
ax.text(13.33, 0.8, "SINGLE-CELL RESOLUTION ANALYSIS", ha="center", **FONT_TINY)

col1_x = 11.6
col2_x = 13.33
col3_x = 15.1

ax.text(col1_x, 4.8, "Cell Barcode", **FONT_LABEL)
ax.text(col2_x, 4.8, "Perturbation", **FONT_LABEL)
ax.text(col3_x, 4.8, "Transcriptome", **FONT_LABEL)

ax.plot([col1_x - 0.7, col1_x + 0.7], [4.6, 4.6], color="#cbd5e1", lw=1.5)
ax.plot([col2_x - 0.7, col2_x + 0.7], [4.6, 4.6], color="#cbd5e1", lw=1.5)
ax.plot([col3_x - 0.7, col3_x + 0.7], [4.6, 4.6], color="#cbd5e1", lw=1.5)

hm = ["#f1f5f9", "#cbd5e1", "#64748b", "#1e293b"]


def draw_row(y, cbc, guide, guide_color, hm_colors):
    # CBC
    ax.add_patch(
        patches.Rectangle((col1_x - 0.45, y - 0.18), 0.9, 0.36, facecolor=C_BC_1)
    )
    ax.text(col1_x, y, cbc, **FONT_MONO)

    # Guide
    ax.add_patch(
        patches.Rectangle((col2_x - 0.5, y - 0.18), 1.0, 0.36, facecolor=guide_color)
    )
    ax.text(col2_x, y, guide, **FONT_MONO)

    # Heatmap
    for i, c in enumerate(hm_colors):
        ax.add_patch(
            patches.Rectangle(
                (col3_x - 0.6 + i * 0.24, y - 0.15), 0.22, 0.3, facecolor=c
            )
        )


draw_row(4.0, "CBC_1", "sg-Red", C_RED, [hm[0], hm[3], hm[1], hm[0], hm[2]])
draw_row(3.3, "CBC_2", "sg-Blue", C_BLUE, [hm[2], hm[0], hm[3], hm[1], hm[0]])
draw_row(2.6, "CBC_3", "sg-Grn", C_GREEN, [hm[1], hm[1], hm[0], hm[3], hm[2]])

ax.text(col1_x, 1.9, "...", fontsize=16, color="#94a3b8", ha="center")
ax.text(col2_x, 1.9, "...", fontsize=16, color="#94a3b8", ha="center")
ax.text(col3_x, 1.9, "...", fontsize=16, color="#94a3b8", ha="center")

# Summary
ax.annotate(
    "",
    xy=(col2_x, 1.2),
    xytext=(col2_x, 1.7),
    arrowprops=dict(arrowstyle="->", color="#cbd5e1", lw=2),
)
ax.annotate(
    "",
    xy=(col3_x, 1.2),
    xytext=(col3_x, 1.7),
    arrowprops=dict(arrowstyle="->", color="#cbd5e1", lw=2),
)

ax.add_patch(
    patches.Rectangle(
        (col2_x - 0.75, 0.5),
        1.5,
        0.6,
        facecolor="white",
        edgecolor=C_RED,
        lw=2.5,
        zorder=5,
    )
)
ax.text(
    col2_x,
    0.8,
    "Genotype",
    fontdict={
        "fontsize": 13,
        "weight": "bold",
        "color": C_RED,
        "ha": "center",
        "va": "center",
    },
    zorder=10,
)

ax.text(
    (col2_x + col3_x) / 2,
    0.8,
    "=",
    fontdict={
        "fontsize": 24,
        "weight": "bold",
        "color": "#94a3b8",
        "ha": "center",
        "va": "center",
    },
)

ax.add_patch(
    patches.Rectangle(
        (col3_x - 0.75, 0.5),
        1.5,
        0.6,
        facecolor="white",
        edgecolor=C_BLUE,
        lw=2.5,
        zorder=5,
    )
)
ax.text(
    col3_x,
    0.8,
    "Phenotype",
    fontdict={
        "fontsize": 13,
        "weight": "bold",
        "color": C_BLUE,
        "ha": "center",
        "va": "center",
    },
    zorder=10,
)

plt.tight_layout()
plt.savefig(
    "slides/perturb-seq-reanalysis/what_is_perturb_seq.pdf",
    format="pdf",
    bbox_inches="tight",
    pad_inches=0.1,
)
