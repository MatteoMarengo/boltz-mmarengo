#!/usr/bin/env python3
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.colors as mcolors
import matplotlib.cm as cm
import sys
import os

def generate_colored_table(csv_file, output_png):
    # Load data
    df = pd.read_csv(csv_file)

    # Prepare display DataFrame (format strings)
    display_df = df.copy()
    display_df['align_rmsd'] = display_df['align_rmsd'].map(lambda x: f"{x:.2f}")
    display_df['us_rmsd']    = display_df['us_rmsd'].map(lambda x: f"{x:.2f}")
    display_df['tm_score']   = display_df['tm_score'].map(lambda x: f"{x:.3f}")

    # Numeric arrays for colormapping
    align_vals = df['align_rmsd'].astype(float).values
    us_vals    = df['us_rmsd'].astype(float).values
    tm_vals    = df['tm_score'].astype(float).values

    # Normalizers
    align_norm = mcolors.Normalize(vmin=np.nanmin(align_vals), vmax=np.nanmax(align_vals))
    us_norm    = mcolors.Normalize(vmin=np.nanmin(us_vals),    vmax=np.nanmax(us_vals))
    tm_norm    = mcolors.Normalize(vmin=np.nanmin(tm_vals),    vmax=np.nanmax(tm_vals))

    # Colormaps
    cmap_rmsd = cm.get_cmap('RdYlGn_r')  # low-green, high-red
    cmap_tm   = cm.get_cmap('RdYlGn')    # low-red, high-green

    # Create figure sized by table dimensions
    nrows, ncols = display_df.shape
    fig, ax = plt.subplots(figsize=(ncols * 1.2, nrows * 0.4))
    ax.axis('off')

    # Create table
    table = ax.table(
        cellText=display_df.values,
        colLabels=display_df.columns,
        loc='center',
        cellLoc='center'
    )
    table.auto_set_font_size(False)
    table.set_fontsize(10)
    table.scale(1, 1.5)

    # Color cells
    for (row, col), cell in table.get_celld().items():
        # Header row
        if row == 0:
            cell.set_facecolor('#40466e')
            cell.set_text_props(color='white', weight='bold')
            continue
        # pdb_id column
        if col == 0:
            cell.set_facecolor('white')
            continue

        col_name = display_df.columns[col]
        if col_name == 'align_rmsd':
            val = float(df.at[row-1, 'align_rmsd'])
            color = cmap_rmsd(align_norm(val))
        elif col_name == 'us_rmsd':
            val = float(df.at[row-1, 'us_rmsd'])
            color = cmap_rmsd(us_norm(val))
        elif col_name == 'tm_score':
            val = float(df.at[row-1, 'tm_score'])
            color = cmap_tm(tm_norm(val))
        else:
            color = 'white'

        cell.set_facecolor(color)

    plt.tight_layout()
    plt.savefig(output_png, dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved colored table to {output_png}")

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print("Usage: python make_colored_table.py <input_csv> <output_png>")
        sys.exit(1)
    csv_file = sys.argv[1]
    out_png  = sys.argv[2]
    # Ensure output directory exists
    os.makedirs(os.path.dirname(out_png) or '.', exist_ok=True)
    generate_colored_table(csv_file, out_png)
