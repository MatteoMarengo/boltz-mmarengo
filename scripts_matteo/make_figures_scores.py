import os
import sys
import pandas as pd
import matplotlib.pyplot as plt

# Usage:
# python plot_metrics.py <input_csv> <output_dir>

def main(csv_path, out_dir):
    # Create output directory
    os.makedirs(out_dir, exist_ok=True)

    # Load data
    df = pd.read_csv(csv_path)

    # Convert columns to numeric
    df['align_rmsd'] = pd.to_numeric(df['align_rmsd'], errors='coerce')
    df['us_rmsd']    = pd.to_numeric(df['us_rmsd'],    errors='coerce')
    df['tm_score']   = pd.to_numeric(df['tm_score'],   errors='coerce')

    # 1. Histogram of sequence-based RMSD
    plt.figure()
    df['align_rmsd'].hist(bins=40)
    plt.xlabel('Sequence-based RMSD (Å)')
    plt.ylabel('Count')
    plt.title('Distribution of Sequence-based RMSD')
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, 'hist_align_rmsd.png'), dpi=300)
    plt.close()

    # 2. Histogram of US-align RMSD
    plt.figure()
    df['us_rmsd'].hist(bins=40)
    plt.xlabel('US-align RMSD (Å)')
    plt.ylabel('Count')
    plt.title('Distribution of US-align RMSD')
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, 'hist_us_rmsd.png'), dpi=300)
    plt.close()

    # 3. Histogram of TM-score
    plt.figure()
    df['tm_score'].hist(bins=40)
    plt.xlabel('TM-score')
    plt.ylabel('Count')
    plt.title('Distribution of TM-score')
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, 'hist_tm_score.png'), dpi=300)
    plt.close()

    # 4. Scatter: align_rmsd vs us_rmsd
    plt.figure()
    plt.scatter(df['align_rmsd'], df['us_rmsd'])
    plt.xlabel('Sequence-based RMSD (Å)')
    plt.ylabel('US-align RMSD (Å)')
    plt.title('Sequence vs US-align RMSD')
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, 'scatter_align_us.png'), dpi=300)
    plt.close()

    # 5. Scatter: us_rmsd vs tm_score
    plt.figure()
    plt.scatter(df['us_rmsd'], df['tm_score'])
    plt.xlabel('US-align RMSD (Å)')
    plt.ylabel('TM-score')
    plt.title('US-align RMSD vs TM-score')
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, 'scatter_us_tm.png'), dpi=300)
    plt.close()

    # 6. Scatter: align_rmsd vs tm_score
    plt.figure()
    plt.scatter(df['align_rmsd'], df['tm_score'])
    plt.xlabel('Sequence-based RMSD (Å)')
    plt.ylabel('TM-score')
    plt.title('Sequence RMSD vs TM-score')
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, 'scatter_align_tm.png'), dpi=300)
    plt.close()

    print(f"Saved figures in {out_dir}")


if __name__ == '__main__':
    if len(sys.argv) != 3:
        print('Usage: python plot_metrics.py <input_csv> <output_dir>')
        sys.exit(1)
    main(sys.argv[1], sys.argv[2])
