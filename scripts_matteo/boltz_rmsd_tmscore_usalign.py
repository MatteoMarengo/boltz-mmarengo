#!/usr/bin/env python3
import os
import sys
import csv
import re
import subprocess
import requests

# Path to the USalign binary. If it's on your PATH, leave as-is.
USALIGN_CMD = "USalign"

BOLTZ_CSV_HEADER = ["pdb_id","align_rmsd","align_atoms"]
OUTPUT_HEADER    = [
    "pdb_id",
    "align_rmsd", "align_atoms",
    "us_rmsd", "us_atoms",
    "tm_score"
]

def download_mmCIF(pdb_id, dest_dir):
    pdb = pdb_id.upper()
    fn = os.path.join(dest_dir, f"{pdb}.cif")
    if not os.path.exists(fn):
        url = f"https://files.rcsb.org/download/{pdb}.cif"   
        print(f"  ↳ downloading reference {pdb}.cif …")
        r = requests.get(url); r.raise_for_status()
        with open(fn, "wb") as f:
            f.write(r.content)
    return fn

def compute_usalign(ref_cif, pred_cif):
    """
    Runs USalign and returns:
        (rmsd: float, aligned_length: int, tm_score: float)
    or (None, None, None) on failure.
    """
    cmd = [USALIGN_CMD, ref_cif, pred_cif, "-outfmt", "2"]  # quiet mode
    try:
        proc = subprocess.run(cmd, capture_output=True, text=True, check=True)
    except FileNotFoundError:
        print("  ⚠️  USalign not found (check USALIGN_CMD or your PATH)")
        return None, None, None
    except subprocess.CalledProcessError as e:
        print("  ⚠️  USalign returned non-zero exit code")
        if e.stdout:
            print("    --- USalign stdout: ---")
            print("\n".join("    " + line for line in e.stdout.splitlines()))
        if e.stderr:
            print("    --- USalign stderr: ---")
            print("\n".join("    " + line for line in e.stderr.splitlines()))
        return None, None, None

    # Skip header line and parse tab-separated values
    lines = proc.stdout.strip().splitlines()
    if len(lines) < 2:
        print("  ⚠️  USalign output has no data lines")
        print("    --- Full USalign output: ---")
        for line in lines:
            print("    " + line)
        return None, None, None

    # Parse second line (first data line)
    cols = lines[1].split()
    if len(cols) < 11:
        print("  ⚠️  unexpected number of columns in USalign output")
        return None, None, None

    try:
        tm_score = float(cols[2])  # TM1 column
        rmsd     = float(cols[4])  # RMSD column
        length   = int(cols[10])   # Lali column (aligned residue count)
    except ValueError as e:
        print(f"  ⚠️  failed to parse numeric values: {e}")
        return None, None, None

    return rmsd, length, tm_score


def main(main_dir, boltz_csv, output_csv):
    # Read the Boltz CSV
    with open(boltz_csv, newline="") as f:
        reader = csv.DictReader(f)
        if reader.fieldnames != BOLTZ_CSV_HEADER:
            raise ValueError(f"Expected header {BOLTZ_CSV_HEADER}, got {reader.fieldnames}")
        boltz_rows = list(reader)

    results = []
    for row in boltz_rows:
        pdb_id = row["pdb_id"]
        subdir = os.path.join(main_dir, pdb_id)
        base   = f"{pdb_id}_polymer_only"
        boltz_out = os.path.join(subdir, f"boltz_results_{base}")
        pred_cif  = os.path.join(
            boltz_out,
            "predictions",
            base,
            f"{base}_model_0.cif"
        )
        if not os.path.exists(pred_cif):
            print(f"[{pdb_id}] ⚠️ predicted CIF not found at {pred_cif}, skipping")
            continue

        # Ensure we have the reference CIF
        ref_cif = os.path.join(subdir, f"{pdb_id}.cif")
        if not os.path.exists(ref_cif):
            ref_cif = download_mmCIF(pdb_id, subdir)

        print(f"[{pdb_id}] running USalign…")
        us_rmsd, us_atoms, tm = compute_usalign(ref_cif, pred_cif)
        if us_rmsd is not None:
            print(f"  → USalign RMSD = {us_rmsd:.3f} Å over {us_atoms} residues")
            print(f"  → TM-score (from USalign) = {tm:.3f}")

        results.append({
            "pdb_id":      pdb_id,
            "align_rmsd":  row["align_rmsd"],
            "align_atoms": row["align_atoms"],
            "us_rmsd":     f"{us_rmsd:.3f}" if us_rmsd is not None else "",
            "us_atoms":    us_atoms or "",
            "tm_score":    f"{tm:.3f}"       if tm is not None else ""
        })

    # Write out the combined CSV
    with open(output_csv, "w", newline="") as outf:
        writer = csv.DictWriter(outf, fieldnames=OUTPUT_HEADER)
        writer.writeheader()
        writer.writerows(results)

    print(f"\nWrote USalign & TM-score results to {output_csv}")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: boltz_rmsd_usalign.py <main_folder> <boltz_csv> <output_csv>")
        sys.exit(1)
    main(sys.argv[1], sys.argv[2], sys.argv[3])