#!/usr/bin/env python3
import csv
import argparse

def filter_pdb(input_csv, output_list, year_cutoff=2024):
    """
    Read `input_csv`, filter rows where deposition_year >= year_cutoff,
    and write pdb IDs one per line into `output_list`.
    """
    with open(input_csv, newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        # make sure the headers exist
        if 'deposition_year' not in reader.fieldnames or 'pdb' not in reader.fieldnames:
            raise ValueError("CSV must contain 'deposition_year' and 'pdb' columns")
        
        # filter and collect
        filtered = []
        for row in reader:
            try:
                year = int(row['deposition_year'])
            except ValueError:
                continue  # skip rows with non-integer year
            if year >= year_cutoff:
                filtered.append(row['pdb'].strip())
    
    # write out one pdb per line
    with open(output_list, 'w') as out:
        for pdb_id in filtered:
            out.write(f"{pdb_id}\n")

    print(f"Wrote {len(filtered)} entries to {output_list}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Filter a CSV of PDB entries by deposition_year and output a .list file"
    )
    parser.add_argument(
        "input_csv", 
        help="Path to input CSV (must have columns 'deposition_year' and 'pdb')"
    )
    parser.add_argument(
        "output_list", 
        help="Path to write filtered PDB IDs (one per line, e.g. filtered.list)"
    )
    parser.add_argument(
        "--year", "-y",
        type=int,
        default=2024,
        help="Minimum deposition_year to include (default: 2024)"
    )
    args = parser.parse_args()

    filter_pdb(args.input_csv, args.output_list, year_cutoff=args.year)
