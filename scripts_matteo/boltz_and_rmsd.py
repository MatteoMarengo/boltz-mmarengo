# #!/usr/bin/env python3
# import os, sys, csv, subprocess, requests
# from Bio.PDB import MMCIFParser, Superimposer
# from collections import defaultdict

# BOLTZ_CMD = "boltz"

# def download_mmCIF(pdb_id, dest_dir):
#     pdb = pdb_id.upper()
#     fn = os.path.join(dest_dir, f"{pdb}.cif")
#     if not os.path.exists(fn):
#         url = f"https://files.rcsb.org/download/{pdb}.cif"
#         print(f"  ↳ downloading reference {pdb}.cif …")
#         r = requests.get(url); r.raise_for_status()
#         with open(fn, "wb") as f: f.write(r.content)
#     return fn

# def extract_polymer_atoms(structure):
#     atoms = {}
#     model = next(structure.get_models())
#     for chain in model:
#         for res in chain:
#             het, resseq, icode = res.get_id()
#             if het == " ":
#                 for atom in res:
#                     atoms[(chain.id, resseq, atom.get_name())] = atom
#     return atoms

# def run_boltz_and_compute_rmsd(main_dir, output_csv):
#     parser = MMCIFParser(QUIET=True)
#     results = []

#     for sub in sorted(os.listdir(main_dir)):
#         subdir = os.path.join(main_dir, sub)
#         if not os.path.isdir(subdir):
#             continue

#         # 1) find the polymer_only YAML and set base = e.g. "8ruh_polymer_only"
#         base = None
#         for fn in os.listdir(subdir):
#             if fn.endswith("_polymer_only.yaml"):
#                 base, _ = os.path.splitext(fn)
#                 yaml_path = os.path.join(subdir, fn)
#                 break
#         if base is None:
#             continue

#         print(f"\n[{sub}] processing {base}:")
#         # 2) run boltz predict
#         boltz_out = os.path.join(subdir, f"boltz_results_{base}")
#         cmd = [
#             BOLTZ_CMD, "predict", yaml_path,
#             "--use_msa_server", "--use_potentials"
#         ]
#         subprocess.run(cmd, cwd=subdir, check=True)

#         # 3) locate predicted model_0.cif
#         pred_cif = os.path.join(
#             boltz_out,
#             "predictions",
#             base,
#             f"{base}_model_0.cif"
#         )
#         if not os.path.exists(pred_cif):
#             print(f"  ⚠️  predicted CIF not found at {pred_cif}")
#             continue

#         # 4) download reference mmCIF
#         ref_cif = download_mmCIF(base.replace("_polymer_only",""), subdir)

#         # 5) parse & extract polymer atoms
#         struct_ref  = parser.get_structure("ref",  ref_cif)
#         struct_pred = parser.get_structure("pred", pred_cif)
#         atoms_ref  = extract_polymer_atoms(struct_ref)
#         atoms_pred = extract_polymer_atoms(struct_pred)

#         # 6) match keys and compute RMSD
#         common = sorted(set(atoms_ref) & set(atoms_pred))
#         if not common:
#             print("  ⚠️  no overlapping atoms—skipping")
#             continue

#         ref_list  = [atoms_ref[k] for k in common]
#         pred_list = [atoms_pred[k] for k in common]
#         sup = Superimposer()
#         sup.set_atoms(ref_list, pred_list)
#         rms = sup.rms
#         print(f"  ✅ RMSD = {rms:.3f} Å over {len(common)} atoms")

#         results.append({"pdb_id": base.replace("_polymer_only",""), "rmsd": f"{rms:.3f}"})

#     # write out CSV
#     with open(output_csv, "w", newline="") as csvf:
#         writer = csv.DictWriter(csvf, fieldnames=["pdb_id","rmsd"])
#         writer.writeheader()
#         writer.writerows(results)
#     print(f"\nWrote RMSD results to {output_csv}")

# if __name__ == "__main__":
#     if len(sys.argv) != 3:
#         print("Usage: run_boltz_and_rmsd.py <main_folder> <output_csv>")
#         sys.exit(1)
#     run_boltz_and_compute_rmsd(sys.argv[1], sys.argv[2])

#!/usr/bin/env python3
import os
import sys
import csv
import subprocess
import requests

from Bio.PDB import MMCIFParser, Superimposer
from collections import defaultdict

# Try to import the CEAligner for CE‐based superposition
try:
    from Bio.PDB.CEAlign import CEAligner
    HAS_CE = True
except ImportError:
    HAS_CE = False

BOLTZ_CMD = "boltz"

def download_mmCIF(pdb_id, dest_dir):
    pdb = pdb_id.upper()
    fn = os.path.join(dest_dir, f"{pdb}.cif")
    if not os.path.exists(fn):
        url = f"https://files.rcsb.org/download/{pdb}.cif"
        print(f"  ↳ downloading reference {pdb}.cif …")
        r = requests.get(url); r.raise_for_status()
        with open(fn, "wb") as f: f.write(r.content)
    return fn

def extract_polymer_atoms(structure):
    atoms = {}
    model = next(structure.get_models())
    for chain in model:
        for res in chain:
            het, resseq, icode = res.get_id()
            if het == " ":
                for atom in res:
                    atoms[(chain.id, resseq, atom.get_name())] = atom
    return atoms

def run_boltz_and_compute_rmsd(main_dir, output_csv):
    parser = MMCIFParser(QUIET=True)
    results = []

    for sub in sorted(os.listdir(main_dir)):
        subdir = os.path.join(main_dir, sub)
        if not os.path.isdir(subdir):
            continue

        # 1) locate the polymer‐only YAML and derive its base name
        base = None
        for fn in os.listdir(subdir):
            if fn.endswith("_polymer_only.yaml"):
                base, _ = os.path.splitext(fn)
                yaml_path = os.path.join(subdir, fn)
                break
        if base is None:
            continue

        print(f"\n[{sub}] processing {base}:")
        # 2) run Boltz
        boltz_out = os.path.join(subdir, f"boltz_results_{base}")
        cmd = [
            BOLTZ_CMD, "predict", yaml_path,
            "--use_msa_server", "--use_potentials"
        ]
        subprocess.run(cmd, cwd=subdir, check=True)

        # 3) find the predicted CIF
        pred_cif = os.path.join(
            boltz_out,
            "predictions",
            base,
            f"{base}_model_0.cif"
        )
        if not os.path.exists(pred_cif):
            print(f"  ⚠️  predicted CIF not found at {pred_cif}")
            continue

        # 4) download reference mmCIF
        ref_id = base.replace("_polymer_only", "")
        ref_cif = download_mmCIF(ref_id, subdir)

        # 5) parse both structures
        struct_ref  = parser.get_structure("ref",  ref_cif)
        struct_pred = parser.get_structure("pred", pred_cif)

        atoms_ref  = extract_polymer_atoms(struct_ref)
        atoms_pred = extract_polymer_atoms(struct_pred)

        # 6) sequence‐based superposition (align)
        common = sorted(set(atoms_ref) & set(atoms_pred))
        align_rms = None
        if common:
            ref_list  = [atoms_ref[k] for k in common]
            pred_list = [atoms_pred[k] for k in common]
            sup = Superimposer()
            sup.set_atoms(ref_list, pred_list)
            align_rms = sup.rms
            print(f"  align RMSD = {align_rms:.3f} Å over {len(common)} atoms")
        else:
            print("  ⚠️  no common atoms for align, skipping")

        # 7) CE‐based superposition (cealign), if available
        ce_rms = None
        ce_count = None
        if HAS_CE:
            ce = CEAligner()
            ce.set_reference(struct_ref)
            ce.align(struct_pred)
            ce_rms = ce.rms
            # CEAligner stores its alignment as a list of residue‐pair fragments:
            # we can sum their lengths to get atom count approximation
            ce_count = sum(len(frag) for frag in ce.aln) if hasattr(ce, "aln") else None
            print(f"  CEalign RMSD = {ce_rms:.3f} Å over {ce_count or '?'} residues")
        else:
            print("  ℹ️  CEAligner not available—install Biopython ≥1.79 for cealign")

        # 8) collect results
        results.append({
            "pdb_id":    ref_id,
            "align_rmsd": f"{align_rms:.3f}" if align_rms is not None else "",
            "align_atoms": len(common),
            "ce_rmsd":   f"{ce_rms:.3f}" if ce_rms is not None else "",
            "ce_atoms":  ce_count or ""
        })

    # 9) write out CSV
    with open(output_csv, "w", newline="") as csvf:
        fieldnames = ["pdb_id","align_rmsd","align_atoms","ce_rmsd","ce_atoms"]
        writer = csv.DictWriter(csvf, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(results)
    print(f"\nWrote RMSD results to {output_csv}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: run_boltz_and_rmsd.py <main_folder> <output_csv>")
        sys.exit(1)
    run_boltz_and_compute_rmsd(sys.argv[1], sys.argv[2])

