#!/usr/bin/env python3
import os, sys, csv, subprocess, requests
from Bio.PDB import MMCIFParser, Superimposer
from collections import defaultdict

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

def run_cealign_rmsd(ref_cif, pred_cif):
    """
    Run PyMOL's cealign via subprocess and return the RMSD.
    """
    cmd = [
        "pymol", "-c", "--",
        "-r", ref_cif,
        "-p", pred_cif,
        "-x", """
        fetch %s, type=cif, async=0
        fetch %s, type=cif, async=0
        cealign %s, %s
        """ % (ref_cif, pred_cif, pred_cif, ref_cif)
    ]
    
    # Replace with actual PDB IDs or paths if needed
    cmd[4] = ref_cif.split("/")[-1].replace(".cif", "")
    cmd[6] = pred_cif.split("/")[-1].replace(".cif", "")

    try:
        result = subprocess.run(
            ["pymol", "-c"],
            input="""
from pymol import cmd
cmd.load("%s", "ref")
cmd.load("%s", "pred")
cmd.cealign("pred", "ref")
            """.strip() % (ref_cif, pred_cif),
            text=True,
            capture_output=True,
            check=True
        )
        
        # Parse RMSD from output
        for line in result.stdout.splitlines():
            if "RMSD" in line:
                rmsd_line = line.strip()
                rmsd = float(rmsd_line.split()[1])
                return rmsd

        print("⚠️ CEalign RMSD not found in output.")
        return None

    except Exception as e:
        print("⚠️ Error running CEalign:", str(e))
        return None

def run_boltz_and_compute_rmsd(main_dir, output_csv):
    parser = MMCIFParser(QUIET=True)
    results = []

    for sub in sorted(os.listdir(main_dir)):
        subdir = os.path.join(main_dir, sub)
        if not os.path.isdir(subdir):
            continue

        base = None
        for fn in os.listdir(subdir):
            if fn.endswith("_polymer_only.yaml"):
                base, _ = os.path.splitext(fn)
                yaml_path = os.path.join(subdir, fn)
                break
        if base is None:
            continue

        print(f"\n[{sub}] processing {base}:")
        boltz_out = os.path.join(subdir, f"boltz_results_{base}")
        cmd = [BOLTZ_CMD, "predict", yaml_path, "--use_msa_server", "--use_potentials"]
        subprocess.run(cmd, cwd=subdir, check=True)

        pred_cif = os.path.join(boltz_out, "predictions", base, f"{base}_model_0.cif")
        if not os.path.exists(pred_cif):
            print(f"  ⚠️  predicted CIF not found at {pred_cif}")
            continue

        ref_cif = download_mmCIF(base.replace("_polymer_only",""), subdir)

        # BioPython Superimposer RMSD
        struct_ref = parser.get_structure("ref", ref_cif)
        struct_pred = parser.get_structure("pred", pred_cif)
        atoms_ref = extract_polymer_atoms(struct_ref)
        atoms_pred = extract_polymer_atoms(struct_pred)
        common = sorted(set(atoms_ref) & set(atoms_pred))

        if not common:
            print("  ⚠️  no overlapping atoms—skipping")
            continue

        ref_list = [atoms_ref[k] for k in common]
        pred_list = [atoms_pred[k] for k in common]
        sup = Superimposer()
        sup.set_atoms(ref_list, pred_list)
        biopython_rms = sup.rms

        # CEAlign RMSD
        cealign_rms = run_cealign_rmsd(ref_cif, pred_cif)

        print(f"  ✅ BioPython RMSD = {biopython_rms:.3f} Å over {len(common)} atoms")
        if cealign_rms is not None:
            print(f"  ✅ CEalign RMSD    = {cealign_rms:.3f} Å")

        results.append({
            "pdb_id": base.replace("_polymer_only",""),
            "rmsd_biopython": f"{biopython_rms:.3f}",
            "rmsd_cealign": f"{cealign_rms:.3f}" if cealign_rms else "ERROR"
        })

    # Write out CSV
    with open(output_csv, "w", newline="") as csvf:
        writer = csv.DictWriter(csvf, fieldnames=["pdb_id", "rmsd_biopython", "rmsd_cealign"])
        writer.writeheader()
        writer.writerows(results)
    print(f"\nWrote RMSD results to {output_csv}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: run_boltz_and_rmsd.py <main_folder> <output_csv>")
        sys.exit(1)
    run_boltz_and_compute_rmsd(sys.argv[1], sys.argv[2])