#!/usr/bin/env python3
import sys
import yaml
import requests
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from collections import defaultdict

# === 3-letter codes ===
RNA_CODES = {
    "A": "A", "G": "G", "C": "C", "U": "U",
    "ADE": "A", "GUA": "G", "CYT": "C", "URA": "U"
}
DNA_CODES = {
    "DA": "A", "DG": "G", "DC": "C", "DT": "T",
    "DAD": "A", "GUA": "G", "CYT": "C", "THY": "T"
}
AA_CODES = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V"
}
STANDARD_RESIDUES = set(RNA_CODES) | set(DNA_CODES) | set(AA_CODES)

# === polyatomic species to ignore ===
IGNORED_POLYATOMIC = {"SO4", "PO4", "NH4"}

def download_cif(pdb_id, out_path=None):
    pdb = pdb_id.upper()
    url = f'https://files.rcsb.org/download/{pdb}.cif'
    print(f"🌐 Downloading {pdb}.cif from {url}")
    resp = requests.get(url)
    resp.raise_for_status()
    fn = out_path or f"{pdb}.cif"
    with open(fn, 'wb') as f:
        f.write(resp.content)
    print(f"✅ Saved mmCIF to {fn}")
    return fn

def extract_sequences_from_cif(cif_path):
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure("structure", cif_path)
    model = next(structure.get_models())
    residues_by_chain = defaultdict(list)

    for chain in model:
        for res in chain:
            het, _, _ = res.get_id()
            if het == " ":
                residues_by_chain[chain.id].append(res)

    sequences = {}
    for chain_id, residues in residues_by_chain.items():
        seq = ""
        moltype = "rna"
        for res in residues:
            r = res.get_resname().strip()
            if r in RNA_CODES:
                seq += RNA_CODES[r]; moltype = "rna"
            elif r in DNA_CODES:
                seq += DNA_CODES[r]; moltype = "dna"
            elif r in AA_CODES:
                seq += AA_CODES[r]; moltype = "protein"
            else:
                seq += "X"
        sequences[chain_id] = (moltype, seq)
    return sequences

def extract_ligands_from_cif(cif_path):
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure("structure", cif_path)
    model = next(structure.get_models())

    ligands = []
    for chain in model:
        for res in chain:
            het, resseq, icode = res.get_id()
            if het != " ":
                mon_id = res.get_resname().strip()
                # skip standard residues
                if mon_id in STANDARD_RESIDUES:
                    continue
                # skip monatomic "ligands" (waters, metal ions, halides...)
                if len(list(res.get_atoms())) == 1:
                    continue
                # skip common polyatomic buffers/ions
                if mon_id in IGNORED_POLYATOMIC:
                    continue
                ligands.append({
                    "ccd": mon_id,
                    "chain": chain.id,
                    "resseq": resseq,
                    "icode": icode
                })
    return ligands

def write_yaml(sequences, ligands, filename, include_ligands=True):
    docs = []
    counter = ord("A")

    # polymers first
    for chain_id, (moltype, seq) in sequences.items():
        docs.append({moltype: {"id": chr(counter), "sequence": seq}})
        counter += 1

    # then each ligand instance
    if include_ligands:
        for lig in ligands:
            docs.append({"ligand": {"id": chr(counter), "ccd": lig["ccd"]}})
            counter += 1

    with open(filename, "w") as f:
        yaml.dump({"sequences": docs}, f, sort_keys=False)
    print(f"✅ Wrote YAML to {filename}")

def main(cif_file):
    print(f"📂 Reading {cif_file}")
    sequences = extract_sequences_from_cif(cif_file)
    ligands = extract_ligands_from_cif(cif_file)
    print(f"🧬 Found {len(sequences)} polymer chains")
    print(f"🔗 Found {len(ligands)} small-molecule ligands: "
          f"{', '.join(l['ccd'] for l in ligands) or 'none'}")

    base = cif_file.rsplit(".", 1)[0]
    write_yaml(sequences, ligands, filename=f"{base}_polymer_only.yaml", include_ligands=False)
    if ligands:
        write_yaml(sequences, ligands, filename=f"{base}_polymer_ligand.yaml", include_ligands=True)
    else:
        print("⚠️ No small-molecule ligands found, skipping ligand YAML")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python cif_to_yaml.py <PDB_ID>")
        sys.exit(1)
    pdb_id = sys.argv[1]
    cif_path = download_cif(pdb_id)
    main(cif_path)
