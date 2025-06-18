#!/usr/bin/env python3
import sys
import os
import argparse
import yaml
import requests
from Bio.PDB.MMCIFParser import MMCIFParser
from collections import defaultdict

# === 3-letter codes for polymer residues ===
RNA_CODES = { "A":"A","G":"G","C":"C","U":"U",
              "ADE":"A","GUA":"G","CYT":"C","URA":"U" }
DNA_CODES = { "DA":"A","DG":"G","DC":"C","DT":"T",
              "DAD":"A","GUA":"G","CYT":"C","THY":"T" }
AA_CODES  = { "ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C",
              "GLN":"Q","GLU":"E","GLY":"G","HIS":"H","ILE":"I",
              "LEU":"L","LYS":"K","MET":"M","PHE":"F","PRO":"P",
              "SER":"S","THR":"T","TRP":"W","TYR":"Y","VAL":"V" }
STANDARD_RESIDUES = set(RNA_CODES)|set(DNA_CODES)|set(AA_CODES)

# === polyatomic species to ignore ===
IGNORED_POLYATOMIC = {
    "1MG","2MA","2MG","3TD","4D4","4OC","4SU",
    "5MC","5MU","6MZ","D2T","G7M","H2U","IAS",
    "MA6","MEQ","MS6","OMC","OMG","OMU","PSU",
    "PUT","RSP","SPD","SPM","T6A","UR3","SO4","PO4","NH4","M2G","SF4","0TD"
}


def download_cif(pdb_id, dest_dir):
    pdb = pdb_id.upper()
    url = f'https://files.rcsb.org/download/{pdb}.cif'
    resp = requests.get(url); resp.raise_for_status()
    fn = os.path.join(dest_dir, f"{pdb}.cif")
    with open(fn,"wb") as f: f.write(resp.content)
    return fn

def extract_sequences(cif_path):
    parser = MMCIFParser(QUIET=True)
    st = parser.get_structure("s", cif_path)
    model = next(st.get_models())
    chains = defaultdict(list)
    for ch in model:
        for res in ch:
            if res.get_id()[0]==" ":
                chains[ch.id].append(res)
    seqs={}
    for cid, residues in chains.items():
        seq=""; mol="rna"
        for r in residues:
            name=r.get_resname().strip()
            if name in RNA_CODES:
                seq+=RNA_CODES[name]; mol="rna"
            elif name in DNA_CODES:
                seq+=DNA_CODES[name]; mol="dna"
            elif name in AA_CODES:
                seq+=AA_CODES[name]; mol="protein"
            else:
                seq+="X"
        seqs[cid]=(mol,seq)
    return seqs

def extract_ligands(cif_path):
    parser = MMCIFParser(QUIET=True)
    st = parser.get_structure("s", cif_path)
    model = next(st.get_models())
    ligs=[]
    for ch in model:
        for res in ch:
            if res.get_id()[0]!=" ":
                code=res.get_resname().strip()
                if code in STANDARD_RESIDUES: continue
                if len(list(res.get_atoms()))==1: continue
                if code in IGNORED_POLYATOMIC: continue
                ligs.append(code)
    return sorted(ligs)

def write_yamls(seqs, ligs, dest_dir, base):
    # polymer-only
    docs=[]; idx=ord("A")
    for _,(mol,seq) in seqs.items():
        docs.append({mol:{"id":chr(idx),"sequence":seq}})
        idx+=1
    fn1=os.path.join(dest_dir,f"{base}_polymer_only.yaml")
    with open(fn1,"w") as f: yaml.dump({"sequences":docs},f,sort_keys=False)
    # polymer+ligand
    fn2=None
    if ligs:
        docs2=docs.copy()
        for l in ligs:
            docs2.append({"ligand":{"id":chr(idx),"ccd":l}})
            idx+=1
        fn2=os.path.join(dest_dir,f"{base}_polymer_ligand.yaml")
        with open(fn2,"w") as f: yaml.dump({"sequences":docs2},f,sort_keys=False)
    return fn1,fn2

def main():
    p = argparse.ArgumentParser(description="Batch PDB→YAML (one folder per PDB).")
    p.add_argument("list_file", help="text file with one PDB ID per line")
    p.add_argument("out_dir",   help="output directory for all subfolders")
    args=p.parse_args()

    with open(args.list_file) as f:
        pdbs=[L.strip() for L in f if L.strip()]

    os.makedirs(args.out_dir, exist_ok=True)
    for pdb in pdbs:
        sub = os.path.join(args.out_dir, pdb.lower())
        os.makedirs(sub, exist_ok=True)
        print(f"Processing {pdb} → {sub}")
        cif = download_cif(pdb, sub)
        seqs = extract_sequences(cif)
        ligs = extract_ligands(cif)
        print(f"  chains={len(seqs)} ligs={len(ligs)} ({', '.join(ligs) or 'none'})")
        base = pdb.lower()
        y1,y2 = write_yamls(seqs, ligs, sub, base)
        print(f"  wrote: {os.path.basename(y1)}", end="")
        if y2: print(f", {os.path.basename(y2)}")
        else:  print("")

if __name__=="__main__":
    main()
