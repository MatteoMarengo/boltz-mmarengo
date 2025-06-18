import requests
import sys
import yaml

def fetch_entry(pdb_id):
    url = f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
    resp = requests.get(url)
    if resp.status_code != 200:
        raise Exception(f"Could not fetch data for {pdb_id}")
    return resp.json()

def parse_from_entry(entry):
    polymer_chains = {}
    ligands = []

    # Parse polymers
    polymers = entry.get("polymer_entities", [])
    for i, entity in enumerate(polymers):
        mol_type = entity["entity_poly"]["type"].lower()
        seq = entity["entity_poly"].get("pdbx_seq_one_letter_code_can")
        if not seq:
            continue
        chains = entity.get("rcsb_polymer_entity_container_identifiers", {}).get("auth_asym_ids", [])
        for chain in chains:
            if "protein" in mol_type:
                type_ = "protein"
            elif "rna" in mol_type:
                type_ = "rna"
            elif "dna" in mol_type:
                type_ = "dna"
            else:
                continue
            polymer_chains[chain] = {"type": type_, "sequence": seq}

    # Parse ligands (non-polymer entities)
    nonpolys = entry.get("nonpolymer_entities", [])
    for entity in nonpolys:
        lig_id = entity.get("nonpolymer_comp", {}).get("chem_comp", {}).get("id")
        if lig_id:
            ligands.append(lig_id)

    return polymer_chains, list(set(ligands))

def generate_yaml(polymer_chains, ligands=None, filename="output.yaml", include_ligands=True):
    if ligands is None:
        ligands = []

    sequences = []
    chain_counter = ord('A')

    for chain_id, info in polymer_chains.items():
        sequences.append({info['type']: {
            'id': chr(chain_counter),
            'sequence': info['sequence']
        }})
        chain_counter += 1

    if include_ligands:
        for lig in ligands:
            sequences.append({'ligand': {
                'id': chr(chain_counter),
                'ccd': lig
            }})
            chain_counter += 1

    with open(filename, 'w') as f:
        yaml.dump({'sequences': sequences}, f, sort_keys=False)

    print(f"✅ Wrote YAML to {filename}")

def main(pdb_id):
    print(f"🔍 Fetching metadata for PDB ID: {pdb_id}")
    entry = fetch_entry(pdb_id)

    print("📦 Parsing polymer chains and ligands...")
    polymer_chains, ligands = parse_from_entry(entry)

    print(f"✅ Found {len(polymer_chains)} polymer chains.")
    print(f"✅ Found {len(ligands)} ligands: {', '.join(ligands) if ligands else 'None'}")

    generate_yaml(polymer_chains, filename=f"{pdb_id}_polymer_only.yaml", include_ligands=False)

    if ligands:
        generate_yaml(polymer_chains, ligands=ligands, filename=f"{pdb_id}_polymer_ligand.yaml", include_ligands=True)
    else:
        print("⚠️ No ligands found. Skipping ligand-containing YAML.")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python produce_yaml_files.py <PDB_ID>")
        sys.exit(1)

    pdb_id = sys.argv[1].strip().lower()
    main(pdb_id)
