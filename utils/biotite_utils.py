from biotite.structure.io.pdb import PDBFile
import biotite.structure.io.pdb as pdb
import biotite.structure.io.pdbx as pdbx
from pathlib import Path
import numpy as np

three_to_one_protein = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D",
    "CYS": "C", "GLN": "Q", "GLU": "E", "GLY": "G",
    "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S",
    "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    "SEC": "U", "PYL": "O", "ASX": "B", "GLX": "Z",
    "XLE": "J", "UNK": "X"
}
three_to_one_dna = {
    "DA": "A", "DC": "C", "DG": "G", "DT": "T",
}
three_to_one_rna = {
    "A": "A", "C": "C", "G": "G", "U": "U"
}

three_to_one = {**three_to_one_protein, **three_to_one_dna, **three_to_one_rna}

def cif2pdb(args):
    cif_path, output_dir = args
    pdb_path = output_dir / (cif_path.stem + ".pdb")
    
    try:
        cif_file = pdbx.CIFFile.read(str(cif_path))
        structure = pdbx.get_structure(cif_file)

        pdb_file = pdb.PDBFile()
        pdb_file.set_structure(structure)
        pdb_file.write(str(pdb_path))

        return True
    except Exception as e:
        return False

def get_intervals(seen):
    nums = sorted(seen)
    intervals = []
    for n in nums:
        if not intervals or n > intervals[-1][1] + 1:
            intervals.append([n, n])
        else:
            intervals[-1][1] = n
    return [tuple(x) for x in intervals]

def extract_entity_sequences(pdb_path: Path):
    if pdb_path.suffix == '.pdb':
        pdb_file = PDBFile.read(pdb_path)
        structure = pdb_file.get_structure()[0]
    elif pdb_path.suffix == '.cif':
        cif_file = pdbx.CIFFile.read(pdb_path)
        structure = pdbx.get_structure(cif_file)[0]
    else:
        raise NotImplementedError

    chain_seqs = {}
    for chain_id in np.unique(structure.chain_id):
        chain_atoms = structure[structure.chain_id == chain_id]
        chain_atoms = chain_atoms[chain_atoms.res_name[:] != "HOH"]
        seq = ""
        seen = set()
        for res_id, res_name in zip(chain_atoms.res_id, chain_atoms.res_name):
            if res_id in seen:
                continue
            seen.add(res_id) 
            one_letter = three_to_one.get(res_name.upper(), "X")
            seq += one_letter
        
        res_intervals = get_intervals(seen)
        if res_name in three_to_one_protein.keys():
            entity_type = "protein"
        elif res_name in three_to_one_dna.keys():
            entity_type = "dna"
        elif res_name in three_to_one_rna.keys():
            entity_type = "rna"
        else:
            print(f"fail to figure out chain {chain_id} entity type with its resname {res_name}")
        
        chain_seqs[str(chain_id)] = {"entity_type": entity_type, 
                                     "intervals": [f"{chain_id}/{f[0]}-{f[1]}" for f in res_intervals], 
                                     "seq": str(seq)}
    return chain_seqs

def extract_chain(pdb_path, output_dir, chain_id):
    pdb_path = Path(pdb_path)
    if pdb_path.suffix == '.pdb':
        pdb_file = PDBFile.read(pdb_path)
        structure = pdb_file.get_structure()[0]
    elif pdb_path.suffix == '.cif':
        cif_file = pdbx.CIFFile.read(pdb_path)
        structure = pdbx.get_structure(cif_file)[0]
    else:
        raise NotImplementedError
    
    structure = structure[structure.chain_id[:] == chain_id]
    
    pdb_path = output_dir / (pdb_path.stem + ".pdb")
    pdb_file = pdb.PDBFile()
    pdb_file.set_structure(structure)
    pdb_file.write(str(pdb_path))
