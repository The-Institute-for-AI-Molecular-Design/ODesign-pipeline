import json
import argparse
import numpy as np
import multiprocessing as mp

from pathlib import Path
from tqdm import tqdm
from biotite.structure.io.pdb import PDBFile, get_structure
import biotite.structure.io.pdbx as pdbx

from utils.biotite_utils import cif2pdb
     
def get_redesign_residue(
    file: str,
    redesign_chain_id_list: list
):
    # read file
    if file.endswith('.pdb'):
        pdb_file = PDBFile.read(file)
        atom_array = get_structure(pdb_file, model=1)
    elif file.endswith('.cif'):
        cif_file = pdbx.PDBxFile.read(file)
        atom_array = pdbx.get_structure(cif_file, model=1)
    else:
        raise NotImplementedError
    
    # get fixed chain for ligandmpnn
    fixed_residues = []
    for redesign_chain_id in redesign_chain_id_list:
        for res_id in np.unique(atom_array.res_id[atom_array.chain_id == redesign_chain_id]):
            icode = atom_array.ins_code[(atom_array.res_id == res_id) & (atom_array.chain_id == redesign_chain_id)][0]
            fixed_residues.append(f'{redesign_chain_id}{res_id}{icode}')
        fixed_residues = " ".join(map(str, fixed_residues))
    return fixed_residues

   
def json_maker(
        input_path: str,
        output_path: str,
        redesign_chain_id: str,
):
    path = Path(input_path)
    all_pdb_files = list(path.rglob('*.pdb'))
    all_cif_files = list(path.rglob('*.cif'))
    all_files = all_pdb_files + all_cif_files
    all_files = [str(f) for f in all_files if f.is_file()]

    json_dict = {}
    redesign_chain_id = [c for c in redesign_chain_id.split(",")]
    for f in all_files:
        json_dict[f] = get_redesign_residue(f, redesign_chain_id_list=redesign_chain_id)
    
    with open(output_path, 'w') as f:
        json.dump(json_dict, f, indent=2)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--input_path',
        type=str,
        required=True
    )
    parser.add_argument(
        '--output_path',
        type=str,
        help='path to save ligandmpnn json file',
        required=True
    )
    parser.add_argument(
        '--redesign_chain_id',
        type=str,
        default="A"
    )
    args = parser.parse_args()
    
    print(f"prepare ligandmpnn input in args: {args}")
    
    # convert backbone cif files to pdb files
    tmp_dir = Path(args.output_path) / "bb_pdb"
    tmp_dir.mkdir(exist_ok=True)
    
    cif2pdb_inputs = [(f, tmp_dir) for f in Path(args.input_path).rglob('*.cif')]
    num_cores = mp.cpu_count()
    num_processes = max(64, int(num_cores * 0.8))
    with mp.Pool(processes=num_processes) as pool:
        all_results = [
            r for r in tqdm(pool.imap(cif2pdb, cif2pdb_inputs), total=len(cif2pdb_inputs), desc="converting cif files to pdb files") 
            if r is not None
        ]
    print(f"Successful process files: {len(all_results)}")
    
    json_maker(
        input_path=tmp_dir,
        output_path=Path(args.output_path) / "mpnn_input.json",
        redesign_chain_id=args.redesign_chain_id,
    )