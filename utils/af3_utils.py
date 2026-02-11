import json
import os.path as osp
import argparse, json, sys
from pathlib import Path
from tqdm import tqdm
from typing import Tuple, Optional

def read_json(filename):
    """Read a JSON file and return the data."""
    with open(filename, 'r') as f:
        data = json.load(f)
    return data

def write_txt(data, filename):
    """Write data to a text file."""
    with open(filename, 'w') as f:
        f.write(data)

def extract_msa_from_json(data_json, msa_dir=None):
    if msa_dir is None:
        msa_dir = osp.dirname(data_json)
    data_prepare = read_json(data_json)
    
    msa_files = {}
    for sequence in data_prepare['sequences']:
        entity = sequence.get('protein', None)
        if entity is not None:
            unpaired_msa = entity['unpairedMsa']
            paired_msa = entity['pairedMsa']
            seq_id = entity['id']
            seq = entity['sequence']
        msa_file = osp.join(msa_dir, f'{data_json.name[:-len("_data.json")]}_{seq_id}_unpaired_msa.a3m')
        pairedmsa_file = osp.join(msa_dir, f'{data_json.name[:-len("_data.json")]}_{seq_id}_paired_msa.a3m')
        write_txt(unpaired_msa, msa_file)
        write_txt(paired_msa, pairedmsa_file)
        print(f'Extracted msa for {seq_id} from {data_json}')
        
        msa_files.update({
            seq: {
            'unpairedMsaPath': msa_file,
            "pairedMsaPath": pairedmsa_file
            }
        })
    
    msa_json = osp.join(msa_dir, "msa.json")
    with open(msa_json, "w", encoding="utf-8") as f:
        json.dump(msa_files, f, ensure_ascii=False, indent=2)

def extract_template_from_json(data_json, template_dir=None):
    if template_dir is None:
        template_dir = osp.dirname(data_json)
    data_prepare = read_json(data_json)
    templates_all = {}
    for sequence in data_prepare['sequences']:
        entity = sequence.get('protein', None)
        
        template_info = []
        if entity is not None:
            templates = entity['templates']
            seq_id = entity['id']
            for template_id, template in enumerate(templates):
                template_file = osp.join(template_dir, f'{data_json.name[:-len("_data.json")]}_{seq_id}_{template_id}.cif')
                
                write_txt(template['mmcif'], template_file)
        
                template_info.append({
                    "mmcifPath": template_file,
                    'queryIndices': template['queryIndices'],
                    'templateIndices': template['templateIndices']
                    })

                print(f'Extracted template {template_file} from {data_json}')
        templates_all.update(
            {
                entity['sequence']: template_info
            }
        )
                
    queries_file = osp.join(template_dir, f'template_queries.json')
    with open(queries_file, "w", encoding="utf-8") as f:
        json.dump(templates_all, f, ensure_ascii=False, indent=2) 
        
if __name__ == "__main__":    
    
    pa = argparse.ArgumentParser()
    pa.add_argument("--root", required=True, help="Root directory containing AF3 outputs")
    pa.add_argument("--presearch_outdir", required=True, help="Dir to store msa and templates")
    pa.add_argument("--recursive", action="store_true", help="Recurse into subdirectories")
    args = pa.parse_args()
    
    Path(args.presearch_outdir).mkdir(parents=True, exist_ok=True)
    
    root = Path(args.root)
    if not root.exists():
        print(f"[ERR] root not found: {root}", file=sys.stderr); sys.exit(2)

    dirs = []
    if args.recursive:
        dirs = [p for p in root.rglob("*") if p.is_dir()]
    else:
        dirs = [root] + [p for p in root.iterdir() if p.is_dir()]

    def find_data_json(dirpath: Path) -> Optional[Tuple[Path, Path, Optional[Path]]]:
        data_json_candidates = list(dirpath.glob("*_data.json"))
        if not data_json_candidates:
            return None
        
        if data_json_candidates[0].exists():
            return data_json_candidates[0]

    for d in tqdm(sorted(dirs)):
        data_json = find_data_json(d)
        if not data_json:
            continue

        extract_msa_from_json(data_json=data_json, msa_dir=args.presearch_outdir)
        extract_template_from_json(data_json=data_json, template_dir=args.presearch_outdir)
    
    index_json = {}

    msa_json, template = Path(args.presearch_outdir) / f"msa.json", Path(args.presearch_outdir) / f"template_queries.json"
    if msa_json is None or template is None:
        raise f"fail to find msa and template under {args.presearch_outdir}"
    
    msa_json = json.load(open(msa_json, 'r'))
    templates = json.load(open(template, 'r'))
    
    for seq in msa_json.keys():
        presearch_dict = {
            seq:
                {
                    "unpairedMsaPath": msa_json[seq]['unpairedMsaPath'],
                    "pairedMsaPath": msa_json[seq]['pairedMsaPath'],
                    "templates": templates[seq]
                }
        }
    
        index_json.update(presearch_dict)
    index_json_path = osp.join(args.presearch_outdir, "index.json")
    with open(index_json_path, "w", encoding="utf-8") as f:
        json.dump(index_json, f, ensure_ascii=False, indent=2)
        
    
    
    
