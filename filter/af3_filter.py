import argparse
import warnings
import shutil
import pandas as pd
from pathlib import Path
from tqdm import tqdm
from utils.filter_utils import get_af3_confidence, check_filter
from utils.biotite_utils import extract_entity_sequences
import json
import multiprocessing as mp

warnings.filterwarnings("ignore", message="Attribute 'auth_.*' not found")

def af3_filter(filter_intput):
    summary_json, binder_chain, target_chain, cif, filter_json = filter_intput
    
    af3_confidence = get_af3_confidence(summary_json, binder_chain)
    
    filter = json.load(open(filter_json, "r"))
    chains = extract_entity_sequences(cif)
    
    sample_score = af3_confidence
    filter_passed, fail = check_filter(sample_score, filter["AF3_filter"])
    fail.update({"af3_pass": filter_passed})

    sample_info = {
        "sample_id": str(summary_json.stem).replace("_summary_confidences", ""),
        "summary_json": summary_json,
        "target_seq": chains[target_chain]['seq'], 
        "binder_seq": chains[binder_chain]['seq'], 
        "af3_pass": filter_passed
    }
    
    return {**sample_score, **sample_info, **fail}

if __name__ == "__main__":

    args = argparse.ArgumentParser()
    args.add_argument("--src_dir", type=str, 
                      default="./examples/prot_binder/af3_out", 
                      help="the af3 output dir")
    args.add_argument("--out_dir", type=str,
                      default="./examples/prot_binder/filter_out", 
                      help="the af3 filter output dir")
    args.add_argument("--filter_json", type=str,
                      default="./filter/score_json/miniprotein_filter.json", 
                      help="score json")
    args.add_argument("--binder_chain", type=str,
                      default="A", help="binder chain id")
    args.add_argument("--target_chain", default="B")
    args.add_argument("--exp_name", default="prot_binder")
    args = args.parse_args()
    
    print(f"AF3 filter args: {args}")
    
    src_dir = Path(args.src_dir)
    outdir = Path(args.out_dir) / "af3_filter"
    outdir.mkdir(parents=True, exist_ok=True)
    
    assert src_dir.is_dir() and len(list(src_dir.glob("*"))) > 0, f"no files under {src_dir}"
    
    # list input    
    summary_json_list = list(src_dir.rglob("*_summary_confidences.json"))
    inputs = []
    for summary_json in summary_json_list:
        cif = Path(str(summary_json).replace("summary_confidences.json", "model.cif"))
        assert cif.is_file(), print(f"fail to get af3 cif using summary json {summary_json}")
        inputs.append((summary_json, args.binder_chain, args.target_chain, cif, args.filter_json))
    
    assert len(inputs) > 0, f"fail to process af3-like outputs in {src_dir}"
    
    # af3 filter
    num_cores = mp.cpu_count()
    num_processes = max(64, int(num_cores * 0.8))
    print(f"af3 output number: {len(inputs)}")
    with mp.Pool(processes=num_processes) as pool:
        all_results = [
            r for r in tqdm(pool.imap(af3_filter, inputs), total=len(inputs), desc="filtering samples using af3 confidence") 
            if r is not None
        ]
    print(f"Successful process files: {len(all_results)}")
    
    pass_results = [r for r in all_results if r['af3_pass']]
    fail_results = [r for r in all_results if not r['af3_pass']]
    pd.DataFrame(fail_results).to_csv(outdir / Path(f"{args.exp_name}_af3_fail.csv"), index=False)

    print(f"af3 filer: pass {len(pass_results)}, fail {len(fail_results)}")
    
    if len(pass_results) > 0:
        out_df = pd.DataFrame(pass_results)
        
        print("coping filtered files...")
        dst_dir = outdir / Path(f"{args.exp_name}_af3_filtered")
        summary_json_list = []
        for i, summary_json in enumerate(out_df["summary_json"]):
            summary_json = Path(summary_json)
            sample_af3_out = summary_json.parent
            assert sample_af3_out.is_dir(), print(f"File not found: {sample_af3_out}")
            assert not (dst_dir / sample_af3_out.name).is_dir(), print(f"{dst_dir / sample_af3_out.name} already exists")
            
            shutil.copytree(sample_af3_out, dst_dir / sample_af3_out.name)
            summary_json_list.append(dst_dir / sample_af3_out.name / summary_json.name )   

        out_df["summary_json"] = summary_json_list
        out_df.to_csv(outdir / Path(f"{args.exp_name}_af3_pass.csv"), index=False)
        
        print(f"finish af3 filter process")
            