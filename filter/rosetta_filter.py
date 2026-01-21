import argparse
import warnings
import shutil
import pandas as pd
from pathlib import Path
from tqdm import tqdm
from utils.biotite_utils import cif2pdb
import json
import multiprocessing as mp
import shlex
import traceback

warnings.filterwarnings("ignore", message="Attribute 'auth_.*' not found")

def filter_df(df: pd.DataFrame, filter_json):
    filter_dict = json.load(open(filter_json, 'r'))["rosetta_filter"]
    df = df[['decoy', 'ddg', 'sap_score', 'contact_molecular_surface']].copy()
    for filter_name, threshold in filter_dict.items():
        score = pd.to_numeric(df[filter_name], errors='coerce')

        if threshold["higher"]:
            df[f"{filter_name}_pass"] = (score > threshold["threshold"])
        else:
            df[f"{filter_name}_pass"] = (score < threshold["threshold"])

        df[f"{filter_name}_pass"] = df[f"{filter_name}_pass"].map({True:"true", False:"false"})
    
    pass_cols = [f"{name}_pass" for name in filter_dict]
    succ_df = df[df[pass_cols].eq("true").all(axis=1)].copy()
    succ_df["sample_id"] = succ_df['decoy'].str.replace('_0001', '', regex=True).replace("_model", "", regex=True)
    
    fail_df = df.loc[~df.index.isin(succ_df.index)]
    
    return succ_df, fail_df

def rosetta_score(inputs):
    xml_path, pdblist, outlist = inputs
    try:
        import pyrosetta
        from pyrosetta import rosetta

        flags = " ".join([
            "-beta_nov16",
            "-in:file:l", shlex.quote(str(pdblist)),
            "-out:file:score_only", shlex.quote(str(outlist)),
            "-out:file:scorefile_format", "json",
            "-mute", "all"
        ])
        pyrosetta.init(flags)

        with open(xml_path, "r") as fh:
            xml = fh.read()
        xml_obj = rosetta.protocols.rosetta_scripts.XmlObjects.create_from_string(xml)
        mover = xml_obj.get_mover("ParsedProtocol")

        jd2 = rosetta.protocols.jd2.JobDistributor.get_instance()
        jd2.go(mover) 
        return "success"
    except Exception as e:
        print(f"Error processing {pdblist}: {e}")
        traceback.print_exc()
        return f"Failed: {e}"

if __name__ == "__main__":

    args = argparse.ArgumentParser()
    args.add_argument("--src_dir", 
                      required=True, 
                      help="af3-like out dir")
    args.add_argument("--filter_outdir", 
                      required=True, 
                      help="rosetta filter output dir")
    args.add_argument("--filter_json", 
                      default="./filter/score_json/miniprotein_filter.json",
                      help="score json")
    args.add_argument("--binder_chain", default="A")
    args.add_argument("--target_chain", default="B")
    args.add_argument("--exp_name", 
                      default="prot_binder",
                      help="custom experiment name for tag")
    args.add_argument("--rosetta_xml", 
                      default="./filter/rosetta_cmds/ppi.xml",
                      help="rosetta xml for complex eval")
    args.add_argument("--merged", 
                      default=True, 
                      help="whether the rosetta and af3 csvs are merged or not")
    args = args.parse_args()
    
    print(f"rosetta filter args: {args}")
    
    src_dir = Path(args.src_dir)
    filter_outdir = Path(args.filter_outdir) / "rosetta_filter"
    filter_outdir.mkdir(parents=True, exist_ok=True)
    
    assert src_dir.is_dir() and len(list(src_dir.glob("*"))) > 0, f"no files under {src_dir}"
    
    cif_list = list(src_dir.glob("*/*.cif"))
    tmp_dir = filter_outdir / "tmp"
    tmp_pdb_dir = tmp_dir / "pdb"
    tmp_pdb_dir.mkdir(parents=True, exist_ok=True)
    
    # convert cif to pdb
    convert_inputs = [(f, tmp_pdb_dir) for f in cif_list]
    
    num_cores = mp.cpu_count()
    num_processes = max(64, int(num_cores * 0.8))
    print(f"{len(convert_inputs)} cifs to convert")
    with mp.Pool(processes=num_processes) as pool:
        all_results = [
            r for r in tqdm(pool.imap(cif2pdb, convert_inputs), total=len(convert_inputs), desc="converting files to pdb") 
            if r
        ]
    print(f"Successful cnvert: {len(all_results)}")
    
    # split for multiprocess
    pdblist = list(tmp_pdb_dir.glob("*.pdb"))
    
    def split(list, n_chunk):
        k = (len(list) + n_chunk - 1) // n_chunk
        for i in range(0, len(list), k):
            yield list[i:i+k]
            
    pdblist_chunk = split(pdblist, num_processes)
    tmp_pdblist_dir = tmp_dir / "pdblist"
    tmp_pdblist_dir.mkdir(parents=True, exist_ok=True)
    
    for i, chunk in enumerate(pdblist_chunk):
        tmp_pdblist = tmp_pdblist_dir / f"pdblist_part{i}"
        with open(tmp_pdblist, "w") as f:
            for pdb_path in chunk:
                f.write(f"{pdb_path.resolve()}\n")
    
    # rosetta score
    tmp_score_dir = tmp_dir / "score"
    tmp_score_dir.mkdir(parents=True, exist_ok=True)
    score_inputs = [(args.rosetta_xml, str(f), tmp_score_dir / f"score_part{i}.sc") for i, f in enumerate(tmp_pdblist_dir.glob("*"))]
    with mp.Pool(processes=num_processes, maxtasksperchild=1) as pool:
        all_results = list(tqdm(pool.imap(rosetta_score, score_inputs), total=len(score_inputs), desc="rosetta scoring..."))

    # merge score
    with open(filter_outdir / "score.sc", "w") as out:
        for tmp_score in tmp_score_dir.glob("score_part*"):
            with open(tmp_score, "r") as f:
                out.write(f.read())
    # filter
    with open(filter_outdir / "score.sc") as f:
        rows = [json.loads(line) for line in f if line.strip()]
    cols = list(rows[0].keys())
    data = pd.DataFrame(rows)[cols].drop_duplicates(subset=['decoy'], keep='first')
    filter_data, fail_df = filter_df(data, args.filter_json)
    fail_df.to_csv(f"{filter_outdir}/{args.exp_name}_rosetta_fail.csv", index=False)
    print(f"before rosetta {len(data)}, after rosetta {len(filter_data)}")
    
    # gather filtered samples
    if len(filter_data) > 0:
        rosetta_out = filter_outdir / f"{args.exp_name}_rosetta_filtered"
        rosetta_out.mkdir(parents=True, exist_ok=True)
        pdb_list = []
        for i, sample_id in enumerate(filter_data["sample_id"]):
            dst_pdb_dir = rosetta_out / f"{sample_id}.pdb"
            pdb_list.append(dst_pdb_dir)
            srcs = sorted(Path(tmp_pdb_dir).glob(f"{sample_id}*.pdb"))
            
            if len(srcs) != 1:
                raise RuntimeError(f"Expected 1 pdb, got {len(srcs)}: {[s.name for s in srcs]}")
            
            shutil.copy2(srcs[0], dst_pdb_dir)
        filter_data['pdb_dir'] = pdb_list
        filter_data.drop(columns=['decoy']).to_csv(f"{filter_outdir}/{args.exp_name}_rosetta_pass.csv", index=False)
    else:
        print("all samples fail to pass rosetta filter")  
    
    shutil.rmtree(tmp_dir)
    
    # merged af3 and rosetta or not
    if args.merged:
        af3_rosetta_df = pd.read_csv(Path(args.filter_outdir) / "af3_filter" / f"{args.exp_name}_af3_pass.csv")
        rosetta_pass_df = pd.read_csv(Path(args.filter_outdir) / "rosetta_filter" / f"{args.exp_name}_rosetta_pass.csv")
        merged_df = pd.merge(af3_rosetta_df, rosetta_pass_df, on='sample_id', how='inner')

        cols = [
            "sample_id","binder_seq",
            "iptm","binder_ptm","complex_ptm","chain_ptm_avg","ipae_min","ipae_avg",
            "ddg","sap_score","contact_molecular_surface",
            "pdb_dir","summary_json","target_seq",
        ]
        merged_df = merged_df.drop(columns=[c for c in merged_df.columns if c.endswith("_pass")])
        merged_df = merged_df.reindex(columns=cols + [c for c in merged_df.columns if c not in cols])
        merged_df.sort_values(by="ddg").drop_duplicates(subset=["binder_seq"]).to_csv(Path(args.filter_outdir) / f"{args.exp_name}_filter_final.csv", index=False)