import json
import random
import os
import shutil
import os.path as osp
import argparse
import subprocess
from pathlib import Path
from collections import deque
from time import sleep
from datetime import datetime

def split_json(input_json, gpu_list, exp_name, tmp_dir):
    with open(input_json, "r") as f:
        data = json.load(f)
    data = list(data.items())
    n = len(data)
    k = len(gpu_list)
    chunk_size = (n + k - 1) // k
    json_paths = []
    for i, gpu_id in enumerate(gpu_list):
        chunk = dict(data[i*chunk_size : (i+1)*chunk_size])
        if not chunk:
            continue
        chunk_path = os.path.join(tmp_dir, f"{exp_name}_part_{gpu_id}.json")
        with open(chunk_path, "w") as f:
            json.dump(chunk, f)
        json_paths.append(chunk_path)
    return json_paths

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="mpnn")
    # inputs
    parser.add_argument("--input_json", type=str, required=True, help="input json file")
    parser.add_argument("--output_dir", type=str, required=True, help="dirs to output mpnn results")
    parser.add_argument('--exp_name', type=str, help='a tag used to differentiate your json, tmp dir and files')
    parser.add_argument("--binder_chain_id", type=str, default="A")
    parser.add_argument("--batch_size", type=int, default=1, help="inverse folding times for one backbone.")
    parser.add_argument("--number_of_batches", type=int, default=1, help="Number of times to design sequence using a chosen batch size.")
    parser.add_argument("--temperature", type=str, default=0.1)
    parser.add_argument("--omit_AA", type=str, default="C")
    parser.add_argument("--mpnn_model_type", type=str, default="protein_mpnn", choices=["protein_mpnn", "ligand_mpnn"])
    parser.add_argument("--protein_mpnn_ckpts", type=str, default="./LigandMPNN/model_params/proteinmpnn_v_48_020.pt")
    parser.add_argument("--ligand_mpnn_ckpts", type=str, default="./LigandMPNN/model_params/ligandmpnn_v_32_010_25.pt")
    parser.add_argument('--gpus',type=str, help="List of GPU IDs to use, e.g. --gpus 1,2,3")
    args = parser.parse_args()
    
    args.gpus = [int(v) for v in args.gpus.split(',')]
    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    print(f"mpnn args: {args}, time_stamp: {ts}")
    
    output_dir = Path(args.output_dir)
    tmp_dir = output_dir / Path(f"{args.exp_name}_tmp_{ts}")
    tmp_json_dir = tmp_dir / Path('input_json') 
    
    output_dir.mkdir(parents=True, exist_ok=True)
    tmp_json_dir.mkdir(parents=True, exist_ok=True)
        
    #===========================================================
    # 1. split mpnn json
    #===========================================================
    
    available_gpus = list(args.gpus)
    if not available_gpus:
        raise ValueError("No GPU IDs provided via --gpus.")

    json_paths = split_json(args.input_json, available_gpus, args.exp_name, tmp_json_dir)
    json_tasks = [(json_path, output_dir) for json_path in json_paths]

    if not json_tasks:
        print("❌ No JSON tasks were created (no PDBs found after grouping/conversion).")
        raise SystemExit(1)

    #===========================================================
    # 2. distribute mpnn
    #===========================================================
    # model type
    if args.mpnn_model_type == 'protein_mpnn':
        mpnn_ckpt = args.protein_mpnn_ckpts
    elif args.mpnn_model_type == 'ligand_mpnn':
        mpnn_ckpt = args.ligand_mpnn_ckpts

    task_queue = deque(json_tasks)
    running = []

    def launch_one(json_path, out_dir, gpu_id):
        # set gpu
        env = os.environ.copy()
        env["CUDA_VISIBLE_DEVICES"] = str(gpu_id)
        cmd = [
            "python", "LigandMPNN/run.py",
            "--model_type", f"{args.mpnn_model_type}",
            "--seed", f"{random.randint(0, 2**31 - 1)}",
            "--pdb_path_multi", f"{json_path}",
            "--redesigned_residues_multi", f"{json_path}",
            "--out_folder", f"{out_dir}",
            "--temperature", f"{args.temperature}",
            "--omit_AA", f"{args.omit_AA}",
            "--parse_atoms_with_zero_occupancy", "1",
            f"--checkpoint_{args.mpnn_model_type}", f"{mpnn_ckpt}",
            "--batch_size", f"{args.batch_size}",
            "--number_of_batches", f"{args.number_of_batches}",
        ]
        desc = f"GPU{gpu_id}: {osp.basename(str(json_path))}"
        print(f"[LAUNCH] {desc}")
        proc = subprocess.Popen(cmd, env=env)
        return {'proc': proc, 'gpu': gpu_id, 'desc': desc}

    for gpu in available_gpus:
        if task_queue:
            json_path, out_dir = task_queue.popleft()
            running.append(launch_one(json_path, out_dir, gpu))

    while running:
        still_running = []
        for item in running:
            ret = item['proc'].poll()
            if ret is None:
                still_running.append(item)
                continue
            status = "OK" if ret == 0 else f"EXIT={ret}"
            print(f"[DONE] {item['desc']} -> {status}")

            if task_queue:
                json_path, out_dir = task_queue.popleft()
                still_running.append(launch_one(json_path, out_dir, item['gpu']))
        running = still_running
        if running:
            sleep(2)

    print("✅ All LigandMPNN jobs finished.")
    
    shutil.rmtree(tmp_dir)
    print("✅ Delete tmp pdb files and jsons.")

