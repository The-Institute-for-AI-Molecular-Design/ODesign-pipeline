#!/bin/bash

#######################################################################
# ODesign Inference Demo Script
# ---------------------------------------------------------------------
# Modify the arguments below, then simply run:
#     bash inference_demo.sh
#######################################################################

# 1. Choose inference model name (REQUIRED)
#     Options:
#        - odesign_base_prot_flex
#        - odesign_base_prot_rigid
#        - odesign_base_ligand_rigid
#        - odesign_base_na_rigid
infer_model_name="odesign_base_prot_flex"

# 2. Design modality
#     Specify "dna" or "rna" if using nucleic acid design model,
#     otherwise you may leave this argument empty
design_modality=""

# 3. Data root directory (default: ./ODesign/data)
data_root_dir="./ODesign/data"

# 4. Checkpoint root directory (default: ./ODesign/ckpt)
ckpt_root_dir="./ODesign/ckpt"

# 5. Path to input JSON file (REQUIRED)
input_json_path="./examples/prot_binder/prot_binder.json"

# 6. Custom experiment name for this inference run (optional)
#     If empty -> auto generate: infer_${infer_model_name}
exp_name="protein_binding_protein_design"

# 7. Seeds used for inference sampling (default: [42])
#     Format: "[42, 123, 777]" if multiple seeds
seeds="[42]"

# 8. Number of generated samples per seed (default: 5)
N_sample=20

# 9. Whether to use MSA (default: false)
#     Set true ONLY if your input JSON contains msa information
use_msa=false

# 10. Number of dataloader workers (default: 4)
num_workers=4

# 11. CUDA Device Setup (default: 0)
gpus="0"

# 12. InverseFolding Topk (default: 8)
invfold_topk=8

# 13. Output directory for ODesign inference results (default: ./examples/prot_binder/odesign_out/)
output_dir=./examples/prot_binder/odesign_out/

#!/bin/bash
while [[ $# -gt 0 ]]; do
    case $1 in
        --infer_model_name) infer_model_name="$2"; shift 2;;
        --design_modality) design_modality="$2"; shift 2;;
        --data_root_dir) data_root_dir="$2"; shift 2;;
        --ckpt_root_dir) ckpt_root_dir="$2"; shift 2;;
        --input_json_path) input_json_path="$2"; shift 2;;
        --exp_name) exp_name="$2"; shift 2;;
        --seeds) seeds="$2"; shift 2;;
        --N_sample) N_sample="$2"; shift 2;;
        --use_msa) use_msa="$2"; shift 2;;
        --num_workers) num_workers="$2"; shift 2;;
        --invfold_topk) invfold_topk="$2"; shift 2;;
        --gpus) gpus="$2"; shift 2;;
        --output_dir) output_dir="$2"; shift 2;;
        *) echo "❌ Error: Unknown option: $1"; exit 1;;
    esac
done


# auto-assign design_modality
if [[ "$infer_model_name" == *"prot"* ]]; then
    design_modality="protein"
elif [[ "$infer_model_name" == *"ligand"* ]]; then
    design_modality="ligand"
elif [[ "$infer_model_name" == *"na"* ]]; then
    if [[ -z "$design_modality" ]]; then
        echo "❌ ERROR: Please specify design_modality as 'dna' or 'rna'."
        exit 1
    fi
fi

# If exp_name is empty, create default experiment name
if [[ -z "$exp_name" ]]; then
    exp_name="infer_${infer_model_name}"
fi

# Set CUDA devices
export CUDA_VISIBLE_DEVICES=$gpus
GPU_COUNT=$(echo $gpus | awk -F',' '{print NF}')

#######################################################################
# Config Summary
#######################################################################

echo "-----------------------------------------------------------"
echo "🚀 Start ODesign Inference"
echo "-----------------------------------------------------------"
echo "Model                   : $infer_model_name"
echo "Design Modality         : $design_modality"
echo "Input JSON              : $input_json_path"
echo "Experiment Name         : $exp_name"
echo "Seeds                   : $seeds"
echo "Samples per Seed        : $N_sample"
echo "Use MSA                 : $use_msa"
echo "Dataloader Workers      : $num_workers"
echo "Data Root Dir           : $data_root_dir"
echo "Checkpoint Root Dir     : $ckpt_root_dir"
echo "CUDA_VISIBLE_DEVICES    : $CUDA_VISIBLE_DEVICES"
echo "InverseFolding Topk     : $invfold_topk"
echo "ODesign Outdir          : $output_dir"
echo "-----------------------------------------------------------"
echo ""


#######################################################################
# Single-GPU Inference Command
# ---------------------------------------------------------------------
# You don't need to modify the following command
#######################################################################

# Launch inference
python ./ODesign/scripts/inference.py \
    exp="train_${infer_model_name}" \
    data_root_dir="$data_root_dir" \
    ckpt_root_dir="$ckpt_root_dir" \
    exp.infer_model_name="$infer_model_name" \
    exp.design_modality="$design_modality" \
    exp.input_json_path="$input_json_path" \
    exp.exp_name="$exp_name" \
    exp.seeds="$seeds" \
    exp.model.sample_diffusion.N_sample="$N_sample" \
    exp.use_msa="$use_msa" \
    exp.num_workers="$num_workers" \
    exp.invfold_topk="$invfold_topk" \
    hydra.run.dir="$output_dir"


#######################################################################
# 🔥 Optional: Distributed Multi-GPU Inference Example
# ---------------------------------------------------------------------
# Uncomment and modify the following block if you want to run inference
# with multiple GPUs using torchrun.
#######################################################################

# NPROC=$GPU_COUNT             # Number of GPUs used
# NODE_RANK=0         # 0 for single-machine multi-GPU
# NNODES=1            # Number of nodes (keep 1 if only one machine)

# # Auto-generate random port
# MASTER_PORT=$((10000 + RANDOM % 50000))
# MASTER_ADDR=$(hostname -I | awk '{print $1}')

# # Launch inference
# torchrun \
#     --nproc_per_node=$NPROC \
#     --nnodes=$NNODES \
#     --node_rank=$NODE_RANK \
#     --master_addr=$MASTER_ADDR \
#     --master_port=$MASTER_PORT \
#     ./ODesign/scripts/inference.py \
#         exp="train_${infer_model_name}" \
#         data_root_dir="$data_root_dir" \
#         ckpt_root_dir="$ckpt_root_dir" \
#         exp.infer_model_name="$infer_model_name" \
#         exp.design_modality="$design_modality" \
#         exp.input_json_path="$input_json_path" \
#         exp.exp_name="$exp_name" \
#         exp.seeds="$seeds" \
#         exp.model.sample_diffusion.N_sample="$N_sample" \
#         exp.use_msa="$use_msa" \
#         exp.num_workers="$num_workers" \
#         exp.invfold_topk="$invfold_topk" \
#         hydra.run.dir="$output_dir"



echo "-----------------------------------------------------------"
echo "🎉 FINISHED: ODesign inference completed successfully!"
echo "-----------------------------------------------------------"

