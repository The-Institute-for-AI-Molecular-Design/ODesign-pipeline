exp_name=pdl1_prot_binder
af3_out=./examples/prot_binder/af3/refold_results
filter_outdir=./examples/prot_binder/20260120_112530/filter_out
filter_json=./filter/score_json/miniprotein_filter_relax.json
binder_chain=A
target_chain=B
rosetta_xml=./filter/rosetta_cmds/ppi.xml

while (($#)); do case "$1" in
  --exp_name)        exp_name="$2";        shift 2 ;;
  --af3_out)         af3_out="$2";         shift 2 ;;
  --filter_outdir)   filter_outdir="$2";   shift 2 ;;
  --filter_json)     filter_json="$2";     shift 2 ;;
  --binder_chain)    binder_chain="$2";    shift 2 ;;
  --target_chain)    target_chain="$2";    shift 2 ;;
  --rosetta_xml)     rosetta_xml="$2";     shift 2 ;;
  *) echo "Error: Unknown argument '$1'"; exit 1 ;;
esac; done


export PYTHONPATH="$(pwd):$PYTHONPATH"

echo === starting af3 filter... ===
python ./filter/af3_filter.py \
    --src_dir "$af3_out" \
    --out_dir "$filter_outdir" \
    --filter_json "$filter_json" \
    --binder_chain "$binder_chain" \
    --target_chain "$target_chain" \
    --exp_name "$exp_name"

echo === starting rosetta filter... ===
python ./filter/rosetta_filter.py \
    --src_dir "$filter_outdir/af3_filter/${exp_name}_af3_filtered" \
    --filter_outdir "$filter_outdir" \
    --filter_json "$filter_json" \
    --binder_chain "$binder_chain" \
    --target_chain "$target_chain" \
    --exp_name "$exp_name" \
    --rosetta_xml "$rosetta_xml" \
    --merged True