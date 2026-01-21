import json

def get_af3_confidence(summary_json, binder_chain):
    try:
        with open(summary_json, 'r') as f:
            data = json.load(f)
    except json.decoder.JSONDecodeError:
        print(f"Warning: JSONDecodeError in file {summary_json}. Skipping...")
        return None
    except FileNotFoundError:
        print(f"File not found: {summary_json}")
        return None
    
    binder_id = ord(binder_chain) - 65
    complex_chains = len(data["chain_pair_pae_min"])
    
    if complex_chains > 2:
        iptm = data["chain_iptm"][binder_id]
    else:
        iptm = data["iptm"]
    complex_ptm = data["ptm"]
    
    chain_ptm_avg = sum(data["chain_ptm"]) / len(data["chain_ptm"])
    binder_ptm = data["chain_ptm"][binder_id]
    
    ipae_list = [data["chain_pair_pae_min"][binder_id][i] for i in range(0, complex_chains) if i != binder_id] + \
        [data["chain_pair_pae_min"][i][binder_id] for i in range(0, complex_chains) if i != binder_id]
        
    ipae_min = min(ipae_list)
    ipae_avg = sum(ipae_list) / len(ipae_list)
    
    return {
        "iptm": iptm,
        "binder_ptm": binder_ptm,
        "complex_ptm": complex_ptm,
        "chain_ptm_avg": f"{chain_ptm_avg:.2f}",
        "ipae_min": ipae_min,
        "ipae_avg": f"{ipae_avg:.2f}"
    }

def check_filter(sample_score_dict: dict, filter_dict: dict):
    
    fail = {}
    check_filter = True
    for filter, threshold in filter_dict.items():
        score = sample_score_dict.get(filter)
        
        if threshold["higher"]:
            if float(score) < float(threshold["threshold"]):
                fail.update({f"{filter}_pass": 'false'})
                check_filter = False
            else:
                fail.update({f"{filter}_pass": 'true'})
        else:
            if float(score) > float(threshold["threshold"]):
                fail.update({f"{filter}_pass": 'false'})
                check_filter = False
            else:
                fail.update({f"{filter}_pass": 'true'})
                
    if check_filter:
        return True, fail
    else:
        return False, fail
    