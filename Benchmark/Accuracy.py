#%%
import pandas as pd
#%%
def get_correct(guessed_LRPs_for_int:list, true_LRPs_for_int:list):
    return set(guessed_LRPs_for_int) & set(true_LRPs_for_int)

def get_false_pos(guessed_LRPs_for_int:list, unreliable_lrps:list):
    return set(guessed_LRPs_for_int).intersection(set(unreliable_lrps))

#%%
#false positive rates calculated for each cell-type interaction
def calc_false_pos(guessed_LRPs:dict, true_LRPs:dict, total_lrps: dict, model_name = ""):
    evt_false_pos = {}
    for cell_type in true_LRPs.keys():
        guessed_ind = model_name + "_" + cell_type

        # Check if the guessed key exists
        if guessed_ind not in guessed_LRPs:
            print(f"Warning: {guessed_ind} not found in guessed_LRPs")
            print(f"Available keys: {list(guessed_LRPs.keys())}")
            evt_false_pos[cell_type] = 0.0
            continue
        
        #false pos rate calculation
        guessed_lrps = guessed_LRPs[guessed_ind]
        all_lrps = total_lrps[cell_type]
        true_lrps = true_LRPs[cell_type]
        unreliable_lrps = set(total_lrps[cell_type]) - set(true_lrps)

        false_pos = get_false_pos(guessed_lrps, unreliable_lrps)
        guessed_correct = get_correct(guessed_lrps, true_lrps)

        # FPR = FP / (FP + TN)
        denominator = len(unreliable_lrps)
        print(f"False Positives: {len(false_pos)}, Denominator: {denominator}, all lrps: {len(all_lrps)}, guessed correct: {len(guessed_correct)}")
        if denominator == 0:
            evt_false_pos[cell_type] = 0.0
        else:
            evt_false_pos[cell_type] = len(false_pos) / denominator
            
    return evt_false_pos 
#%%
def get_interactions(other_model_results):
    # Get all unique interactions
    all_ints = list(set(["_".join(i.split("_")[1:]) for i in other_model_results.columns]))
    lrps_per_int = {}
    
    for interaction in all_ints:
        int_cols = other_model_results.loc[:,other_model_results.columns.str.contains(interaction)]
        true_rows = int_cols[int_cols.any(axis=1)]
        # Count unique events (row indices) that were predicted by at least one model for this interaction
        lrps_per_int[interaction] = list(set(true_rows.index))
        
    return lrps_per_int
#%%
# get totals per cell-type interaction-- i.e. for each interaction get CCC events marked as "true" by at least one model
def get_totals_per_interaction(other_model_results):
    lrps_per_int = get_interactions(other_model_results)
    totals_per_int = {key: len(value) for key, value in lrps_per_int.items()}
    
    return totals_per_int
#%%
#calculate false positive rates for all models
def calc_false_pos_all(other_model_results:pd.DataFrame, true_LRPs, agreed_LRPs=None):
    models = list(set([column.split("_")[0] for column in other_model_results.columns]))
    all_models_false_pos = {}
    lrps_per_int = get_interactions(other_model_results)
   
    for model in models:
        model_cols = other_model_results.loc[:,other_model_results.columns.str.contains(model, regex=False)]
        model_results = model_cols[model_cols.any(axis=1)]

        model_results = model_results.to_dict()
        print(model_results.keys())        
        false_pos = calc_false_pos(model_results, true_LRPs, lrps_per_int, model_name=model)
        all_models_false_pos[model] = false_pos
    
    return pd.DataFrame(all_models_false_pos)
# %%
