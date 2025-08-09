#%%
import pandas as pd
#%%
def get_num_correct(guessed_LRPs_for_int:list, true_LRPs_for_int:list):
    return len(set(guessed_LRPs_for_int) & set(true_LRPs_for_int))

def get_num_false_pos(guessed_LRPs_for_int:list, true_LRPs_for_int:list):
    return len(set(guessed_LRPs_for_int) - set(true_LRPs_for_int))

#%%
#false positive rates calculated for each cell-type interaction
def calc_false_pos(guessed_LRPs:dict, true_LRPs:dict, interaction_total_dict: dict, model_name = ""):
    evt_false_pos = {}
    for cell_type in true_LRPs.keys():
        total = interaction_total_dict[cell_type]
        guessed_ind = model_name + "_" + cell_type
        
        # Check if the guessed key exists
        if guessed_ind not in guessed_LRPs:
            print(f"Warning: {guessed_ind} not found in guessed_LRPs")
            print(f"Available keys: {list(guessed_LRPs.keys())}")
            evt_false_pos[cell_type] = 0.0
            continue
        
        #false pos rate calculation
        unreliable = get_num_false_pos(guessed_LRPs[guessed_ind], true_LRPs[cell_type])
        true_positives = len(true_LRPs[cell_type])
        pred_positives = len(guessed_LRPs[guessed_ind])
        print(f"Unreliable: {unreliable}, True positives: {true_positives}, Predicted positives: {pred_positives}")

        # FPR = FP / (FP + TN)
        denominator = pred_positives - true_positives
        if denominator == 0:
            evt_false_pos[cell_type] = 0.0
        else:
            evt_false_pos[cell_type] = unreliable / denominator
            
    return evt_false_pos 

# get totals per cell-type interaction-- i.e. for each interaction get CCC events marked as "true" by at least one model
def get_totals_per_interaction(other_model_results):
    # Get all unique interactions
    all_ints = list(set(["_".join(i.split("_")[1:]) for i in other_model_results.columns]))
    totals_per_int = {}
    
    # Get total unique events across all interactions
    all_true_rows = other_model_results[other_model_results.any(axis=1)]
    total_unique_events = len(all_true_rows)
    print(f"Debug: Total unique events across all interactions: {total_unique_events}")
    
    for interaction in all_ints:
        int_cols = other_model_results.loc[:,other_model_results.columns.str.contains(interaction)]
        true_rows = int_cols[int_cols.any(axis=1)]
        # Count unique events (row indices) that were predicted by at least one model for this interaction
        totals_per_int[interaction] = len(true_rows)
    
    print(f"Debug: Totals per interaction: {totals_per_int}")
    
    return totals_per_int

#calculate false positive rates for all models
def calc_false_pos_all(other_model_results:pd.DataFrame, true_LRPs):
    models = list(set([column.split("_")[0] for column in other_model_results.columns]))
    all_models_false_pos = {}
    totals_per_int = get_totals_per_interaction(other_model_results)
    
    print(f"Debug: DataFrame shape: {other_model_results.shape}")
    print(f"Debug: DataFrame columns: {list(other_model_results.columns[:10])}...")  # First 10 columns
    print(f"Debug: Models found: {models}")
    print(f"Debug: True_LRPs keys: {list(true_LRPs.keys())}")
   
    for model in models:
        model_cols = other_model_results.loc[:,other_model_results.columns.str.contains(model, regex=False)]
        model_results = model_cols[model_cols.any(axis=1)]
        model_results = model_results.to_dict()
        
        # Check a few sample values to see what the data looks like
        if len(model_results) > 0:
            sample_key = list(model_results.keys())[0]
            sample_values = list(model_results[sample_key])[:5]  # First 5 values
            print(f"  - Sample values from {sample_key}: {sample_values}")

        false_pos = calc_false_pos(model_results, true_LRPs, totals_per_int, model_name=model)
        all_models_false_pos[model] = false_pos
    
    return pd.DataFrame(all_models_false_pos)

# %%
