#%%
import pandas as pd
import numpy as np

#%% 
def add_agreed_ints(other_model_results: pd.DataFrame, agreed_LRPs):
    other_model_results = other_model_results.copy()
    
    for interaction, lrps in agreed_LRPs.items():
        print(f"{interaction}: {lrps}")
        other_model_results.loc[lrps, f"3-or-More-Models_{interaction}"] = True
        other_model_results[f"3-or-More-Models_{interaction}"] = other_model_results[f"3-or-More-Models_{interaction}"].fillna(False)
    return other_model_results

#%%
def get_events_per_interaction(other_model_results: pd.DataFrame):
    all_ints = ["_".join(i.split("_")[1:]) for i in other_model_results.columns]
    evts_per_int = {}
    for interaction in all_ints:
        int_cols = other_model_results.loc[:,other_model_results.columns.str.contains(interaction)]
        true_rows = int_cols[int_cols.any(axis=1)]
        evts_per_int[interaction] = true_rows
    return evts_per_int

def get_totals_per_interaction(totals_per_int):
    return {key: sum(val) for key, val in totals_per_int}
    
#%%
def get_agreed_LRPs(other_model_results, min_agreement=3):
    interactions = ["_".join(i.split("_")[1:]) for i in other_model_results.columns]
    agreed_LRPs = {} # agreed CCC events by cell type
    # Pre-extract the method name from each column, once:
    col2method = {col: col.split("_")[0] for col in other_model_results.columns}

    #TODO: this is really inefficient probably
    for cell_interaction in interactions:
        print(f"getting agreed LRPs for {cell_interaction}")
        int_df = other_model_results.loc[:,other_model_results.columns.str.contains(cell_interaction)]
        agreed_LRPs[cell_interaction] = []
        for lrp in int_df.itertuples():
            true_cols = int_df.columns[int_df.loc[lrp.Index]]
            methods_true = {col2method[c] for c in true_cols}
            if len(methods_true) >= min_agreement:
                agreed_LRPs[cell_interaction].append(lrp.Index)
    return agreed_LRPs
