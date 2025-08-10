
#%%
import pandas as pd
import numpy as np
#%%

#false positive rates calculated for each cell-type interaction
def calc_false_pos(guessed_LRPs:dict, true_LRPs:dict, interaction_total_dict: dict, model_name = ""):
    evt_false_pos = {}
    for cell_type in true_LRPs.keys():
        total = interaction_total_dict[cell_type]
        guessed_ind = model_name + "_" + cell_type
        unreliable = sum(guessed_LRPs[guessed_ind].values()) - len(true_LRPs[cell_type])
        evt_false_pos[cell_type] = unreliable/(total - len(true_LRPs[cell_type]))
    return evt_false_pos 

# get totals per cell-type interaction-- i.e. for each interaction get CCC events marked as "true" by at least one model
def get_totals_per_interaction(other_model_results):
    all_ints = ["_".join(i.split("_")[1:]) for i in other_model_results.columns]
    totals_per_int = {}
    for interaction in all_ints:
        int_cols = other_model_results.loc[:,other_model_results.columns.str.contains(interaction)]
        true_rows = int_cols[int_cols.any(axis=1)]
        totals_per_int[interaction] = len(true_rows)
    return totals_per_int

#calculate false positive rates for all models
def calc_false_pos_all(other_model_results:pd.DataFrame, true_LRPs):
    models = list(set([column.split("_")[0] for column in other_model_results.columns]))
    all_models_false_pos = {}
    totals_per_int = get_totals_per_interaction(other_model_results)

    for model in models:
        model_cols = other_model_results.loc[:,other_model_results.columns.str.contains(model, regex=False)]
        model_results = model_cols[model_cols.any(axis=1)]
        model_results = model_results.to_dict()

        false_pos = calc_false_pos(model_results, true_LRPs, totals_per_int, model_name=model)
        all_models_false_pos[model] = false_pos
    
    return pd.DataFrame(all_models_false_pos)




#%%
e = pd.read_csv("../GSE208253/S1/info/expression_matrix.csv", index_col=0)
c = pd.read_csv("../GSE208253/S1/info/coordinates.csv", index_col=0)
#%%
lrp = get_LRP_from_name("CXCL12_CXCR4", LRP_database, complexes, aliases, e)
print(lrp)
lrp = subset_spots(e, lrp)
radius = calc_diffusion_radius(c)
#%%
lrp = get_interacting(lrp, c, radius)
#%%
lrp.overlap_rate = calc_overlap_rate(lrp, c, radius)
visualize_overlap(lrp, c, radius, plot_interacting=True)
print(lrp.overlap_rate)

#%%
unorg = list(set(chain.from_iterable(list(agreed_LRPs.values())))) # just all LRPs
print(unorg[252])
# %%
print(unpack_complex("ACVR1_ACVR2A", complexes))

# %%
LRP_database.shape
# %%
