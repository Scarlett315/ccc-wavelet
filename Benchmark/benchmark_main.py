#%%
from tkinter import E
import pandas as pd
import numpy as np
import json
import seaborn as sns
from Agreement import *
from GetTrueCCC import *
from Accuracy import *
import matplotlib.pyplot as plt
from Visualization import *
#%%
# reading files
other_model_results = pd.read_csv("Human_BRCA/fig_1b.csv", index_col=0)
coords = pd.read_csv("../Data/human_breast_cancer/info/coordinates.csv", index_col=0)
expression = pd.read_csv("../Data/human_breast_cancer/info/expression_unfiltered.csv", index_col=0)
LRP_database = pd.read_csv("CellChatDB/interaction_input_CellChatDB.csv", index_col=0)
aliases = pd.read_csv("Aliases_HGNC.tsv", sep="\t", index_col=0)
complexes = pd.read_csv("CellChatDB/complex_input_CellChatDB.csv", index_col=0)
#%%
with open("Human_BRCA/reliable_CCC.json", "r") as f:
    accurate_CCC = json.load(f)
with open("Human_BRCA/agreed_LRPs.json", "r") as f:
    agreed_LRPs = json.load(f)
# %%
other_model_results.head()
# %%
other_model_results = remove_models(other_model_results, ["stLearn", "SpaTalk", "Giotto", "NichNet"])
other_model_results_agree = add_agreed_ints(other_model_results, agreed_LRPs)
other_model_results_agree.head(50)
#%%
false_pos = calc_false_pos_all(other_model_results_agree, accurate_CCC)
false_pos = false_pos.T
row_averages = false_pos.mean(axis=1)
sorted_false_pos = false_pos.loc[row_averages.sort_values().index]
#%%
fig, ax = plt.subplots(figsize=(10, 3)) 
sns.heatmap(sorted_false_pos, ax=ax, vmin=0, vmax=1)
# %%
print(get_totals_per_interaction(other_model_results))

# %%
fig, ax = plt.subplots()
plt.xticks(rotation=90)
sorted_false_pos.T.boxplot(rot=90)
#%%
e = pd.read_csv("../GSE208253/S1/info/expression_matrix.csv", index_col=0)
c = pd.read_csv("../GSE208253/S1/info/coordinates.csv", index_col=0)

lrp = get_LRP_from_name("CXCL12_CXCR4", LRP_database, complexes, aliases, e)
print(lrp)
lrp = subset_spots(e, lrp)
radius = calc_diffusion_radius(c)
lrp = get_interacting(lrp, c, radius)
visualize_overlap(lrp, c, radius, plot_interacting=True)
overlap = calc_overlap_rate(lrp, c, radius)
print(overlap)

#%%
unorg = list(set(chain.from_iterable(list(agreed_LRPs.values())))) # just all LRPs
print(unorg[252])
# %%
print(unpack_complex("ACVR1_ACVR2A", complexes))

# %%
LRP_database.shape
# %%
