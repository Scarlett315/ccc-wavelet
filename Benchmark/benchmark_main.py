#%%
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
LRP_database = pd.read_csv("CellChatDB/interaction_input_CellChatDB.csv", index_col=0)
aliases = pd.read_csv("Aliases_HGNC.tsv", sep="\t", index_col=0)
complexes = pd.read_csv("CellChatDB/complex_input_CellChatDB.csv", index_col=0)

#%%
# reading files
other_model_results = pd.read_csv("Human_BRCA/fig_1b.csv", index_col=0)
coords = pd.read_csv("../Data/human_breast_cancer/info/coordinates.csv", index_col=0)
expression = pd.read_csv("../Data/human_breast_cancer/info/expression_unfiltered.csv", index_col=0)

#%%
with open("Human_BRCA/reliable_CCC.json", "r") as f:
    accurate_CCC = json.load(f)
with open("Human_BRCA/agreed_LRPs.json", "r") as f:
    agreed_LRPs = json.load(f)

# %%
other_model_results = remove_models(other_model_results, ["stLearn", "SpaTalk", "Giotto", "NichNet"])
agreed_LRPs_rm = get_agreed_LRPs(other_model_results)
#%%
with open("Human_BRCA/agreed_LRPs_rm.json", "w") as f:
    json.dump(agreed_LRPs_rm, f)
#reliable_CCC = get_accurate_CCC_events(agreed_LRPs_rm, coords, expression, LRP_database, aliases, complexes)
#%%
other_model_results_agree = add_agreed_ints(other_model_results, agreed_LRPs_rm)
other_model_results_agree.head(50)
#%%
false_pos = calc_false_pos_all(other_model_results_agree, accurate_CCC)
false_pos = false_pos.T

print(false_pos.head())

# Sort by row averages (models) and column averages (interactions)
# Calculate row averages (average false positive rate per model)
row_means = false_pos.mean(axis=1)
# Calculate column averages (average false positive rate per interaction)
col_means = false_pos.mean(axis=0)
    
false_pos_sorted = false_pos.loc[row_means.sort_values(ascending=True).index]
#false_pos_sorted = false_pos_sorted[col_means.sort_values(ascending=True).index]

#%%
fig, ax = plt.subplots(figsize=(10, 3)) 
sns.heatmap(false_pos_sorted, ax=ax, vmin=0, vmax=1, cmap="seismic_r")

# %%
fig, ax = plt.subplots()
plt.xticks(rotation=90)
false_pos_sorted.T.boxplot(rot=90)

# %%
