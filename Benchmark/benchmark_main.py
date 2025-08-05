#%%
#from Accuracy import *
import pandas as pd
import numpy as np
import json
import seaborn as sns
from Agreement import *
from GetTrueCCC import *
from Visualization import *
#%%
# reading files
other_model_results = pd.read_csv("extended_data_fig1/fig_1b.csv", index_col=0)
coords = pd.read_csv("../Data/human_breast_cancer/info/coordinates.csv", index_col=0)
expression = pd.read_csv("../Data/human_breast_cancer/info/expression_unfiltered.csv", index_col=0)
LRP_database = pd.read_csv("CellChatDB/interaction_input_CellChatDB.csv", index_col=0)
aliases = pd.read_csv("Aliases_HGNC.tsv", sep="\t", index_col=0)
#%%
lrp = get_LRP_from_row("ANGPTL2_ITGA5_ITGB1", LRP_database, "CellChatDB")
print(lrp)
S, R = subset_spots(expression, lrp)
radius = calc_diffusion_radius(coords)
visualize_overlap(S, R, coords, radius, plot_interacting=True, secreted=True)
overlap = calc_overlap_rate(lrp, coords, expression, LRP_db=LRP_database)
print(overlap)
#get_interacting(S, R, coords, radius)
#%%
# Open the JSON file in read mode
with open('reliable_CCC_events.json', 'r') as file:
    accurate_CCC = json.load(file)
#%%
false_pos = calc_false_pos_all(other_model_results, accurate_CCC)
false_pos = false_pos.T
row_averages = false_pos.mean(axis=1)
sorted_false_pos = false_pos.loc[row_averages.sort_values().index]

#%%
fig, ax = plt.subplots(figsize=(10, 3)) 
sns.heatmap(sorted_false_pos, ax=ax)
# %%
print(get_totals_per_interaction(other_model_results))

# %%
fig, ax = plt.subplots()
plt.xticks(rotation=90)
sorted_false_pos.T.boxplot(rot=90)


# %%
other_model_results.head()
agreed_LRPs = get_agreed_LRPs(other_model_results)
# %%
other_model_results_agree = add_agreed_ints(other_model_results, agreed_LRPs)
other_model_results_agree.head(50)
# %%
aliases = pd.read_csv("Aliases_HGNC.tsv", sep="\t", index_col=0)
aliases.head()
# %%
exp_r, name = standardize_name("ANGPTL2_PIRB", aliases, LRP_database, expression)
print(name)

# %%
print(aliases.loc["HGNC:7", "Alias symbols"].split(", "))
# %%
print(get_LRP_name("ANGPTL2_LILRB1", aliases, LRP_database))

# %%
print(get_gene_name("PIRB", aliases, LRP_database))
# %%
