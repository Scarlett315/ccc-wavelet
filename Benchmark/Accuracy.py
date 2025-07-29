#%%
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.neighbors import NearestNeighbors

#%%
# reading files
data = pd.read_csv("/Users/scarlett/Programming/ccc-wavelet/Benchmark/extended_data_fig1/fig_1c.csv", index_col=0)
data_name ="Human Bronchus Dataset (1c)"
data.head()
data_cols = data.columns

coords = pd.read_csv("../GSE208253/S1/info/coordinates.csv", index_col=0)

#%%
def get_methods(data):
    interactions = {}
    columns = data.columns
    for col in columns:
        if data[col].any():
            col = col.split("_")
            method = col[0]
            interaction = "_".join(col[1:])

            if interaction not in interactions:
                interactions[interaction] = [method]
            else:
                if method not in interactions[interaction]:
                    interactions[interaction].append(method)
    return interactions 

def create_method_interaction_df(data):
    interactions = get_methods(data)
    method_interaction_df = pd.DataFrame(columns=interactions.keys())

    for interaction in interactions.keys():
        for method in interactions[interaction]:
            if method not in method_interaction_df.index:
                method_interaction_df.loc[method] = False
            method_interaction_df.loc[method, interaction] = True
    return method_interaction_df

#%% 
def get_agreement(method_interaction_df, add_to=False):
    method_interaction_df = method_interaction_df.copy()
    method_agreement = pd.DataFrame(columns=method_interaction_df.columns)

    true_pos = []
    for name, series in method_interaction_df.items():
        num_agreed = series.value_counts()[True]
        if num_agreed >= 3:
            true_pos.append(name)
            if add_to:
                match num_agreed:
                    case s if s >= 5:
                        method_agreement.loc[f"more than 5 agreed", name]=True
                    case s if s >= 4:
                        method_agreement.loc[f"more than 4 agreed", name]=True
                    case s if s >= 3:
                        method_agreement.loc[f"more than 3 agreed", name]=True
                    
    method_agreement = method_agreement.fillna(False)
    method_interaction_df = pd.concat([method_interaction_df, method_agreement], axis=0)
    if add_to:
        return method_interaction_df, true_pos
    else:
        return true_pos

#%%
# "diffusion rate was determined by using the coordinates
# of each cell/spot in S (cell/spot set expressing ligands) as the center of a circle with a specific radius"
# radius = avg dist. between all cell/spots and their 6 nearest neighbors
# non-secretory LRP radius = radius
# secretory LRP radius = radius * 2

# be careful with units
def calc_diffusion_radius(S_coords):
    S_coords_arr = np.array(S_coords)[:,:2]
    nbrs = NearestNeighbors(n_neighbors=7).fit(S_coords_arr) # self + 6 neighbors
    distances, _ = nbrs.kneighbors(S_coords_arr)
    #print(distances[:,1:])
    averaged_distances = np.mean(distances[:,1:], axis=1)

    radius_ns = np.mean(averaged_distances)
    radius_s = radius_ns * 2
    return radius_ns, radius_s

#def calc_overlap_rate(radii, )

#def filter_agreement(true_pos, data):
    

#%%
# create heatmap
method_interaction_df = create_method_interaction_df(data)
method_interaction_df, true_pos = get_agreement(method_interaction_df, add_to=True)
print(method_interaction_df)
#%%
sns.heatmap(method_interaction_df, cmap="YlGnBu")
plt.title(f"{data_name} - Method Agreement")
plt.show()

#%%
print(calc_diffusion_radius(coords))

# %%
