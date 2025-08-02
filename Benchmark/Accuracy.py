#%%
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
import seaborn as sns
from sklearn.neighbors import NearestNeighbors
from scipy.spatial import cKDTree
import matplotlib.patches as mpatches

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
def get_agreed_interactions(method_interaction_df, add_to=False):
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
    

def get_agreed_LRPs(other_method_results, min_agreement=3):
    agreed_ints = get_agreed_interactions(create_method_interaction_df(other_method_results))
    agreed_LRPs = {} # agreed CCC events by cell type
    # Pre-extract the method name from each column, once:
    col2method = {col: col.split("_")[0] for col in other_method_results.columns}

    for cell_interaction in agreed_ints:
        int_df = other_model_results.loc[:,other_model_results.columns.str.contains(cell_interaction)]
        agreed_LRPs[cell_interaction] = []
        for lrp in int_df.itertuples():
            true_cols = int_df.columns[int_df.loc[lrp.Index]]
            methods_true = {col2method[c] for c in true_cols}
            if len(methods_true) >= min_agreement:
                agreed_LRPs[cell_interaction].append(lrp.Index)

    return agreed_LRPs

#%%
# gets LRP from the row name in the other model results table
def get_LRP_from_row(name, LRP_db, separator="_"):
    parts = name.split(separator)
    LRP = [None, None]
    # deal with complexes-- there can be 1-2 pairs
    # if there's ambiguity hopefully it's their problem :)
    num_genes = len(parts)
    if num_genes >=3:
        parts_left = parts

        possible_ligand_complex = "COMPLEX:" + "_".join([parts[0], parts[1]])
        possible_receptor_complex = "COMPLEX:" + "_".join([parts[-2], parts[-1]])

        if possible_ligand_complex in LRP_db["ligand"].tolist():
            LRP[0] = possible_ligand_complex
            del parts_left[0:1]
        if possible_receptor_complex in LRP_db["receptor"].tolist():
            LRP[1] = possible_receptor_complex
            del parts_left[-1]
            del parts_left[-1]
        
        #fill the other value
        if LRP[0] == None:
            LRP[0] = "_".join(parts_left)
        elif LRP[1] == None:
            LRP[1] = "_".join(parts_left)
    else:
        LRP = parts
    return LRP
#%%
# radius = avg dist. between all cell/spots and their 6 nearest neighbors
# non-secretory LRP radius = radius
# secretory LRP radius = radius * 2

# be careful with units
#TODO: convert into a proper unit of measurement before calculation? idk if it's needed though since
# it's the number of sets but might just be good to do
def calc_diffusion_radius(coords):
    coords_arr = np.array(coords)[:,:2]
    nbrs = NearestNeighbors(n_neighbors=7).fit(coords_arr) # self + 6 neighbors
    distances, _ = nbrs.kneighbors(coords_arr)
    #print(distances[:,1:])
    averaged_distances = np.mean(distances[:,1:], axis=1)

    radius = np.mean(averaged_distances)
    return radius

#%%
# get spots/barcodes expressing ligands and receptors
def subset_spots(expression, LRP, threshold=0.1):
    # TODO: remove when you fix the freaking expression thing to actually export correctly 
    expression = expression.drop("gene", axis=1) 

    ligand = LRP[0]
    receptor = LRP[1]

    ligand_expression = expression.loc[ligand]
    receptor_expression = expression.loc[receptor]

    S = ligand_expression[ligand_expression > threshold]
    R = receptor_expression[receptor_expression > threshold]

    return S, R
#%%
# determine which spots/barcodes are interacting based on spatial distance
# a lot of this is kind of unnecessary but whatever
def get_interacting(S, R, coords, radius):
    '''
    S: subset of cell/spots expressing ligands
    R: subset of cell/spots expressing receptors
    '''
    S_interacting = set()
    R_interacting = set()

    ligand_coords = coords.loc[S.index, ["x", "y"]] # ligand coordinates ("center")
    tree = cKDTree(coords.loc[R.index, ["x", "y"]]) # receptor coordinates
    for barcode, row in ligand_coords.iterrows():
        if barcode in coords.index:
            indices = tree.query_ball_point((row.x, row.y), radius, p=2) #p=2 -> Euclidean
            in_radius = R.iloc[indices]
            
            if len(in_radius) > 0:
                S_interacting.add(barcode)
                R_interacting.update(in_radius.index)
                #print(f"Barcode {barcode} in S interacting with {len(in_radius)} spots in R")
        else:
            print(f"Barcode {barcode} not found in coords")
    return list(S_interacting), list(R_interacting)

#%%
def visualize_interacting(S, R, coords, radius, secreted, plot_interacting = False, down_factor = 1):
    ligand_coords = coords.loc[S.index] # ligand coordinates ("center")
    receptor_coords = coords.loc[R.index]

    ligand_coords["x"], ligand_coords["y"] = ligand_coords["y"], ligand_coords["x"]
    receptor_coords["x"], receptor_coords["y"] = receptor_coords["y"], receptor_coords["x"]

    fig, ax = plt.subplots(layout="constrained")
    ax.set_aspect("equal")
    ax.invert_yaxis()
    plt.title(f"{S.name} --> {R.name}, scaled: 1/{down_factor}")

    if plot_interacting:
        S_int, R_int = get_interacting(S, R, coords, radius)
    if secreted:
        radius *= 2

    for spot in ligand_coords.itertuples():
        color = "#0021ff"
        if (plot_interacting != False) and (spot.Index in S_int):
            color = "#00e7ff"
        ax.plot(spot.x / down_factor, spot.y / down_factor, 'o', color=color, markersize=0.5)
        circle = Circle((spot.x / down_factor, spot.y / down_factor), radius / down_factor,
                        edgecolor='none', facecolor=color, alpha=0.1, label='Radius')
        ax.add_patch(circle)

    for spot in receptor_coords.itertuples():
        color = "#FF0006"
        if (plot_interacting != False) and (spot.Index in R_int):
            color = "#FF9B00"

        ax.plot(spot.x / down_factor, spot.y / down_factor, 'o', color=color, markersize=0.5)


    # add a legend
    if plot_interacting != False:
        colors = {'ligands (interacting)': "#00e7ff", 'ligands (not interacting)': "#0021ff", 'receptors (interacting)': "#FF9B00", "receptors (not interacting)": "#FF0006"}
    else:
        colors = {'ligands': "#0021ff", "receptors": "#FF0006"}

    handles = [mpatches.Patch(color=color, label=label) for label, color in colors.items()]

    plt.legend(handles=handles, loc='lower left', title="Key", bbox_to_anchor = (0, -0.5))
    plt.show()


#%%
def calc_overlap_rate(LRP, coords, expression, threshold=0, LRP_db=None):
    S, R = subset_spots(expression, LRP, threshold)
    sum_SR = len(S) + len(R) 

    radius = calc_diffusion_radius(coords) 

    # secreted versus non-secreted
    if LRP_db is not None:
        ligand_type = LRP_db.at[LRP[0] + "|" + LRP[1], "ligand_annotation"]
        #print(ligand_type)
        if ligand_type == "Secreted": # secreted = non-secreted * 2
            radius *=2
            print(f"{LRP} - secreted")

    S_interacting, R_interacting = get_interacting(S, R, coords, radius) 

    sum_SR_interacting = len(S_interacting) + len(R_interacting)
    
    return sum_SR_interacting / sum_SR 

#%%
def determine_CCC_accuracy(LRP, coords, expression, threshold=0.95, LRP_db=None):
    '''
    determine whether a CCC event is accurate or not
    '''
    overlap_rate = calc_overlap_rate(LRP, coords, expression, LRP_db=LRP_db)
    if overlap_rate >= threshold:
        return True
    else:
        return False

def get_accurate_CCC_events(other_method_results, coords, expression, LRP_db, overlap_threshold=0.95):
    agreed_LRPs = get_agreed_LRPs(other_method_results)
    true_LRPs = {}

    for cell_type in agreed_LRPs.keys(): # each cell-type interaction (multiple CCC events)
        true_LRPs[cell_type] = []
        for event in agreed_LRPs[cell_type]: # each CCC event
            LRP = get_LRP_from_row(event, LRP_db)
            print(LRP)
            reliable = determine_CCC_accuracy(LRP, coords, expression, threshold=overlap_threshold, LRP_db=LRP_db)
            print(f"{event}: {reliable}")
            if reliable:
                true_LRPs[cell_type].append(LRP)
    
    return true_LRPs

#%%
#false positive rates calculated for each cell-type interaction
def calc_false_pos(guessed_LRPs:dict, true_LRPs:dict):
    evt_false_pos = {}
    for cell_type in true_LRPs.keys():
        unreliable = len(guessed_LRPs[cell_type]) - len(true_LRPs[cell_type])
        evt_false_pos[cell_type] = unreliable/len(true_LRPs[cell_type])
    return evt_false_pos 

#calculate false positive rates for all models
def calc_false_pos_all(other_model_results, coords, expression,  LRP_db, overlap_threshold=0.95):
    true_LRPs = get_accurate_CCC_events(other_model_results, coords, expression, LRP_db)
    models = list(set([column.split("_")[0] for column in other_model_results.columns]))
    all_models_false_pos = {}

    for model in models:
        print(model)
        model_results = other_model_results.loc[:,other_model_results.columns.str.contains(model, regex=False)]
        model_results = model_results.to_dict()
        false_pos = calc_false_pos(model_results, true_LRPs)
        all_models_false_pos[model] = false_pos
    
    return pd.DataFrame(all_models_false_pos)


#%%
# reading files
other_model_results = pd.read_csv("extended_data_fig1/fig_1c.csv", index_col=0)
coords = pd.read_csv("../GSE208253/S1/info/coordinates.csv", index_col=0)
expression = pd.read_csv("../GSE208253/S1/info/expression_filtered.csv", index_col=0)
LRP_database = pd.read_csv("STCaseDB_Human.csv", index_col=0)
#%%
lrp = ["CTHRC1", "FZD6"]
S, R = subset_spots(expression, lrp)
radius = calc_diffusion_radius(coords)
visualize_interacting(S, R, coords, radius, True, plot_interacting=True, down_factor = 10)
#print(f"num ligands {len(S)}")
#print(f"num receptor{len(R)}")
#print(calc_overlap_rate(lrp, coords, expression, 0.01))

#%%
calc_false_pos_all(other_model_results, coords, expression, LRP_database)

# %%
