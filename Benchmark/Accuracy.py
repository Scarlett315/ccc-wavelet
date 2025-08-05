#%%
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from sklearn.neighbors import NearestNeighbors
from scipy.spatial import cKDTree
import matplotlib.patches as mpatches
from itertools import chain

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
    

def get_agreed_LRPs(other_model_results, min_agreement=3):
    agreed_ints = get_agreed_interactions(create_method_interaction_df(other_model_results))
    agreed_LRPs = {} # agreed CCC events by cell type
    # Pre-extract the method name from each column, once:
    col2method = {col: col.split("_")[0] for col in other_model_results.columns}

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
def get_LRP_from_row_STCaseDB(name, LRP_db, separator="_"):
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

def get_LRP_from_row_CellChatDB(name, LRP_db, separator = "_"):
    lrp = [None, None]
    try:
        lrp[0] = LRP_db.loc[name, "ligand"]
        lrp[1] = LRP_db.loc[name, "receptor"]
    except KeyError: # complexes can be A_B or B_A  >:(
        parts = name.split(separator)
        pattern = "".join([f"(?=.*{part})" for part in parts])
        r = LRP_db.loc[LRP_db.index.str.contains(pattern, regex=True, case=False)]

        if r.empty or len(r) > 1: 
            raise KeyError("LRP, as named, is not anywhere in LRP database or is ambiguous. Throwing error...")
        r = r.squeeze()
        if len(r.name.split("_")) != len(parts):
            raise KeyError("LRP, as named, is not anywhere in LRP database or is ambiguous. Throwing error...")

        lrp[0] = r.at["ligand"]
        lrp[1] = r.at["receptor"]            
    return lrp

def get_LRP_from_row(name, LRP_db, LRP_db_name, separator="_"):
    if LRP_db_name == "CellChatDB":
        return get_LRP_from_row_CellChatDB(name, LRP_db)
    elif LRP_db_name == "STCaseDB":
        return get_LRP_from_row_STCaseDB(name, LRP_db)
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
    ligand = LRP[0].split("_")
    receptor = LRP[1].split("_")

    #try:
    ligand_expressed = expression.loc[ligand]
    receptor_expressed = expression.loc[receptor]
        
    #except KeyError:
        #print(f"Couldn't find one of {ligand} (ligand) or {receptor} (receptor) in the filtered expression matrix")
        #return (None, None)

    # Check if ALL genes in the complex are expressed above the threshold at each spot
    ligand_spots = (ligand_expressed > threshold).all(axis=0)
    receptor_spots = (receptor_expressed > threshold).all(axis=0)
   
    S = ligand_expressed.loc[:,ligand_spots]
    R = receptor_expressed.loc[:,receptor_spots]
    S = pd.Series({"barcodes": list(S.columns)}, name=LRP[0])
    R = pd.Series({"barcodes": list(R.columns)}, name=LRP[1])
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

    ligand_coords = coords.loc[S["barcodes"], ["x", "y"]] # ligand coordinates ("center")
    tree = cKDTree(coords.loc[R["barcodes"], ["x", "y"]]) # receptor coordinates
    for barcode, row in ligand_coords.iterrows():
        if barcode in coords.index:
            indices = tree.query_ball_point((row.x, row.y), radius, p=2) #p=2 -> Euclidean
            in_radius = np.array(R["barcodes"])[indices]
            
            if len(in_radius) > 0:
                S_interacting.add(barcode)
                R_interacting.update(in_radius)
                #print(f"Barcode {barcode} in S interacting with {len(in_radius)} spots in R")
        else:
            print(f"Barcode {barcode} not found in coords")
    return list(S_interacting), list(R_interacting)

#%%
def visualize_overlap(S, R, coords, radius, secreted, plot_interacting = False, down_factor = 1):
    ligand_coords = coords.loc[S.barcodes] # ligand coordinates ("center")
    receptor_coords = coords.loc[R.barcodes]

    ligand_coords["x"], ligand_coords["y"] = ligand_coords["y"], ligand_coords["x"]
    receptor_coords["x"], receptor_coords["y"] = receptor_coords["y"], receptor_coords["x"]

    fig, ax = plt.subplots(layout="constrained")
    ax.set_aspect("equal")
    ax.invert_yaxis()
    plt.title(f"{S.name} --> {R.name}, scaled: 1/{down_factor}")

    if secreted:
        radius *= 2
    if plot_interacting:
        S_int, R_int = get_interacting(S, R, coords, radius)
    
    for spot in ligand_coords.itertuples():
        color = "#FF9B00"
        if (plot_interacting) and (spot.Index in S_int):
            color = "#FF0006"
        ax.plot(spot.x / down_factor, spot.y / down_factor, 'o', color=color, markersize=0.5)
        circle = Circle((spot.x / down_factor, spot.y / down_factor), radius / down_factor,
                        edgecolor='none', facecolor=color, alpha=0.1, label='Radius')
        ax.add_patch(circle)

    for spot in receptor_coords.itertuples():
        color = "#00e7ff"
        if (plot_interacting) and (spot.Index in R_int):
            color = "#0021ff"

        ax.plot(spot.x / down_factor, spot.y / down_factor, 'o', color=color, markersize=0.5)


    # add a legend
    if plot_interacting:
        colors = {'ligands (interacting)': "#FF0006", 'ligands (not interacting)': "#FF9B00", 'receptors (interacting)': "#0021ff", "receptors (not interacting)": "#00e7ff"}
    else:
        colors = {'ligands': "#FF0006", "receptors": "#0021ff"}

    handles = [mpatches.Patch(color=color, label=label) for label, color in colors.items()]

    plt.legend(handles=handles, loc='lower left', title="Key", bbox_to_anchor = (0, -0.5))
    plt.show()


#%%
def calc_overlap_rate(LRP, coords, expression, threshold=0.1, LRP_db=None):
    sub = subset_spots(expression, LRP, threshold)
    if sub is None:
        return None
    else:
        S, R = sub

    sum_SR = len(S["barcodes"]) + len(R["barcodes"])

    radius = calc_diffusion_radius(coords) 

    # secreted versus non-secreted
    if LRP_db is not None:
        ligand_type = LRP_db.at[LRP[0] + "_" + LRP[1], "annotation"]
        #print(ligand_type)
        #CellChatDB
        if ligand_type == "Secreted Signaling" or ligand_type == "ECM-Receptor": # secreted = non-secreted * 2
            radius *=2
            #print(f"{LRP} - secreted")

    S_interacting, R_interacting = get_interacting(S, R, coords, radius) 

    sum_SR_interacting = len(S_interacting) + len(R_interacting)
    
    return sum_SR_interacting / sum_SR 

#%%
def determine_CCC_accuracy(LRP, coords, expression, threshold=0.95, LRP_db=None):
    '''
    determine whether a CCC event is accurate or not
    '''
    if LRP == None:
        return False
    overlap_rate = calc_overlap_rate(LRP, coords, expression, LRP_db=LRP_db)

    if overlap_rate >= threshold:
        return True
    else:
        return False

def get_accurate_CCC_events(other_method_results, coords, expression, LRP_db, LRP_db_name="CellChatDB", overlap_threshold=0.95):
    agreed_LRPs = get_agreed_LRPs(other_method_results) # organized by cell-type interaction
    agreed_LRPs_unorganized = list(set(chain.from_iterable(list(agreed_LRPs.values())))) # just all LRPs
    true_LRPs = {}
    not_found = []

    print(f"Getting accurate CCC events out of {len(agreed_LRPs_unorganized)} agreed events and \n\
          {len(agreed_LRPs.keys())} cell-type interactions")
    for i, event in enumerate(agreed_LRPs_unorganized, start=1): # each CCC event
        try:
            LRP = get_LRP_from_row(event, LRP_db, LRP_db_name)
            if LRP == None:
                reliable = False
                overlap_rate = -1
            overlap_rate = calc_overlap_rate(LRP, coords, expression, LRP_db=LRP_db)
            reliable = (overlap_rate >= overlap_threshold)
        except KeyError:
            print(f"Couldn't find CCC event {event}; trying alternate names...")
            #TODO: try alternate names!!!
            print(f"Couldn't find CCC event {event} somewhere in the process. Check for alternate naming...")
            not_found.append(event)
            continue

        print(f"{i}. {event}: {reliable}, overlap = {overlap_rate}")
        if reliable:
            cell_ints = [key for key, val in agreed_LRPs.items() if (event in val)] # all cell-type interactions this CCC evt is part of 
            for key in cell_ints:
                if key in true_LRPs:
                    true_LRPs[key].append(event)
                else:
                    true_LRPs[key] = [event]

    print(f"Couldn't find these {len(not_found)} LRPs: ")
    print(*not_found, sep="\n  - ")

    print(f"{len(list(chain.from_iterable(list(true_LRPs.values()))))} true LRPs out of {len(agreed_LRPs_unorganized)} total")
    return true_LRPs

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