#%%
from numba import none
import pandas as pd
import numpy as np
from sklearn.neighbors import NearestNeighbors
from scipy.spatial import cKDTree
from itertools import chain
from Agreement import *
#%%
class LRP:
    def __init__(self, name:str, ligand:str, receptor:str, secreted:bool):
        self.name = name # name in LRP database
        self.ligand = ligand # name in expression/results
        self.receptor = receptor # name in expression/results
        self.secreted = secreted
        self.S = None
        self.R = None
        self.S_int = set()
        self.R_int = set()
        self.overlap_rate = None
    
    def add_subsets(self, S, R):
        self.S = S
        self.R = R

    def add_interacting(self, S_int, R_int):
        self.S_int = S_int
        self.R_int = R_int
    
    def add_overlap_rate(self, overlap_rate):
        self.overlap_rate = overlap_rate

    def __str__(self):
        return f"{self.name} ({self.ligand} + {self.receptor}); secreted: {self.secreted}"

#%%
def get_LRP_from_name(name, LRP_db, complexes, aliases, expression,separator = "_"):
    def compare_num_parts(row, parts):
        return len(row.index.str.split(separator)) == len(parts)
    # LRP (as named in expression/results) is in LRP_db
    if name in LRP_db.index:
        r = LRP_db.loc[name]
    else:
        parts = name.split(separator)

        pattern = "".join([f"(?=.*{part})" for part in parts])
        r = LRP_db.loc[LRP_db.index.str.contains(pattern, regex=True, case=False)]
        
        #filter out rows that don't have the exact number of parts
        r = r[r.apply(compare_num_parts, parts=parts, axis=1)]

        # coudn't find or multiple matches
        if r.empty or len(r) > 1: 
            raise KeyError("LRP, as named, is not anywhere in LRP database or is ambiguous. Throwing error...")
        r = r.iloc[0]

    secreted = r.at["annotation"] == "Secreted Signaling" or r.at["annotation"] == "ECM-Receptor"
    ligand_in_db = unpack_complex(r.at["ligand"], complexes)
    receptor_in_db = unpack_complex(r.at["receptor"], complexes)

    ligand = [get_gene_name_in_exp(gene, aliases, expression) for gene in ligand_in_db]
    receptor = [get_gene_name_in_exp(gene, aliases, expression) for gene in receptor_in_db]

    ligand = "_".join(ligand)
    receptor = "_".join(receptor)

    return LRP(name, ligand, receptor, secreted)
#%%
def unpack_complex(name, complexes):
    if name in complexes.index:
        # Get all components (all elements in the row, not skipping the first)
        components = list(complexes.loc[name])
        
        # Filter out NaN values
        components = [comp for comp in components if (comp is not None and not pd.isna(comp))]
        return components
    else:
        # Check if name is a string before splitting
        if isinstance(name, str):
            return list(name.split("_"))
        else:
            print(f"Warning: name is not a string, returning as single item: {name}")
            return [str(name)]

#%%
def find_gene_row(symbol, aliases):
    # First check if it's an approved symbol
    approved_match = aliases[aliases["Approved symbol"] == symbol]
    if not approved_match.empty:
        return approved_match.iloc[0]
    
    # Then check if it's in any of the alias symbols
    def symbol_in_aliases(row):
        if pd.isna(row["Alias symbols"]):
            return False
        alias_list = row["Alias symbols"].split(", ")
        return symbol in alias_list
    
    alias_match = aliases[aliases.apply(symbol_in_aliases, axis=1)]
    if not alias_match.empty:
        return alias_match.iloc[0]
    
    return None

#%%
# get the name/alias of the gene that is used in the LRP database
def get_gene_name_in_exp(name, aliases, expression):
    if name in expression.index:
        return name
    else:
        # Use the new helper function to find the row
        gene_row = find_gene_row(name, aliases)
        if gene_row is None:
            print(f"Couldn't find {name} in aliases.")
            return None
        else:
            # Get all possible symbols (approved + aliases)
            approved_symbol = gene_row["Approved symbol"]
            alias_symbols = gene_row["Alias symbols"].split(", ") if pd.notna(gene_row["Alias symbols"]) else []
            all_symbols = [approved_symbol] + alias_symbols
            
            common = list(set(all_symbols) & set(expression.index))
            if len(common) != 1:
                print(f"{name} is ambiguous or not in sexpression.")
                return None
            else:
                return list(common)[0]
#%%
# radius = avg dist. between all cell/spots and their 6 nearest neighbors
# non-secretory LRP radius = radius
# secretory LRP radius = radius * 2

# be careful with units
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
def subset_spots(expression, lrp, threshold=0.1):
    ligand = lrp.ligand.split("_")
    receptor = lrp.receptor.split("_")

    ligand_expressed = expression.loc[ligand]
    receptor_expressed = expression.loc[receptor]

    # Check if ALL genes in the complex are expressed above the threshold at each spot
    ligand_spots = (ligand_expressed > threshold).all(axis=0)
    receptor_spots = (receptor_expressed > threshold).all(axis=0)
   
    S = ligand_expressed.loc[:,ligand_spots]
    R = receptor_expressed.loc[:,receptor_spots]
    S = pd.Series({"barcodes": list(S.columns)}, name=lrp.ligand)
    R = pd.Series({"barcodes": list(R.columns)}, name=lrp.receptor)
    lrp.add_subsets(S, R)
    return lrp
#%%
# determine which spots/barcodes are interacting based on spatial distance
# a lot of this is kind of unnecessary but whatever
def get_interacting(lrp, coords, radius):
    '''
    S: subset of cell/spots expressing ligands
    R: subset of cell/spots expressing receptors
    '''
    S, R = lrp.S, lrp.R
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
    lrp.add_interacting(list(S_interacting), list(R_interacting))
    return lrp
#%%
def calc_overlap_rate(lrp, coords, radius_unsecreted):
    S, R = lrp.S, lrp.R
    sum_SR = len(S["barcodes"]) + len(R["barcodes"])

    # secreted versus non-secreted
    if lrp.secreted:
        radius = radius_unsecreted * 2
    else:
        radius = radius_unsecreted

    lrp = get_interacting(lrp, coords, radius) 
    sum_SR_interacting = len(lrp.S_int) + len(lrp.R_int)
    
    return sum_SR_interacting / sum_SR 
#%%
def get_accurate_CCC_events(agreed_LRPs, coords, expression, LRP_db, aliases, complexes, overlap_threshold=0.95):
    def is_reliable(event, radius):
        lrp = get_LRP_from_name(event, LRP_db, complexes, aliases, expression)
        lrp = subset_spots(expression, lrp)
        lrp.overlap_rate = calc_overlap_rate(lrp, coords, radius)
        lrp.reliable = lrp.overlap_rate >= overlap_threshold
        return lrp
    
    agreed_LRPs_unorganized = list(set(chain.from_iterable(list(agreed_LRPs.values())))) # just all LRPs
    true_LRPs = {}
    not_found = []

    print(f"Getting accurate CCC events out of {len(agreed_LRPs_unorganized)} agreed events and \n{len(agreed_LRPs.keys())} cell-type interactions")
    radius_unsec = calc_diffusion_radius(coords)
    for i, event in enumerate(agreed_LRPs_unorganized, start=1): # each CCC event
        try:
            lrp = is_reliable(event, radius_unsec)
        except KeyError:
            print(f"Couldn't find CCC event {event} somewhere in the process.")
            not_found.append(event)
            continue

        print(f"{i}. {event}: {lrp.reliable}, overlap = {lrp.overlap_rate}")
        if lrp.reliable:
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
# %%
