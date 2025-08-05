#%%
import pandas as pd
import numpy as np
from sklearn.neighbors import NearestNeighbors
from scipy.spatial import cKDTree
from itertools import chain
from Agreement import *
#%%
class lrp:
    def __init__(self, name, ligand, receptor, secreted):
        self.name = name
        self.ligand = ligand
        self.receptor = receptor
        self.secreted = secreted

    def __str__(self):
        return f"{self.name} ({self.ligand} + {self.receptor})"
#%%
# gets LRP from the row name in the other model results table
# this doesn't work ignore it
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

def get_LRP_from_row_CellChatDB(name, LRP_db, aliases, separator = "_"):
    lrp = [None, None]
    try:
        if name in LRP_db.index:
            lrp[0] = LRP_db.loc[name, "ligand"]
            lrp[1] = LRP_db.loc[name, "receptor"]
        else:
            name_alt = get_alt_name(name, aliases, LRP_db)
            lrp[0] = name_alt[0]
            lrp[1] = name_alt[1]

    except KeyError: # complexes can be A_B or B_A  >:(
        parts = name.split(separator)
        pattern = "".join([f"(?=.*{part})" for part in parts])
        r = LRP_db.loc[LRP_db.index.str.contains(pattern, regex=True, case=False)]

        # coudn't find or multiple matches
        if r.empty: 
            raise KeyError("LRP, as named, is not anywhere in LRP database or is ambiguous. Throwing error...")
        r = r.squeeze()

        # found all parts but had an extra gene
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
    
def gene_in_expression(gene, expression):
    print(expression.index)
    if gene in expression.index:
        return True
    else:
        return False
    
def get_genes_in_LRP_db(LRP_db:pd.DataFrame):
    ligand_genes = list(LRP_db["ligand"])
    receptor_genes = list(LRP_db["receptor"])
    combined = list(set(ligand_genes + receptor_genes))
    return combined

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
def get_gene_name(name, aliases, LRP_db):
    genes_in_LRP_db = get_genes_in_LRP_db(LRP_db)
    if name in genes_in_LRP_db:
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
            
            common = list(set(all_symbols) & set(genes_in_LRP_db))
            print(f"Common: {common}")
            if len(common) != 1:
                print(f"Couldn't find {name} in aliases.")
                return None
            else:
                return list(common)[0]
#%%
def get_LRP_name(name, aliases, LRP_db):
    if name in LRP_db.index:
        return name
    else:
        final_name = []
        genes = name.split("_")
        for gene in genes:
            gene_in_lrp = get_gene_name(gene, aliases, LRP_db)
            final_name.append(gene_in_lrp)

        # make sure that new LRP is in the database
        new_LRP = "_".join(final_name)
        if new_LRP in LRP_db.index:
            return new_LRP
        else:
            print(f"Couldn't find {name}.")
            return None
#%%
def standardize_name(event: str, aliases: pd.DataFrame, LRP_db: pd.DataFrame, expression: pd.DataFrame):
    def split_to_list(s):
        return s.split(", ") if pd.notna(s) else []
    aliases = aliases.copy()
    aliases["Alias symbols"] = aliases["Alias symbols"].apply(split_to_list)

    parts = event.split("_")
    genes_in_LRP_db = get_genes_in_LRP_db(LRP_db)
    expression = expression.copy()
    new_parts = []

    for gene in parts:
        gene_row = aliases[
            aliases["Approved symbol"].eq(gene) |
            aliases["Alias symbols"].fillna("").apply(lambda x: gene in x)
            ]
        print(gene)

        if gene_row.empty:
            print(f"Alias for {gene} not found.")
            new_parts.append(gene)
            continue

        #getting all possible names
        approved = gene_row["Approved symbol"].values[0]
        aliases_list = gene_row["Alias symbols"]
        possible_names = set(aliases_list + [approved])

        print(f"{gene} possible names: {possible_names}")

        common_exp = list(possible_names & set(expression.index)) # what the gene is called in expression
        common_lrp = list(possible_names & set(genes_in_LRP_db)) # what the gene is called in LRP database
        
        if len(common_exp) != 1 or len(common_lrp) != 1:
            print(f"{gene} has multiple or no matches.")
            continue
        else:
            common_exp = str(common_exp[0])
            common_lrp = str(common_lrp[0])
        
        if common_exp != common_lrp:
            expression.loc[common_exp, "gene"] = common_lrp
        new_parts.append(common_lrp)

    new_name = "_".join(new_parts)
    return expression, new_name
    
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
def subset_spots(expression, LRP, threshold=0.1):
    ligand = LRP[0].split("_")
    receptor = LRP[1].split("_")

    ligand_expressed = expression.loc[ligand]
    receptor_expressed = expression.loc[receptor]

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
def calc_overlap_rate(LRP, coords, radius_unsecreted, sub:tuple, LRP_db=None):
    if sub is None:
        return None
    else:
        S, R = sub

    sum_SR = len(S["barcodes"]) + len(R["barcodes"])

    # secreted versus non-secreted
    if LRP_db is not None:
        ligand_type = LRP_db.at[LRP[0] + "_" + LRP[1], "annotation"]
        if ligand_type == "Secreted Signaling" or ligand_type == "ECM-Receptor": # secreted = non-secreted * 2
            radius = radius_unsecreted * 2
        else:
            radius = radius_unsecreted

    S_interacting, R_interacting = get_interacting(S, R, coords, radius) 

    sum_SR_interacting = len(S_interacting) + len(R_interacting)
    
    return sum_SR_interacting / sum_SR 
#%%
def get_accurate_CCC_events(agreed_LRPs, coords, expression, LRP_db, aliases:pd.DataFrame=None, LRP_db_name="CellChatDB", overlap_threshold=0.95):
    def is_reliable(event, radius):
        LRP = get_LRP_from_row(event, LRP_db, LRP_db_name)
        if LRP == None:
            reliable = False
            overlap_rate = -1
        sub = subset_spots(expression, LRP)
        overlap_rate = calc_overlap_rate(LRP, coords, radius, sub, LRP_db=LRP_db)
        reliable = (overlap_rate >= overlap_threshold)
        return reliable, overlap_rate
    
    def try_alt_names(event):
        if aliases is None:
            return None
        parts = event.split("_")
        missing = [part for part in parts if (not name_in_LRP_db(part, LRP_db))]
        
    agreed_LRPs_unorganized = list(set(chain.from_iterable(list(agreed_LRPs.values())))) # just all LRPs
    true_LRPs = {}
    not_found = []

    print(f"Getting accurate CCC events out of {len(agreed_LRPs_unorganized)} agreed events and \n{len(agreed_LRPs.keys())} cell-type interactions")
    radius_unsec = calc_diffusion_radius(coords)
    for i, event in enumerate(agreed_LRPs_unorganized, start=1): # each CCC event
        try:
            reliable, overlap_rate = is_reliable(event, radius_unsec)
        except KeyError:
            print(f"Couldn't find CCC event {event}; trying alternate names...")
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