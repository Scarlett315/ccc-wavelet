#%%
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
import matplotlib.patches as mpatches
from GetTrueCCC import *
#%%
def visualize_overlap(lrp, coords, radius, plot_interacting = False):
    S, R = lrp.S, lrp.R
    coords = coords.copy()
    coords["x"], coords["y"] = coords["y"], coords["x"]
    ligand_coords = coords.loc[S["barcodes"]] # ligand coordinates ("center")
    receptor_coords = coords.loc[R["barcodes"]]
    autocrine_coords = coords.loc[lrp.autocrine]

    fig, ax = plt.subplots(layout="constrained")
    ax.set_aspect("equal")
    ax.invert_yaxis()
    plt.title(f"{S.name} --> {R.name}, overlap: {lrp.overlap_rate}")

    # get_interacting now handles radius doubling, but we need it for visualization
    if lrp.secreted:
        radius *= 2
    if plot_interacting:
        S_int, R_int = set(lrp.S_int), set(lrp.R_int)
        

    for spot in ligand_coords.itertuples():
        color = "#FF9B00"
        if (plot_interacting) and (spot.Index in S_int):
            color = "#FF0006"
        ax.plot(spot.x, spot.y, 'o', color=color, markersize=0.5)
        circle = Circle((spot.x, spot.y), radius,
                        edgecolor='none', facecolor=color, alpha=0.1, label='Radius')
        ax.add_patch(circle)

    for spot in receptor_coords.itertuples():
        color = "#00e7ff"
        is_interacting = (plot_interacting) and (spot.Index in R_int)
        if is_interacting:
            color = "#0021ff"


        ax.plot(spot.x, spot.y, 'o', color=color, markersize=0.5)
    
    for spot in autocrine_coords.itertuples():
        color = "#FF00FF"
        ax.plot(spot.x, spot.y, 'o', color=color, markersize=0.5)

    # add a legend
    if plot_interacting:
        colors = {'ligands (interacting)': "#FF0006", 'ligands (not interacting)': "#FF9B00", 'receptors (interacting)': "#0021ff", "receptors (not interacting)": "#00e7ff", "autocrine": "#FF00FF"}
    else:
        colors = {'ligands': "#FF0006", "receptors": "#0021ff", "autocrine": "#FF00FF"}

    handles = [mpatches.Patch(color=color, label=label) for label, color in colors.items()]

    plt.legend(handles=handles, loc='lower left', title="Key", bbox_to_anchor = (0, -0.6))
    plt.show()

# %%

