#%%
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
import matplotlib.patches as mpatches
from GetTrueCCC import *
#%%
def visualize_overlap(lrp, coords, radius, plot_interacting = False, down_factor = 1):
    S, R = lrp.S, lrp.R
    ligand_coords = coords.loc[S.barcodes] # ligand coordinates ("center")
    receptor_coords = coords.loc[R.barcodes]

    ligand_coords["x"], ligand_coords["y"] = ligand_coords["y"], ligand_coords["x"]
    receptor_coords["x"], receptor_coords["y"] = receptor_coords["y"], receptor_coords["x"]

    fig, ax = plt.subplots(layout="constrained")
    ax.set_aspect("equal")
    ax.invert_yaxis()
    plt.title(f"{S.name} --> {R.name}, scaled: 1/{down_factor}")

    if lrp.secreted:
        radius *= 2
    if plot_interacting:
        S_int, R_int = lrp.S_int, lrp.R_int

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

# %%

