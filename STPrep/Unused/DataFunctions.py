# %%
#imports
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy.sparse as scisparse
# %%
# Read files & reconstruct matrix
def readExpFromNpz(filePath):
    load = np.load(filePath, allow_pickle = True)
    # Extract keys
    row = load["row"]
    col = load["col"]
    data = load["data"]
    shape = tuple(load["shape"])  

    sparse = scisparse.coo_matrix((data, (row, col)), shape=shape)
    dense = sparse.toarray()
    return pd.DataFrame(dense)

#%%
def plotExpMatrix(geneMatrix, geneName):
    nonzero_y, nonzero_x = np.nonzero(geneMatrix)
    values = geneMatrix.values[nonzero_y, nonzero_x]
    # Nonzero points only (change later?)
    plt.figure(figsize=(10, 10))
    sc = plt.scatter(
        nonzero_x,
        nonzero_y,
        c=values,
        cmap="plasma",
        s=15,             
        edgecolors='none',
        vmin=0.01          
    )
    plt.colorbar(sc, label="Expression Level")
    plt.axis("equal")

    #inverting x-axis to match picture
    plt.ylim(plt.ylim()[::-1])

    plt.title(f"Spatial Gene Expression (for gene {geneName})")
    plt.show()

# %%

# %%
# example
#DCN = readExpFromNpz("DCN_expression_matrix.npz")
#DCN.head()
#plotExpMatrix(DCN, "DCN")