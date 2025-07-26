
# %%
#imports
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy.sparse as scisparse
from DataFunctions import plotExpMatrix
from ResampleMatrix import plotResampledMatrix

# %%
# reading data
expression = pd.read_csv("Her2ST_data/count-matrices/A1.tsv", sep='\t', index_col=0)
coordinates = pd.read_csv("Her2ST_data/spot-selections/A1_selection.tsv", sep='\t')
# %% I used this to make "coordinates_fixed" because the format was different between the files
print("Expression")
print(expression.info)
print("Coordinates")
print(coordinates.info)
# %% 
# getting 1 gene for testing
selectedGene = "SAMD11"
#print(expression.columns)
mask = (expression.columns.str.contains(selectedGene))
selected = expression.loc[:, mask]

#set barcodes as index
#selectedRow = selectedRow.squeeze()
#selectedRow.index.name = None
selected.head(10)
#%% 
# create spatial expression matrix for the gene
genes = np.zeros((30, 30), dtype=int)

## change so that you use the expression matrix only
for barcode in coordinates.itertuples():
    x = barcode.x
    y = barcode.y
    spot_id = f"{x}x{y}"
    if spot_id in selected.index:
        genes[x][y] = selected.loc[spot_id].values[0]
    else:
        print(f"{barcode} not found in selectedRow.")

# %%
genes = pd.DataFrame(genes)
genes.head(30)
# %%
plotResampledMatrix(genes, selectedGene)
# %%
# Conversion to sparse matrix and export
#geneMatrix.to_csv('DCN_expression_matrix.csv', index=False)
sparseGeneMatrix = genes.astype(pd.SparseDtype("int", 0))

scisparse.save_npz(f"{selectedGene}_expression_matrix", sparseGeneMatrix.sparse.to_coo() , compressed=True)