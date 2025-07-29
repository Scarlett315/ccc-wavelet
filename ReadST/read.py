#%%
import anndata as ad
import scanpy as sc
from scipy.io import mmwrite


# Load the .h5ad file
adata = sc.read_h5ad("../Data/human_lung/human_lung.h5ad")
# %%
mmwrite("../Data/human_lung/matrix.mtx", adata.X) # counts matrix
adata.obs.to_csv("../Data/human_lung/obs.csv") # cell metadata
adata.var.to_csv("../Data/human_lung/var.csv") # gene metadata
# %%
print(adata.var)
# %%
