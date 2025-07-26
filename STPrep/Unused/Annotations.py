# %% 
# imports
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.pyplot as plt
import matplotlib.cm as cm

# %% 
a = pd.read_csv("GSE208253/export/pathologist_annotations.csv", index_col=0)
s = pd.read_csv("GSE208253/export/coordinates.csv", index_col=0)
print(a.head())
# %%
def annotationMatrix(annotations, spatial):
    annSpatial = np.zeros((spatial["y"].max(), spatial["x"].max()))
    annSpatial = annSpatial.astype(object)
    print(f"{spatial["x"].max()}, {spatial["y"].max()}") # 20871, 19884
    i = 0
    for row in annotations.itertuples():
        barcode = row.Index
        spatRow = spatial.loc[barcode]
        #print(row.pathologist_anno_x)

        annSpatial[spatRow["y"]][spatRow["x"]] = row.pathologist_anno_x
    return annSpatial
# %%
def concatAnnotations(annotations, spatial):
    return pd.concat([spatial, annotations], axis=1)

#%%
#test
a = a.rename(columns=lambda x: x.replace('.', '_'))
ann = annotationMatrix(a, s)# %%
concatAnnotations(a, s).to_csv("GSE208253/export/combined_annotations_spatial.csv")
# %%
