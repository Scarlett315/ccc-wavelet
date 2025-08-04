#%% 
# imports
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from PIL import Image
import json
from IPython.display import display
from ImagePrep import cropImage, overlayST
from ResampleMatrix import *
from WaveletCoeffProcessing import readWaveletCoeffs, plotWavelets, resampleCoeffs, exportCoeffs
from ImageSegmentationMasks import createPolygons, drawMask, concatAnnotations
# %%
# reading files
# I have data structured as GSE208253 / 
#     - sample (one for each)
#         - info (expression + coordinate matrices, annotations)
#         - process (wavelet coefficients)
#         - raw_data (scalefactor, high-res image, other raw data processed in R/Seurat)
#     - general (annotation color key, cropped images and masks)

sample = "S1"
geneTemp = "A2M"
dataPath = f"../../Computational Biology Group Project/GSE208253/" #replace with your data path
expression = pd.read_csv(f"{dataPath}/{sample}/info/expression_filtered.csv")
coordinates = pd.read_csv(f"{dataPath}/{sample}/info/coordinates.csv", index_col=0)
imagePath = f"{dataPath}/{sample}/raw_data/spatial/tissue_hires_image.png"
image = Image.open(imagePath)
annotations = concatAnnotations(pd.read_csv(f"{dataPath}/{sample}/info/pathologist_annotations.csv", index_col=0), coordinates)
annotations = annotations.rename(columns=lambda x: x.replace('.', '_'))
colorFile = f"{dataPath}/general/annotation_colors.json"
#path = f"GSE208253/export/{sample}"
bounds = getBounds(coordinates)


scalePath = f"{dataPath}/{sample}/raw_data/spatial/scalefactors_json.json"
with open(scalePath, 'r') as f:
    scales = json.load(f) 
scaleF = scales["tissue_hires_scalef"]

#%%
p = createPolygons(annotations, scaleF, label="cluster_annotations")
print(p)
colorMask = drawMask(p, image, colorFile, noLabel=True)
#colorMask.save(f"../GSE208253/{sample}/info/annotationPlot.png")
annotations.head()
overlayST(image, colorMask)
# %%
# Cropping Image
print("Bounds: ", bounds)
imgCropped = cropImage(image, coordinates, scaleF)
#maskCropped = cropImage(colorMask, coordinates, scaleF, False)
#imgCropped.save(f"../GSE208253/general/images/S{sample}_cropped.png")
display(imgCropped)
print(imgCropped.size)

# %%
# Resampling Expression
selectedGene = "A2M"
selectedRow = expression.loc[expression['Unnamed: 0'] == selectedGene]
testRow = pd.DataFrame(columns=expression.columns) #visualize all barcodes for testing
testRow.loc[0] = 3
print(testRow)
r = resampleEfficient(testRow, coordinates, 32, bounds)
plotResampledMatrix(r, "A2M", D=32)
#np.array_equal(r1, r2)

#r1.equals(r2)

#%%
u = upsampleToImage(r, coordinates, 128, scaleF, wv=4)
print(u.shape)
plotResampledMatrix(u, "Upsampled test")
img_as_arr = np.array(imgCropped)
print(f"Upsampled Matrix Shape: {u.shape}")
print(f"Cropped Image Size: {imgCropped.size}")
print(f"Cropped Image As Array Shape: {img_as_arr.shape}")
# %%
# Wavelet Coefficients
coeffs = readWaveletCoeffs(f"../GSE208253/S1/process/S1_wv22L2/A2M_resampled_128/")
resampledCoeffs = resampleCoeffs(coeffs, coordinates, 128, scaleF)
plotWavelets(resampledCoeffs, "A2M")
#drawMatrix(u).show()
#%%
#plotResampledMatrix(coeffs['L2_B00'], "Wavelet A2M")
up, map = upsampleToImage(coeffs['L2_B00'], coordinates, 128, scaleF, 4, exportMapping=True)
#plotResampledMatrix(up, "Upsampled A2M")
print(up.shape)
for y, x in zip(np.nonzero(up)[0], np.nonzero(up)[1]):
    print(f"x: {x}, y: {y}")
    mask = map[(map["x"] == x) & (map["y"] == y)]
    print(mask)

print(up[196, 595])
print(up.shape)

#ap.to_csv(f"test.csv")
#%%
print(f"Cropped Image Size: {imgCropped.size}")
imgCoeff = drawMatrix(resampledCoeffs['L2_B00'])
print(f"Upsampled Coefficient Shape: {resampledCoeffs['L2_B01'].shape}")
#%%
print(resampledCoeffs['L2_B00'][62, 504])
# %%
m = mapBarcodes(coordinates, 1545, 1475)
m.to_csv(f"../GSE208253/{sample}/info/coordinates_mapped_test.csv")

# %%
print(f"Cropped Image Size: {imgCropped.size}")
print(f"upsampled shape: {u.shape}")
overlayST(imgCropped, drawMatrix(u))

#%%
e = pd.read_csv("/Users/scarlett/Programming/ccc-wavelet/Data/human_breast_cancer/info/expression_filtered.csv", index_col=0)
c = pd.read_csv("/Users/scarlett/Programming/ccc-wavelet/Data/human_breast_cancer/info/coordinates.csv", index_col=0)
#a = resampleEfficient(e.loc["BGN"], c, 128, getBounds(c))
a = pd.read_csv("/Users/scarlett/Programming/ccc-wavelet/Data/human_breast_cancer/process/resampled_128/BGN_resampled_128.csv", index_col=False)
print(a.head())

a = a.apply(pd.to_numeric, errors='coerce')
a = a.to_numpy()
print(a.shape)
plotResampledMatrix(a, "hi", D=128)

# %%
