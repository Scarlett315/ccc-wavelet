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
sample = "S1"
geneTemp = "A2M"
expression = pd.read_csv(f"../GSE208253/{sample}/info/expression_filtered.csv")
coordinates = pd.read_csv(f"../GSE208253/{sample}/info/coordinates.csv", index_col=0)
imagePath = f"../GSE208253/{sample}/raw_data/spatial/tissue_hires_image.png"
image = Image.open(imagePath)
annotations = concatAnnotations(pd.read_csv(f"../GSE208253/{sample}/info/pathologist_annotations.csv", index_col=0), coordinates)
annotations = annotations.rename(columns=lambda x: x.replace('.', '_'))
colorFile = "../GSE208253/general/annotation_colors.json"
#path = f"GSE208253/export/{sample}"
bounds = getBounds(coordinates)


scalePath = f"../GSE208253/{sample}/raw_data/spatial/scalefactors_json.json"
with open(scalePath, 'r') as f:
    scales = json.load(f) 
scaleF = scales["tissue_hires_scalef"]

#%%
p = createPolygons(annotations, scaleF)
colorMask = drawMask(p, image, colorFile)
#colorMask.save(f"../GSE208253/{sample}/info/annotationPlot.png")
annotations.head()
overlayST(image, colorMask)
# %%
# Cropping Image
print("Bounds: ", bounds)
imgCropped = cropImage(image, coordinates, scaleF, False)
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
r = resampleEfficient(testRow, coordinates, 128, bounds)
plotResampledMatrix(r, "A2M", D=128)
#np.array_equal(r1, r2)

#r1.equals(r2)

#%%
u = upsampleToImage(r, coordinates, 128, scaleF, 1)
print(u.shape)
plotResampledMatrix(u, "Upsampled test")
# %%
# Wavelet Coefficients
coeffs = readWaveletCoeffs(f"../GSE208253/S1/process/S1_wv22L2/A2M_resampled_128/")
resampledCoeffs = resampleCoeffs(coeffs, coordinates, 128, scaleF)
print("Original Coefficients:")
for key, value in coeffs.items():
    print(f"{key}: {value.shape}")

#print("Resampled Coefficients:")
for key, value in resampledCoeffs.items():
    print(f"{key}: {value.shape}")
plotWavelets(resampledCoeffs, "A2M")
#drawMatrix(u).show()
#%%
np.set_printoptions(threshold=np.inf)
plotResampledMatrix(coeffs['L2_B00'], "Wavelet A2M")
test = upsampleToImage(coeffs['L2_B00'], coordinates, 128, scaleF, 4)
plotResampledMatrix(test, "Upsampled A2M")

#%%
print(f"Cropped Image Size: {imgCropped.size}")
imgCoeff = drawMatrix(resampledCoeffs['L2_B00'])
print(f"Upsampled Coefficient Shape: {resampledCoeffs['L2_B01'].shape}")
overlayST(imgCropped, imgCoeff)
#%%
print(resampledCoeffs['L2_B00'][62, 504])
# %%
m = mapBarcodes(coordinates, 1545, 1475)
m.to_csv(f"../GSE208253/{sample}/info/coordinates_mapped_test.csv")

# %%
