#%% 
# imports
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from PIL import Image
import json
from ImagePrep import cropImage
from ResampleMatrix import getBounds, resampleMatrix, plotResampledMatrix, upsampleToImage, drawMatrix, resampleAllGenes
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
path = f"../GSE208253/{sample}/process"
bounds = getBounds(coordinates)


scalePath = f"../GSE208253/{sample}/raw_data/spatial/scalefactors_json.json"
with open(scalePath, 'r') as f:
    scales = json.load(f) 
scaleF = scales["tissue_hires_scalef"]
#%% Single Gene Downsampling (test)
r = expression.loc[expression['Unnamed: 0'] == "A2M"]
r = resampleMatrix(r, coordinates, 128, bounds)
r.to_csv(f"{path}/A2M_resampled_128.csv")
plotResampledMatrix(r, "A2M")
#%%
# Gene Downsampling
D = 128
resampleAllGenes(expression, scaleF, D, bounds, f"{path}/expression_downsampled_{D}")
#%%
# Wavelet Coefficients Upsampling
#%%
# Image Cropping & Test Alignment
