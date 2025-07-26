#%%
#Imports
import numpy as np
import pandas as pd
import json
from ImagePrep import *
from ResampleMatrix import *
from WaveletCoeffProcessing import *

#%%
sample = "S1"
geneTemp = "A2M"
expression = pd.read_csv(f"../GSE208253/{sample}/info/expression_filtered.csv", index_col=0)
coordinates = pd.read_csv(f"../GSE208253/{sample}/info/coordinates.csv", index_col=0)
colorFile = "../GSE208253/general/annotation_colors.json"
#path = f"GSE208253/export/{sample}"

scalePath = f"../GSE208253/{sample}/raw_data/spatial/scalefactors_json.json"
with open(scalePath, 'r') as f:
    scales = json.load(f) 
scaleF = scales["tissue_hires_scalef"]

imagePath = f"../GSE208253/{sample}/raw_data/spatial/tissue_hires_image.png"
image = Image.open(imagePath)
imgCropped = cropImage(image, coordinates, scaleF, False)

#%%
coordsInPixels = coordinates.copy()
coordsInPixels["x"] = coordsInPixels["x"] * scaleF
coordsInPixels["y"] = coordsInPixels["y"] * scaleF
coordsInPixels.head()
# %%
# %%
# Wavelet Coefficients
coeffs = readWaveletCoeffs(f"../GSE208253/S1/process/S1_wv22L2/A2M_resampled_128/")
resampledCoeffs = resampleCoeffs(coeffs, coordinates, 128, scaleF)
print(resampledCoeffs.keys())
#drawMatrix(u).show()

# %%
testCoeff = upsampleToImage(coeffs["L2_B01"], coordinates, 128, scaleF, 4)
#%%
def calculateNewCoords(x, y, xMin, xMax, yMin, yMax, imgLength, imgHeight):
    x_rel = (x - xMin) / (xMax - xMin)
    y_rel = (y - yMin) / (yMax - yMin)
    x_img = int(np.clip(round(x_rel * (imgLength-1)), 0, imgLength-1))
    y_img = int(np.clip(round(y_rel * (imgHeight-1)), 0, imgHeight-1))
    return x_img, y_img
#%%
count = 0
coords_new = pd.DataFrame(columns=["barcode", "x", "y"])
for barcode in coordinates1.itertuples():
    x_img, y_img = calculateNewCoords(barcode.x, barcode.y, xMin, xMax, yMin, yMax, resLength, resHeight)
    value = testCoeff[y_img, x_img]
    count += 1

    newRow = pd.DataFrame([{"barcode": barcode.Index, "x": x_img, "y": y_img}])
    coords_new = pd.concat([coords_new, newRow], ignore_index=True)

    print(f"Barcode: {barcode.Index}, x: {barcode.x}, y: {barcode.y}, Matrix indices: ({x_img}, {y_img}), Value: {value}")
print(f"Number of barcodes: {count}")
coords_new.to_csv(f"../GSE208253/{sample}/info/coordinates_img.csv")

# %%
coords_new = mapBarcodes(coordinates, resLength, resHeight)
coords_new.to_csv(f"../GSE208253/{sample}/info/coordinates_img_2.csv")
# %%