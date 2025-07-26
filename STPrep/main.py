import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from PIL import Image
import json
from ImagePrep import cropImage
from WaveletCoeffProcessing import readWaveletCoeffs, plotWavelets, resampleCoeffs, exportCoeffs
from ImageSegmentationMasks import createPolygons, drawMask, concatAnnotations
from ResampleMatrix import resampleEfficient, getBounds
import os

def resampleAllGenesMain(Y, S, D, bounds, pathExport):
    print("Rows in Y: ", len(Y))
    for index, r in Y.iterrows():
        name = Y["Unnamed: 0"][index]
        print(f"{index}: {name}")

        row = Y.loc[Y['Unnamed: 0'] == name]
        r = resampleEfficient(row, S, D, bounds)
        os.makedirs(pathExport, exist_ok=True)
        r.to_csv(f"{pathExport}/{name}_resampled_{str(D)}.csv")

def upsampleAllCoeffsMain(path, S, D, scaleFactor):
    num = 1
    for folder in os.listdir(path):
        folder_path = os.path.join(path, folder)
        geneName = folder.split("_")[0]
        if not os.path.isdir(folder_path) or folder.split("_")[1] != "resampled":
            continue  
        print(f"{num}: {geneName}")
        num += 1
        coeffs = readWaveletCoeffs(f"{path}/{folder}")
        resampled = resampleCoeffs(coeffs, S, D, scaleFactor)
        exportCoeffs(resampled, f"/Volumes/Samsung USB/UpsampledCoeffs/S1_wv22L2", geneName)
        
def main():
    sample = "S1"
    expression = pd.read_csv(f"../GSE208253/{sample}/info/expression_filtered.csv")
    coordinates = pd.read_csv(f"../GSE208253/{sample}/info/coordinates.csv", index_col=0)
    bounds = getBounds(coordinates)
    pathExport = f"../GSE208253/{sample}/process"

    scalePath = f"../GSE208253/{sample}/raw_data/spatial/scalefactors_json.json"
    with open(scalePath, 'r') as f:
        scales = json.load(f) 
    scaleFactor = scales["tissue_hires_scalef"]
    #resampleAllGenesMain(expression, coordinates, 128, bounds, f"{pathExport}/resampled")
    upsampleAllCoeffsMain(f"../GSE208253/{sample}/process/S1_wv22L2/", coordinates, 128, scaleFactor)

if __name__ == "__main__":
    main()