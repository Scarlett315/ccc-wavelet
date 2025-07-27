import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from PIL import Image
import json
from ImagePrep import cropImage
from WaveletCoeffProcessing import readWaveletCoeffs, plotWavelets, resampleCoeffs, exportCoeffs
from ImageSegmentationMasks import createPolygons, drawMask, concatAnnotations
from ResampleMatrix import resampleEfficient, getBounds, upsampleToImage
import os

#I have data structured as GSE208253 / 
#     - sample (one for each)
#         - info (expression + coordinate matrices, annotations)
#         - process (wavelet coefficients)
#         - raw_data (scalefactor, high-res image, other raw data processed in R/Seurat)
#     - general (annotation color key, cropped images and masks)

def readFiles(sample, dataPath):
    expression = pd.read_csv(f"{dataPath}/{sample}/info/expression_filtered.csv")
    coordinates = pd.read_csv(f"{dataPath}/{sample}/info/coordinates.csv", index_col=0)
    imagePath = f"{dataPath}/{sample}/raw_data/spatial/tissue_hires_image.png"
    image = Image.open(imagePath)
    annotations = concatAnnotations(pd.read_csv(f"{dataPath}/{sample}/info/pathologist_annotations.csv", index_col=0), coordinates)
    annotations = annotations.rename(columns=lambda x: x.replace('.', '_'))
    colorFile = f"{dataPath}/general/annotation_colors.json"
    bounds = getBounds(coordinates)

    scalePath = f"{dataPath}/{sample}/raw_data/spatial/scalefactors_json.json"
    with open(scalePath, 'r') as f:
        scales = json.load(f) 
    scaleF = scales["tissue_hires_scalef"]

    return expression, coordinates, image, annotations, colorFile, scaleF


def resampleAllGenesMain(Y, S, D, bounds, pathExport):
    print("Rows in Y: ", len(Y))
    for index, r in Y.iterrows():
        name = Y["Unnamed: 0"][index]
        print(f"{index}: {name}")

        row = Y.loc[Y['Unnamed: 0'] == name]
        r = resampleEfficient(row, S, D, bounds)
        os.makedirs(pathExport, exist_ok=True)
        r.to_csv(f"{pathExport}/{name}_resampled_{str(D)}.csv")

def upsampleAllCoeffsMain(pathToCoeffs, S, D, scaleFactor, exportPath):
    num = 1
    for folder in os.listdir(pathToCoeffs):
        folder_path = os.path.join(pathToCoeffs, folder)
        geneName = folder.split("_")[0]
        if not os.path.isdir(folder_path) or folder.split("_")[1] != "resampled":
            continue  
        print(f"{num}: {geneName}")
        
        coeffs = readWaveletCoeffs(f"{pathToCoeffs}/{folder}")
        resampled = resampleCoeffs(coeffs, S, D, scaleFactor)
        if num == 1:  #redoing upsampling for a random coefficient to get the mapping (checked to be equal across genes & coeffs)
            print(f"Reading from: {pathToCoeffs}")
            print(f"Saving to: {exportPath}")
            _, mapping = upsampleToImage(coeffs['L2_B00'], S, D, scaleFactor, wv=4, exportMapping=True)
            mapping.to_csv(f"{exportPath}/coordinates_img.csv")

        exportCoeffs(resampled, exportPath, geneName)
        num += 1
        
def main():
    sample = "S1"
    dataPath = f"../Computational Biology Group Project/GSE208253" #replace with your data path

    expression, coordinates, image, annotations, colorFile, scaleFactor = readFiles(sample, dataPath)
    bounds = getBounds(coordinates)
    pathExport = f"/Volumes/Samsung USB/UpsampledCoeffs/S1_wv22L2"
    #resampleAllGenesMain(expression, coordinates, 128, bounds, f"{pathExport}/resampled")
    #upsampleAllCoeffsMain(f"{dataPath}/{sample}/process/S1_wv22L2/", coordinates, 128, scaleFactor, pathExport)
    test_orig = np.load(f"{dataPath}/{sample}/process/S1_wv22L2/PPL_resampled_128/L1_B10.npy")
    test_up = pd.read_csv(f"{pathExport}/PPL_upsampled_coeffs/PPL_L1_B10.csv")
    print(np.count_nonzero(test_up))
    print(test_up.shape)
    print(np.count_nonzero(test_orig))
    print(test_orig.shape)

if __name__ == "__main__":
    main()