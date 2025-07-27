# %%imports
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from PIL import Image
import cv2
from ImagePrep import cropImage, overlayST
from IPython.display import display
from ResampleMatrix import getBounds
import os
import json
#%%
root = "../../Computational Biology Group Project/GSE208253/"
# reading files
def readFiles(num):
    sample = f"S{num}"
    #expression = pd.read_csv(f"../GSE208253/{sample}/info/expression_filtered.csv")
    coordinates = pd.read_csv(f"{root}/{sample}/info/coordinates.csv", index_col=0)
    imagePath = f"{root}/{sample}/raw_data/spatial/tissue_hires_image.png"
    image = Image.open(imagePath)

    scalePath = f"{root}/{sample}/raw_data/spatial/scalefactors_json.json"
    with open(scalePath, 'r') as f:
        scales = json.load(f) 
    scaleFactor = scales["tissue_hires_scalef"]
    return coordinates, image, scaleFactor
#%%
path = f"{root}/general/masks"
for file in os.listdir(path):
    if file.endswith(".png"):
        sample = str(file.split("_")[1].split(".")[0])
        coords, img, scaleF= readFiles(sample)
        print(file)
        mask = Image.open(os.path.join(path, file))
        croppedMask = cropImage(mask, coords, scaleF)
        croppedImg = cropImage(img, coords, scaleF)
        
        display(overlayST(croppedImg, croppedMask))
        croppedMask.save(f"{root}/general/masks/cropped/{file.replace('.png', '.png')}")
        croppedImg.save(f"{root}/general/images/S{sample}_cropped.png")
# %%
