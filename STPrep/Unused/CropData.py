#%%
# imports
import pandas as pd
from PIL import Image
import json
#%% 
# read files
#coordinates = pd.read_csv("coordinates.csv")
#imagePath = "Sampl#eSTData/spatial/tissue_hires_image.png"
#image = Image.open(imagePath)

# %%
# cropping image
cropInPixels = [397, 415, 1730, 1700] # you have to manually do this
    # top left x, top left y, bottom right x, bottom right y
#croppedImg = image.crop((cropInPixels[0], cropInPixels[1], cropInPixels[2], cropInPixels[3])) 
#croppedImg.show()
# %%
# function for cropping ST
def cropSTData(spatial, cropPixels, scaleFactor):
    cropSpots = [x / scaleFactor for x in cropPixels]

    croppedST = spatial[
    (spatial["x"] >= cropSpots[0]) & (spatial["x"] <= cropSpots[2]) &
    (spatial["y"] >= cropSpots[1]) & (spatial["y"] <= cropSpots[3])]
    return croppedST
#%% 
#test
#scalePath = "SampleSTData/spatial/scalefactors_json.json"
#with open(scalePath, 'r') as f:
    #scales = json.load(f) 
#print(scales["tissue_hires_scalef"])
#coordinatesCropped = cropSTData(coordinates, cropInPixels, scales["tissue_hires_scalef"])
#coordinatesCropped.to_csv("coordinates_cropped.csv")

# %%
