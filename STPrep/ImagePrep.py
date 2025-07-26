#%% 
# imports
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from PIL import Image
import cv2
from ResampleMatrix import getBounds
from PIL import ImageDraw
# %%
def imgAsMat(path):
    img = Image.open(path)
    matrix = np.asarray(img)

    # separate channels
    red = pd.DataFrame(matrix[:,:,0])
    green = pd.DataFrame(matrix[:,:,1])
    blue = pd.DataFrame(matrix[:,:,2])
    return red.to_numpy(), green.to_numpy(), blue.to_numpy() 

def matToImg(red, green, blue):
    matrix = np.dstack((red, green, blue))
    return Image.fromarray(matrix)

def showChannels(red, green, blue): #probably not very useful but fun to look at
    matrix = np.vstack((red, green, blue))
    return Image.fromarray(matrix)

def normalizeChannels(red, green, blue):
    nRed = cv2.normalize(red, None, alpha=0, beta=255, norm_type=cv2.NORM_MINMAX)
    nGreen = cv2.normalize(green, None, alpha=0, beta=255, norm_type=cv2.NORM_MINMAX)
    nBlue = cv2.normalize(blue, None, alpha=0, beta=255, norm_type=cv2.NORM_MINMAX)
    return red, green, blue

#%%
def cropImage(img, S, scaleFactor, flip=False, padding=0):
    S = S.copy()
    def swapAxes(spatial):
        spatial['x'], spatial['y'] = spatial['y'], spatial['x']
        return spatial

    def toPixels(spot):
        pixel = int(spot * scaleFactor)
        return pixel
    
    S = swapAxes(S)
    S['x'] = S['x'].apply(toPixels)
    S['y'] = S['y'].apply(toPixels)
    xMin, xMax, yMin, yMax = getBounds(S)

    #print(xMin, xMax, yMin, yMax)
    shape = [xMin - padding, yMin - padding, xMax + padding, yMax + padding]
    print(getBounds(S))
    cropped = img
 
    if flip:
        cropped = img.transpose(Image.FLIP_LEFT_RIGHT)
        cropped = cropped.crop(shape) 
        cropped = cropped.transpose(Image.FLIP_LEFT_RIGHT)
    else:
        cropped = cropped.crop(shape)
    
    return cropped
#%%
def overlayST(img, STImg, origin = (0,0)):
    img = img.copy()
    STImg = STImg.copy()
    img = img.convert("RGBA")
    STImg.putalpha(150)
    img.paste(STImg, origin, STImg)
    return img
# %%
# example
#r, g, b = imgAsMat("GSE208253/export/tissue_hires_image_cropped.png")

#nR, nG, nB = normalizeChannels(r, g, b)
#Image.fromarray(nR).save("GSE208253/export/image/normalized_red.png")
#Image.fromarray(nG).save("GSE208253/export/image/normalized_green.png")
#Image.fromarray(nB).save("GSE208253/export/image/normalized_blue.png")

#matToImg(nR, nG, nB).save("GSE208253/export/image/normalized_all.png")
#matToImg(nR, nG, nB).show()
#img = showChannels(nR, nG, nB)
#img.show()