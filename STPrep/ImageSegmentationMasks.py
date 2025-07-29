#%%import cv2
import numpy as np
from PIL import Image, ImageDraw
import pandas as pd
from IPython.display import display
import json

#%%
# read color file
def readColorFile(path):
    with open(path, 'r') as f:
        color_dict = json.load(f)
    labels = color_dict["labels"]
    colors = color_dict["colors"]
    for key, value in colors.items():
        colors[key] = tuple(value)
    return labels, colors

#%% creating polygons from files
def createPolygons(a, scaleFactor, label=None):
    if label is None:
        labels = a["pathologist_anno_x"].unique()
        col_name = "pathologist_anno_x"
    else:
        labels = a[label].unique()
        col_name = label
    print(labels)
    polygons = {}
    for l in range(len(labels)):
        #print(labels[l])
        polygons[labels[l]] = []

    for annotation in a.itertuples():
        polygons[getattr(annotation, col_name)].append((annotation.x * scaleFactor, annotation.y * scaleFactor))
    return polygons
#%%
def drawMask(polygons, img, colorFile, radius=2, noLabel = False):
    labels, colors = readColorFile(colorFile)
    if noLabel:
        labels = {label: label for label in polygons.keys()}
        print(labels)

    width, height = img.size
    imgSize = (int(width), int(height))
    print(imgSize)

    # create mask
    colorMask = Image.new("RGB", imgSize, (0, 0, 0))
    labelMask = Image.new("L", imgSize, 0) # blank label mask (ML)
    
    # Drawing setup
    draw_color = ImageDraw.Draw(colorMask)
    draw_label = ImageDraw.Draw(labelMask)

    labelNum = 0
    pointsPlotted = 0
    for label, points in polygons.items():
        if noLabel:
            color = list(colors)[labelNum]
        else:
            color = labels[label][1]
        print(f"Drawing label '{label}' with {len(points)} points, color {color}")
        for x, y in points:
            bbox = (x - radius, y - radius, x + radius, y + radius)
            draw_color.ellipse(bbox, fill=colors[color])
            #draw_label.ellipse(bbox, fill=colors[color])
            pointsPlotted += 1
        labelNum += 1
    print(f"Plotted {pointsPlotted} points")


    colorMask = colorMask.transpose(Image.FLIP_LEFT_RIGHT)
    colorMask = colorMask.rotate(90)

    display(colorMask)
    return colorMask
    #display(labelMask)

# %%
def concatAnnotations(annotations, spatial):
    def strip_suffix(barcode):
        return str(barcode).rsplit('-', 1)[0]
    spatial = spatial.copy()
    annotations = annotations.copy()
    print(spatial.index)
    spatial.index = spatial.index.map(strip_suffix)
    annotations.index = annotations.index.map(strip_suffix)
    return replaceNaN(pd.concat([spatial, annotations], axis=1))

# %%
def replaceNaN(annotations):
    annotations = annotations.copy()
    annotations = annotations.fillna("NaN")
    return annotations

#%% 
#test
#p = createPolygons(annotations, 0.087)
#drawMask(p, image)
# %%
