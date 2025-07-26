# %% 
# imports
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.ndimage import zoom
from PIL import Image, ImageDraw

# %%
def getBounds(S):
    '''
    Gets the "rectangle" that can cover all locations
    S: the matrix mapping barcodes to locations with columns in the form [header, "x", "y", "cell"]
    '''
    xMin = np.min(S["x"])
    xMax = np.max(S["x"]) + 1
    yMin = np.min(S["y"])
    yMax = np.max(S["y"]) + 1
    return xMin, xMax, yMin, yMax
# %% 
def partitionMatrix(D, bounds):
    xMin, xMax, yMin, yMax = bounds[0], bounds[1], bounds[2], bounds[3]
    xPartLength = (xMax - xMin) / D
    yPartLength = (yMax - yMin) / D
    return np.zeros(shape=(D, D)), xPartLength, yPartLength

#%% 
def resampleMatrix(gene, S, D, bounds):
    '''
    Resamples the matrix to create an evenly-spaced representation
    in matrix form for the wavelet transform.

    Arguments:
        S: spatial matrix 
        D (int): size of gene expression matrix
        bounds [int list]: [xMin, xMax, yMin, yMax]
    '''
    xMin, yMin = bounds[0],  bounds[2]
    resampled, xPartLength, yPartLength = partitionMatrix(D, bounds)
    filled = 0
    for x in range(D): # computing local averages for each spot
        for y in range(D):
            xLow = xMin + x * xPartLength
            yLow = yMin + y * yPartLength
            xHigh = xMin + (x + 1) * xPartLength
            yHigh = yMin + (y + 1) * yPartLength

            boxMask = ((S["x"] >= xLow) & (S["x"] < xHigh) 
                    & (S["y"] >= yLow) & (S["y"] < yHigh))
            coordsInBox = S[boxMask]

            #print(coordsInBox["cell"])
            
            barcodesInBox = coordsInBox["cell"].values

            # getting expression values
            expVals = []
            for b in barcodesInBox:
                if b in gene.columns:
                    expVals.append(gene[b].values[0])  # 1 row
            if expVals:
                resampled[y, x] = np.mean(expVals)
                filled += 1
    print(f"Filled {filled} out of {D*D} cells ({(filled / (D*D)) * 100:.2f}%)")

    # Transform the final matrix to match image orientation
    # Rotate 90° counterclockwise and flip over x-axis
    #transformed = np.flipud(resampled.T)
    #transformed = Image.fromarray(resampled)
    #transformed = transformed.transpose(Image.FLIP_LEFT_RIGHT)
    #transformed = transformed.rotate(90)
    #transformed = np.asarray(transformed)

    print(resampled.shape)
    return pd.DataFrame(resampled)
#%%
def resampleEfficient(gene, S, D, bounds):
    xMin, xMax, yMin, yMax = bounds[0], bounds[1], bounds[2] , bounds[3]
    sums = np.zeros(shape=(D, D))
    nums = np.zeros(shape=(D, D))
    for row in S.itertuples():
        x = row.x
        y = row.y
        expression = gene[row.cell].iloc[0]
        # calculate the downsampled coordinates of x and y
        #res_x = int(D * (x - xMin) / (xMax - xMin))
        #res_y = int(D * (y - yMin) / (yMax - yMin))
        res_x, res_y = calcDownsampledCoords(x, y, D, bounds)
        

        sums[res_y, res_x] += expression
        nums[res_y, res_x] += 1

    nums[nums == 0] = 1
    resampled = sums / nums

    # Transform the final matrix to match image orientation
    # Rotate 90° counterclockwise and flip over x-axis
    #transformed = np.flipud(resampled.T)
    #transformed = Image.fromarray(resampled)
    #transformed = transformed.transpose(Image.FLIP_LEFT_RIGHT)
    #transformed = transformed.rotate(90)
    #transformed = np.asarray(transformed)

    return pd.DataFrame(resampled)

# %%
def resampleAllGenes(Y, S, D, bounds, pathExport):
    for index, r in Y.iterrows():
        name = Y["Unnamed: 0"][index]
        print(name)

        row = Y.loc[Y['Unnamed: 0'] == name]
        r = resampleMatrix(row, S, D, bounds)
    
        r.to_csv(f"{pathExport}/{name}_resampled_{str(D)}.csv")
    
# %%
def plotResampledMatrix(geneMatrix, geneName, ax=None, title=None, vmin=None, vmax=None, s=50, show_colorbar=True, D = None):
    nonzero_y, nonzero_x = np.nonzero(geneMatrix)
    values = np.asarray(geneMatrix)[nonzero_y, nonzero_x]
    if ax is None:
        fig, ax = plt.subplots(figsize=(10, 10))
        show_plot = True
    else:
        show_plot = False
    sc = ax.scatter(nonzero_x, nonzero_y, c=values, cmap='Blues', s=s, vmin=vmin, vmax=vmax)
    if show_colorbar:
        plt.colorbar(sc, label='Value', ax=ax)
    ax.set_title(title if title else f'Resampled Gene Expression Matrix for {geneName}')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.xaxis.tick_top()
    ax.xaxis.set_label_position('top')
    if D:
        ax.set_xlim(0, D) 
        ax.set_ylim(D, 0)
    else:
        ax.set_ylim(bottom=0)
        ax.set_xlim(left=0)
        ax.invert_yaxis()

    ax.set_aspect('equal')
    if show_plot:
        plt.tight_layout()
        plt.show()
    return sc
#%% 
def drawMatrix(geneMatrix, radius=3):
    nonzero_y, nonzero_x = np.nonzero(geneMatrix)
    imgSize = (int(geneMatrix.shape[1]), int(geneMatrix.shape[0]))

    img = Image.new("L", imgSize, 100)
    draw = ImageDraw.Draw(img)
    
    for x, y in zip(nonzero_x, nonzero_y):
        bbox = [x - radius, y - radius, x + radius + 1, y + radius + 1]
        draw.ellipse(bbox, fill=255)
    #img = img.transpose(Image.FLIP_LEFT_RIGHT)
    #img = img.rotate(180)
    return img
#%%
def upsampleToImage(geneMatrix, S, D, scaleFactor, wv=1):
    '''
    Upsamples gene matrix to the scale of the image and removes barcodes
    that aren't there in the original
    
    Args:
        geneMatrix: The resampled gene matrix (DxD)
        S: Spatial coordinates dataframe
        D: Downsampled matrix dimension
        scaleFactor: Scale factor for upsampling
        bounds: returned by getBounds(S)
    '''
    geneMatrix = np.asarray(geneMatrix)
    print("Nonzero values before mapping: ", np.count_nonzero(geneMatrix))
    def swapAxes(spatial):
        spatial['x'], spatial['y'] = spatial['y'], spatial['x']
        return spatial

    def toPixels(spot):
        pixel = int(spot * scaleFactor)
        return pixel

    S_orig = S.copy()
    S = swapAxes(S.copy())
    S["x"] = S["x"].apply(toPixels)
    S["y"] = S["y"].apply(toPixels)

    bounds = getBounds(S)
    print(f"Bounds: {bounds}")

    _, xPartLength, yPartLength = partitionMatrix(D, bounds)    
    imgLength = int(xPartLength * D) # length of orig in pixels
    imgHeight = int(yPartLength * D) # height of orig in pixels

    origLength = imgHeight
    origHeight = imgLength

    print(f"length, height: {imgLength}, {imgHeight}")

    upsampled = np.zeros((imgHeight, imgLength))


    origBounds = getBounds(S_orig)
    print("BOUNDS: ", origBounds)

    for _, row in S_orig.iterrows():
        print(f"Barcode x, y: {row['x']}, {row['y']}")
        x_downsampled, y_downsampled = calcDownsampledCoords(row["x"], row["y"], D, origBounds, wv=wv)
        x_img, y_img = calculateImgCoords(row["x"], row["y"], origBounds, origLength, origHeight)
        print(f"x_img, y_img: {x_img}, {y_img}")
        print(f"x_downsampled, y_downsampled: {x_downsampled}, {y_downsampled}")
        print(f"geneMatrix value: {geneMatrix[y_downsampled, x_downsampled]}")
        upsampled[y_img, x_img] = geneMatrix[y_downsampled, x_downsampled]

        print(f"name: {row['cell']}; Barcode x, y: {row['x']}, {row['y']}; Upsampled x, y: {x_img}, {y_img}")
    print(f"Nonzero values after mapping: {np.count_nonzero(upsampled)}")

    #plotResampledMatrix(upsampled, "hi")
    return upsampled
    
#%%
def calculateImgCoords(x, y, bounds, imgLength, imgHeight):
    xMin, xMax, yMin, yMax = bounds[0], bounds[1], bounds[2], bounds[3]
    x_rel = (x - xMin) / (xMax - xMin)
    y_rel = (y - yMin) / (yMax - yMin)
    #x_img = int(np.clip(round(x_rel * (imgLength-1)), 0, imgLength-1))
    #y_img = int(np.clip(round(y_rel * (imgHeight-1)), 0, imgHeight-1))
    x_img = int(x_rel * (imgLength-1))
    y_img = int(y_rel * (imgHeight-1))

    #x_img, y_img = -y_img + imgHeight - 1, -x_img + imgLength - 1
    x_img, y_img = y_img, x_img
    return x_img, y_img

#%%
def calcDownsampledCoords(x, y, D, bounds, wv=1):
    xMin, xMax, yMin, yMax = bounds[0], bounds[1], bounds[2], bounds[3]
     # calculate the downsampled coordinates of x and y
    res_x = int(D * (x - xMin) / (wv * (xMax - xMin)))
    res_y = int(D * (y - yMin) / (wv * (yMax - yMin)))

    # Transform: rotate 90° counterclockwise and flip over x-axis
    # Combined transformation: (x, y) → (y, x)
    x_transformed = res_y
    y_transformed = res_x
    return x_transformed, y_transformed
#%%
def mapBarcodes(coordinates, imgLength, imgHeight):
    coords_new = pd.DataFrame(columns=["barcode", "x", "y"])
    coordinates = coordinates.copy()
    coordinates['x'], coordinates['y'] = coordinates['y'], coordinates['x']
    xMin, xMax, yMin, yMax = getBounds(coordinates)
    for barcode in coordinates.itertuples():
        x_img, y_img = calculateImgCoords(barcode.x, barcode.y, xMin, xMax, yMin, yMax, imgLength, imgHeight)
        newRow = pd.DataFrame([{"barcode": barcode.Index, "x": x_img, "y": y_img}])
        coords_new = pd.concat([coords_new, newRow], ignore_index=True)
    coords_new.set_index("barcode", inplace=True)
    return coords_new