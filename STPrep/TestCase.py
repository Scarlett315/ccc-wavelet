#%%
from ResampleMatrix import *
from ImagePrep import *
from PIL import Image, ImageDraw
from IPython.display import display
import random
#%%
gene = pd.DataFrame(columns=["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P", "Q", "R", "S", "T", "U", "V", "W", "X", "Y", "Z"])
S = pd.DataFrame(columns=["barcode", "x", "y", "cell"])

alpha = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
def maketest(low, high, n):
    for i in range(n):
        S.loc[i] = [alpha[i], random.randint(low, high), random.randint(low, high), alpha[i]]
        gene.loc[i] = [random.randint(0, 5) for _ in range(26)]

maketest(30, 200, 24)
S.loc[25] = ["Z", 199, 199, "Z"]
gene["Z"] = 10
S.set_index("barcode", inplace=True)
print(S)


#%%
bounds = getBounds(S)
print(bounds)
resampled = resampleEfficient(gene, S, 100, bounds)
# %%
print(resampled)
print(resampled.shape)
plotResampledMatrix(resampled, "test")

# %%
#%% make an "image"
img = Image.new("RGB", (146, 865), (128, 128, 128))
draw = ImageDraw.Draw(img)
display(drawMatrix(resampled, 1))
display(img)
cropped = cropImage(img, S, 1)
display(cropped)
# %%
