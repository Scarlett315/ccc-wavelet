#%%
# they didn't provide tissue_lowres_image so i just had to make them 
# but it doeesn't actually correspond with the scalefactor so it's kinda useless
# it's literally just b/c Seurat doesn't accept the folder if it doesn't have the img for some reason
from PIL import Image
import os
#%%
def compress_image(input_path, output_path, optimize=True):
    """
    Compresses an image using Pillow.

    Args:
        input_path (str): Path to the original image file.
        output_path (str): Path to ssave the compressed image.
        quality (int): JPEG quality setting (0-100). Lower values mean more compression.
        optimize (bool): If True, perform an extra optimization pass.
    """
    try:
        img = Image.open(input_path)
    
        if img.mode == 'RGBA' and output_path.lower().endswith('.jpg'):
            # Convert RGBA to RGB for JPEG saving to avoid OSError
            img = img.convert('RGB')

        img = img.convert("P", palette=Image.ADAPTIVE, colors=10)
        img.save(f'{output_path}/tissue_lowres_image.png', optimize=optimize)
        print(f"Image compressed and saved to: {output_path}")

    except FileNotFoundError:
        print(f"Error: Input file not found at {input_path}")
    except Exception as e:
        print(f"An error occurred: {e}")
#%%
sample = "S12"
path = f"../GSE208253/{sample}/raw_data/spatial"
compress_image(f'{path}/tissue_hires_image.png',path)
# %%
