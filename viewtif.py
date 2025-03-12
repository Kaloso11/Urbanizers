import rasterio
import numpy as np
import matplotlib.pyplot as plt



# Load the TIFF image
raster_path = "3318DC_16_2016_1143_RGB_RECT.tif"

with rasterio.open(raster_path) as tif:
    img = tif.read([1, 2, 3])  # Reading the first three bands (RGB)
    img = np.moveaxis(img, 0, -1)  # Reorder dimensions from (bands, height, width) to (height, width, bands)

# Convert image to uint8 if necessary
if img.dtype != np.uint8:
    img = (255 * (img - img.min()) / (img.max() - img.min())).astype(np.uint8)

# Show the image
plt.figure(figsize=(10, 10))
plt.imshow(img[:6000, :6000])
plt.xticks([]), plt.yticks([])  # Hide tick values on X and Y axis
plt.show()


