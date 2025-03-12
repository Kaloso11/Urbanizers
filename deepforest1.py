from deepforest import main
from deepforest import get_data
import matplotlib.pyplot as plt
import rasterio
from deepforest import utilities

# load pretained model from deepforest
model = main.deepforest()
model.use_release()

# define data file
raster_path = "3318DC_16_2016_1143_RGB_RECT.tif"

# open image with raster
# img = rasterio.open(raster_path)

try:
    with rasterio.open(raster_path) as dataset:
        print("Raster opened successfully!")
        print("Raster shape:", dataset.shape)
except Exception as e:
    print("Error opening raster:", str(e))
    

# predicting the image with deepforest model
predicted_raster = model.predict_tile(raster_path, patch_size=400, patch_overlap=0.1)

# annotates the rater of the predictions into a shape file
gdf = utilities.annotations_to_shapefile(predicted_raster, transform=dataset.transform, crs=dataset.crs)

# saves the shapefile of the predictions as a geojson file
gdf.to_file('preds.geojson',driver='GeoJSON')