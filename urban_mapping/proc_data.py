import os
import json
import geopandas as gpd
import pystac_client
import stackstac
import planetary_computer
from shapely import to_geojson, bounds
import rasterio
import xarray as xr
import rioxarray as rxr
import numpy as np
from sqlalchemy import create_engine
from geocube.api.core import make_geocube
from xrpatcher import XRDAPatcher

import setup_env
setup_env.add_envs()

print('Procs Loaded')

config = dict(os.environ)

def db_scale(x):
    return 10 * np.log10(x)

def min_max_scaler(x):
    return (x-x.min())/(x.max()-x.min())

def get_lc_masks():
    print('Loading data masks')
    connection_string = "postgresql+psycopg2://{user}:{pwd}@{host}:{port}/{dbname}".format(
        user=config['_DB_USERNAME'],
        pwd=config['_DB_PWD'],
        host=config['_DB_ADDRESS'],
        port=config['_DB_PORT'],
        dbname=config['_DATABASE']
    )
    conn = create_engine(connection_string).connect()
    sql = 'SELECT * FROM landcover.lc_masks WHERE lc_class=6'
    masks = gpd.GeoDataFrame.from_postgis(sql=sql, con=conn, geom_col='geometry', crs='EPSG:4326')
    masks['lc_class'] = masks['lc_class'].replace(6,1)
    
    return masks

def get_lc_masks_native():
    print('Loading data masks')
    connection_string = "postgresql+psycopg2://{user}:{pwd}@{host}:{port}/{dbname}".format(
        user=config['_DB_USERNAME'],
        pwd=config['_DB_PWD'],
        host=config['_DB_ADDRESS'],
        port=config['_DB_PORT'],
        dbname=config['_DATABASE']
    )
    conn = create_engine(connection_string).connect()
    sql = 'SELECT * FROM landcover.lc_masks'
    masks = gpd.GeoDataFrame.from_postgis(sql=sql, con=conn, geom_col='geometry', crs='EPSG:4326')
    print(masks[['lc_id','lc_name','lc_class']].head(10))

    return masks

def mask_to_tile(masks, sub_tiles, sub_id):
    print('Converting clipping mask to tile')
    
    sub_tile_info = sub_tiles[sub_tiles.id == sub_id].to_crs('EPSG:32630')

    return masks.clip(sub_tile_info)

def mask_to_tile(masks, sub_tiles, sub_id):
    print('Converting clipping mask to tile')
    
    sub_tile_info = sub_tiles[sub_tiles.id == sub_id].to_crs('EPSG:32630')

    return masks.clip(sub_tile_info)

def mask_to_xr(masks, composites):
    print('Converting mask to raster')
    composites = composites.rio.write_crs('EPSG:32630')
    sub_mask = masks.reset_index().to_crs('EPSG:32630')
    out_grid= make_geocube(vector_data=sub_mask, measurements=["lc_class"], like=composites['DEM'], fill=0.0)
    # sub_tile = out_grid.rio.clip([sub_aoi], all_touched=True)
    
    # class_ts_arr = np.stack([sub_tile["lc_class"] for a in range(len(composites.time))])
    # sub_y, sub_x = sub_tile["lc_class"].shape[0],sub_tile["lc_class"].shape[1]
    
    
    return out_grid

def mask_to_xr_stac_ts(masks, composites):
    print('Converting mask to raster')
    composites = composites.rio.write_crs('EPSG:32630')
    sub_mask = masks.reset_index().to_crs('EPSG:32630')
    out_grid= make_geocube(vector_data=sub_mask, measurements=["lc_class"], like=composites['DEM'], fill=0.0)
    
    
    class_ts_arr = np.stack([out_grid["lc_class"].values for a in range(len(composites.time))])
    
    class_ts =  xr.DataArray(
            class_ts_arr,
            coords={
                "time": composites.time,
                "y": composites.y,
                "x": composites.x,
            },
            dims=["time","y","x"],
            name = 'LC_CLASS')
    
    return class_ts

def mask_to_xr_old(masks, tile_df, tile_id, composites):
    print('Converting mask to raster')
    composites = composites.rio.write_crs('EPSG:32630')
    base_tile = tile_df[tile_df.id == tile_id]['geometry'].to_crs('EPSG:32630')
    print(base_tile)
    base_aoi = json.loads(to_geojson(base_tile.values[0]))
    # sub_aoi = json.loads(to_geojson(sub_tile_info['geometry'].values[0]))
    print(masks[masks.tile_id==tile_id][['lc_id','lc_name','lc_class']])
    out_grid= make_geocube(vector_data=masks[masks.tile_id==tile_id], measurements=["lc_class"], like=composites['DEM'], fill=0.0)
    # sub_tile = out_grid.rio.clip([sub_aoi], all_touched=True)
    
    # class_ts_arr = np.stack([sub_tile["lc_class"] for a in range(len(composites.time))])
    # sub_y, sub_x = sub_tile["lc_class"].shape[0],sub_tile["lc_class"].shape[1]
    
    
    return out_grid
    
    # return xr.DataArray(
    #         class_ts_arr,
    #         coords={
    #             "time": composites.time,
    #             "y": composites.y[:sub_y],
    #             "x": composites.x[:sub_x],
    #         },
    #         dims=["time","y","x"],
    #         name = 'LC_CLASS')

def get_satellite_items(aoi, date_range, bbox, cloud_cover=5):
    print('Searching for satellite images')
    catalog = pystac_client.Client.open(
        "https://planetarycomputer.microsoft.com/api/stac/v1",
        modifier=planetary_computer.sign_inplace
    )
    # 2020-01-01/2020-12-31
    search_s1 = catalog.search(
        collections=["sentinel-1-rtc"],
        intersects=aoi,
        datetime=date_range,
        query={'sat:orbit_state':{"eq":"ascending"}}
    )

    search_s2 = catalog.search(
        collections=["sentinel-2-l2a"],
        intersects=aoi,
        datetime=date_range,
        query={"eo:cloud_cover": {"lt": 5}},
    )

    search_dem = catalog.search(
        collections=["cop-dem-glo-30"],
        intersects=aoi,
    )

    # Check how many items were returned
    s1_items = search_s1.item_collection()
    s2_items = search_s2.item_collection()
    dem_items = search_dem.item_collection()
    
    print(f"Returned S1:{len(s1_items)} S2:{len(s2_items)} DEM:{len(dem_items)} Items")
    
    use_bands = ['B02', 'B03', 'B04', 'B05', 'B06', 'B07', 'B08', 'B09','B11', 'B8A', 'SCL']
    s1_ds = stackstac.stack(s1_items, bounds_latlon=bbox, epsg=32630, resolution=10)
    s2_ds = stackstac.stack(s2_items, assets=use_bands, bounds_latlon=bbox, epsg=32630, resolution=10)
    dem_ds = stackstac.stack(dem_items, bounds_latlon=bbox, epsg=32630, resolution=10)
    
    return s1_ds, s2_ds, dem_ds

    
def get_sar_composites(da):
    ## Process SAR
    print('Preparing SAR')
    sar_xr = xr.Dataset({'VH_DB':(('time','y','x'), min_max_scaler(db_scale(da.sel(band='vh'))).data),
                        'VV_DB':(('time','y','x'), min_max_scaler(db_scale(da.sel(band='vv'))).data)})
    sar_xr.coords['time'] = da.time
    sar_xr.coords['y'] = da.y
    sar_xr.coords['x'] = da.x

    # sar_xr['RVI'] = (4*da.sel(band='vh'))/(da.sel(band='vv')+da.sel(band='vh'))
    return sar_xr.resample(time='ME').median('time')

def get_dem_composites(da, composites):
    ## Process DEM
    print('Preparing DEM')
    dem_dat = da.sel(band='data')
    dem_dat = dem_dat.groupby('time').mean()
    dem_dat_sc = xr.apply_ufunc(min_max_scaler, dem_dat.data.compute(), input_core_dims=[['time', 'lat', 'lon']])

    ## Stack DEM for each Timestamp
    dem_ts_arr = np.stack([dem_dat_sc[0] for a in range(len(composites.time))])
    return xr.DataArray(
            dem_ts_arr,
            coords={
                "time": composites.time,
                "y": composites.y,
                "x": composites.x,
            },
            dims=["time","y","x"],
            name = 'DEM')


def get_opt_composites(da):
    ## Process Optical Sentinel 2
    print('Preparing S2')
    das = {'B02':(('time','y','x'), da.sel(band='B02').data),
        'B03':(('time','y','x'), da.sel(band='B03').data),
        'B04':(('time','y','x'), da.sel(band='B04').data),
        'B08':(('time','y','x'), da.sel(band='B08').data),
        'B11':(('time','y','x'), da.sel(band='B11').data),
        'SCL':(('time','y','x'), da.sel(band='SCL').data)}

    data1 = xr.Dataset(das)
    data1.coords['time'] = da.time
    data1.coords['y'] = da.y
    data1.coords['x'] = da.x

    # Calc Indices
    ndvi = (data1['B08']-data1['B04'])/(data1['B08']+data1['B04'])
    evi2 = 2.4*(data1['B08']-data1['B04'])/(data1['B08']+data1['B04']+1)
    ndmi = (data1['B08']-data1['B11'])/(data1['B08']+data1['B11'])
    ndwi = (data1['B03']-data1['B08'])/(data1['B03']+data1['B08'])
    ndbi = (data1['B11']-data1['B08'])/(data1['B11']+data1['B08'])

    data1['NDVI'] = (('time','y','x'), ndvi.data)
    data1['EVI2'] = (('time','y','x'), evi2.data)
    data1['NDMI'] = (('time','y','x'), ndmi.data)
    data1['NDWI'] = (('time','y','x'), ndwi.data)
    data1['NDBI'] = (('time','y','x'), ndbi.data)

    return data1.resample(time='ME').median('time')

# Process Masks
def make_data_patch(da, patch_size, stride):
    print('Preparing Patches')
    patches = {"x": patch_size, "y": patch_size}
    strides = {"x": stride, "y": stride}
    
    composite_batches = XRDAPatcher(
        da=da,
        patches=patches,
        strides=strides
    )
    
    return composite_batches