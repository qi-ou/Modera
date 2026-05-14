import geopandas as gpd
import shapely.speedups
gpd.io.file.fiona.drvsupport.supported_drivers['KML'] = 'rw'
shapely.speedups.enable()
import numpy as np
import matplotlib.pyplot as plt
from shapely.geometry import Point
import os
from osgeo import gdal
import warnings
warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=RuntimeWarning)
warnings.filterwarnings("ignore", category=Warning)

class OpenTif:
    """ a Class that stores the band array and metadata of a Gtiff file."""
    def __init__(self, filename, sigfile=None, incidence=None, heading=None, N=None, E=None, U=None):
        self.ds = gdal.Open(filename)
        self.basename = os.path.splitext(os.path.basename(filename))[0]
        self.band = self.ds.GetRasterBand(1)
        self.data = self.band.ReadAsArray()
        self.xsize = self.ds.RasterXSize
        self.ysize = self.ds.RasterYSize
        self.left = self.ds.GetGeoTransform()[0]
        self.top = self.ds.GetGeoTransform()[3]
        self.xres = self.ds.GetGeoTransform()[1]
        self.yres = self.ds.GetGeoTransform()[5]
        self.right = self.left + self.xsize * self.xres
        self.bottom = self.top + self.ysize * self.yres
        self.projection = self.ds.GetProjection()
        pix_lin, pix_col = np.indices((self.ds.RasterYSize, self.ds.RasterXSize))
        self.lat, self.lon = self.top + self.yres*pix_lin, self.left+self.xres*pix_col
    def extract_pixel_value(self, lon, lat, max_width=200):
        x = int((lon - self.left)/self.xres + 0.5)
        y = int((lat - self.top) / self.yres + 0.5)
        # increase window size in steps of 2 until there are non-nan values in the window
        # starting from 2 with 5x5 window because if 1x1 window, stdev will be zero
        # if we use the std of values instead of the corresponding sigma files as stdev
        for n in np.arange(2, max_width+1, 2):
            pixel_values = self.data[y - n: y + n + 1, x - n: x + n + 1]
            index = np.nonzero(~np.isnan(pixel_values))
            if len(index[0]) > 10:
                # print(n, pixel_values)
                break
        pixel_value = np.nanmean(pixel_values)
        # weighted mean with sigma
        # pixel_sigma = self.sigma[y - n: y + n + 1, x - n: x + n + 1]
        # print(pixel_sigma)
        # pixel_value = np.nansum(pixel_values[index] * pixel_sigma[index]) / np.nansum(pixel_sigma[index])   # takes care of the case where sigma is nan but data is not nan
        stdev = np.nanstd(pixel_values)  # by using nanstd(pixel_values), we are not taking into account the quality of the pixels here.
        return pixel_value, stdev

def load_chris():
    """ Load GPS data into a geopandas dataframe"""
    gps_gdf = gpd.GeoDataFrame(crs="EPSG:4326")
    gps_gdf['geometry'] = None
    index = 0
    gps2d_file = "../gps/chris_tianshan_7_july_2023/tianshan_tol2.1_minocc2.5_dist1_edges_2D_7July2023_weighted.dat"
    gps3d_file = "../gps/chris_tianshan_7_july_2023/tianshan_tol2.1_minocc2.5_dist1_edges_3D_7July2023_weighted.dat"

    fl = open(gps2d_file, "r").readlines()
    for line in fl:
        lon, lat, Ve, Vn, dVe, dVn, Cen, sta = line.split()
        gps_gdf.loc[index, 'geometry'] = Point(float(lon), float(lat))
        gps_gdf.loc[index, 've'] = float(Ve)  # eastern velocity in mm/yr in fixed eurasia reference frame
        gps_gdf.loc[index, 'vn'] = float(Vn)  # northern velocity in mm/yr in fixed eurasia reference frame
        gps_gdf.loc[index, 'vu'] = 0  # 2d gps has no vertical velocity
        gps_gdf.loc[index, 'se'] = float(dVe)  # sigma ve
        gps_gdf.loc[index, 'sn'] = float(dVn)  # sigma vn
        gps_gdf.loc[index, 'su'] = 1  # random value for sigma vu for 2d gps, will be used for data culling
        index += 1

    fl = open(gps3d_file, "r").readlines()
    for line in fl:
        lon, lat, Ve, Vn, Vu, dVe, dVn, dVu, Cen, Ceu, Cnu, sta = line.split()
        gps_gdf.loc[index, 'geometry'] = Point(float(lon), float(lat))
        gps_gdf.loc[index, 've'] = float(Ve)  # eastern velocity in mm/yr in fixed eurasia reference frame
        gps_gdf.loc[index, 'vn'] = float(Vn)  # northern velocity in mm/yr in fixed eurasia reference frame
        gps_gdf.loc[index, 'vu'] = float(Vu)  # northern velocity in mm/yr in fixed eurasia reference frame
        gps_gdf.loc[index, 'se'] = float(dVe)  # sigma ve
        gps_gdf.loc[index, 'sn'] = float(dVn)  # sigma vn
        gps_gdf.loc[index, 'su'] = float(dVu)  # sigma vn
        index += 1

    return gps_gdf
