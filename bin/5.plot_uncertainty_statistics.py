
# from merge_tif import OpenTif
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

if __name__ == "__main__":

#     # define map extent with longitude and latitude in degrees
#     west = 63
#     east = 98
#     south = 36
#     north = 52

#     # fault
#     fault_file = "/exports/geos.ed.ac.uk/comet/qou/datasets/vectors/faults/gem-global-active-faults-master/kml/gem_active_faults_harmonized.kml"
#     faults = gpd.read_file(fault_file, driver='KML')
#     fig, ax = plt.subplots(4, 2, sharex='all', sharey='all', figsize=(8, 12))
#     for x in ax.flatten():
#         faults.plot(ax=x, linewidth=0.5, color='grey')
#         x.set_xlim((west, east))
#         x.set_ylim((south, north))

#     # load and compare InSAR and GPS Ve
#     ve = OpenTif("../los_full/decompose/ve_masked_3.tif")
#     se = OpenTif("../los_full/decompose/ve_sig_3.tif")

#     gps = load_chris()
#     # gps = gps[gps["se"] < 1]   # only keep gps with ve > -8 mm/yr to avoid outliers, which are likely due to bad gps data or bad reference frame transformation
#     print(len(gps))
#     gps.loc[:, 'InSAR_ve'] = [ve.extract_pixel_value(point.x, point.y)[0] for point in gps['geometry']]
#     gps.loc[:, 'InSAR_se'] = [se.extract_pixel_value(point.x, point.y)[0] for point in gps['geometry']]
#     gps["InSAR-GNSS Ve"] = gps['InSAR_ve'] -gps["ve"]
#     gps.dropna(inplace=True)
#     print(len(gps))

#     # plot Ve
#     gps.plot("ve", ax=ax[0,0], vmin=-10, vmax=10, cmap='bwr', markersize=5, legend=True)
#     ax[0,0].set_title("GNSS Ve")
#     gps.plot("InSAR_ve", ax=ax[1,0], vmin=-10, vmax=10, cmap='bwr', markersize=5, legend=True)
#     ax[1,0].set_title("InSAR Ve")
#     gps.plot("InSAR-GNSS Ve", ax=ax[2,0], vmin=-5, vmax=5, cmap='bwr', markersize=5, legend=True)
#     ax[2,0].set_title("InSAR-GNSS Ve")
#     gps.plot("se", ax=ax[3,0], vmin=0, vmax=1, cmap='viridis', markersize=5, legend=True)
#     ax[3,0].set_title("GNSS $\sigma$(Ve)")

#     # load and compare InSAR and GPS Vu
#     vu = OpenTif("../los_full/decompose/vu_3.tif")
#     su = OpenTif("../los_full/decompose/vu_sig_3.tif")

#     gps3d = gps[gps["vu"] != 0]  # only keep gps with non-zero vertical velocity, which are gps3d
#     # gps3d = gps3d[gps3d["su"] < 1.5]  # only keep gps with vu > -8 mm/yr to avoid outliers, which are likely due to bad gps data or bad reference frame transformation
#     print(len(gps3d))
#     gps3d.loc[:, 'InSAR_vu'] = [vu.extract_pixel_value(point.x, point.y)[0] for point in gps3d['geometry']]
#     gps3d.loc[:, 'InSAR_su'] = [su.extract_pixel_value(point.x, point.y)[0] for point in gps3d['geometry']]
#     gps3d["InSAR-GNSS Vu"] = gps3d['InSAR_vu'] -gps3d["vu"]
#     # gps3d = gps3d[gps3d["InSAR_vu"] > -15]  # only keep gps with vu > -8 mm/yr to avoid outliers, which are likely due to bad gps data or bad reference frame transformation
#     gps3d.dropna(inplace=True)
#     print(len(gps3d))

#     # plot Vu
#     gps3d.plot("vu", ax=ax[0,1], vmin=-5, vmax=5, cmap='bwr', markersize=5, legend=True)
#     ax[0,1].set_title("GNSS Vu")
#     gps3d.plot("InSAR_vu", ax=ax[1,1], vmin=-5, vmax=5, cmap='bwr', markersize=5, legend=True)
#     ax[1,1].set_title("InSAR Vu")
#     gps3d.plot("InSAR-GNSS Vu", ax=ax[2,1], vmin=-5, vmax=5, cmap='bwr', markersize=5, legend=True)
#     ax[2,1].set_title("InSAR-GNSS Vu")
#     gps3d.plot("su", ax=ax[3,1], vmin=0, vmax=1, cmap='viridis', markersize=5, legend=True)
#     ax[3,1].set_title("GNSS $\sigma$(Vu)")

#     plt.tight_layout()
#     fig.savefig("../los_full/uncertainty/compare_gps_mask3.png", format='PNG', dpi=300, bbox_inches='tight', transparent=True)


    ve = OpenTif("../los_full/decompose/ve_masked_3.tif")
    vu = OpenTif("../los_full/decompose/vu_3.tif")

    gps = load_chris()
    gps.loc[:, 'InSAR_ve'] = [ve.extract_pixel_value(point.x, point.y)[0] for point in gps['geometry']]
    gps["InSAR-GNSS Ve"] = gps['InSAR_ve'] -gps["ve"]
    gps.dropna(inplace=True)

    gps3d = gps[gps["vu"] != 0]  # only keep gps with non-zero vertical velocity, which are gps3d
    gps3d.loc[:, 'InSAR_vu'] = [vu.extract_pixel_value(point.x, point.y)[0] for point in gps3d['geometry']]
    gps3d["InSAR-GNSS Vu"] = gps3d['InSAR_vu'] -gps3d["vu"]
    gps3d.dropna(inplace=True)

    #####################
    # plot 2 histogram for Ve and Vu differences between InSAR and GPS, and GPS uncertainty
    fig, ax = plt.subplots(2, 1, figsize=(3, 5))

    ax[0].hist(gps["InSAR-GNSS Ve"], bins=np.arange(-10, 10, 0.5), color="C0", )
    ax[0].annotate("RMS: %.1f " % (np.nanstd(gps["InSAR-GNSS Ve"])), xy=(0.03, 0.93), xycoords='axes fraction', ha="left", va="top", fontsize=12)
    ax[0].set_xlabel("mm/yr")
    ax[0].set_title("InSAR-GNSS Ve")
    ax[0].set_xlim(-10, 10)
    ax[0].set_ylabel("Count")
    ax[0].yaxis.tick_right()
    ax[0].yaxis.set_label_position("right")


    ax[1].hist(gps3d["InSAR-GNSS Vu"], bins=np.arange(-10, 10, 0.5), color="C0", )
    ax[1].annotate("RMS: %.1f " % (np.nanstd(gps3d["InSAR-GNSS Vu"])), xy=(0.03, 0.93), xycoords='axes fraction', ha="left", va="top", fontsize=12)
    ax[1].set_xlabel("mm/yr")
    ax[1].set_title("InSAR-GNSS Vu")
    ax[1].set_xlim(-10, 10)
    ax[1].set_ylabel("Count")
    ax[1].yaxis.tick_right()
    ax[1].yaxis.set_label_position("right")

    plt.tight_layout()
    fig.savefig("../los_full/uncertainty/compare_gps_insar_ve_vu_scatter_hist_slides.png", format='PNG', dpi=300, bbox_inches='tight', transparent=True)



#     fig, ax = plt.subplots(3,2, sharex='all')

#     # plot histogram for Ve
#     ax[0,0].hist(gps["InSAR-GNSS Ve"], bins=np.arange(-10, 10, 0.5), color="C0")
#     ax[0,0].set_title("InSAR-GNSS Ve")
#     ax[0,0].annotate("mean: %.1f \n std: %.1f " % (np.nanmean(gps["InSAR-GNSS Ve"]), np.nanstd(gps["InSAR-GNSS Ve"])), xy=(0.03, 0.93), xycoords='axes fraction', ha="left", va="top", fontsize=12)
#     ax[0,0].set_ylabel("Count")


#     # plot histogram for Vu
#     ax[0,1].hist(gps3d["InSAR-GNSS Vu"], bins=np.arange(-10, 10, 0.5), color="C0")
#     ax[0,1].set_title("InSAR-GNSS Vu")
#     ax[0,1].annotate("mean: %.1f \n std: %.1f " % (np.nanmean(gps3d["InSAR-GNSS Vu"]), np.nanstd(gps3d["InSAR-GNSS Vu"])), xy=(0.03, 0.93), xycoords='axes fraction', ha="left", va="top", fontsize=12)


#     # plot histogram for GPS Se
#     ax[1,0].hist(gps["se"], bins=np.arange(-10, 10, 0.5), color="C0")
#     ax[1,0].set_title(r'GNSS $\sigma$(Ve)')
#     ax[1,0].annotate("mean: {:.1f}".format(np.nanmean(gps["se"])), xy=(0.03, 0.93), xycoords='axes fraction', ha="left", va="top", fontsize=12)
#     ax[1,0].set_xlabel("mm/yr")
#     ax[1,0].set_ylabel("Count")


#     # plot histogram for GPS Su
#     ax[1,1].hist(gps3d["su"], bins=np.arange(-10, 10, 0.5), color="C0")
#     ax[1,1].set_title(r'GNSS $\sigma$(Vu)')
#     ax[1,1].annotate("mean: {:.1f}".format(np.nanmean(gps3d["su"])), xy=(0.03, 0.93), xycoords='axes fraction', ha="left", va="top", fontsize=12)
#     ax[1,1].set_xlabel("mm/yr")
#     # ax[1,1].set_ylabel("Count")

#     ax[2,0].scatter(gps["InSAR_ve"], gps["ve"], c=gps["se"], vmin=0, vmax=1, cmap='viridis', s=5)
#     ax[2,0].plot([-10, 10], [-10, 10], c='k', linewidth=1)
#     ax[2,0].set_title("Ve")
#     ax[2,0].set_xlabel("InSAR Ve (mm/yr)")
#     ax[2,0].set_ylabel("GNSS Ve (mm/yr)")

#     ax[2,1].scatter(gps3d["InSAR_vu"], gps3d["vu"], c=gps3d["su"], vmin=0, vmax=1, cmap='viridis', s=5)
#     ax[2,1].plot([-10, 10], [-10, 10], c='k', linewidth=1)
#     ax[2,1].set_title("Vu")
#     ax[2,1].set_xlabel("InSAR Vu (mm/yr)")
#     ax[2,1].set_ylabel("GNSS Vu (mm/yr)")

#     plt.tight_layout()
#     fig.savefig("../los_full/uncertainty/compare_gps_mask3_hist.png", format='PNG', dpi=300, bbox_inches='tight', transparent=True)

#     ######################
#     # plot 2x2 scatter and histogram for Ve and Vu differences between InSAR and GPS, and GPS uncertainty
#     fig, ax = plt.subplots(2, 2, figsize=(8, 7))
#     ax[0,0].errorbar(gps["InSAR_ve"], gps["ve"], xerr=gps["InSAR_se"], yerr=gps["se"], fmt='o', ecolor='lightgray', elinewidth=1, capsize=2, markersize=5, color="C0", alpha=0.7)
#     ax[0,0].plot([-10, 10], [-10, 10], c='k', linewidth=1)
#     ax[0,0].set_xlabel("InSAR Ve (mm/yr)")
#     ax[0,0].set_ylabel("GNSS Ve (mm/yr)")
#     ax[0,0].set_xlim(-10, 10)
#     ax[0,0].set_ylim(-10, 10)
#     ax[0,0].set_aspect('equal')

#     ax[0,1].errorbar(gps3d["InSAR_vu"], gps3d["vu"], xerr=gps3d["InSAR_su"], yerr=gps3d["su"], fmt='o', ecolor='lightgray', elinewidth=1, capsize=2, markersize=5, color="C0", alpha=0.7)
#     ax[0,1].plot([-10, 10], [-10, 10], c='k', linewidth=1)
#     ax[0,1].set_xlabel("InSAR Vu (mm/yr)")
#     ax[0,1].set_ylabel("GNSS Vu (mm/yr)")
#     ax[0,1].set_xlim(-10, 10)
#     ax[0,1].set_ylim(-10, 10)
#     ax[0,1].set_aspect('equal')

#     ax[1,0].hist(gps["InSAR-GNSS Ve"], bins=np.arange(-10, 10, 0.5), color="C0", )
#     ax[1,0].annotate("RMS: %.1f " % (np.nanstd(gps["InSAR-GNSS Ve"])), xy=(0.03, 0.93), xycoords='axes fraction', ha="left", va="top", fontsize=12)
#     ax[1,0].set_xlabel("InSAR-GNSS Ve (mm/yr)")
#     ax[1,0].set_ylabel("Count")
#     ax[1,0].set_xlim(-10, 10)

#     ax[1,1].hist(gps3d["InSAR-GNSS Vu"], bins=np.arange(-10, 10, 0.5), color="C0", )
#     ax[1,1].annotate("RMS: %.1f " % (np.nanstd(gps3d["InSAR-GNSS Vu"])), xy=(0.03, 0.93), xycoords='axes fraction', ha="left", va="top", fontsize=12)
#     ax[1,1].set_xlabel("InSAR-GNSS Vu (mm/yr)")
#     ax[1,1].set_ylabel("Count")
#     ax[1,1].set_xlim(-10, 10)

#     plt.tight_layout()
#     fig.savefig("../los_full/uncertainty/compare_gps_insar_ve_vu_scatter_hist.png", format='PNG', dpi=300, bbox_inches='tight', transparent=True)




# a1_diff=OpenTif("../los_full/uncertainty/a1_diff_TianShan.tif")
# a2_diff=OpenTif("../los_full/uncertainty/a2_diff_TianShan.tif")
# d1_diff=OpenTif("../los_full/uncertainty/d1_diff_TianShan.tif")
# d2_diff=OpenTif("../los_full/uncertainty/d2_diff_TianShan.tif")
# ve_diff_A1D1_A2D2=OpenTif("../los_full/uncertainty/ve_diff_a1d1_a2d2.tif")
# vu_diff_A1D1_A2D2=OpenTif("../los_full/uncertainty/vu_diff_a1d1_a2d2.tif")
# ve_diff_A1D2_A2D1=OpenTif("../los_full/uncertainty/ve_diff_a1d2_a2d1.tif")
# vu_diff_A1D2_A2D1=OpenTif("../los_full/uncertainty/vu_diff_a1d2_a2d1.tif")
# sig_ve=OpenTif("../los_full/uncertainty/ve_sig_TianShan_asig1.2_dsig1.9.tif")
# sig_vu=OpenTif("../los_full/uncertainty/vu_sig_TianShan_asig1.2_dsig1.9.tif")

# import numpy as np
# import matplotlib.pyplot as plt
# from matplotlib.ticker import FormatStrFormatter

# # Helper function to calculate RMS
# def calc_rms(data_array):
#     return np.sqrt(np.nanmean(data_array**2))


# lon_min = sig_ve.left
# lon_max = sig_ve.right
# lat_min = sig_ve.bottom
# lat_max = sig_ve.top
# spatial_extent = (lon_min, lon_max, lat_min, lat_max)

# # ==========================================
# # 2. CALCULATE ZOOM LIMITS (Based on non-NaN data)
# # ==========================================
# # Create a mask of valid data (using a1_diff as the reference for all plots)
# valid_mask = ~np.isnan(sig_ve.data)
# valid_rows = np.where(valid_mask.any(axis=1))[0]
# valid_cols = np.where(valid_mask.any(axis=0))[0]

# # Calculate what percentage of the image the data occupies to map pixels to Lat/Lon
# total_rows, total_cols = sig_ve.data.shape
# lon_range = lon_max - lon_min
# lat_range = lat_max - lat_min

# # Convert pixel bounding box to Lat/Lon bounding box for zooming
# # Note: lat_max is at row 0 (top of image), lat_min is at the last row (bottom)
# zoom_lon_min = lon_min + (valid_cols[0] / total_cols) * lon_range
# zoom_lon_max = lon_min + (valid_cols[-1] / total_cols) * lon_range
# zoom_lat_max = lat_max - (valid_rows[0] / total_rows) * lat_range
# zoom_lat_min = lat_max - (valid_rows[-1] / total_rows) * lat_range


# # ==========================================
# # 3. SETUP FIGURE & PLOT
# # ==========================================
# fig = plt.figure(figsize=(10, 12), layout='constrained')
# gs = fig.add_gridspec(5, 3, width_ratios=[1, 1, 0.02], wspace=0.05)

# ax = np.empty((5, 2), dtype=object)
# for i in range(5):
#     for j in range(2):
#         if i == 0 and j == 0:
#             ax[i, j] = fig.add_subplot(gs[i, j])
#         else:
#             ax[i, j] = fig.add_subplot(gs[i, j], sharex=ax[0,0], sharey=ax[0,0])

# # --- DIFFERENCE PLOTS (Rows 0-3) ---
# # Notice we added `extent=spatial_extent` to every imshow call!
# im_diff = ax[2,0].imshow(a1_diff.data, vmin=-2.5, vmax=2.5, cmap='bwr', extent=spatial_extent)
# ax[2,0].set_title(f"A1 diff (RMS: {calc_rms(a1_diff.data):.1f})")

# ax[2,1].imshow(a2_diff.data, vmin=-2.5, vmax=2.5, cmap='bwr', extent=spatial_extent)
# ax[2,1].set_title(f"A2 diff (RMS: {calc_rms(a2_diff.data):.1f})")

# ax[3,0].imshow(d1_diff.data, vmin=-2.5, vmax=2.5, cmap='bwr', extent=spatial_extent)
# ax[3,0].set_title(f"D1 diff (RMS: {calc_rms(d1_diff.data):.1f})")

# ax[3,1].imshow(d2_diff.data, vmin=-2.5, vmax=2.5, cmap='bwr', extent=spatial_extent)
# ax[3,1].set_title(f"D2 diff (RMS: {calc_rms(d2_diff.data):.1f})")

# ax[0,0].imshow(ve_diff_A1D1_A2D2.data, vmin=-2.5, vmax=2.5, cmap='bwr', extent=spatial_extent)
# ax[0,0].set_title(f"Ve diff A1D1-A2D2 (RMS: {calc_rms(ve_diff_A1D1_A2D2.data):.1f})")

# ax[0,1].imshow(vu_diff_A1D1_A2D2.data, vmin=-2.5, vmax=2.5, cmap='bwr', extent=spatial_extent)
# ax[0,1].set_title(f"Vu diff A1D1-A2D2 (RMS: {calc_rms(vu_diff_A1D1_A2D2.data):.1f})")

# ax[1,0].imshow(ve_diff_A1D2_A2D1.data, vmin=-2.5, vmax=2.5, cmap='bwr', extent=spatial_extent)
# ax[1,0].set_title(f"Ve diff A1D2-A2D1 (RMS: {calc_rms(ve_diff_A1D2_A2D1.data):.1f})")

# ax[1,1].imshow(vu_diff_A1D2_A2D1.data, vmin=-2.5, vmax=2.5, cmap='bwr', extent=spatial_extent)
# ax[1,1].set_title(f"Vu diff A1D2-A2D1 (RMS: {calc_rms(vu_diff_A1D2_A2D1.data):.1f})")


# # --- SIGMA PLOTS (Row 4) ---
# im_sig = ax[4,0].imshow(sig_ve.data, vmin=0, vmax=2, cmap='viridis', extent=spatial_extent)
# ax[4,0].set_title(f"Sigma Ve (RMS: {calc_rms(sig_ve.data):.1f})")

# ax[4,1].imshow(sig_vu.data, vmin=0, vmax=2, cmap='viridis', extent=spatial_extent)
# ax[4,1].set_title(f"Sigma Vu (RMS: {calc_rms(sig_vu.data):.1f})")


# # --- AXES FORMATTING & ZOOMING ---
# # Apply the zoom calculated in step 2 (Only need to do this to one axis since they are shared)
# ax[0,0].set_xlim(zoom_lon_min, zoom_lon_max)
# ax[0,0].set_ylim(zoom_lat_min, zoom_lat_max)

# # Format axes and hide inner labels
# for a in ax.flat:
#     # Set labels for all
#     # a.set_xlabel("Longitude")
#     # a.set_ylabel("Latitude")

#     # # Optional: Format the numbers to look like coordinates with degree symbols
#     # a.xaxis.set_major_formatter(FormatStrFormatter('%.2f°'))
#     # a.yaxis.set_major_formatter(FormatStrFormatter('%.2f°'))

#     # This is the magic command: it hides the labels for inner plots automatically!
#     a.label_outer()


# # --- ADD COLORBARS ---
# cax_diff = fig.add_subplot(gs[0:4, 2])
# cbar_diff = fig.colorbar(im_diff, cax=cax_diff)
# cbar_diff.set_label("Differences, mm/yr")

# cax_sig = fig.add_subplot(gs[4, 2])
# cbar_sig = fig.colorbar(im_sig, cax=cax_sig)
# cbar_sig.set_label("Uncertainty, mm/yr")

# fig.savefig("../los_full/uncertainty/uncertainty_estimates.png", format='PNG', dpi=300, bbox_inches='tight', transparent=True)
