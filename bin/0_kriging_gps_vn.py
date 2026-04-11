#####################
# Applying the PyKrige package to interpolate point data into surface data.
# Qi Ou, University of Leeds
# 3 Sep 2022
# If you use this script for your study, please cite:
# Ou, Q., Daout, S., Weiss, J. R., Shen, L., Lazecký, M., Wright, T. J., & Parsons, B. E. (2022). Large-scale interseismic strain mapping of the NE Tibetan Plateau from Sentinel-1 interferometry. Journal of Geophysical Research: Solid Earth, 127, e2022JB024176. https://doi.org/10.1029/2022JB024176
##############

import geopandas as gpd
from pykrige.uk import UniversalKriging
from pykrige.kriging_tools import write_asc_grid
import numpy as np
import itertools
import matplotlib.pyplot as plt
from shapely.geometry import Point, Polygon




def polyfit2d(x, y, z, order=3):
    """
    Fit a (default 3rd order) polynomial model predictions to GPS points
    @param x: lon of GPS points
    @param y: lat of GPS points
    @param z: velocity of GPS points
    @return: m: polynomial model
    """
    ncols = (order + 1)**2
    G = np.zeros((x.size, ncols))
    ij = itertools.product(range(order+1), range(order+1))
    for k, (i,j) in enumerate(ij):
        G[:,k] = x**i * y**j
    m, _, _, _ = np.linalg.lstsq(G, z, rcond=None)
    return m


def polyval2d(x, y, m):
    """
    Evaluate polynomial model predictions at given locations based on the polynomial model
    @param x: x array
    @param y: y array
    @param m: polynomial model
    @return:  z array
    """
    order = int(np.sqrt(len(m))) - 1
    ij = itertools.product(range(order+1), range(order+1))
    z = np.zeros_like(x)
    for a, (i,j) in zip(m, ij):
        z += a * x**i * y**j
    return z


def grid_gps_universal_spherical(gps_df, gridx, gridy, comp, zz):
    """ 
    Universal kriging
    @param gps_df: gps dataframe
    @param gridx: x locations for interpolation
    @param gridy: y locations for interpolation
    @param comp: which component (column name) in the GPS dataframe you want to interpolate
    @param zz: best-fit polynomial surface of the same dimension as the interpolated map
    @param m: polynomial model obtained from polyfit2d
    @return: interpolated velocity field and uncertainty field
    """
    # Kriging
    model = "spherical"
    UK = UniversalKriging(
        gps_df.geometry.x,
        gps_df.geometry.y,
        gps_df[comp],
        variogram_model=model,
        weight=True,
        nlags=20,
        enable_plotting=True,
    )

    vel_interpolated, var_interpolated = UK.execute("grid", gridx, gridy)

    return vel_interpolated, np.sqrt(var_interpolated)


def grid_gps(gps_df, gridx, gridy, comp):
    """
    Kriging process
    @param gps_df: GPS dataframe
    @param gridx: x locations where you want to derive an interpolated value
    @param gridy: y locations where you want to derive an interpolated value
    @param comp: which component (column name) in the GPS dataframe you want to interpolate
    @return: interpolated velocity field and uncertainty field
    """
    # Fit a 3rd order, 2d polynomial
    m = polyfit2d(gps_df.geometry.x, gps_df.geometry.y, gps_df[comp])
    xx, yy = np.meshgrid(gridx, gridy)
    zz = polyval2d(xx, yy, m)

    # Plot polynomial fit to data points for universal kriging
    fit, ax = plt.subplots()
    im = ax.imshow(zz, extent=(west, east, south, north), origin='lower', vmin=vel_vmin, vmax=vel_vmax)
    gps_df.plot(comp, ax=ax, vmin=vel_vmin, vmax=vel_vmax, edgecolor="w")
    plt.colorbar(im, ax=ax)
    ax.set_xlim((west, east))
    ax.set_ylim((south, north))
    ax.set_title("Polynomial Surface Fit as Drift")
    plt.show()

    vel_interpolated, sig_interpolated = grid_gps_universal_spherical(gps_df, gridx, gridy, comp, zz)

    return vel_interpolated, sig_interpolated


def plot_interpolation(vel_interpolated, sig_interpolated):
    """ Plot the interpolated velocity and uncertainty maps """
    # plot vel_interpolated
    fig, ax = plt.subplots(2, figsize=(4.4, 6), sharex='all')
    im = ax[0].imshow(vel_interpolated, extent=(west, east, south, north), origin='lower', vmin=vel_vmin, vmax=vel_vmax)  
    plt.colorbar(im, ax=ax[0], label="mm/yr")
    ax[0].set_title("V{}".format(component))
    ax[0].set_xlim((west, east))
    ax[0].set_ylim((south, north))

    # plot sig_interpolated
    sig_vmin = np.nanpercentile(sig_interpolated, 0.5)
    sig_vmax = np.nanpercentile(sig_interpolated, 99.5)
    im = ax[1].imshow(sig_interpolated, extent=(west, east, south, north), origin='lower', vmin=sig_vmin, vmax=sig_vmax)
    # gps.plot("sn", ax=ax[0, 1], vmin=vmin, vmax=vmax, edgecolor="w") # to see where the GNSS are
    plt.colorbar(im, ax=ax[1], label="mm/yr")
    ax[1].set_title("V{}_sig".format(component))
    ax[1].set_xlim((west, east))
    ax[1].set_ylim((south, north))

    plt.show()


def prepare_clean_gps(west, east, south, north, culling_thresh=0.7):
    """ Perform data culling and plot the culling outcome to
    check if there are GPS data with large sigma mixed in with good GPS.
    Remove bad GPS with an adjustable threshold sigma called culling_thresh  """

    # load gps in the box extent for kriging, needs to be larger than InSAR coverage
    all_gps = load_wang_liang('_Eurasia')
    # all_gps = load_chris()
    map_box = Polygon([(west, south), (east, south), (east, north), (west, north)])
    gps_in_box = all_gps[all_gps.within(map_box)]

    # data culling, get rid of contradictory clustering values by removing data with high uncertainties
    gps = gps_in_box[gps_in_box[sig] <= culling_thresh]

    # define vmin vmax
    sig_vmin = 0
    sig_vmax = 2
    vel_vmin = np.nanpercentile(gps[vel], 1)
    vel_vmax = np.nanpercentile(gps[vel], 99)

    # plotting data culling
    fig, ax = plt.subplots()
    gps.sort_values(sig).plot(sig, ax=ax, legend=True, vmin=sig_vmin, vmax=sig_vmax)
    gps[gps[sig] > culling_thresh].plot(sig, edgecolor='r', ax=ax, vmin=sig_vmin, vmax=sig_vmax)
    ax.set_title('GNSS Sigma V{} (mm/yr) , red=({}>{})'.format(component, sig, culling_thresh))
    plt.show()

    return gps, sig_vmin, sig_vmax, vel_vmin, vel_vmax


def load_wang_liang(reference):
    """ Load GPS data into a geopandas dataframe
    - I needed to combine GPS velocities from two publications
    - this section could be much simpler if using GPS from just one publication
    - only need lon, lat, vel and sig for the purpose of interpolation"""
    gps_gdf = gpd.GeoDataFrame(crs="EPSG:4326")
    gps_gdf['geometry'] = None
    index = 0
    wang = '/Users/qi/OneDrive/WFH/insar/GPS/Wang_Shen_2020_GPS_JGR/Table.S4_results'
    fl = open(wang, "r").readlines()
    for line in fl:
        if not line.startswith(('*', "Table")):
            sta, lon, lat, Ve_I, Ve_E, dVe, Vn_I, Vn_E, dVn, Cen = line.split()
            gps_gdf.loc[index, 'geometry'] = Point(float(lon), float(lat))
            gps_gdf.loc[index, 'station'] = sta[:4]
            if reference == '_Eurasia':
                gps_gdf.loc[index, 've'] = float(Ve_E)  # eastern velocity in mm/yr in fixed eurasia reference frame
                gps_gdf.loc[index, 'vn'] = float(Vn_E)  # northern velocity in mm/yr in fixed eurasia reference frame
            elif reference == '_ITRF2008':
                gps_gdf.loc[index, 've'] = float(Ve_I)  # eastern velocity in mm/yr in ITRF20082008 reference frame
                gps_gdf.loc[index, 'vn'] = float(Vn_I)  # northern velocity in mm/yr in ITRF20082008 reference frame
            else:
                raise ValueError('reference has to be _Eurasia or _ITRF2008')
            gps_gdf.loc[index, 'se'] = float(dVe)  # sigma ve
            gps_gdf.loc[index, 'sn'] = float(dVn)  # sigma vn
            gps_gdf.loc[index, 'cen'] = float(Cen)  # correlation between east and northern velocities
            gps_gdf.loc[index, 'vu'] = float(0)  # up velocity in mm/yr in ITRF20082008 reference frame
            gps_gdf.loc[index, 'su'] = float(0)  # sigma vu in ITRF20082008 reference frame
            index += 1

    liang = '/Users/qi/OneDrive/WFH/insar/GPS/Liang_etal_2013_supporting/GPS_Liang_etal_2013.csv'
    fl = open(liang, "r").readlines()
    for line in fl:
        if not line.startswith(('*', "Sta.", ",", " Eurasia", " ITRF2008")):
            sta, lon, lat, VN, σvN, VE, σvE, CNE, VU, σvU, Vu, σvu = line.split(',')
            if VU != '-':  # only input vertical if Liang as a vertical value
                if sta in gps_gdf['station'].values:  # only input vertical if Wang has the same station
                    index = gps_gdf.index[gps_gdf['station'] == sta].tolist()[0]  # get index in dataframe based on station
                    if reference == '_Eurasia':
                        gps_gdf.loc[index, 'vu'] = float(Vu)  # up velocity in mm/yr relative to stable northern neighbour
                        gps_gdf.loc[index, 'su'] = float(σvu)  # sigma vu relative to stable northern neighbour
                    elif reference == '_ITRF2008':
                        gps_gdf.loc[index, 'vu'] = float(VU)  # up velocity in mm/yr in ITRF2008 reference frame
                        gps_gdf.loc[index, 'su'] = float(σvU)  # sigma vu in ITRF20082008 reference frame
    return gps_gdf


if __name__ == "__main__":
    ##########################
    #      Define inputs     #
    ##########################

    # define sampling resolution in degree
    res = 0.2

    # define map extent with longitude and latitude in degrees
    west = 60
    east = 96
    south = 34
    north = 53

    # define component name
    component = 'n'

    # define data culling threshold
    cleaning_threshold = 1  # gps sigma Vn in mm/yr => only Vn with lower uncertainties are used for kriging.

    # define output directory
    out_dir = "../gps/"

    ##########################
    #      Loading data      #
    ##########################

    # load gps and perform data culling
    print("load gps and perform data culling")
    sig = 's{}'.format(component)  # only used for data culling, not for weighting interpolation
    vel = 'v{}'.format(component)  # used for kriging interpolation
    # the function prepare_clean_gps calls the function load_wang_liang for loading gps from two sources,
    # feel free to modify the loading commands to suit your use case
    gps, sig_vmin, sig_vmax, vel_vmin, vel_vmax = prepare_clean_gps(west, east, south, north, cleaning_threshold)

    ##########################
    #         Kriging        #
    ##########################

    # Kriging
    print("Kriging starts ... ")
    gridx = np.arange(west, east+res, res)
    gridy = np.arange(south, north+res, res)
    vel_interpolated, sig_interpolated = grid_gps(gps, gridx, gridy, vel)

    # plot kriging results
    print("Plot Kriging Results ... ")
    plot_interpolation(vel_interpolated, sig_interpolated)

    # write data to file
    print("Saving Kriging Results ... ")
    write_asc_grid(gridx, gridy, vel_interpolated, filename=out_dir+"{}_interpolated.grd".format(vel))
    write_asc_grid(gridx, gridy, sig_interpolated, filename=out_dir+"{}_interpolated.grd".format(sig))

