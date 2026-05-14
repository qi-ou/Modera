from plot_lib import *
if __name__ == "__main__":

    # define map extent with longitude and latitude in degrees
    west = 63
    east = 98
    south = 36
    north = 52

    # fault
    fault_file = "/exports/geos.ed.ac.uk/comet/qou/datasets/vectors/faults/gem-global-active-faults-master/kml/gem_active_faults_harmonized.kml"
    faults = gpd.read_file(fault_file, driver='KML')
    fig, ax = plt.subplots(4, 2, sharex='all', sharey='all', figsize=(8, 12))
    for x in ax.flatten():
        faults.plot(ax=x, linewidth=0.5, color='grey')
        x.set_xlim((west, east))
        x.set_ylim((south, north))

    # load and compare InSAR and GPS Ve
    ve = OpenTif("../los_full/decompose/ve_masked_3.tif")
    se = OpenTif("../los_full/decompose/ve_sig_3.tif")

    gps = load_chris()
    # gps = gps[gps["se"] < 1]   # only keep gps with ve > -8 mm/yr to avoid outliers, which are likely due to bad gps data or bad reference frame transformation
    print(len(gps))
    gps.loc[:, 'InSAR_ve'] = [ve.extract_pixel_value(point.x, point.y)[0] for point in gps['geometry']]
    gps.loc[:, 'InSAR_se'] = [se.extract_pixel_value(point.x, point.y)[0] for point in gps['geometry']]
    gps["InSAR-GNSS Ve"] = gps['InSAR_ve'] -gps["ve"]
    gps.dropna(inplace=True)
    print(len(gps))

    # plot Ve
    gps.plot("ve", ax=ax[0,0], vmin=-10, vmax=10, cmap='bwr', markersize=5, legend=True)
    ax[0,0].set_title("GNSS Ve")
    gps.plot("InSAR_ve", ax=ax[1,0], vmin=-10, vmax=10, cmap='bwr', markersize=5, legend=True)
    ax[1,0].set_title("InSAR Ve")
    gps.plot("InSAR-GNSS Ve", ax=ax[2,0], vmin=-5, vmax=5, cmap='bwr', markersize=5, legend=True)
    ax[2,0].set_title("InSAR-GNSS Ve")
    gps.plot("se", ax=ax[3,0], vmin=0, vmax=1, cmap='viridis', markersize=5, legend=True)
    ax[3,0].set_title("GNSS $\sigma$(Ve)")

    # load and compare InSAR and GPS Vu
    vu = OpenTif("../los_full/decompose/vu_3.tif")
    su = OpenTif("../los_full/decompose/vu_sig_3.tif")

    gps3d = gps[gps["vu"] != 0]  # only keep gps with non-zero vertical velocity, which are gps3d
    # gps3d = gps3d[gps3d["su"] < 1.5]  # only keep gps with vu > -8 mm/yr to avoid outliers, which are likely due to bad gps data or bad reference frame transformation
    print(len(gps3d))
    gps3d.loc[:, 'InSAR_vu'] = [vu.extract_pixel_value(point.x, point.y)[0] for point in gps3d['geometry']]
    gps3d.loc[:, 'InSAR_su'] = [su.extract_pixel_value(point.x, point.y)[0] for point in gps3d['geometry']]
    gps3d["InSAR-GNSS Vu"] = gps3d['InSAR_vu'] -gps3d["vu"]
    # gps3d = gps3d[gps3d["InSAR_vu"] > -15]  # only keep gps with vu > -8 mm/yr to avoid outliers, which are likely due to bad gps data or bad reference frame transformation
    gps3d.dropna(inplace=True)
    print(len(gps3d))

    # plot Vu
    gps3d.plot("vu", ax=ax[0,1], vmin=-5, vmax=5, cmap='bwr', markersize=5, legend=True)
    ax[0,1].set_title("GNSS Vu")
    gps3d.plot("InSAR_vu", ax=ax[1,1], vmin=-5, vmax=5, cmap='bwr', markersize=5, legend=True)
    ax[1,1].set_title("InSAR Vu")
    gps3d.plot("InSAR-GNSS Vu", ax=ax[2,1], vmin=-5, vmax=5, cmap='bwr', markersize=5, legend=True)
    ax[2,1].set_title("InSAR-GNSS Vu")
    gps3d.plot("su", ax=ax[3,1], vmin=0, vmax=1, cmap='viridis', markersize=5, legend=True)
    ax[3,1].set_title("GNSS $\sigma$(Vu)")

    plt.tight_layout()
    fig.savefig("../los_full/uncertainty/compare_gps_mask3.png", format='PNG', dpi=300, bbox_inches='tight', transparent=True)
