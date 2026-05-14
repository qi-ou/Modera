from plot_lib import *
if __name__ == "__main__":

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

    ######################
    # plot 2x2 scatter and histogram for Ve and Vu differences between InSAR and GPS, and GPS uncertainty
    fig, ax = plt.subplots(2, 2, figsize=(8, 7))
    ax[0,0].errorbar(gps["InSAR_ve"], gps["ve"], xerr=gps["InSAR_se"], yerr=gps["se"], fmt='o', ecolor='lightgray', elinewidth=1, capsize=2, markersize=5, color="C0", alpha=0.7)
    ax[0,0].plot([-10, 10], [-10, 10], c='k', linewidth=1)
    ax[0,0].set_xlabel("InSAR Ve (mm/yr)")
    ax[0,0].set_ylabel("GNSS Ve (mm/yr)")
    ax[0,0].set_xlim(-10, 10)
    ax[0,0].set_ylim(-10, 10)
    ax[0,0].set_aspect('equal')

    ax[0,1].errorbar(gps3d["InSAR_vu"], gps3d["vu"], xerr=gps3d["InSAR_su"], yerr=gps3d["su"], fmt='o', ecolor='lightgray', elinewidth=1, capsize=2, markersize=5, color="C0", alpha=0.7)
    ax[0,1].plot([-10, 10], [-10, 10], c='k', linewidth=1)
    ax[0,1].set_xlabel("InSAR Vu (mm/yr)")
    ax[0,1].set_ylabel("GNSS Vu (mm/yr)")
    ax[0,1].set_xlim(-10, 10)
    ax[0,1].set_ylim(-10, 10)
    ax[0,1].set_aspect('equal')

    ax[1,0].hist(gps["InSAR-GNSS Ve"], bins=np.arange(-10, 10, 0.5), color="C0", )
    ax[1,0].annotate("RMS: %.1f " % (np.nanstd(gps["InSAR-GNSS Ve"])), xy=(0.03, 0.93), xycoords='axes fraction', ha="left", va="top", fontsize=12)
    ax[1,0].set_xlabel("InSAR-GNSS Ve (mm/yr)")
    ax[1,0].set_ylabel("Count")
    ax[1,0].set_xlim(-10, 10)

    ax[1,1].hist(gps3d["InSAR-GNSS Vu"], bins=np.arange(-10, 10, 0.5), color="C0", )
    ax[1,1].annotate("RMS: %.1f " % (np.nanstd(gps3d["InSAR-GNSS Vu"])), xy=(0.03, 0.93), xycoords='axes fraction', ha="left", va="top", fontsize=12)
    ax[1,1].set_xlabel("InSAR-GNSS Vu (mm/yr)")
    ax[1,1].set_ylabel("Count")
    ax[1,1].set_xlim(-10, 10)

    plt.tight_layout()
    fig.savefig("../los_full/uncertainty/compare_gps_insar_ve_vu_scatter_hist.png", format='PNG', dpi=300, bbox_inches='tight', transparent=True)
