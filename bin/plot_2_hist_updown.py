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
