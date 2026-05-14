from plot_lib import *
if __name__ == "__main__":

    a1_diff=OpenTif("../los_full/uncertainty/a1_diff_TianShan.tif")
    a2_diff=OpenTif("../los_full/uncertainty/a2_diff_TianShan.tif")
    d1_diff=OpenTif("../los_full/uncertainty/d1_diff_TianShan.tif")
    d2_diff=OpenTif("../los_full/uncertainty/d2_diff_TianShan.tif")
    ve_diff_A1D1_A2D2=OpenTif("../los_full/uncertainty/ve_diff_a1d1_a2d2.tif")
    vu_diff_A1D1_A2D2=OpenTif("../los_full/uncertainty/vu_diff_a1d1_a2d2.tif")
    ve_diff_A1D2_A2D1=OpenTif("../los_full/uncertainty/ve_diff_a1d2_a2d1.tif")
    vu_diff_A1D2_A2D1=OpenTif("../los_full/uncertainty/vu_diff_a1d2_a2d1.tif")
    sig_ve=OpenTif("../los_full/uncertainty/ve_sig_TianShan_asig1.2_dsig1.9.tif")
    sig_vu=OpenTif("../los_full/uncertainty/vu_sig_TianShan_asig1.2_dsig1.9.tif")

    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib.ticker import FormatStrFormatter

    # Helper function to calculate RMS
    def calc_rms(data_array):
        return np.sqrt(np.nanmean(data_array**2))


    lon_min = sig_ve.left
    lon_max = sig_ve.right
    lat_min = sig_ve.bottom
    lat_max = sig_ve.top
    spatial_extent = (lon_min, lon_max, lat_min, lat_max)

    # ==========================================
    # 2. CALCULATE ZOOM LIMITS (Based on non-NaN data)
    # ==========================================
    # Create a mask of valid data (using a1_diff as the reference for all plots)
    valid_mask = ~np.isnan(sig_ve.data)
    valid_rows = np.where(valid_mask.any(axis=1))[0]
    valid_cols = np.where(valid_mask.any(axis=0))[0]

    # Calculate what percentage of the image the data occupies to map pixels to Lat/Lon
    total_rows, total_cols = sig_ve.data.shape
    lon_range = lon_max - lon_min
    lat_range = lat_max - lat_min

    # Convert pixel bounding box to Lat/Lon bounding box for zooming
    # Note: lat_max is at row 0 (top of image), lat_min is at the last row (bottom)
    zoom_lon_min = lon_min + (valid_cols[0] / total_cols) * lon_range
    zoom_lon_max = lon_min + (valid_cols[-1] / total_cols) * lon_range
    zoom_lat_max = lat_max - (valid_rows[0] / total_rows) * lat_range
    zoom_lat_min = lat_max - (valid_rows[-1] / total_rows) * lat_range


    # ==========================================
    # 3. SETUP FIGURE & PLOT
    # ==========================================
    fig = plt.figure(figsize=(10, 12), layout='constrained')
    gs = fig.add_gridspec(5, 3, width_ratios=[1, 1, 0.02], wspace=0.05)

    ax = np.empty((5, 2), dtype=object)
    for i in range(5):
        for j in range(2):
            if i == 0 and j == 0:
                ax[i, j] = fig.add_subplot(gs[i, j])
            else:
                ax[i, j] = fig.add_subplot(gs[i, j], sharex=ax[0,0], sharey=ax[0,0])

    # --- DIFFERENCE PLOTS (Rows 0-3) ---
    # Notice we added `extent=spatial_extent` to every imshow call!
    im_diff = ax[2,0].imshow(a1_diff.data, vmin=-2.5, vmax=2.5, cmap='bwr', extent=spatial_extent)
    ax[2,0].set_title(f"A1 diff (RMS: {calc_rms(a1_diff.data):.1f})")

    ax[2,1].imshow(a2_diff.data, vmin=-2.5, vmax=2.5, cmap='bwr', extent=spatial_extent)
    ax[2,1].set_title(f"A2 diff (RMS: {calc_rms(a2_diff.data):.1f})")

    ax[3,0].imshow(d1_diff.data, vmin=-2.5, vmax=2.5, cmap='bwr', extent=spatial_extent)
    ax[3,0].set_title(f"D1 diff (RMS: {calc_rms(d1_diff.data):.1f})")

    ax[3,1].imshow(d2_diff.data, vmin=-2.5, vmax=2.5, cmap='bwr', extent=spatial_extent)
    ax[3,1].set_title(f"D2 diff (RMS: {calc_rms(d2_diff.data):.1f})")

    ax[0,0].imshow(ve_diff_A1D1_A2D2.data, vmin=-2.5, vmax=2.5, cmap='bwr', extent=spatial_extent)
    ax[0,0].set_title(f"Ve diff A1D1-A2D2 (RMS: {calc_rms(ve_diff_A1D1_A2D2.data):.1f})")

    ax[0,1].imshow(vu_diff_A1D1_A2D2.data, vmin=-2.5, vmax=2.5, cmap='bwr', extent=spatial_extent)
    ax[0,1].set_title(f"Vu diff A1D1-A2D2 (RMS: {calc_rms(vu_diff_A1D1_A2D2.data):.1f})")

    ax[1,0].imshow(ve_diff_A1D2_A2D1.data, vmin=-2.5, vmax=2.5, cmap='bwr', extent=spatial_extent)
    ax[1,0].set_title(f"Ve diff A1D2-A2D1 (RMS: {calc_rms(ve_diff_A1D2_A2D1.data):.1f})")

    ax[1,1].imshow(vu_diff_A1D2_A2D1.data, vmin=-2.5, vmax=2.5, cmap='bwr', extent=spatial_extent)
    ax[1,1].set_title(f"Vu diff A1D2-A2D1 (RMS: {calc_rms(vu_diff_A1D2_A2D1.data):.1f})")


    # --- SIGMA PLOTS (Row 4) ---
    im_sig = ax[4,0].imshow(sig_ve.data, vmin=0, vmax=2, cmap='viridis', extent=spatial_extent)
    ax[4,0].set_title(f"Sigma Ve (RMS: {calc_rms(sig_ve.data):.1f})")

    ax[4,1].imshow(sig_vu.data, vmin=0, vmax=2, cmap='viridis', extent=spatial_extent)
    ax[4,1].set_title(f"Sigma Vu (RMS: {calc_rms(sig_vu.data):.1f})")


    # --- AXES FORMATTING & ZOOMING ---
    # Apply the zoom calculated in step 2 (Only need to do this to one axis since they are shared)
    ax[0,0].set_xlim(zoom_lon_min, zoom_lon_max)
    ax[0,0].set_ylim(zoom_lat_min, zoom_lat_max)

    # Format axes and hide inner labels
    for a in ax.flat:
        # Set labels for all
        # a.set_xlabel("Longitude")
        # a.set_ylabel("Latitude")

        # # Optional: Format the numbers to look like coordinates with degree symbols
        # a.xaxis.set_major_formatter(FormatStrFormatter('%.2f°'))
        # a.yaxis.set_major_formatter(FormatStrFormatter('%.2f°'))

        # This is the magic command: it hides the labels for inner plots automatically!
        a.label_outer()


    # --- ADD COLORBARS ---
    cax_diff = fig.add_subplot(gs[0:4, 2])
    cbar_diff = fig.colorbar(im_diff, cax=cax_diff)
    cbar_diff.set_label("Differences, mm/yr")

    cax_sig = fig.add_subplot(gs[4, 2])
    cbar_sig = fig.colorbar(im_sig, cax=cax_sig)
    cbar_sig.set_label("Uncertainty, mm/yr")

    fig.savefig("../los_full/uncertainty/uncertainty_estimates.png", format='PNG', dpi=300, bbox_inches='tight', transparent=True)
