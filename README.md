# README

This MODERA repository contains a list of scripts for **Mosaicking**, **Decomposing** and **Error Analysing** InSAR Line-of-Sight velocities derived from InSAR time series analysis tools such as LiCSBAS (https://github.com/comet-licsar/LiCSBAS).

If you use these scripts, please cite:
- Ou, Q., Daout, S., Weiss, J. R., Shen, L., Lazecký, M., Wright, T. J., & Parsons, B. E. (2022). Large-scale interseismic strain mapping of the NE Tibetan Plateau from Sentinel-1 interferometry. Journal of Geophysical Research: Solid Earth, 127, e2022JB024176. https://doi.org/10.1029/2022JB024176
- Ou, Q., Elliott, J., Maghsoudi, Y., Rollins, C., Lazecky, M., & Wright, T. (2025). Extension of Tian Shan along a nascent shear zone. https://doi.org/10.21203/RS.3.RS-7529996/V1

## merge_tif.py
This script contains 3 classes: **OpenTif**, **Overlap** and **Merge**.
- **OpenTif** opens a GeoTiff file (optionally with uncertainty layer and line-of-sight information) and stores geographical metadata. 
- **Overlap** finds the overlapping pixels between two **OpenTif** classes. 
- **Merge** stitches two or more **OpenTif** classes by doing nothing, or solving for a constant (style=_mode), or a ramp along range (style="_range") or along azimuth (style="_azimuth").

This script is used to stitch line-of-sight geometry files (N/E/U, or incidence/heading angles, and corrected vstd) along track. 

## 0_kriging_gps_vn.py
This script can be used to interpolate point data (e.g., GNSS velocities) through a chosen Kriging algorithm implemented in the pykrige pacakge.

## 1_plate_motion_correction.py
This script applies plate motion correction to line-of-sight velocities decomposed from InSAR times series. The plate motion model rasters are interpolated from an array of latitude and longitude plate motion values calculated using the UNAVCO plate motion calculator (https://www.unavco.org/software/geodetic-utilities/plate-motion-calculator/plate-motion-calculator.html)


## 2_remove_reference_effect_from_vstd.py
This script fits a semivariogram function to the uncertainty profiles away from the reference point, and scales the uncertainty map based on the best-fit semivariogram function, hence remove the reference effect that shows the arteficial near-zero uncertainty around the reference point.

## 3_invert_gps_overlap_with_dummy.py
This script mosaics InSAR line-of-sight velocities along track by solving for the best-fit planar ramp for each frame, exploiting the 2D and 3D GNSS velocities in the track and the frame overlap between neighbouring InSAR frames along track. This script provides the option of including dummy GNSS points sampled from interpolated GNSS velocity fields and forward-calculated point velocities from rigid body rotation to better constrain ramps for frames with inadequate GNSS distribution. The results are tracks of InSAR line-of-sight velocities in a common reference frame provided by the GNSS velocities.

## 4_decompose_into_ve_vu_with_uncertainty_analysis.py
This script decomposes ascending and descending line-of-sight velocity tracks in a common reference frame into east and vertical velocities by first removing the contribution from the interpolated north velocities. The overlapping tracks in the same line-of-sight are arranged into two staggering layers (**A0** with all odd numbered ascending tracks, **A1** with all even numbered ascending tracks, **D0** with all odd numbered descending tracks, and **D1** with all even number descending tracks). 

- With **predict_los_False**, **half_geometry=False** and **constant_los_sig=False**, this script performs decompsition with all tracks and the vstd (with reference effect corrected and then stitched along track). 

- With **predict_los=True**, it removes one of the **A0/A1/D0/D1** geometries at a time, to predict the missing geometry with the remaining three geometries, and calculate the RMS differences between the removed (**observed**) and **predicted** (from east and vertical velocities decomposed from single ascending and single descending observations) line-of-sight velocities. 

- With **half_gemetry=True**, it decomposes **A0/D0** and **A1/D1** into east and vertical velocities, respectively, which overlap in diamonds covered by four tracks. The RMS differences between two east and between two vertical velocities represent variance of east and vertical velocities decomposed from single ascending with single descending observations. The same is done for **A0/D1** and **A1/D0**. The uncertainties of the east and vertical velocities can be propagated through the line-of-sight vector to yield the uncertainties of the **predicted** line-of-sight velocity in the above step. Removing the variance of the predicted layer from the RMS difference gives the variance of the observed line-of-sight velocities.

- With **constant_los_sig=True**, one can enforce a fixed uncertainty for east velocity and for vertical velocity, derived from the assessment above, to forward calculate the pixel-wise uncertainties for the decomposed east and vertical velocities. These uncertainties are generally larger than that derived from time series analysis because of the added uncertainties e.g., due to the mosaiking step. 