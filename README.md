# Climate_CornYield

To run the project and generate the figures:

1. Start with running read_macametdata.R and read_macametdata_hind.R: 
Read the source MACA-METDATA netcdf projections and hindcasts and extract the subset county level data, save the extracted tmax,tmin,pr,rhmax,rhmin as large matrices. These two scripts are time consuming due to reading large datasets.

*read_macametdata_par.R is a parallel version using foreach as comparison.

2. Run METrawdata_dataframe.R: 
Read the source METDATA (observation data), save the extracted data as a matrix, then calculate the weather variables (GDD,EDD,VPD) and save as a dataframe for further use.

3. Run macametmodel_dataframe.R: 
Process the MACA-METDATA matrices and save as dataframes for further use.

4. Run meanclimateprojection.R: 
Calculate the mean climate projection, save as matrix (together with MACA-METDATA hindcast).

5. Run macamet_linearshifted_dataframe.R:
Process the above matrix: shift the observational METDATA based on the change between hindcast and mean projection (shift temperature based on difference and shift precipitation
/relative humidity based on ratio). We focus on two time windows: 2020-2049 and 2070-2099, so we get two shifted climate projections. Save these projections as dataframes for further use.

6. Run Parametric_sampling.R:
Sample all three uncertainty sources. This is the most time consuming script due to the sampling of parametric uncertainty. Save the calculated best fits and samples of yield projection as matrices.

7. Run Parametric_sampling_linearshiftedclimate.R
Same as 6, but for the two 30-year shifted climate projection. Also save the best fits and samples as matrices.

8. Run 30yravg_data.R:
For each of the 30-year window: calculate the 30-year average data.

9. Run Yieldproj_macamodels.R:
Plot the time series plot (Figure 2) and the marginal distribution of yield projections in each of the two 30-year time window (Figure 3).

10. Run cumulative_sens.R:
Calculate the cumulative sensitiviy at each stage and make the plot (Figure 4).

Some notes below:
Figure 1 is a flow diagram, not plotted by these scripts.
The required source MACA-METDATA data are in "/gpfs/group/kzk10/default/private/data_archive/MACAv2-METDATA/raw" folder, METDATA in "/gpfs/group/kzk10/default/public/METDATA/raw" folder. The required yield data and growing area data are "yielddata.csv" and "harvest_area.csv".

Author: 
Haochen Ye
hxy46@psu.edu
