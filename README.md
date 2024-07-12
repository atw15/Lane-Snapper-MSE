# Lane Snapper MSE

## Table of Contents
- [Objectives](#objectives)
- [Methodology](#methodology)
- [Operating Models](#operating-models)
- [MSE Simulations](#mse-simulations)
- [Performance Metrics](#performance-metrics)
- [References](#references)
- [Acknowledgments](#acknowledgments)

## Objectives

The objectives of this MSE analysis are:
1. To examine the implications of the updates to the DLM Toolkit model set-up that would explain some of the differences between this analysis and the SEDAR results.
2. To determine whether the revised growth estimates influence the selection of the best Management Procedure (MP) for setting the quota.
3. To evaluate whether size limits should be adjusted considering this new information.

## Methodology

To tackle the objectives, we first compare differences among MSEs conducted using the 2017 packages and dependencies to MSEs conducted using the current versions of the DLM Toolkit when using the same or different operating models (OM). Using the updated versions of the DLM Toolkit, we examine 13 MPs, comparing the OM from SEDAR-49 and the OM with updated growth parameters. Simulations were conducted in the R environment using the Data Limited Methods toolkit (Carruthers et al., 2020b; R Core Team, 2022).

## Operating Models

The OM is the main component of an MSE, describing all the characteristics of the fishery. Within an OM, there are four main parts (Carruthers et al., 2020a):
1. **Stock** - biological parameters describing the stock.
2. **Fleet** - fishing parameters describing the fleet.
3. **Observation (Obs)** - parameters describing how the fishery data is generated from the simulated data, mimicking existing sampling programs.
4. **Implementation (Imp)** - management implementation parameters that describe the effectiveness of management procedures and regulations.

The first OM used parameters and model settings as close as possible to the NOAA SEDAR 49 Update report (Sagarese et al., 2017), while the second OM used updated growth parameters (L∞ = 360.61, K = 0.42, t0 = -0.73) estimated from samples collected from 2015-2022 (n=974).

For detailed parameters of the NOAA OM, refer to the SEDAR report:
NOAA. (2016). Sedar 49 Stock Assessment Report Gulf of Mexico Data-Limited Species: Red Drum, Lane Snapper, Wenchman, Yellowmouth Grouper, Speckled Hind, Snowy Grouper, Almaco Jack, Lesser Amberjack. Retrieved from [SEDAR 49 SAR report](http://sedarweb.org/docs/sar/SEDAR_49_SAR_report.pdf)

## MSE Simulations

1. **MSE-1:** SEDAR 49 replica using their OM and the 2017 versions of R, DLM Toolkit, and dependencies.
2. **MSE-2:** Uses OM with updated life history parameters and the 2017 versions of R, DLM Toolkit, and dependencies.
3. **MSE-3:** Uses OM with updated life history parameters current versions of R, DLM Toolkit, and dependencies.
4. **MSE-4:** SEDAR 49 replica using their OM and current versions of R, DLM Toolkit, and dependencies.

MSE-1 and MSE-2 were conducted using an “renv” project environment with R-version 3.5.1 and DLM Toolkit version 3.2.2. MSE-3 and MSE-4 used R version 4.2.2, the updated DLM Toolkit version 6.0.6, and MSEtool version 3.6.1. 

## Performance Metrics

Performance metrics were designed so that higher values are more desirable. For most performance metrics, a probability of >50% is considered to meet the criteria and ≤50% does not. Each MP was evaluated to determine which, if any, are feasible and meet the performance metrics criterion. Of the feasible methods, the best-performing method was recommended for management. The Long-Term Yield (LTY) was calculated for all NOAA MPs, with criteria as follows: Good (67-100%), Fair (33-67.0%), and Poor (0-33%).

## References

- Carruthers, T. et al. (2020a). [Title of the Reference].
- Carruthers, T. et al. (2020b). [Title of the Reference].
- R Core Team. (2022). R: A language and environment for statistical computing.
- Sagarese, S., Isely, J., & Smith, M. (2016). Sedar 49-Aw-05: Review of Operating Model Parameters for Sedar 49: Lane Snapper.
- Sagarese, S., Isely, J. J., & Smith, M. (2017). Sedar 49: Comprehensive Review of Lane Snapper Results and Data Triage Results.
- NOAA. (2016). Sedar 49 Stock Assessment Report Gulf of Mexico Data-Limited Species: Red Drum, Lane Snapper, Wenchman, Yellowmouth Grouper, Speckled Hind, Snowy Grouper, Almaco Jack, Lesser Amberjack. Retrieved from [SEDAR 49 SAR report](http://sedarweb.org/docs/sar/SEDAR_49_SAR_report.pdf)
- `renv` for package management. See the `renv` documentation for details.

## Acknowledgments

This project was developed under the guidance of Dr. Elizabeth Babcock. Her expertise and support were invaluable in completing this analysis.