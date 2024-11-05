# Benchmarking NOOA's Tsunami Runup on a Conical Island Using Thetis

## conical_island_benchmarking
This script attempts to simulate one of NOAA's benchmark cases which can be found here: https://nctr.pmel.noaa.gov/benchmark/Laboratory/Laboratory_ConicalIsland/index.html. This is a simplified application of model verification and it should be noted that greater success can be found by using more advanced techniques like that of the adjoint. This script essentially mimmics the bathymetry of that used in the experiment, and loops through different combinations of viscosity and drag coefficients, outputting velocity and elevation readings at specified locations, set up to be the same locations at which elevation reading are taken in the experiment.

## conical_island_plots
This script is rather clunky, however displays the range of elevation results obtained by the model. Due to the number of output files, it would be recommended to plot the legend separately to each plot.
