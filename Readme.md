## Information on Data ##

Biomass and fluorescence growth data is in the files "biomass_ecoli_new.dat" and "Gemini_ecoli_new.dat" respectively. These data can be loaded into R. The wells are as follows: 
AA01-AA08: 100:0 (tagged:untagged)
AB01-AB08: 90:10
AC01-AC08: 50:50
AD01-AD08: 30:70
AE01-AE08: 70:30
AF01-AF08: 0:100

## Running a model ##
This code was run on the R platform (version 4.3.2) and depends on the following R libraries: 
* rlang
* dplyr
* cmdstanr
* posterior
* bayesplot

STAN (version 2.36.0) needs to be installed (see https://mc-stan.org/install/ for help).

Due to a slight but clear anomaly in the early data for AE04, the code below excludes this. The code also excludes all **04 data from the analysis to maintain a balanced design of 7 replicates per initial condition.

All output data is returned as .csv files. Headings are common for STAN. Hierachical parameters have [#] next to their names to associate them to a particular well. Hyperparameters are called [parameter]_mean, or [parameter]_sd as in the .STAN code. 

### Single Species Models ###
The MCMC code for single-species models is under Single_species/mc_single_all. Run this as a standard R script and it will automatically cycle through all tagged and untagged data for each model.
This has an approximate runtime of 2h, depending on your system. It will output files entitled 'draws_[strain_name]_[model_name].csv' (e.g. 'draws_untagged_richards.csv') which contain all of the data from the post-burnin chains. 

### Two Species Models ### 
Two species MCMC code is in the "Two_species" directory. There is one code for each model considered with model names: 
* Lotka-Volterra = LV
* Gilpin-Ayala = Ayala
* Ram et al. = Ram
* Hybrid Richards-Gompertz = GR

The relevant R code is then run mc_[model_name].R. 
Runtime for all two-species models is upwards of 16 hours, often taking more than 24h to complete.





