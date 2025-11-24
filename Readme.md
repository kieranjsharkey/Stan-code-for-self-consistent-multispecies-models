### Information on Data ###

Wells are as follows: 
AA01-AA08: 100:0 (tagged:untagged)
AB01-AB08: 90:10
AC01-AC08: 50:50
AD01-AD08: 30:70
AE01-AE08: 70:30
AF01-AF08: 0:100

### Running a model ###
Please ensure that you have all appopriate R libraries installed first. This code depends upon: 
* rlang
* dplyr
* cmdstanr
* posterior
* bayesplot

Please also ensure that you have STAN correctly installed on your system as well (see https://mc-stan.org/install/ for help)

## Single Species Models ### 
The MCMC for single-species models is under Single_species/mc_single_all. Run this as a standard R script and it will automatically cycle through all tagged and untagged data for each model.
This has an approximate runtime of 2h, depending on your system. It will output files entitled 'draws_[strain_name]_[model_name].csv' (e.g. 'draws_untagged_richards.csv') which contain all of the data from the post-burnin chains. 

## Two Species Models ### 
To run a two species model, look in the Two_species. We built four models: 
* Lotka-Volterra = LV
* Gilpin-Ayala = Ayala
* Ram et al. = Ram
* Two-species Richards (hybrid Richards-Gompertz) = GR

To run a MCMC model, run mc_[model_name].R as you would for the single species above. 
Runtime for all two-species models is upwards of 16 hours, often taking more than 24h to complete. It is recommended that you ensure that you are obtaining standard outputs from the single-species models before attempting to run two-species.

## Data ###
All data is returned as .csv files. Headings are common for STAN. Hierachical parameters have [#] next to their names to associate them to a particular well. Hyperparameters are called [parameter]_mean, or [parameter]_sd as in the .STAN code. 