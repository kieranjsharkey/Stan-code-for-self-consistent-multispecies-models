####Load Libraries####
library(rlang)
library(dplyr)
library(cmdstanr)
library(posterior)
library(bayesplot)

#### Select Wells ####
# If you would like to fit to different wells, please read readme.md first to 
# see what the initial conditions were in each well. 
# It is important that each list represents a single initial condition. The 
# code may break if this is not true (not tested). 
# Due to abnormal data in column 4, we did not fit to this column and do not recommend doing so, but the data is still provided. 

wells_to_fit <- c(list("AB01","AB02","AB03","AB05","AB06","AB07","AB08"),
                  list("AF01","AF02","AF03","AF05","AF06","AF07","AF08"),
                  list("AE01","AE02","AE03","AE05","AE06","AE07","AE08"))

n_wells <- length(wells_to_fit)

initial_conditions_A <- c(rep(0.9,7),rep(0,7),rep(0.7,7))
initial_conditions_B <- 1 - initial_conditions_A
y0 <- array(c(initial_conditions_A,initial_conditions_B),dim = c(n_wells,2))


#### Pre-set parameters, change at your own risk. ####

t_star = 9 #The time point that rho_1 switches
rho_2 = 3159.1 #Phase 2 scaling
rho_1 = 1034.8 #Phase 1 scaling

#### Data cutoff point ####
# Can vary to consider more/less data. 
length <- 44 #The amount of data we want to take before cutoff.  1 = 0.25h

#### Load Data ####
current_wd <- getwd()

load('../Gemini_ecoli_new.dat')
load('../Biomass_ecoli_new.dat')

Data_Biomass <- data.frame(Data_Biomass)
Data_gemini <- data.frame(Data_gemini)

Time <- data.frame(time = Data_gemini$`Time`[1:length])
model_T <- data.frame(time = c(Data_gemini$`Time`[1:length],Data_Biomass$`Time`[1:length]))
model_Time <- model_T[order(model_T$time),]

gemini_time <- is.element(data.frame(model_Time)$model_Time, Data_gemini$`Time`[1:length])
bio_time <- is.element(data.frame(model_Time)$model_Time, Data_Biomass$`Time`[1:length])

tg <- model_Time[gemini_time]
tb <- model_Time[bio_time]

fitting_biomass <- list()
fitting_gemini <- list()
yA <- data.frame(tg)
yB <- data.frame(tb)
for (well in wells_to_fit) {
  print(well)
  fitting_biomass <- append(fitting_biomass, list(Data_Biomass[1:length,well[[1]]]))
  fitting_gemini  <- append(fitting_gemini, list(Data_gemini[1:length,well[[1]]]))
  yA <- cbind(yA,as.data.frame(Data_gemini[1:length,well[[1]]]))
  yB <- cbind(yB,as.data.frame(Data_Biomass[1:length,well[[1]]]))
}

times <- model_Time
  

names(yA) <- c("time", paste0("well_",seq(1,n_wells,1)))
names(yB) <- c("time", paste0("well_",seq(1,n_wells,1)))
y_combined <- cbind(yA[,-1],yB[,-1])
y_array <- array(as.numeric(unlist(y_combined)),dim = c(length,n_wells,2))

# Need y0 to be a 2D list of vectors

data_list <- list(
  T = length,
  n_wells = n_wells,
  ts = model_Time,
  t_index_g = which(gemini_time),
  t_index_b = which(bio_time),
  y = y_array,
  y0 = y0,
  rho_1 = rho_1,
  rho_2 = rho_2,
  t_star = t_star
)

file <- file.path(paste0(current_wd,"/stan_version_GR.stan"))
mod <- cmdstan_model(file)
mod$print()


fit <- mod$sample(
  data = data_list,
  seed = 1,
  chains = 5,
  parallel_chains = 5,
  iter_warmup = 2500,
  iter_sampling = 10000,            # total number of iterations per chain
  refresh = 50 # print update every 50 iters
)



draws <- fit$draws() %>%
  posterior::as_draws_df()

write.csv(draws, file = "draws_full_model_GR.csv")