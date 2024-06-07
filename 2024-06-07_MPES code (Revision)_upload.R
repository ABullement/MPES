# Start run time calcs # ----
start_time <- Sys.time()
# Load libraries (install them if missing) ----
# install.packages("survextrap", repos=c('https://chjackson.r-universe.dev',
#                                        'https://cloud.r-project.org'))
pacman::p_load(flexsurv, tidyr, rstan, dplyr, survminer, openxlsx,
               expertsurv, rjags, R2jags, survtools, R2WinBUGS, survextrap,
               survHE, INLA, pracma)

# Set working directory and seed ----
setwd("XXXX")
set.seed(99999)

# General settings = ----
plot_x_limit <- 50 # true upper limit for the model fitting
output_plot_x_limit <- 30 # upper limit for the outputted graph

# Reading in data # ----
gp.data.raw <-read.xlsx(xlsxFile = "data/nationallifetables3yruk.xlsx", sheet = "2018-2020", startRow = 6, colNames = TRUE, skipEmptyCols = TRUE)

# Select case study # ----
# Choose case study by commenting out all options except one:
# case_study <- "SCCHN"
case_study <- "Melanoma"
# case_study <- "NSCLC"

# For NSCLC scenario - reduce survival estimates (Figure 3, Panel C)?
reduce_ext <- FALSE
reduce_ext_perc <- 0.9

# For Jackson's MPES, enable or disable relative survival?
# (switch located here for ease of running)
jackson_relsurv <- TRUE

# Data tidying # ----

if(case_study == "SCCHN"){
  # SCCHN (head and neck)
  trial.data.raw <- read.csv("data/bonner_data.csv")
  trial.data.raw.long <- read.csv("data/bonner_data_2010.csv")
  ext.data.discrete.raw <- read.csv("data/seer_data.csv")
  ta.comp <- read.csv("data/NICE TA_SCCHN.csv")
} else if(case_study == "Melanoma"){
  # Melanoma
  trial.data.raw <- read.csv("data/ca184024_data_orig.csv")
  trial.data.raw.long <- read.csv("data/ca184024_data.csv")
  ext.data.discrete.raw <- read.csv("data/ajcc_data.csv")
  ta.comp <- read.csv("data/NICE TA_Melanoma.csv")
} else {
  # NSCLC (lung)
  trial.data.raw <- read.csv("data/reck_data_2016.csv")
  trial.data.raw.long <- read.csv("data/reck_data_2021.csv")
  ext.data.discrete.raw <- read.csv("data/herbst_data_v2.csv")
  ta.comp <- read.csv("data/NICE TA_NSCLC.csv")
}

# Manual override of case study - use to include different data sets
# trial.data.raw <- read.csv("data/XXXX.csv")
# trial.data.raw.long <- read.csv("data/XXXX.csv")
# ext.data.discrete.raw <- read.csv("data/XXXX.csv")

if (case_study == "NSCLC"){
  if (reduce_ext == TRUE){
    ext.data.discrete.raw$r <- floor(ext.data.discrete.raw$r * reduce_ext_perc)
  }
}

# Tidy up the data
trial.data <- subset(trial.data.raw, treat == 1)
trial.data$treat <- as.factor(trial.data$treat)
rownames(trial.data) <- NULL #need this to avoid gen gamma code for informative prior breaking, as relies on ID column
km_function <- survfit(Surv(t,dth)~1,data=trial.data,type="kaplan-meier")
km_function.long <- survfit(Surv(t,dth)~1,data=subset(trial.data.raw.long, treat == 1),type="kaplan-meier")

trial.data_c <- subset(trial.data.raw, treat == 0)
trial.data_c$treat <- as.factor(trial.data_c$treat)
rownames(trial.data_c) <- NULL #need this to avoid gen gamma code for informative prior breaking, as relies on ID column
km_function_c <- survfit(Surv(t,dth)~1,data=trial.data_c,type="kaplan-meier")
km_function.long_c <- survfit(Surv(t,dth)~1,data=subset(trial.data.raw.long, treat == 0),type="kaplan-meier")

# Guyot MPES function # ----
run_mpes <- function(single_arm_mode = FALSE, #use to denote running for one arm only, or two arms
                     gen_pop_adjust = TRUE, #use to include or exclude general population mortality adjustment
                     hr_adjust = TRUE, #use to include or exclude HR prior for treatment effect
                     max_surv = plot_x_limit, #upper limit for extracting survival probabilities
                     iter_spec = 1000, #number of iterations for Bayesian analysis
                     gp_specs = c(34, 30196, 23080), #general population mortality inputs - year, starting survivors, ending survivors
                     hr_specs = c(1, 0.1, 6, 28), #HR inputs - mean, SD, start time, end time
                     knot_loc = FALSE, #default FALSE means function will pick 2 internal knots
                     trial_source = trial.data.raw, #source data for trial
                     ext_source = ext.data.discrete.raw, #source data for external
                     gp_choice = TRUE){
  ## Read in data
  trial.dat <- trial_source # Patient-level data from the trial
  ext.dat <- ext_source # Conditional survival data 
  
  n_pred_feed <- max_surv*12
  warmup_spec <- 0.1*iter_spec #note: this can be changed as needed!
  
  if (single_arm_mode == TRUE){
    arm_select <- 0 #choose to subset the trial.dat by the flag used by treatment arm
    trial.dat <- subset(trial.dat, treat == arm_select)
    trial.dat$treat <- 1 #Need to set this back to 1 for the stan code to function as intended!
  }
  
  if (gen_pop_adjust == TRUE){
    genpop.data <- data.frame(t_n = gp_specs[1]*12, t_r = (gp_specs[1]+1)*12, n = gp_specs[2], r = gp_specs[3]) # general population data
    if (gp_choice == TRUE){
      gp_choice_int <- 1
    } else {
      gp_choice_int <- 2
    }
  }
  
  if (hr_adjust == TRUE && single_arm_mode == FALSE){
    N_hr_feed <- hr_specs[4] - hr_specs[3] + 1
  }
  
  logt <- log(trial.dat$t) # Create a column for log times - this is used later in the code
  
  ## Bootstrap data for range of plausible starting values
  b_dat <- list(b1 = sample(1:nrow(trial.dat), 1000, replace = TRUE),
                b2 = sample(1:nrow(trial.dat), 1000, replace = TRUE),
                b3 = sample(1:nrow(trial.dat), 1000, replace = TRUE)) # 3 chains used
  
  ## Set knot positions (and calculate lambdas) ----
  k <- knot_loc
  if (length(k) == 1){
    k <- c(0, log(median(trial.dat$t)), log(median(ext.dat$t_n)), log(max_surv*12))
  }
  lambda <- (k[length(k)] - k[2:(length(k)-1)]) / (k[length(k)] - k[1])
  
  ## Fit spline models to intercept only model to get starting parameters that the model can fit
  if(single_arm_mode==TRUE){
    spl_dat <- list(
      r1 = flexsurvspline(Surv(t, dth) ~ 1,
                          data = trial.dat[b_dat[[1]], ],
                          bknots = c(min(ifelse(trial.dat$dth[b_dat[[1]]], log(trial.dat$t[b_dat[[1]]]), NA), na.rm = TRUE),
                                     max(ifelse(trial.dat$dth[b_dat[[1]]], log(trial.dat$t[b_dat[[1]]]), NA), na.rm = TRUE))),
      r2 = flexsurvspline(Surv(t, dth) ~ 1,
                          data = trial.dat[b_dat[[2]], ],
                          bknots = c(min(ifelse(trial.dat$dth[b_dat[[2]]], log(trial.dat$t[b_dat[[2]]]), NA), na.rm = TRUE),
                                     max(ifelse(trial.dat$dth[b_dat[[2]]], log(trial.dat$t[b_dat[[2]]]), NA), na.rm = TRUE))),
      r3 = flexsurvspline(Surv(t, dth) ~ 1,
                          data = trial.dat[b_dat[[3]], ],
                          bknots = c(min(ifelse(trial.dat$dth[b_dat[[3]]], log(trial.dat$t[b_dat[[3]]]), NA), na.rm = TRUE),
                                     max(ifelse(trial.dat$dth[b_dat[[3]]], log(trial.dat$t[b_dat[[3]]]), NA), na.rm = TRUE)))
    )
  }
  else {
    spl_dat <- list(
      r1 = flexsurvspline(Surv(t, dth) ~ gamma1(treat),
                          data = trial.dat[b_dat[[1]], ],
                          bknots = c(min(ifelse(trial.dat$dth[b_dat[[1]]], log(trial.dat$t[b_dat[[1]]]), NA), na.rm = TRUE),
                                     max(ifelse(trial.dat$dth[b_dat[[1]]], log(trial.dat$t[b_dat[[1]]]), NA), na.rm = TRUE))),
      r2 = flexsurvspline(Surv(t, dth) ~ gamma1(treat),
                          data = trial.dat[b_dat[[2]], ],
                          bknots = c(min(ifelse(trial.dat$dth[b_dat[[2]]], log(trial.dat$t[b_dat[[2]]]), NA), na.rm = TRUE),
                                     max(ifelse(trial.dat$dth[b_dat[[2]]], log(trial.dat$t[b_dat[[2]]]), NA), na.rm = TRUE))),
      r3 = flexsurvspline(Surv(t, dth) ~ gamma1(treat),
                          data = trial.dat[b_dat[[3]], ],
                          bknots = c(min(ifelse(trial.dat$dth[b_dat[[3]]], log(trial.dat$t[b_dat[[3]]]), NA), na.rm = TRUE),
                                     max(ifelse(trial.dat$dth[b_dat[[3]]], log(trial.dat$t[b_dat[[3]]]), NA), na.rm = TRUE)))
    )
  }
  
  ## Set up data list for Stan program (short version)
  data_feed_base <-
    list(N = nrow(trial.dat),
         m = length(k)-1,
         is_censored = 1 - trial.dat$dth,
         log_times = logt,
         N_pred = n_pred_feed,
         log_pred_t = log(1:n_pred_feed),
         N_ext = nrow(ext.dat),
         n_ext = ext.dat$n,
         r_ext = ext.dat$r,
         idx1_ext = ext.dat$t_n,
         idx2_ext = ext.dat$t_r
    )
  
    #next section depends on choice of number of knots
  #selects relevant part based on length of knots vector k - order is 1, 2, and 3 internal knots
  
  if((length(k)-1)==2){
    ## Data needed for two-arm version
    data_feed_twoarm <-
      list(basis_evals = rbind(
        (1 - trial.dat$treat) * logt,
        (1 - trial.dat$treat) * (pmax(0.0, logt - k[2])^3 - lambda * pmax(0.0, logt - k[1])^3 - (1 - lambda) * pmax(0.0, logt - k[3])^3),
        trial.dat$treat * logt,
        trial.dat$treat * (pmax(0.0, logt - k[2])^3 - lambda * pmax(0.0, logt - k[1])^3 - (1 - lambda) * pmax(0.0, logt - k[3])^3)
      ),
      deriv_basis_evals = rbind(
        (1 - trial.dat$treat),
        (1 - trial.dat$treat) * 3 * (pmax(0.0, logt - k[2])^2 - lambda * pmax(0.0, logt - k[1])^2 - (1 - lambda) * pmax(0.0, logt - k[3])^2),
        trial.dat$treat,
        trial.dat$treat * 3 * (pmax(0.0, logt - k[2])^2 - lambda * pmax(0.0, logt - k[1])^2 - (1 - lambda) * pmax(0.0, logt - k[3])^2)
      ),
      N_pred = n_pred_feed,
      log_pred_t = log(1:n_pred_feed),
      basis_c = rbind(
        log(1:n_pred_feed),
        (pmax(0.0, log(1:n_pred_feed) - k[2])^3 - lambda * pmax(0.0, log(1:n_pred_feed) - k[1])^3 - (1 - lambda) * pmax(0.0, log(1:n_pred_feed) - k[3])^3),
        rep(0.0, n_pred_feed),
        rep(0.0, n_pred_feed)
      ),
      basis_t = rbind(
        rep(0.0, n_pred_feed),
        rep(0.0, n_pred_feed),
        log(1:n_pred_feed),
        (pmax(0.0, log(1:n_pred_feed) - k[2])^3 - lambda * pmax(0.0, log(1:n_pred_feed) - k[1])^3 - (1 - lambda) * pmax(0.0, log(1:n_pred_feed) - k[3])^3)
      ),
      deriv_basis_c = rbind(
        rep(1, n_pred_feed),
        3 * (pmax(0.0, log(1:n_pred_feed) - k[2])^2 - lambda * pmax(0.0, log(1:n_pred_feed) - k[1])^2 - (1 - lambda) * pmax(0.0, log(1:n_pred_feed) - k[3])^2),
        rep(0.0, n_pred_feed),
        rep(0.0, n_pred_feed)
      ),
      deriv_basis_t = rbind(
        rep(0.0, n_pred_feed),
        rep(0.0, n_pred_feed),
        rep(1, n_pred_feed),
        3 * (pmax(0.0, log(1:n_pred_feed) - k[2])^2 - lambda * pmax(0.0, log(1:n_pred_feed) - k[1])^2 - (1 - lambda) * pmax(0.0, log(1:n_pred_feed) - k[3])^2)
      )
      )
    
    ## Data needed for one-arm version
    data_feed_onearm <-
      list(basis_evals = rbind(
        trial.dat$treat * logt,
        trial.dat$treat * (pmax(0.0, logt - k[2])^3 - lambda * pmax(0.0, logt - k[1])^3 - (1 - lambda) * pmax(0.0, logt - k[3])^3)
      ),
      deriv_basis_evals = rbind(
        trial.dat$treat,
        trial.dat$treat * 3 * (pmax(0.0, logt - k[2])^2 - lambda * pmax(0.0, logt - k[1])^2 - (1 - lambda) * pmax(0.0, logt - k[3])^2)
      ),
      N_pred = n_pred_feed,
      log_pred_t = log(1:n_pred_feed),
      basis_t = rbind(
        log(1:n_pred_feed),
        (pmax(0.0, log(1:n_pred_feed) - k[2])^3 - lambda * pmax(0.0, log(1:n_pred_feed) - k[1])^3 - (1 - lambda) * pmax(0.0, log(1:n_pred_feed) - k[3])^3)
      ),
      deriv_basis_t = rbind(
        rep(1, n_pred_feed),
        3 * (pmax(0.0, log(1:n_pred_feed) - k[2])^2 - lambda * pmax(0.0, log(1:n_pred_feed) - k[1])^2 - (1 - lambda) * pmax(0.0, log(1:n_pred_feed) - k[3])^2)
      )
      )
    
  }
  else if((length(k)-1)==3){
    
    ## Data needed for two-arm version
    data_feed_twoarm <-
      list(basis_evals = rbind(
        (1 - trial.dat$treat) * logt,
        (1 - trial.dat$treat) * (pmax(0.0, logt - k[2])^3 - lambda[1] * pmax(0.0, logt - k[1])^3 - (1 - lambda[1]) * pmax(0.0, logt - k[4])^3),
        (1 - trial.dat$treat) * (pmax(0.0, logt - k[3])^3 - lambda[2] * pmax(0.0, logt - k[1])^3 - (1 - lambda[2]) * pmax(0.0, logt - k[4])^3),
        trial.dat$treat * logt,
        trial.dat$treat * (pmax(0.0, logt - k[2])^3 - lambda[1] * pmax(0.0, logt - k[1])^3 - (1 - lambda[1]) * pmax(0.0, logt - k[4])^3),
        trial.dat$treat * (pmax(0.0, logt - k[3])^3 - lambda[2] * pmax(0.0, logt - k[1])^3 - (1 - lambda[2]) * pmax(0.0, logt - k[4])^3)
      ),
      deriv_basis_evals = rbind(
        (1 - trial.dat$treat),
        (1 - trial.dat$treat) * 3 * (pmax(0.0, logt - k[2])^2 - lambda[1] * pmax(0.0, logt - k[1])^2 - (1 - lambda[1]) * pmax(0.0, logt - k[4])^2),
        (1 - trial.dat$treat) * 3 * (pmax(0.0, logt - k[3])^2 - lambda[2] * pmax(0.0, logt - k[1])^2 - (1 - lambda[2]) * pmax(0.0, logt - k[4])^2),
        trial.dat$treat,
        trial.dat$treat * 3 * (pmax(0.0, logt - k[2])^2 - lambda[1] * pmax(0.0, logt - k[1])^2 - (1 - lambda[1]) * pmax(0.0, logt - k[4])^2),
        trial.dat$treat * 3 * (pmax(0.0, logt - k[3])^2 - lambda[2] * pmax(0.0, logt - k[1])^2 - (1 - lambda[2]) * pmax(0.0, logt - k[4])^2)
      ),
      N_pred = n_pred_feed,
      log_pred_t = log(1:n_pred_feed),
      basis_c = rbind(
        log(1:n_pred_feed),
        (pmax(0.0, log(1:n_pred_feed) - k[2])^3 - lambda[1] * pmax(0.0, log(1:n_pred_feed) - k[1])^3 - (1 - lambda[1]) * pmax(0.0, log(1:n_pred_feed) - k[4])^3),
        (pmax(0.0, log(1:n_pred_feed) - k[3])^3 - lambda[2] * pmax(0.0, log(1:n_pred_feed) - k[1])^3 - (1 - lambda[2]) * pmax(0.0, log(1:n_pred_feed) - k[4])^3),
        rep(0.0, n_pred_feed),
        rep(0.0, n_pred_feed),
        rep(0.0, n_pred_feed)
      ),
      basis_t = rbind(
        rep(0.0, n_pred_feed),
        rep(0.0, n_pred_feed),
        rep(0.0, n_pred_feed),
        log(1:n_pred_feed),
        (pmax(0.0, log(1:n_pred_feed) - k[2])^3 - lambda[1] * pmax(0.0, log(1:n_pred_feed) - k[1])^3 - (1 - lambda[1]) * pmax(0.0, log(1:n_pred_feed) - k[4])^3),
        (pmax(0.0, log(1:n_pred_feed) - k[3])^3 - lambda[2] * pmax(0.0, log(1:n_pred_feed) - k[1])^3 - (1 - lambda[2]) * pmax(0.0, log(1:n_pred_feed) - k[4])^3)
      ),
      deriv_basis_c = rbind(
        rep(1, n_pred_feed),
        3 * (pmax(0.0, log(1:n_pred_feed) - k[2])^2 - lambda[1] * pmax(0.0, log(1:n_pred_feed) - k[1])^2 - (1 - lambda[1]) * pmax(0.0, log(1:n_pred_feed) - k[4])^2),
        3 * (pmax(0.0, log(1:n_pred_feed) - k[3])^2 - lambda[2] * pmax(0.0, log(1:n_pred_feed) - k[1])^2 - (1 - lambda[2]) * pmax(0.0, log(1:n_pred_feed) - k[4])^2),
        rep(0.0, n_pred_feed),
        rep(0.0, n_pred_feed),
        rep(0.0, n_pred_feed)
      ),
      deriv_basis_t = rbind(
        rep(0.0, n_pred_feed),
        rep(0.0, n_pred_feed),
        rep(0.0, n_pred_feed),
        rep(1, n_pred_feed),
        3 * (pmax(0.0, log(1:n_pred_feed) - k[2])^2 - lambda[1] * pmax(0.0, log(1:n_pred_feed) - k[1])^2 - (1 - lambda[1]) * pmax(0.0, log(1:n_pred_feed) - k[4])^2),
        3 * (pmax(0.0, log(1:n_pred_feed) - k[3])^2 - lambda[2] * pmax(0.0, log(1:n_pred_feed) - k[1])^2 - (1 - lambda[2]) * pmax(0.0, log(1:n_pred_feed) - k[4])^2)
      )
      )
    
    ## Data needed for one-arm version
    data_feed_onearm <-
      list(basis_evals = rbind(
        trial.dat$treat * logt,
        trial.dat$treat * (pmax(0.0, logt - k[2])^3 - lambda[1] * pmax(0.0, logt - k[1])^3 - (1 - lambda[1]) * pmax(0.0, logt - k[4])^3),
        trial.dat$treat * (pmax(0.0, logt - k[3])^3 - lambda[2] * pmax(0.0, logt - k[1])^3 - (1 - lambda[2]) * pmax(0.0, logt - k[4])^3)
      ),
      deriv_basis_evals = rbind(
        trial.dat$treat,
        trial.dat$treat * 3 * (pmax(0.0, logt - k[2])^2 - lambda[1] * pmax(0.0, logt - k[1])^2 - (1 - lambda[1]) * pmax(0.0, logt - k[4])^2),
        trial.dat$treat * 3 * (pmax(0.0, logt - k[3])^2 - lambda[2] * pmax(0.0, logt - k[1])^2 - (1 - lambda[2]) * pmax(0.0, logt - k[4])^2)
      ),
      N_pred = n_pred_feed,
      log_pred_t = log(1:n_pred_feed),
      basis_t = rbind(
        log(1:n_pred_feed),
        (pmax(0.0, log(1:n_pred_feed) - k[2])^3 - lambda[1] * pmax(0.0, log(1:n_pred_feed) - k[1])^3 - (1 - lambda[1]) * pmax(0.0, log(1:n_pred_feed) - k[4])^3),
        (pmax(0.0, log(1:n_pred_feed) - k[3])^3 - lambda[2] * pmax(0.0, log(1:n_pred_feed) - k[1])^3 - (1 - lambda[2]) * pmax(0.0, log(1:n_pred_feed) - k[4])^3)
      ),
      deriv_basis_t = rbind(
        rep(1, n_pred_feed),
        3 * (pmax(0.0, log(1:n_pred_feed) - k[2])^2 - lambda[1] * pmax(0.0, log(1:n_pred_feed) - k[1])^2 - (1 - lambda[1]) * pmax(0.0, log(1:n_pred_feed) - k[4])^2),
        3 * (pmax(0.0, log(1:n_pred_feed) - k[3])^2 - lambda[2] * pmax(0.0, log(1:n_pred_feed) - k[1])^2 - (1 - lambda[2]) * pmax(0.0, log(1:n_pred_feed) - k[4])^2)
      )
      )
    
  }
  else {
    ## Data needed for two-arm version
    data_feed_twoarm <-
      list(basis_evals = rbind(
        (1 - trial.dat$treat) * logt,
        (1 - trial.dat$treat) * (pmax(0.0, logt - k[2])^3 - lambda[1] * pmax(0.0, logt - k[1])^3 - (1 - lambda[1]) * pmax(0.0, logt - k[5])^3),
        (1 - trial.dat$treat) * (pmax(0.0, logt - k[3])^3 - lambda[2] * pmax(0.0, logt - k[1])^3 - (1 - lambda[2]) * pmax(0.0, logt - k[5])^3),
        (1 - trial.dat$treat) * (pmax(0.0, logt - k[4])^3 - lambda[3] * pmax(0.0, logt - k[1])^3 - (1 - lambda[3]) * pmax(0.0, logt - k[5])^3),
        trial.dat$treat * logt,
        trial.dat$treat * (pmax(0.0, logt - k[2])^3 - lambda[1] * pmax(0.0, logt - k[1])^3 - (1 - lambda[1]) * pmax(0.0, logt - k[5])^3),
        trial.dat$treat * (pmax(0.0, logt - k[3])^3 - lambda[2] * pmax(0.0, logt - k[1])^3 - (1 - lambda[2]) * pmax(0.0, logt - k[5])^3),
        trial.dat$treat * (pmax(0.0, logt - k[4])^3 - lambda[3] * pmax(0.0, logt - k[1])^3 - (1 - lambda[3]) * pmax(0.0, logt - k[5])^3)
      ),
      deriv_basis_evals = rbind(
        (1 - trial.dat$treat),
        (1 - trial.dat$treat) * 3 * (pmax(0.0, logt - k[2])^2 - lambda[1] * pmax(0.0, logt - k[1])^2 - (1 - lambda[1]) * pmax(0.0, logt - k[5])^2),
        (1 - trial.dat$treat) * 3 * (pmax(0.0, logt - k[3])^2 - lambda[2] * pmax(0.0, logt - k[1])^2 - (1 - lambda[2]) * pmax(0.0, logt - k[5])^2),
        (1 - trial.dat$treat) * 3 * (pmax(0.0, logt - k[4])^2 - lambda[3] * pmax(0.0, logt - k[1])^2 - (1 - lambda[3]) * pmax(0.0, logt - k[5])^2),
        trial.dat$treat,
        trial.dat$treat * 3 * (pmax(0.0, logt - k[2])^2 - lambda[1] * pmax(0.0, logt - k[1])^2 - (1 - lambda[1]) * pmax(0.0, logt - k[5])^2),
        trial.dat$treat * 3 * (pmax(0.0, logt - k[3])^2 - lambda[2] * pmax(0.0, logt - k[1])^2 - (1 - lambda[2]) * pmax(0.0, logt - k[5])^2),
        trial.dat$treat * 3 * (pmax(0.0, logt - k[4])^2 - lambda[3] * pmax(0.0, logt - k[1])^2 - (1 - lambda[3]) * pmax(0.0, logt - k[5])^2)
      ),
      N_pred = n_pred_feed,
      log_pred_t = log(1:n_pred_feed),
      basis_c = rbind(
        log(1:n_pred_feed),
        (pmax(0.0, log(1:n_pred_feed) - k[2])^3 - lambda[1] * pmax(0.0, log(1:n_pred_feed) - k[1])^3 - (1 - lambda[1]) * pmax(0.0, log(1:n_pred_feed) - k[5])^3),
        (pmax(0.0, log(1:n_pred_feed) - k[3])^3 - lambda[2] * pmax(0.0, log(1:n_pred_feed) - k[1])^3 - (1 - lambda[2]) * pmax(0.0, log(1:n_pred_feed) - k[5])^3),
        (pmax(0.0, log(1:n_pred_feed) - k[4])^3 - lambda[3] * pmax(0.0, log(1:n_pred_feed) - k[1])^3 - (1 - lambda[3]) * pmax(0.0, log(1:n_pred_feed) - k[5])^3),
        rep(0.0, n_pred_feed),
        rep(0.0, n_pred_feed),
        rep(0.0, n_pred_feed),
        rep(0.0, n_pred_feed)
      ),
      basis_t = rbind(
        rep(0.0, n_pred_feed),
        rep(0.0, n_pred_feed),
        rep(0.0, n_pred_feed),
        rep(0.0, n_pred_feed),
        log(1:n_pred_feed),
        (pmax(0.0, log(1:n_pred_feed) - k[2])^3 - lambda[1] * pmax(0.0, log(1:n_pred_feed) - k[1])^3 - (1 - lambda[1]) * pmax(0.0, log(1:n_pred_feed) - k[5])^3),
        (pmax(0.0, log(1:n_pred_feed) - k[3])^3 - lambda[2] * pmax(0.0, log(1:n_pred_feed) - k[1])^3 - (1 - lambda[2]) * pmax(0.0, log(1:n_pred_feed) - k[5])^3),
        (pmax(0.0, log(1:n_pred_feed) - k[4])^3 - lambda[3] * pmax(0.0, log(1:n_pred_feed) - k[1])^3 - (1 - lambda[3]) * pmax(0.0, log(1:n_pred_feed) - k[5])^3)
      ),
      deriv_basis_c = rbind(
        rep(1, n_pred_feed),
        3 * (pmax(0.0, log(1:n_pred_feed) - k[2])^2 - lambda[1] * pmax(0.0, log(1:n_pred_feed) - k[1])^2 - (1 - lambda[1]) * pmax(0.0, log(1:n_pred_feed) - k[5])^2),
        3 * (pmax(0.0, log(1:n_pred_feed) - k[3])^2 - lambda[2] * pmax(0.0, log(1:n_pred_feed) - k[1])^2 - (1 - lambda[2]) * pmax(0.0, log(1:n_pred_feed) - k[5])^2),
        3 * (pmax(0.0, log(1:n_pred_feed) - k[4])^2 - lambda[3] * pmax(0.0, log(1:n_pred_feed) - k[1])^2 - (1 - lambda[3]) * pmax(0.0, log(1:n_pred_feed) - k[5])^2),
        rep(0.0, n_pred_feed),
        rep(0.0, n_pred_feed),
        rep(0.0, n_pred_feed),
        rep(0.0, n_pred_feed)
      ),
      deriv_basis_t = rbind(
        rep(0.0, n_pred_feed),
        rep(0.0, n_pred_feed),
        rep(0.0, n_pred_feed),
        rep(0.0, n_pred_feed),
        rep(1, n_pred_feed),
        3 * (pmax(0.0, log(1:n_pred_feed) - k[2])^2 - lambda[1] * pmax(0.0, log(1:n_pred_feed) - k[1])^2 - (1 - lambda[1]) * pmax(0.0, log(1:n_pred_feed) - k[5])^2),
        3 * (pmax(0.0, log(1:n_pred_feed) - k[3])^2 - lambda[2] * pmax(0.0, log(1:n_pred_feed) - k[1])^2 - (1 - lambda[2]) * pmax(0.0, log(1:n_pred_feed) - k[5])^2),
        3 * (pmax(0.0, log(1:n_pred_feed) - k[4])^2 - lambda[3] * pmax(0.0, log(1:n_pred_feed) - k[1])^2 - (1 - lambda[3]) * pmax(0.0, log(1:n_pred_feed) - k[5])^2)
      )
      )
    
    ## Data needed for one-arm version
    data_feed_onearm <-
      list(basis_evals = rbind(
        trial.dat$treat * logt,
        trial.dat$treat * (pmax(0.0, logt - k[2])^3 - lambda[1] * pmax(0.0, logt - k[1])^3 - (1 - lambda[1]) * pmax(0.0, logt - k[5])^3),
        trial.dat$treat * (pmax(0.0, logt - k[3])^3 - lambda[2] * pmax(0.0, logt - k[1])^3 - (1 - lambda[2]) * pmax(0.0, logt - k[5])^3),
        trial.dat$treat * (pmax(0.0, logt - k[4])^3 - lambda[3] * pmax(0.0, logt - k[1])^3 - (1 - lambda[3]) * pmax(0.0, logt - k[5])^3)
      ),
      deriv_basis_evals = rbind(
        trial.dat$treat,
        trial.dat$treat * 3 * (pmax(0.0, logt - k[2])^2 - lambda[1] * pmax(0.0, logt - k[1])^2 - (1 - lambda[1]) * pmax(0.0, logt - k[5])^2),
        trial.dat$treat * 3 * (pmax(0.0, logt - k[3])^2 - lambda[2] * pmax(0.0, logt - k[1])^2 - (1 - lambda[2]) * pmax(0.0, logt - k[5])^2),
        trial.dat$treat * 3 * (pmax(0.0, logt - k[4])^2 - lambda[3] * pmax(0.0, logt - k[1])^2 - (1 - lambda[3]) * pmax(0.0, logt - k[5])^2)
      ),
      N_pred = n_pred_feed,
      log_pred_t = log(1:n_pred_feed),
      basis_t = rbind(
        log(1:n_pred_feed),
        (pmax(0.0, log(1:n_pred_feed) - k[2])^3 - lambda[1] * pmax(0.0, log(1:n_pred_feed) - k[1])^3 - (1 - lambda[1]) * pmax(0.0, log(1:n_pred_feed) - k[5])^3),
        (pmax(0.0, log(1:n_pred_feed) - k[3])^3 - lambda[2] * pmax(0.0, log(1:n_pred_feed) - k[1])^3 - (1 - lambda[2]) * pmax(0.0, log(1:n_pred_feed) - k[5])^3),
        (pmax(0.0, log(1:n_pred_feed) - k[4])^3 - lambda[3] * pmax(0.0, log(1:n_pred_feed) - k[1])^3 - (1 - lambda[3]) * pmax(0.0, log(1:n_pred_feed) - k[5])^3)
      ),
      deriv_basis_t = rbind(
        rep(1, n_pred_feed),
        3 * (pmax(0.0, log(1:n_pred_feed) - k[2])^2 - lambda[1] * pmax(0.0, log(1:n_pred_feed) - k[1])^2 - (1 - lambda[1]) * pmax(0.0, log(1:n_pred_feed) - k[5])^2),
        3 * (pmax(0.0, log(1:n_pred_feed) - k[3])^2 - lambda[2] * pmax(0.0, log(1:n_pred_feed) - k[1])^2 - (1 - lambda[2]) * pmax(0.0, log(1:n_pred_feed) - k[5])^2),
        3 * (pmax(0.0, log(1:n_pred_feed) - k[4])^2 - lambda[3] * pmax(0.0, log(1:n_pred_feed) - k[1])^2 - (1 - lambda[3]) * pmax(0.0, log(1:n_pred_feed) - k[5])^2)
      )
      )
    
  }
  
  if(gen_pop_adjust == TRUE){
    data_feed_genpop <-
      list(n_GP = genpop.data$n,
           r_GP = genpop.data$r,
           r_GP_t = genpop.data$r,
           idx1_GP = genpop.data$t_n,
           idx2_GP = genpop.data$t_r,
           GP_choice = gp_choice_int)
  }
  
  if(hr_adjust == TRUE && single_arm_mode == FALSE){
    data_feed_hr <-
      list(N_hr = N_hr_feed,
           hr_mean = rep(hr_specs[1], N_hr_feed),
           hr_sd = rep(hr_specs[2], N_hr_feed),
           hr_idx = (hr_specs[3]:hr_specs[4])*12)
  }
  
  if(single_arm_mode==TRUE){
    data_feed <- append(data_feed_base,data_feed_onearm)
  }
  else {
    data_feed <- append(data_feed_base,data_feed_twoarm)
  }
  
  if(gen_pop_adjust==TRUE){
    data_feed <- append(data_feed,data_feed_genpop)
  }
  
  if(hr_adjust == TRUE && single_arm_mode==FALSE){
    data_feed <- append(data_feed,data_feed_hr)
  }
  
  ## Set up initialisation list for three chains in Stan
  if(single_arm_mode==TRUE){
    init_feed <- list(list(gammas = c(coef(spl_dat[[1]])[2],rep(0,length(k)-2)),
                           gamma_intercept = coef(spl_dat[[1]])[1]),
                      list(gammas = c(coef(spl_dat[[2]])[2],rep(0,length(k)-2)),
                           gamma_intercept = coef(spl_dat[[2]])[1]),
                      list(gammas = c(coef(spl_dat[[3]])[2],rep(0,length(k)-2)),
                           gamma_intercept = coef(spl_dat[[3]])[1]))
  }
  else {
    init_feed <- list(list(gammas = c(coef(spl_dat[[1]])[2],rep(0,length(k)-2),
                                      sum(coef(spl_dat[[1]])[2:3]),rep(0,length(k)-2)),
                           gamma_intercept = coef(spl_dat[[1]])[1]),
                      list(gammas = c(coef(spl_dat[[2]])[2],rep(0,length(k)-2),
                                      sum(coef(spl_dat[[2]])[2:3]),rep(0,length(k)-2)),
                           gamma_intercept = coef(spl_dat[[2]])[1]),
                      list(gammas = c(coef(spl_dat[[3]])[2],rep(0,length(k)-2),
                                      sum(coef(spl_dat[[3]])[2:3]),rep(0,length(k)-2)),
                           gamma_intercept = coef(spl_dat[[3]])[1]))
  }
  
  # Run Stan model
  ## Determine the correct Stan code
  if(single_arm_mode == TRUE & gen_pop_adjust == TRUE){
    guyot.tran <- stanc("MPES Stan - single arm.stan")
  } else if(single_arm_mode == TRUE & gen_pop_adjust == FALSE){
    guyot.tran <- stanc("MPES Stan - single arm no gp.stan")
  } else if(hr_adjust == TRUE & gen_pop_adjust == TRUE){
    guyot.tran <- stanc("MPES Stan - two arm.stan")
  } else if(hr_adjust == TRUE & gen_pop_adjust == FALSE){
    guyot.tran <- stanc("MPES Stan - two arm no gp.stan")
  } else if(hr_adjust == FALSE & gen_pop_adjust == TRUE){
    guyot.tran <- stanc("MPES Stan - two arm no hr.stan")
  } else {
    guyot.tran <- stanc("MPES Stan - two arm no gp no hr.stan")
  }
  guyot.comp <- stan_model(stanc = guyot.tran)
  
  ## Sample from Stan model
  pars_feed <- if(single_arm_mode==TRUE){
    c("gammas", "gamma_intercept","S_t", "SC_t","log_lik")
  } else {
    c("gammas", "gamma_intercept","S_c", "S_t", "hr", "SC_c", "SC_t","log_lik")
  }
  
  ## Execute Stan model
  mpes_output <- sampling(guyot.comp, # this line executes the model (i.e., this is the line that takes a long time to run!)
                          data = data_feed, pars = pars_feed, init = init_feed,
                          warmup = warmup_spec, iter = iter_spec,
                          chains = 3, cores = 3)
  
  ## Extract relevant parameters in a data frame
  mpes_summary <- summary(mpes_output, pars = c("gammas", "gamma_intercept", "log_lik"))$summary #summary results
  n_pred_feed <- nrow(as.data.frame(summary(mpes_output, pars = c("S_c"))$summary))
  fit_surv <- bind_rows(mutate(as.data.frame(summary(mpes_output, pars = c("S_c"))$summary),
                               Strata = "treat=0", time = 1:n_pred_feed), #
                        mutate(as.data.frame(summary(mpes_output, pars = c("S_t"))$summary),
                               Strata = "treat=1", time = 1:n_pred_feed)) %>% #n_pred_feed
    mutate(Strata = factor(Strata)) %>%
    rename(lcl = "2.5%", ucl = "97.5%") %>%
    select(Strata, time, mean, lcl, ucl)
  
  ## Return output
  return(fit_surv)
  
} #choose whether to apply the general population contraint to the active arm (TRUE) or not (FALSE)


# Execute Guyot MPES # ----
# Code below for running the base-case models. Comment out and use the manual override code instead if wanting to run a custom model

guyot_start_time <- Sys.time()

if(case_study == "SCCHN"){
  # Base case for SCCHN
  mpes_output <- run_mpes(hr_adjust = TRUE, gen_pop_adjust = TRUE, gp_specs = c(34, 30196, 23080),
                          hr_specs = c(1, 0.1, 6, 35), knot_loc = c(0, 3.39, 5.23, 6.58), iter_spec = 1000, gp_choice = FALSE)
} else if(case_study == "Melanoma"){
  # Panel A
  # mpes_output <- run_mpes(hr_adjust = FALSE, gen_pop_adjust = FALSE, knot_loc = c(0, log(36), log(72), log(20*12)), iter_spec = 1000, gp_choice = FALSE)
  # Panel B
  # mpes_output <- run_mpes(hr_adjust = TRUE, gen_pop_adjust = FALSE, hr_specs = c(1, 0.1, 5, 20), knot_loc = c(0, log(36), log(72), log(20*12)), iter_spec = 1000, gp_choice = FALSE)
  # Panel C
  mpes_output <- run_mpes(hr_adjust = TRUE, gen_pop_adjust = TRUE, gp_specs = c(30, 45791, 41875), hr_specs = c(1, 0.1, 5, 20), knot_loc = c(0, log(36), log(72), log(20*12)), iter_spec = 1000, gp_choice = FALSE)
} else {
  # Panel A (check reduce_ext switch)
  mpes_output <- run_mpes(hr_adjust = TRUE, gen_pop_adjust = TRUE, gp_specs = c(30, 45791, 41875),
                          hr_specs = c(1, 0.1, 10, 20), knot_loc = c(0, log(24), log(72), log(20*12)), iter_spec = 1000, gp_choice = FALSE)
  # Panels B and C (check reduce_ext switch)
  # mpes_output <- run_mpes(hr_adjust = FALSE, gen_pop_adjust = TRUE, gp_specs = c(30, 45791, 41875), hr_specs = c(1, 0.1, 10, 20), knot_loc = c(0, log(24), log(72), log(20*12)), iter_spec = 1000, gp_choice = TRUE)
}

# Manual override
# mpes_output <- run_mpes()

guyot_end_time <- Sys.time()

# Extract Guyot MPES outputs # ----

fit_surv_extract_a <- subset(mpes_output, Strata == "treat=1") #intervention only
fit_surv_extract_c <- subset(mpes_output, Strata == "treat=0") #comparator only
fit_surv_extract <- mpes_output#both arms
mod <- data.frame(time = fit_surv_extract_a$time,
                  surv = fit_surv_extract_a$mean,
                  lcl = fit_surv_extract_a$lcl,
                  ucl = fit_surv_extract_a$ucl) #for later use
mod$lbl <- "Fitted model (active)"

mod_c <- data.frame(time = fit_surv_extract_c$time,
                    surv = fit_surv_extract_c$mean,
                    lcl = fit_surv_extract_c$lcl,
                    ucl = fit_surv_extract_c$ucl) #for later use
mod_c$lbl <- "Fitted model (control)"

km_function_line <- broom::tidy(km_function) %>% mutate(lbl = "Kaplan-Meier")
km_function.long_line <- broom::tidy(km_function.long) %>% mutate(lbl = "Kaplan-Meier")
km_function_line$lbl <- "Original KM (active)"
km_function.long_line$lbl <- "Updated KM (active)"

km_function_line_c <- broom::tidy(km_function_c) %>% mutate(lbl = "Kaplan-Meier")
km_function.long_line_c <- broom::tidy(km_function.long_c) %>% mutate(lbl = "Kaplan-Meier")
km_function_line_c$lbl <- "Original KM (control)"
km_function.long_line_c$lbl <- "Updated KM (control)"

col_list <- c(
  "Fitted model (active)" = "tomato4",
  "Fitted model (control)" = "royalblue4",
  "Original KM (active)" = "tomato3",
  "Updated KM (active)" = "tomato4",
  "Original KM (control)" = "royalblue3",
  "Updated KM (control)" = "royalblue4")

lines_list <- c(
  "Fitted model (active)" = "solid",
  "Fitted model (control)" = "dashed",
  "Original KM (active)" = "solid",
  "Updated KM (active)" = "solid",
  "Original KM (control)" = "dashed",
  "Updated KM (control)" = "dashed")

fill_list <- c(
  "Fitted model (active)" = "tomato1",
  "Fitted model (control)" = "royalblue1",
  "Original KM (active)" = "#FFFFFF",
  "Updated KM (active)" = "#FFFFFF",
  "Original KM (control)" = "#FFFFFF",
  "Updated KM (control)" = "#FFFFFF")

legend_order <- c(
  "Fitted model (active)",
  "Fitted model (control)",
  "Original KM (active)",
  "Updated KM (active)",
  "Original KM (control)",
  "Updated KM (control)")

survival_plot <- ggplot() +
  xlab('Time (years)') + ylab('Survival probability (%)') +
  theme_bw() +
  theme(axis.title.y = element_text(size = 13, angle = 90, face = "bold")) +
  theme(axis.title.x = element_text(size = 13, angle = 0, face = "bold")) +
  theme(axis.text.y = element_text(size = 13)) +
  theme(axis.text.x = element_text(size = 13)) +
  theme(legend.text = element_text(size = 13)) +
  theme(panel.grid = element_blank()) +
  theme(panel.grid.major.x = element_line(color = "grey90", linetype = "solid", size = 0.5)) +
  theme(panel.grid.major.y = element_line(color = "grey90", linetype = "solid", size = 0.5)) +
  theme(plot.margin = margin(1,1,1,1, "cm"),
        legend.position = "top",
        axis.line = element_line(colour = "black"),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  scale_y_continuous(expand = c(0,0), breaks = scales::breaks_extended(n = 11), limits = c(0,100)) +
  scale_x_continuous(expand = c(0,0), breaks = seq(0,output_plot_x_limit*12,12), limits = c(0,output_plot_x_limit*12), labels = function(x) as.integer(x / 12)) +
  guides(col = guide_legend(title = "", nrow = 2),
         linetype = guide_legend(title = "", nrow = 2),
         fill = guide_legend(title = "", nrow = 2)) +
  scale_colour_manual(values = col_list, breaks = legend_order) +
  scale_linetype_manual(values = lines_list, breaks = legend_order) +
  scale_fill_manual(values = fill_list, breaks = legend_order) +
  geom_ribbon(data = km_function.long_line, aes(x = time, y = estimate*100, ymin = conf.low*100, ymax = conf.high*100, fill = lbl), alpha = 0.2) +
  geom_ribbon(data = km_function_line, aes(x = time, y = estimate*100, ymin = conf.low*100, ymax = conf.high*100, fill = lbl), alpha = 0.2) +
  geom_ribbon(data = km_function.long_line_c, aes(x = time, y = estimate*100, ymin = conf.low*100, ymax = conf.high*100, fill = lbl), alpha = 0.2) +
  geom_ribbon(data = km_function_line_c, aes(x = time, y = estimate*100, ymin = conf.low*100, ymax = conf.high*100, fill = lbl), alpha = 0.2) +
  geom_ribbon(data = mod, aes(x = time, y = surv*100, ymin = lcl*100, ymax = ucl*100, fill = lbl), alpha = 0.2) +
  geom_ribbon(data = mod_c, aes(x = time, y = surv*100, ymin = lcl*100, ymax = ucl*100, fill = lbl), alpha = 0.2) +
  geom_step(data = km_function.long_line, aes(x = time, y = estimate*100, col = lbl, lty = lbl), linewidth = 1) +
  geom_step(data = km_function_line, aes(x = time, y = estimate*100, col = lbl, lty = lbl), linewidth = 1) +
  geom_step(data = km_function.long_line_c, aes(x = time, y = estimate*100, col = lbl, lty = lbl), linewidth = 1) +
  geom_step(data = km_function_line_c, aes(x = time, y = estimate*100, col = lbl, lty = lbl), linewidth = 1) +
  geom_line(data = mod, aes(x = time, y = surv*100, col = lbl, lty = lbl), linewidth = 1) +
  geom_line(data = mod_c, aes(x = time, y = surv*100, col = lbl, lty = lbl), linewidth = 1)

# Shows Guyot MPES results # ----
survival_plot

# Execute Jackson MPES # ----
library("survextrap")

run_gp <- function(gp_source = "uk",
                   start_age = 55,
                   prop_fem = 26/60){
  if(gp_source == "uk"){
    life_table <- read.xlsx(xlsxFile = "data/nationallifetables3yruk.xlsx",
                            sheet = "2018-2020", startRow = 6, colNames = TRUE, skipEmptyCols = TRUE) #this is the default data set for uk general population
  } else {
    life_table <- read.csv(gp_source)
  }
  
  life_table$t_plot <- life_table$age-start_age
  life_table_reduced <- subset(life_table,t_plot>=0)
  
  
  #extend so it continues beyond age 100 by carrying forward the last row
  life_table_reduced_add <- life_table_reduced[rep(nrow(life_table_reduced), 31),]
  for(i in 2:nrow(life_table_reduced_add)){
    life_table_reduced_add$age[i] <- life_table_reduced_add$age[i] + i-1
    life_table_reduced_add$t_plot[i] <- life_table_reduced_add$t_plot[i] + i-1
    life_table_reduced_add$lx[i] <- life_table_reduced_add$lx[i-1]*(1-life_table_reduced_add$mx[i-1])
    life_table_reduced_add$'lx.1'[i] <- life_table_reduced_add$'lx.1'[i-1]*(1-life_table_reduced_add$'mx.1'[i-1])
  }
  
  life_table_reduced_add <- subset(life_table_reduced_add,age!=100) #to avoid duplication
  
  life_table_reduced <- rbind(life_table_reduced,life_table_reduced_add)
  
  for(i in 1:nrow(life_table_reduced)){
    if(i==1) { #first row
      life_table_reduced$prop_fem[i] <- prop_fem
      life_table_reduced$lx_both[i] <- life_table_reduced[i,4] * (1 - prop_fem) + life_table_reduced[i,10] * prop_fem #positions 4 and 10 are where mx rates for male and females are presented
    } else { #later rows
      life_table_reduced$prop_fem[i] <- (1 - life_table_reduced[i,8]) * (life_table_reduced$prop_fem[i-1])/((1 - life_table_reduced[i,2]) * (1 - life_table_reduced$prop_fem[i-1]) + (1 - life_table_reduced[i,8]) * (life_table_reduced$prop_fem[i-1]))
      life_table_reduced$lx_both[i] <- life_table_reduced[i,4] * (1 - life_table_reduced$prop_fem[i]) + life_table_reduced[i,10] * life_table_reduced$prop_fem[i] #positions 4 and 10 are where mx rates for male and females are presented
    }
  }
  
  life_table_reduced$S_plot <- life_table_reduced$lx_both/life_table_reduced$lx_both[1]
  life_table_reduced$Strata = "zz_gen pop"
  life_table_reduced$haz <- NA
  
  for(i in c(1:nrow(life_table_reduced)-1)){
    life_table_reduced$haz[i] <- 1-(life_table_reduced$S_plot[i+1]/life_table_reduced$S_plot[i])
  }
  
  life_table_reduced$haz[nrow(life_table_reduced)] <- 1
  
  #due to rounding, it is possible for the extended period to produce a lower hazard in a later year versus a previous year
  #so, we cap the hazard so it has to at least be as large as the previous year:
  
  for(i in 2:nrow(life_table_reduced)){
    if(life_table_reduced$haz[i]<life_table_reduced$haz[i-1]){
      life_table_reduced$haz[i] <- life_table_reduced$haz[i-1]
    }
  }
  
  return(life_table_reduced)
}

extdat <- data.frame(start = ext.data.discrete.raw$t_n,
                     stop = ext.data.discrete.raw$t_r,
                     n = ext.data.discrete.raw$n,
                     r = ext.data.discrete.raw$r,
                     treat = 0)

mspline_data_in_gp <- run_gp()
mspline_data_in <- data.frame(time = mspline_data_in_gp$t_plot*12,
                              hazard = 1 - (1 - mspline_data_in_gp$haz)^(1/12))

jackson_start_time <- Sys.time()

if(jackson_relsurv == TRUE){
  # with relative survival
  nde_mod <- survextrap(Surv(t, dth) ~ treat, data=trial.data.raw, chains=3, external = extdat, backhaz = mspline_data_in)
} else {
  # without relative survival
  nde_mod <- survextrap(Surv(t, dth) ~ treat, data=trial.data.raw, chains=3, external = extdat)
}

jackson_end_time <- Sys.time()

mod_jackson <- data.frame(time = survival(nde_mod, tmax = plot_x_limit*12, newdata = data.frame(treat = 1))$t,
                          surv = survival(nde_mod, tmax = plot_x_limit*12, newdata = data.frame(treat = 1))$median,
                          lcl = survival(nde_mod, tmax = plot_x_limit*12, newdata = data.frame(treat = 1))$lower,
                          ucl = survival(nde_mod, tmax = plot_x_limit*12, newdata = data.frame(treat = 1))$upper) #for later use

mod_jackson_c <- data.frame(time = survival(nde_mod, tmax = plot_x_limit*12, newdata = data.frame(treat = 0))$t,
                            surv = survival(nde_mod, tmax = plot_x_limit*12, newdata = data.frame(treat = 0))$median,
                            lcl = survival(nde_mod, tmax = plot_x_limit*12, newdata = data.frame(treat = 0))$lower,
                            ucl = survival(nde_mod, tmax = plot_x_limit*12, newdata = data.frame(treat = 0))$upper) #for later use

mod_jackson$lbl <- "Fitted model (active)"
mod_jackson_c$lbl <- "Fitted model (control)"

col_list_jackson <- c(
  "Fitted model (active)" = "tomato4",
  "Fitted model (control)" = "royalblue4",
  "Original KM (active)" = "tomato3",
  "Updated KM (active)" = "tomato4",
  "Original KM (control)" = "royalblue3",
  "Updated KM (control)" = "royalblue4")

lines_list_jackson <- c(
  "Fitted model (active)" = "solid",
  "Fitted model (control)" = "dashed",
  "Original KM (active)" = "solid",
  "Updated KM (active)" = "solid",
  "Original KM (control)" = "dashed",
  "Updated KM (control)" = "dashed")

fill_list_jackson <- c(
  "Fitted model (active)" = "tomato1",
  "Fitted model (control)" = "royalblue1",
  "Original KM (active)" = "#FFFFFF",
  "Updated KM (active)" = "#FFFFFF",
  "Original KM (control)" = "#FFFFFF",
  "Updated KM (control)" = "#FFFFFF")

legend_order_jackson <- c(
  "Fitted model (active)",
  "Fitted model (control)",
  "Original KM (active)",
  "Updated KM (active)",
  "Original KM (control)",
  "Updated KM (control)")

survival_plot_jackson <- ggplot() +
  xlab('Time (years)') + ylab('Survival probability (%)') +
  theme_bw() +
  theme(axis.title.y = element_text(size = 13, angle = 90, face = "bold")) +
  theme(axis.title.x = element_text(size = 13, angle = 0, face = "bold")) +
  theme(axis.text.y = element_text(size = 13)) +
  theme(axis.text.x = element_text(size = 13)) +
  theme(legend.text = element_text(size = 13)) +
  theme(panel.grid = element_blank()) +
  theme(panel.grid.major.x = element_line(color = "grey90", linetype = "solid", size = 0.5)) +
  theme(panel.grid.major.y = element_line(color = "grey90", linetype = "solid", size = 0.5)) +
  theme(plot.margin = margin(1,1,1,1, "cm"),
        legend.position = "top",
        axis.line = element_line(colour = "black"),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  scale_y_continuous(expand = c(0,0), breaks = scales::breaks_extended(n = 11), limits = c(0,100)) +
  scale_x_continuous(expand = c(0,0), breaks = seq(0,output_plot_x_limit*12,12), limits = c(0,output_plot_x_limit*12), labels = function(x) as.integer(x / 12)) +
  guides(col = guide_legend(title = "", nrow = 2),
         linetype = guide_legend(title = "", nrow = 2),
         fill = guide_legend(title = "", nrow = 2)) +
  scale_colour_manual(values = col_list_jackson, breaks = legend_order_jackson) +
  scale_linetype_manual(values = lines_list_jackson, breaks = legend_order_jackson) +
  scale_fill_manual(values = fill_list_jackson, breaks = legend_order_jackson) +
  geom_ribbon(data = km_function.long_line, aes(x = time, y = estimate*100, ymin = conf.low*100, ymax = conf.high*100, fill = lbl), alpha = 0.2) +
  geom_ribbon(data = km_function_line, aes(x = time, y = estimate*100, ymin = conf.low*100, ymax = conf.high*100, fill = lbl), alpha = 0.2) +
  geom_ribbon(data = km_function.long_line_c, aes(x = time, y = estimate*100, ymin = conf.low*100, ymax = conf.high*100, fill = lbl), alpha = 0.2) +
  geom_ribbon(data = km_function_line_c, aes(x = time, y = estimate*100, ymin = conf.low*100, ymax = conf.high*100, fill = lbl), alpha = 0.2) +
  geom_ribbon(data = mod_jackson, aes(x = time, y = surv*100, ymin = lcl*100, ymax = ucl*100, fill = lbl), alpha = 0.2) +
  geom_ribbon(data = mod_jackson_c, aes(x = time, y = surv*100, ymin = lcl*100, ymax = ucl*100, fill = lbl), alpha = 0.2) +
  geom_step(data = km_function.long_line, aes(x = time, y = estimate*100, col = lbl, lty = lbl), linewidth = 1) +
  geom_step(data = km_function_line, aes(x = time, y = estimate*100, col = lbl, lty = lbl), linewidth = 1) +
  geom_step(data = km_function.long_line_c, aes(x = time, y = estimate*100, col = lbl, lty = lbl), linewidth = 1) +
  geom_step(data = km_function_line_c, aes(x = time, y = estimate*100, col = lbl, lty = lbl), linewidth = 1) +
  geom_line(data = mod_jackson, aes(x = time, y = surv*100, col = lbl, lty = lbl), linewidth = 1) +
  geom_line(data = mod_jackson_c, aes(x = time, y = surv*100, col = lbl, lty = lbl), linewidth = 1)

# Shows base-case results for Jackson # ----
survival_plot_jackson

# Compare the models # ----

mod_guyot_comp <- mod
mod_guyot_c_comp <- mod_c
mod_jackson_comp <- mod_jackson
mod_jackson_c_comp <- mod_jackson_c

mod_guyot_comp$lbl <- "Fitted model (active, Guyot)"
mod_guyot_c_comp$lbl <- "Fitted model (control, Guyot)"
mod_jackson_comp$lbl <- "Fitted model (active, Jackson)"
mod_jackson_c_comp$lbl <- "Fitted model (control, Jackson)"

col_list_comp <- c(
  "Fitted model (active, Guyot)" = "tomato4",
  "Fitted model (control, Guyot)" = "royalblue4",
  "Fitted model (active, Jackson)" = "seagreen4",
  "Fitted model (control, Jackson)" = "olivedrab4",
  "Original KM (active)" = "black",
  "Updated KM (active)" = "black",
  "Original KM (control)" = "grey50",
  "Updated KM (control)" = "grey50")

lines_list_comp <- c(
  "Fitted model (active, Guyot)" = "solid",
  "Fitted model (control, Guyot)" = "dashed",
  "Fitted model (active, Jackson)" = "solid",
  "Fitted model (control, Jackson)" = "dashed",
  "Original KM (active)" = "solid",
  "Updated KM (active)" = "solid",
  "Original KM (control)" = "dashed",
  "Updated KM (control)" = "dashed")

fill_list_comp <- c(
  "Fitted model (active, Guyot)" = "#FFFFFF",
  "Fitted model (control, Guyot)" = "#FFFFFF",
  "Fitted model (active, Jackson)" = "#FFFFFF",
  "Fitted model (control, Jackson)" = "#FFFFFF",
  "Original KM (active)" = "#FFFFFF",
  "Updated KM (active)" = "#FFFFFF",
  "Original KM (control)" = "#FFFFFF",
  "Updated KM (control)" = "#FFFFFF")

legend_order_comp <- c(
  "Fitted model (active, Guyot)",
  "Fitted model (control, Guyot)",
  "Fitted model (active, Jackson)",
  "Fitted model (control, Jackson)",
  "Original KM (active)",
  "Updated KM (active)",
  "Original KM (control)",
  "Updated KM (control)")


survival_plot_compare <- ggplot() +
  xlab('Time (years)') + ylab('Survival probability (%)') +
  theme_bw() +
  theme(axis.title.y = element_text(size = 13, angle = 90, face = "bold")) +
  theme(axis.title.x = element_text(size = 13, angle = 0, face = "bold")) +
  theme(axis.text.y = element_text(size = 13)) +
  theme(axis.text.x = element_text(size = 13)) +
  theme(legend.text = element_text(size = 13)) +
  theme(panel.grid = element_blank()) +
  theme(panel.grid.major.x = element_line(color = "grey90", linetype = "solid", size = 0.5)) +
  theme(panel.grid.major.y = element_line(color = "grey90", linetype = "solid", size = 0.5)) +
  theme(plot.margin = margin(1,1,1,1, "cm"),
        legend.position = "top",
        axis.line = element_line(colour = "black"),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  scale_y_continuous(expand = c(0,0), breaks = scales::breaks_extended(n = 11), limits = c(0,100)) +
  scale_x_continuous(expand = c(0,0), breaks = seq(0,output_plot_x_limit*12,12), limits = c(0,output_plot_x_limit*12), labels = function(x) as.integer(x / 12)) +
  guides(col = guide_legend(title = "", nrow = 2),
         linetype = guide_legend(title = "", nrow = 2),
         fill = guide_legend(title = "", nrow = 2)) +
  scale_colour_manual(values = col_list_comp, breaks = legend_order_comp) +
  scale_linetype_manual(values = lines_list_comp, breaks = legend_order_comp) +
  scale_fill_manual(values = fill_list_comp, breaks = legend_order_comp) +
  geom_ribbon(data = km_function.long_line, aes(x = time, y = estimate*100, ymin = conf.low*100, ymax = conf.high*100, fill = lbl), alpha = 0) +
  geom_ribbon(data = km_function_line, aes(x = time, y = estimate*100, ymin = conf.low*100, ymax = conf.high*100, fill = lbl), alpha = 0) +
  geom_ribbon(data = km_function.long_line_c, aes(x = time, y = estimate*100, ymin = conf.low*100, ymax = conf.high*100, fill = lbl), alpha = 0) +
  geom_ribbon(data = km_function_line_c, aes(x = time, y = estimate*100, ymin = conf.low*100, ymax = conf.high*100, fill = lbl), alpha = 0) +
  geom_ribbon(data = mod_guyot_comp, aes(x = time, y = surv*100, ymin = lcl*100, ymax = ucl*100, fill = lbl), alpha = 0) +
  geom_ribbon(data = mod_guyot_c_comp, aes(x = time, y = surv*100, ymin = lcl*100, ymax = ucl*100, fill = lbl), alpha = 0) +
  geom_ribbon(data = mod_jackson_comp, aes(x = time, y = surv*100, ymin = lcl*100, ymax = ucl*100, fill = lbl), alpha = 0) +
  geom_ribbon(data = mod_jackson_c_comp, aes(x = time, y = surv*100, ymin = lcl*100, ymax = ucl*100, fill = lbl), alpha = 0) +
  geom_step(data = km_function.long_line, aes(x = time, y = estimate*100, col = lbl, lty = lbl), linewidth = 1) +
  geom_step(data = km_function_line, aes(x = time, y = estimate*100, col = lbl, lty = lbl), linewidth = 1) +
  geom_step(data = km_function.long_line_c, aes(x = time, y = estimate*100, col = lbl, lty = lbl), linewidth = 1) +
  geom_step(data = km_function_line_c, aes(x = time, y = estimate*100, col = lbl, lty = lbl), linewidth = 1) +
  geom_line(data = mod_guyot_comp, aes(x = time, y = surv*100, col = lbl, lty = lbl), linewidth = 1) +
  geom_line(data = mod_guyot_c_comp, aes(x = time, y = surv*100, col = lbl, lty = lbl), linewidth = 1) +
  geom_line(data = mod_jackson_comp, aes(x = time, y = surv*100, col = lbl, lty = lbl), linewidth = 1) +
  geom_line(data = mod_jackson_c_comp, aes(x = time, y = surv*100, col = lbl, lty = lbl), linewidth = 1)

# Shows comparison of Guyot and Jackson # ----
survival_plot_compare

# Same plot but without the legend for the tutorial paper
survival_plot_compare_nolegend <- survival_plot_compare + theme(legend.position = "none")
survival_plot_compare_nolegend

# NICE STA comparisons # ----

if(case_study == "SCCHN"){
  # SCCHN TA145
  ta.comp_df <- data.frame("time" = ta.comp$t_scchn,
                           "surv" = ta.comp$s_scchn,
                           "lcl" = ta.comp$s_scchn-0.01, #set as surv because not available
                           "ucl" = ta.comp$s_scchn+0.01) #set as surv because not available
  ta.comp_df$lbl <- "NICE TA145 model (active)"
  col_list_app <- c("NICE TA145 model (active)" = "black")
  lines_list_app <- c("NICE TA145 model (active)" = "dashed")
  fill_list_app <- c("NICE TA145 model (active)" = "#FFFFFF")
  
} else if(case_study == "Melanoma"){
  # Melanoma TA319
  ta.comp_df <- data.frame("time" = ta.comp$t_mel,
                           "surv" = ta.comp$s_mel,
                           "lcl" = ta.comp$s_mel-0.01, #set as surv because not available
                           "ucl" = ta.comp$s_mel+0.01) #set as surv because not available
  ta.comp_df$lbl <- "NICE TA319 model (active)"
  col_list_app <- c("NICE TA319 model (active)" = "black")
  lines_list_app <- c("NICE TA319 model (active)" = "dashed")
  fill_list_app <- c("NICE TA319 model (active)" = "#FFFFFF")
  
} else {
  # NSCLC TA531
  ta.comp_df <- data.frame("time" = ta.comp$t_nsclc,
                           "surv" = ta.comp$s_nsclc,
                           "lcl" = ta.comp$s_nsclc-0.01, #set as surv because not available
                           "ucl" = ta.comp$s_nsclc+0.01) #set as surv because not available
  ta.comp_df$lbl <- "NICE TA531 model (active)"
  col_list_app <- c("NICE TA531 model (active)" = "black")
  lines_list_app <- c("NICE TA531 model (active)" = "dashed")
  fill_list_app <- c("NICE TA531 model (active)" = "#FFFFFF")
}

legend_order_ta.comp <- legend_order_comp
legend_order_ta.comp <- append(legend_order_ta.comp,ta.comp_df$lbl[1])

col_list_ta.comp <- col_list_comp
col_list_ta.comp <- append(col_list_ta.comp, col_list_app)

lines_list_ta.comp <- lines_list_comp
lines_list_ta.comp <- append(lines_list_ta.comp, lines_list_app)

fill_list_ta.comp <- fill_list_comp
fill_list_ta.comp <- append(fill_list_ta.comp, fill_list_app)

survival_plot_ta.comp <- ggplot() +
  xlab('Time (years)') + ylab('Survival probability (%)') +
  theme_bw() +
  theme(axis.title.y = element_text(size = 13, angle = 90, face = "bold")) +
  theme(axis.title.x = element_text(size = 13, angle = 0, face = "bold")) +
  theme(axis.text.y = element_text(size = 13)) +
  theme(axis.text.x = element_text(size = 13)) +
  theme(legend.text = element_text(size = 13)) +
  theme(panel.grid = element_blank()) +
  theme(panel.grid.major.x = element_line(color = "grey90", linetype = "solid", size = 0.5)) +
  theme(panel.grid.major.y = element_line(color = "grey90", linetype = "solid", size = 0.5)) +
  theme(plot.margin = margin(1,1,1,1, "cm"),
        legend.position = "top",
        axis.line = element_line(colour = "black"),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  scale_y_continuous(expand = c(0,0), breaks = scales::breaks_extended(n = 11), limits = c(0,100)) +
  scale_x_continuous(expand = c(0,0), breaks = seq(0,output_plot_x_limit*12,12), limits = c(0,output_plot_x_limit*12), labels = function(x) as.integer(x / 12)) +
  guides(col = guide_legend(title = "", nrow = 2),
         linetype = guide_legend(title = "", nrow = 2),
         fill = guide_legend(title = "", nrow = 2)) +
  scale_colour_manual(values = col_list_ta.comp, breaks = legend_order_ta.comp) +
  scale_linetype_manual(values = lines_list_ta.comp, breaks = legend_order_ta.comp) +
  scale_fill_manual(values = fill_list_ta.comp, breaks = legend_order_ta.comp) +
  geom_ribbon(data = km_function.long_line, aes(x = time, y = estimate*100, ymin = conf.low*100, ymax = conf.high*100, fill = lbl), alpha = 0.2) +
  geom_ribbon(data = km_function_line, aes(x = time, y = estimate*100, ymin = conf.low*100, ymax = conf.high*100, fill = lbl), alpha = 0.2) +
  geom_ribbon(data = mod_guyot_comp, aes(x = time, y = surv*100, ymin = lcl*100, ymax = ucl*100, fill = lbl), alpha = 0.2) +
  geom_ribbon(data = mod_jackson_comp, aes(x = time, y = surv*100, ymin = lcl*100, ymax = ucl*100, fill = lbl), alpha = 0.2) +
  geom_ribbon(data = ta.comp_df, aes(x = time, y = surv*100, ymin = lcl*100, ymax = ucl*100, fill = lbl), alpha = 0.2) +
  geom_step(data = km_function.long_line, aes(x = time, y = estimate*100, col = lbl, lty = lbl), linewidth = 1) +
  geom_step(data = km_function_line, aes(x = time, y = estimate*100, col = lbl, lty = lbl), linewidth = 1) +
  geom_line(data = mod_guyot_comp, aes(x = time, y = surv*100, col = lbl, lty = lbl), linewidth = 1) +
  geom_line(data = mod_jackson_comp, aes(x = time, y = surv*100, col = lbl, lty = lbl), linewidth = 1) +
  geom_line(data = ta.comp_df, aes(x = time, y = surv*100, col = lbl, lty = lbl), linewidth = 1)

# Shows NICE STA comparison plot # ----
survival_plot_ta.comp

# Run time calcs # ----
end_time <- Sys.time()

#Guyot run time
guyot_end_time - guyot_start_time

#Jackson tun time
jackson_end_time - jackson_start_time

#Total code run time
end_time - start_time
