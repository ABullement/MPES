data {
    int<lower=1> N;                                                 // number of data points
    int<lower=1> m;                                                 // number of basis splines (also number of knots)
    int<lower=0, upper=1> is_censored[N];                           // 0 = event; 1 = censored
    vector[N] log_times;                                            // log(t)
    matrix[2*m,N] basis_evals;                                      // basis functions 1:m = control; (m+1):2m = treatment
    matrix[2*m,N] deriv_basis_evals;                                // derivative of basis functions 1:m = control; (m+1):2m = treatment

    int N_pred;                                                     // number of prediction time points
    vector[N_pred] log_pred_t;                                      // log(t) at each prediction point
    matrix[2*m, N_pred] basis_c;                                    // basis_evals for control arm at each prediction point
    matrix[2*m, N_pred] basis_t;                                    // basis_evals for treatment arm at each prediction point
    matrix[2*m, N_pred] deriv_basis_c;                              // deriv_basis_evals for control arm at each prediction point
    matrix[2*m, N_pred] deriv_basis_t;                              // deriv_basis_evals for treatment arm at each prediction point

    int N_ext;                                                      // number of ext time points
    int n_ext[N_ext];                                               // ext: number alive at each time point
    int r_ext[N_ext];                                               // ext: number alive one year after each time point
    int idx1_ext[N_ext];                                            // ext: prediction index values at each time point
    int idx2_ext[N_ext];                                            // ext: prediction index values one year after each time point
    int n_GP;                                                       // general population: number alive at reference time
    int r_GP;                                                       // general population: number alive one year after reference time (control arm)
    int r_GP_t;                                                     // general population: number alive one year after reference time (treatment arm)
    int idx1_GP;                                                    // general population: prediction index value at reference time
    int idx2_GP;                                                    // general population: prediction index value one year after reference time
    int<lower = 1> GP_choice;                                       // general population: choice about whether to fix treatment arm (1) or not (2)
}

parameters {
    row_vector[2*m] gammas;                                         // regression coefficients for splines
    real gamma_intercept;                                           // intercept gamma
}

transformed parameters {
    vector[N] etas;                                                 // linear predictors for trial data
    vector[N_pred] eta_c;                                           // control linear predictors for prediction
    vector[N_pred] eta_t;                                           // treatment linear predictors for prediction
    vector[N_pred] S_c;                                             // predicted control survival function
    vector[N_pred] S_t;                                             // predicted treatment survival function
    vector[N_pred] h_c;                                             // predicted control hazard function
    vector[N_pred] h_t;                                             // predicted treatment hazard function
    vector[N_pred] hr;                                              // predicted hazard ratio (treatment/control)
    vector[N_pred-12] SC_c;                                         // predicted control conditional survival
    vector[N_pred-12] SC_t;                                         // predicted treatment conditional survival
    vector[N_ext] pred_ext_SC;                                      // predicted control conditional survival at ext time points
    real pred_GP_SC;                                                // predicted control conditional survival at general population reference point
    real pred_GP_SC_t;                                              // predicted treatment conditional survival at general population reference point

    etas = (gammas*basis_evals)' + gamma_intercept;

    eta_c = (gammas*basis_c)' + gamma_intercept;
    eta_t = (gammas*basis_t)' + gamma_intercept;

    S_c = exp(-exp(eta_c));
    S_t = exp(-exp(eta_t));

    h_c = exp(eta_c) .* ((gammas*deriv_basis_c)');
    h_t = exp(eta_t) .* ((gammas*deriv_basis_t)');
    hr = h_t ./ h_c;

    SC_c = S_c[13:N_pred] ./ S_c[1:(N_pred-12)];
    SC_t = S_t[13:N_pred] ./ S_t[1:(N_pred-12)];

    pred_ext_SC = S_c[idx2_ext] ./ S_c[idx1_ext];
    pred_ext_SC = fmax(fmin(pred_ext_SC, 0.9999999),0.0000001);
    pred_GP_SC = S_c[idx2_GP] ./ S_c[idx1_GP];
    pred_GP_SC = fmax(fmin(pred_GP_SC, 0.9999999),0.0000001);
    pred_GP_SC_t = S_t[idx2_GP] ./ S_t[idx1_GP];
    pred_GP_SC_t = fmax(fmin(pred_GP_SC_t, 0.9999999),0.0000001);
}

model {
                                                                    // priors for gamma and gamma_intercept
    gammas ~ normal(0, 10);
    gamma_intercept   ~ normal(0,10);
                                                                    // likelihood for spline model of trial data
    for (i in 1:N) {
        if (is_censored[i] == 0) {
          target +=  etas[i] - exp(etas[i]) - log_times[i] + log(gammas*deriv_basis_evals[,i]);
        }
        else {
           target += -exp(etas[i]);
        }
    }
                                                                    // contraints:
    r_ext ~ binomial(n_ext, pred_ext_SC);                           //  - link control arm to ext conditional survival at time points
    r_GP ~ binomial(n_GP, pred_GP_SC);                              //  - link control arm to general population conditional survival
    if (GP_choice == 1) {
      r_GP_t ~ binomial(n_GP, pred_GP_SC_t);                        //  - link treatment arm to general population conditional survival
    }
}

generated quantities {
    real RMST15_c;                                                  // restricted mean survival time for control at 15 years
    real RMST15_t;                                                  // restricted mean survival time for treatment at 15 years
    real RMST15_diff;                                               // difference in restricted mean survival time at 15 years
    real RMST20_c;                                                  // restricted mean survival time for control at 20 years
    real RMST20_t;                                                  // restricted mean survival time for treatment at 20 years
    real RMST20_diff;                                               // difference in restricted mean survival time at 20 years
    real log_lik;                                                   // log-likelihood

                                                                    // all estimated using trapezoid rule
    RMST15_c = (0.5*1.0) + sum(S_c[1:(15*12-1)]) + (0.5*S_c[15*12]);
    RMST15_t = (0.5*1.0) + sum(S_t[1:(15*12-1)]) + (0.5*S_t[15*12]);
    RMST15_diff = RMST15_t - RMST15_c;
    RMST20_c = (0.5*1.0) + sum(S_c[1:(20*12-1)]) + (0.5*S_c[20*12]);
    RMST20_t = (0.5*1.0) + sum(S_t[1:(20*12-1)]) + (0.5*S_t[20*12]);
    RMST20_diff = RMST20_t - RMST20_c;

    log_lik = sum(-exp(etas) + (1.0 - to_vector(is_censored)) .* (etas - log_times + to_vector(log(gammas*deriv_basis_evals))));
}
