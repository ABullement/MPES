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
    int N_hr;                                                       // hazard ratio: number of reference points
    vector[N_hr] hr_mean;                                           // hazard ratio: expected HR at reference points
    vector[N_hr] hr_sd;                                             // hazard ratio: estimated HR standard deviation at reference points
    int hr_idx[N_hr];                                               // hazard ratio: prediction index values
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
    vector[N_hr] pred_hr;                                           // predicted hazard ratio at hazard ratio reference points

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
    for (i in 1:N_ext) {
        pred_ext_SC[i] = fmax(fmin(pred_ext_SC[i], 0.9999999),0.0000001);
    }

    pred_hr = hr[hr_idx];
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
    hr_mean ~ normal(pred_hr, hr_sd);                               //  - link longer run hazard ratio at reference points
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
