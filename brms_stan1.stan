// generated with brms 2.16.1
functions {
 /* compute correlated group-level effects
  * in the presence of a within-group covariance matrix
  * Args: 
  *   z: matrix of unscaled group-level effects
  *   SD: vector of standard deviation parameters
  *   L: cholesky factor correlation matrix
  *   Lcov: cholesky factor of within-group correlation matrix
  * Returns: 
  *   matrix of scaled group-level effects
  */ 
  matrix scale_r_cor_cov(matrix z, vector SD, matrix L, matrix Lcov) {
    vector[num_elements(z)] z_flat = to_vector(z);
    vector[num_elements(z)] r = rep_vector(0, num_elements(z));
    matrix[rows(L), cols(L)] LC = diag_pre_multiply(SD, L);
    int rows_z = rows(z);
    int rows_L = rows(L);
    // kronecker product of cholesky factors times a vector
    for (icov in 1:rows(Lcov)) {
      for (jcov in 1:icov) {
        if (Lcov[icov, jcov] > 1e-10) { 
          // avoid calculating products between unrelated individuals
          for (i in 1:rows_L) {
            for (j in 1:i) {
              // incremented element of the output vector
              int k = (rows_L * (icov - 1)) + i;
              // applied element of the input vector
              int l = (rows_L * (jcov - 1)) + j;
              r[k] = r[k] + Lcov[icov, jcov] * LC[i, j] * z_flat[l];
            }
          }
        }
      }
    }
    // r is returned in another dimension order than z
    return to_matrix(r, cols(z), rows(z), 0);
  }
  /* multinomial-logit log-PMF
   * Args: 
   *   y: array of integer response values
   *   mu: vector of category logit probabilities
   * Returns:  
   *   a scalar to be added to the log posterior 
   */ 
   real multinomial_logit2_lpmf(int[] y, vector mu) {
     return multinomial_lpmf(y | softmax(mu));
   }
  /* compute monotonic effects
   * Args:
   *   scale: a simplex parameter
   *   i: index to sum over the simplex
   * Returns:
   *   a scalar between 0 and 1
   */
  real mo(vector scale, int i) {
    if (i == 0) {
      return 0;
    } else {
      return rows(scale) * sum(scale[1:i]);
    }
  }
}
data {
  int<lower=1> N;  // total number of observations
  int<lower=1> N_SOcounts;  // number of observations
  int<lower=2> ncat_SOcounts;  // number of categories
  int Y_SOcounts[N_SOcounts, ncat_SOcounts];  // response array
  int trials_SOcounts[N_SOcounts];  // number of trials
  int<lower=1> K_muMF_SOcounts;  // number of population-level effects
  matrix[N_SOcounts, K_muMF_SOcounts] X_muMF_SOcounts;  // population-level design matrix
  int<lower=1> Ksp_muMF_SOcounts;  // number of special effects terms
  int<lower=1> Imo_muMF_SOcounts;  // number of monotonic variables
  int<lower=1> Jmo_muMF_SOcounts[Imo_muMF_SOcounts];  // length of simplexes
  int Xmo_muMF_SOcounts_1[N_SOcounts];  // monotonic variable
  vector[Jmo_muMF_SOcounts[1]] con_simo_muMF_SOcounts_1;  // prior concentration of monotonic simplex
  int<lower=1> K_muMFF_SOcounts;  // number of population-level effects
  matrix[N_SOcounts, K_muMFF_SOcounts] X_muMFF_SOcounts;  // population-level design matrix
  int<lower=1> Ksp_muMFF_SOcounts;  // number of special effects terms
  int<lower=1> Imo_muMFF_SOcounts;  // number of monotonic variables
  int<lower=1> Jmo_muMFF_SOcounts[Imo_muMFF_SOcounts];  // length of simplexes
  int Xmo_muMFF_SOcounts_1[N_SOcounts];  // monotonic variable
  vector[Jmo_muMFF_SOcounts[1]] con_simo_muMFF_SOcounts_1;  // prior concentration of monotonic simplex
  int<lower=1> K_muFMM_SOcounts;  // number of population-level effects
  matrix[N_SOcounts, K_muFMM_SOcounts] X_muFMM_SOcounts;  // population-level design matrix
  int<lower=1> Ksp_muFMM_SOcounts;  // number of special effects terms
  int<lower=1> Imo_muFMM_SOcounts;  // number of monotonic variables
  int<lower=1> Jmo_muFMM_SOcounts[Imo_muFMM_SOcounts];  // length of simplexes
  int Xmo_muFMM_SOcounts_1[N_SOcounts];  // monotonic variable
  vector[Jmo_muFMM_SOcounts[1]] con_simo_muFMM_SOcounts_1;  // prior concentration of monotonic simplex
  int<lower=1> K_muFFMM_SOcounts;  // number of population-level effects
  matrix[N_SOcounts, K_muFFMM_SOcounts] X_muFFMM_SOcounts;  // population-level design matrix
  int<lower=1> Ksp_muFFMM_SOcounts;  // number of special effects terms
  int<lower=1> Imo_muFFMM_SOcounts;  // number of monotonic variables
  int<lower=1> Jmo_muFFMM_SOcounts[Imo_muFFMM_SOcounts];  // length of simplexes
  int Xmo_muFFMM_SOcounts_1[N_SOcounts];  // monotonic variable
  vector[Jmo_muFFMM_SOcounts[1]] con_simo_muFFMM_SOcounts_1;  // prior concentration of monotonic simplex
  int<lower=1> N_IVSOint;  // number of observations
  int Y_IVSOint[N_IVSOint];  // response variable
  int trials_IVSOint[N_IVSOint];  // number of trials
  int<lower=1> K_IVSOint;  // number of population-level effects
  matrix[N_IVSOint, K_IVSOint] X_IVSOint;  // population-level design matrix
  int<lower=1> Ksp_IVSOint;  // number of special effects terms
  int<lower=1> Imo_IVSOint;  // number of monotonic variables
  int<lower=1> Jmo_IVSOint[Imo_IVSOint];  // length of simplexes
  int Xmo_IVSOint_1[N_IVSOint];  // monotonic variable
  vector[Jmo_IVSOint[1]] con_simo_IVSOint_1;  // prior concentration of monotonic simplex
  // data for group-level effects of ID 1
  int<lower=1> N_1;  // number of grouping levels
  int<lower=1> M_1;  // number of coefficients per level
  int<lower=1> J_1_SOcounts[N_SOcounts];  // grouping indicator per observation
  // group-level predictor values
  vector[N_SOcounts] Z_1_muMF_SOcounts_1;
  // data for group-level effects of ID 2
  int<lower=1> N_2;  // number of grouping levels
  int<lower=1> M_2;  // number of coefficients per level
  int<lower=1> J_2_SOcounts[N_SOcounts];  // grouping indicator per observation
  // group-level predictor values
  vector[N_SOcounts] Z_2_muMF_SOcounts_1;
  // data for group-level effects of ID 3
  int<lower=1> N_3;  // number of grouping levels
  int<lower=1> M_3;  // number of coefficients per level
  int<lower=1> J_3_SOcounts[N_SOcounts];  // grouping indicator per observation
  int<lower=1> J_3_IVSOint[N_IVSOint];  // grouping indicator per observation
  matrix[N_3, N_3] Lcov_3;  // cholesky factor of known covariance matrix
  // group-level predictor values
  vector[N_SOcounts] Z_3_muMF_SOcounts_1;
  vector[N_SOcounts] Z_3_muMFF_SOcounts_2;
  vector[N_SOcounts] Z_3_muFMM_SOcounts_3;
  vector[N_SOcounts] Z_3_muFFMM_SOcounts_4;
  vector[N_IVSOint] Z_3_IVSOint_5;
  int<lower=1> NC_3;  // number of group-level correlations
  // data for group-level effects of ID 4
  int<lower=1> N_4;  // number of grouping levels
  int<lower=1> M_4;  // number of coefficients per level
  int<lower=1> J_4_SOcounts[N_SOcounts];  // grouping indicator per observation
  // group-level predictor values
  vector[N_SOcounts] Z_4_muMF_SOcounts_1;
  // data for group-level effects of ID 5
  int<lower=1> N_5;  // number of grouping levels
  int<lower=1> M_5;  // number of coefficients per level
  int<lower=1> J_5_SOcounts[N_SOcounts];  // grouping indicator per observation
  // group-level predictor values
  vector[N_SOcounts] Z_5_muMFF_SOcounts_1;
  // data for group-level effects of ID 6
  int<lower=1> N_6;  // number of grouping levels
  int<lower=1> M_6;  // number of coefficients per level
  int<lower=1> J_6_SOcounts[N_SOcounts];  // grouping indicator per observation
  // group-level predictor values
  vector[N_SOcounts] Z_6_muMFF_SOcounts_1;
  // data for group-level effects of ID 7
  int<lower=1> N_7;  // number of grouping levels
  int<lower=1> M_7;  // number of coefficients per level
  int<lower=1> J_7_SOcounts[N_SOcounts];  // grouping indicator per observation
  // group-level predictor values
  vector[N_SOcounts] Z_7_muMFF_SOcounts_1;
  // data for group-level effects of ID 8
  int<lower=1> N_8;  // number of grouping levels
  int<lower=1> M_8;  // number of coefficients per level
  int<lower=1> J_8_SOcounts[N_SOcounts];  // grouping indicator per observation
  // group-level predictor values
  vector[N_SOcounts] Z_8_muFMM_SOcounts_1;
  // data for group-level effects of ID 9
  int<lower=1> N_9;  // number of grouping levels
  int<lower=1> M_9;  // number of coefficients per level
  int<lower=1> J_9_SOcounts[N_SOcounts];  // grouping indicator per observation
  // group-level predictor values
  vector[N_SOcounts] Z_9_muFMM_SOcounts_1;
  // data for group-level effects of ID 10
  int<lower=1> N_10;  // number of grouping levels
  int<lower=1> M_10;  // number of coefficients per level
  int<lower=1> J_10_SOcounts[N_SOcounts];  // grouping indicator per observation
  // group-level predictor values
  vector[N_SOcounts] Z_10_muFMM_SOcounts_1;
  // data for group-level effects of ID 11
  int<lower=1> N_11;  // number of grouping levels
  int<lower=1> M_11;  // number of coefficients per level
  int<lower=1> J_11_SOcounts[N_SOcounts];  // grouping indicator per observation
  // group-level predictor values
  vector[N_SOcounts] Z_11_muFFMM_SOcounts_1;
  // data for group-level effects of ID 12
  int<lower=1> N_12;  // number of grouping levels
  int<lower=1> M_12;  // number of coefficients per level
  int<lower=1> J_12_SOcounts[N_SOcounts];  // grouping indicator per observation
  // group-level predictor values
  vector[N_SOcounts] Z_12_muFFMM_SOcounts_1;
  // data for group-level effects of ID 13
  int<lower=1> N_13;  // number of grouping levels
  int<lower=1> M_13;  // number of coefficients per level
  int<lower=1> J_13_SOcounts[N_SOcounts];  // grouping indicator per observation
  // group-level predictor values
  vector[N_SOcounts] Z_13_muFFMM_SOcounts_1;
  // data for group-level effects of ID 14
  int<lower=1> N_14;  // number of grouping levels
  int<lower=1> M_14;  // number of coefficients per level
  int<lower=1> J_14_IVSOint[N_IVSOint];  // grouping indicator per observation
  // group-level predictor values
  vector[N_IVSOint] Z_14_IVSOint_1;
  // data for group-level effects of ID 15
  int<lower=1> N_15;  // number of grouping levels
  int<lower=1> M_15;  // number of coefficients per level
  int<lower=1> J_15_IVSOint[N_IVSOint];  // grouping indicator per observation
  // group-level predictor values
  vector[N_IVSOint] Z_15_IVSOint_1;
  // data for group-level effects of ID 16
  int<lower=1> N_16;  // number of grouping levels
  int<lower=1> M_16;  // number of coefficients per level
  int<lower=1> J_16_IVSOint[N_IVSOint];  // grouping indicator per observation
  // group-level predictor values
  vector[N_IVSOint] Z_16_IVSOint_1;
  int prior_only;  // should the likelihood be ignored?
}
transformed data {
  int Kc_muMF_SOcounts = K_muMF_SOcounts - 1;
  matrix[N_SOcounts, Kc_muMF_SOcounts] Xc_muMF_SOcounts;  // centered version of X_muMF_SOcounts without an intercept
  vector[Kc_muMF_SOcounts] means_X_muMF_SOcounts;  // column means of X_muMF_SOcounts before centering
  int Kc_muMFF_SOcounts = K_muMFF_SOcounts - 1;
  matrix[N_SOcounts, Kc_muMFF_SOcounts] Xc_muMFF_SOcounts;  // centered version of X_muMFF_SOcounts without an intercept
  vector[Kc_muMFF_SOcounts] means_X_muMFF_SOcounts;  // column means of X_muMFF_SOcounts before centering
  int Kc_muFMM_SOcounts = K_muFMM_SOcounts - 1;
  matrix[N_SOcounts, Kc_muFMM_SOcounts] Xc_muFMM_SOcounts;  // centered version of X_muFMM_SOcounts without an intercept
  vector[Kc_muFMM_SOcounts] means_X_muFMM_SOcounts;  // column means of X_muFMM_SOcounts before centering
  int Kc_muFFMM_SOcounts = K_muFFMM_SOcounts - 1;
  matrix[N_SOcounts, Kc_muFFMM_SOcounts] Xc_muFFMM_SOcounts;  // centered version of X_muFFMM_SOcounts without an intercept
  vector[Kc_muFFMM_SOcounts] means_X_muFFMM_SOcounts;  // column means of X_muFFMM_SOcounts before centering
  int Kc_IVSOint = K_IVSOint - 1;
  matrix[N_IVSOint, Kc_IVSOint] Xc_IVSOint;  // centered version of X_IVSOint without an intercept
  vector[Kc_IVSOint] means_X_IVSOint;  // column means of X_IVSOint before centering
  for (i in 2:K_muMF_SOcounts) {
    means_X_muMF_SOcounts[i - 1] = mean(X_muMF_SOcounts[, i]);
    Xc_muMF_SOcounts[, i - 1] = X_muMF_SOcounts[, i] - means_X_muMF_SOcounts[i - 1];
  }
  for (i in 2:K_muMFF_SOcounts) {
    means_X_muMFF_SOcounts[i - 1] = mean(X_muMFF_SOcounts[, i]);
    Xc_muMFF_SOcounts[, i - 1] = X_muMFF_SOcounts[, i] - means_X_muMFF_SOcounts[i - 1];
  }
  for (i in 2:K_muFMM_SOcounts) {
    means_X_muFMM_SOcounts[i - 1] = mean(X_muFMM_SOcounts[, i]);
    Xc_muFMM_SOcounts[, i - 1] = X_muFMM_SOcounts[, i] - means_X_muFMM_SOcounts[i - 1];
  }
  for (i in 2:K_muFFMM_SOcounts) {
    means_X_muFFMM_SOcounts[i - 1] = mean(X_muFFMM_SOcounts[, i]);
    Xc_muFFMM_SOcounts[, i - 1] = X_muFFMM_SOcounts[, i] - means_X_muFFMM_SOcounts[i - 1];
  }
  for (i in 2:K_IVSOint) {
    means_X_IVSOint[i - 1] = mean(X_IVSOint[, i]);
    Xc_IVSOint[, i - 1] = X_IVSOint[, i] - means_X_IVSOint[i - 1];
  }
}
parameters {
  vector[Kc_muMF_SOcounts] b_muMF_SOcounts;  // population-level effects
  real Intercept_muMF_SOcounts;  // temporary intercept for centered predictors
  simplex[Jmo_muMF_SOcounts[1]] simo_muMF_SOcounts_1;  // monotonic simplex
  vector[Ksp_muMF_SOcounts] bsp_muMF_SOcounts;  // special effects coefficients
  vector[Kc_muMFF_SOcounts] b_muMFF_SOcounts;  // population-level effects
  real Intercept_muMFF_SOcounts;  // temporary intercept for centered predictors
  simplex[Jmo_muMFF_SOcounts[1]] simo_muMFF_SOcounts_1;  // monotonic simplex
  vector[Ksp_muMFF_SOcounts] bsp_muMFF_SOcounts;  // special effects coefficients
  vector[Kc_muFMM_SOcounts] b_muFMM_SOcounts;  // population-level effects
  real Intercept_muFMM_SOcounts;  // temporary intercept for centered predictors
  simplex[Jmo_muFMM_SOcounts[1]] simo_muFMM_SOcounts_1;  // monotonic simplex
  vector[Ksp_muFMM_SOcounts] bsp_muFMM_SOcounts;  // special effects coefficients
  vector[Kc_muFFMM_SOcounts] b_muFFMM_SOcounts;  // population-level effects
  real Intercept_muFFMM_SOcounts;  // temporary intercept for centered predictors
  simplex[Jmo_muFFMM_SOcounts[1]] simo_muFFMM_SOcounts_1;  // monotonic simplex
  vector[Ksp_muFFMM_SOcounts] bsp_muFFMM_SOcounts;  // special effects coefficients
  vector[Kc_IVSOint] b_IVSOint;  // population-level effects
  real Intercept_IVSOint;  // temporary intercept for centered predictors
  simplex[Jmo_IVSOint[1]] simo_IVSOint_1;  // monotonic simplex
  vector[Ksp_IVSOint] bsp_IVSOint;  // special effects coefficients
  vector<lower=0>[M_1] sd_1;  // group-level standard deviations
  vector[N_1] z_1[M_1];  // standardized group-level effects
  vector<lower=0>[M_2] sd_2;  // group-level standard deviations
  vector[N_2] z_2[M_2];  // standardized group-level effects
  vector<lower=0>[M_3] sd_3;  // group-level standard deviations
  matrix[M_3, N_3] z_3;  // standardized group-level effects
  cholesky_factor_corr[M_3] L_3;  // cholesky factor of correlation matrix
  vector<lower=0>[M_4] sd_4;  // group-level standard deviations
  vector[N_4] z_4[M_4];  // standardized group-level effects
  vector<lower=0>[M_5] sd_5;  // group-level standard deviations
  vector[N_5] z_5[M_5];  // standardized group-level effects
  vector<lower=0>[M_6] sd_6;  // group-level standard deviations
  vector[N_6] z_6[M_6];  // standardized group-level effects
  vector<lower=0>[M_7] sd_7;  // group-level standard deviations
  vector[N_7] z_7[M_7];  // standardized group-level effects
  vector<lower=0>[M_8] sd_8;  // group-level standard deviations
  vector[N_8] z_8[M_8];  // standardized group-level effects
  vector<lower=0>[M_9] sd_9;  // group-level standard deviations
  vector[N_9] z_9[M_9];  // standardized group-level effects
  vector<lower=0>[M_10] sd_10;  // group-level standard deviations
  vector[N_10] z_10[M_10];  // standardized group-level effects
  vector<lower=0>[M_11] sd_11;  // group-level standard deviations
  vector[N_11] z_11[M_11];  // standardized group-level effects
  vector<lower=0>[M_12] sd_12;  // group-level standard deviations
  vector[N_12] z_12[M_12];  // standardized group-level effects
  vector<lower=0>[M_13] sd_13;  // group-level standard deviations
  vector[N_13] z_13[M_13];  // standardized group-level effects
  vector<lower=0>[M_14] sd_14;  // group-level standard deviations
  vector[N_14] z_14[M_14];  // standardized group-level effects
  vector<lower=0>[M_15] sd_15;  // group-level standard deviations
  vector[N_15] z_15[M_15];  // standardized group-level effects
  vector<lower=0>[M_16] sd_16;  // group-level standard deviations
  vector[N_16] z_16[M_16];  // standardized group-level effects
}
transformed parameters {
  vector[N_1] r_1_muMF_SOcounts_1;  // actual group-level effects
  vector[N_2] r_2_muMF_SOcounts_1;  // actual group-level effects
  matrix[N_3, M_3] r_3;  // actual group-level effects
  // using vectors speeds up indexing in loops
  vector[N_3] r_3_muMF_SOcounts_1;
  vector[N_3] r_3_muMFF_SOcounts_2;
  vector[N_3] r_3_muFMM_SOcounts_3;
  vector[N_3] r_3_muFFMM_SOcounts_4;
  vector[N_3] r_3_IVSOint_5;
  vector[N_4] r_4_muMF_SOcounts_1;  // actual group-level effects
  vector[N_5] r_5_muMFF_SOcounts_1;  // actual group-level effects
  vector[N_6] r_6_muMFF_SOcounts_1;  // actual group-level effects
  vector[N_7] r_7_muMFF_SOcounts_1;  // actual group-level effects
  vector[N_8] r_8_muFMM_SOcounts_1;  // actual group-level effects
  vector[N_9] r_9_muFMM_SOcounts_1;  // actual group-level effects
  vector[N_10] r_10_muFMM_SOcounts_1;  // actual group-level effects
  vector[N_11] r_11_muFFMM_SOcounts_1;  // actual group-level effects
  vector[N_12] r_12_muFFMM_SOcounts_1;  // actual group-level effects
  vector[N_13] r_13_muFFMM_SOcounts_1;  // actual group-level effects
  vector[N_14] r_14_IVSOint_1;  // actual group-level effects
  vector[N_15] r_15_IVSOint_1;  // actual group-level effects
  vector[N_16] r_16_IVSOint_1;  // actual group-level effects
  r_1_muMF_SOcounts_1 = (sd_1[1] * (z_1[1]));
  r_2_muMF_SOcounts_1 = (sd_2[1] * (z_2[1]));
  // compute actual group-level effects
  r_3 = scale_r_cor_cov(z_3, sd_3, L_3, Lcov_3);
  r_3_muMF_SOcounts_1 = r_3[, 1];
  r_3_muMFF_SOcounts_2 = r_3[, 2];
  r_3_muFMM_SOcounts_3 = r_3[, 3];
  r_3_muFFMM_SOcounts_4 = r_3[, 4];
  r_3_IVSOint_5 = r_3[, 5];
  r_4_muMF_SOcounts_1 = (sd_4[1] * (z_4[1]));
  r_5_muMFF_SOcounts_1 = (sd_5[1] * (z_5[1]));
  r_6_muMFF_SOcounts_1 = (sd_6[1] * (z_6[1]));
  r_7_muMFF_SOcounts_1 = (sd_7[1] * (z_7[1]));
  r_8_muFMM_SOcounts_1 = (sd_8[1] * (z_8[1]));
  r_9_muFMM_SOcounts_1 = (sd_9[1] * (z_9[1]));
  r_10_muFMM_SOcounts_1 = (sd_10[1] * (z_10[1]));
  r_11_muFFMM_SOcounts_1 = (sd_11[1] * (z_11[1]));
  r_12_muFFMM_SOcounts_1 = (sd_12[1] * (z_12[1]));
  r_13_muFFMM_SOcounts_1 = (sd_13[1] * (z_13[1]));
  r_14_IVSOint_1 = (sd_14[1] * (z_14[1]));
  r_15_IVSOint_1 = (sd_15[1] * (z_15[1]));
  r_16_IVSOint_1 = (sd_16[1] * (z_16[1]));
}
model {
  // likelihood including constants
  if (!prior_only) {
    // initialize linear predictor term
    vector[N_SOcounts] muMF_SOcounts = Intercept_muMF_SOcounts + Xc_muMF_SOcounts * b_muMF_SOcounts;
    // initialize linear predictor term
    vector[N_SOcounts] muMFF_SOcounts = Intercept_muMFF_SOcounts + Xc_muMFF_SOcounts * b_muMFF_SOcounts;
    // initialize linear predictor term
    vector[N_SOcounts] muFMM_SOcounts = Intercept_muFMM_SOcounts + Xc_muFMM_SOcounts * b_muFMM_SOcounts;
    // initialize linear predictor term
    vector[N_SOcounts] muFFMM_SOcounts = Intercept_muFFMM_SOcounts + Xc_muFFMM_SOcounts * b_muFFMM_SOcounts;
    // linear predictor matrix
    vector[ncat_SOcounts] mu_SOcounts[N_SOcounts];
    // initialize linear predictor term
    vector[N_IVSOint] mu_IVSOint = Intercept_IVSOint + Xc_IVSOint * b_IVSOint;
    for (n in 1:N_SOcounts) {
      // add more terms to the linear predictor
      muMF_SOcounts[n] += (bsp_muMF_SOcounts[1]) * mo(simo_muMF_SOcounts_1, Xmo_muMF_SOcounts_1[n]) + r_1_muMF_SOcounts_1[J_1_SOcounts[n]] * Z_1_muMF_SOcounts_1[n] + r_2_muMF_SOcounts_1[J_2_SOcounts[n]] * Z_2_muMF_SOcounts_1[n] + r_3_muMF_SOcounts_1[J_3_SOcounts[n]] * Z_3_muMF_SOcounts_1[n] + r_4_muMF_SOcounts_1[J_4_SOcounts[n]] * Z_4_muMF_SOcounts_1[n];
    }
    for (n in 1:N_SOcounts) {
      // add more terms to the linear predictor
      muMFF_SOcounts[n] += (bsp_muMFF_SOcounts[1]) * mo(simo_muMFF_SOcounts_1, Xmo_muMFF_SOcounts_1[n]) + r_3_muMFF_SOcounts_2[J_3_SOcounts[n]] * Z_3_muMFF_SOcounts_2[n] + r_5_muMFF_SOcounts_1[J_5_SOcounts[n]] * Z_5_muMFF_SOcounts_1[n] + r_6_muMFF_SOcounts_1[J_6_SOcounts[n]] * Z_6_muMFF_SOcounts_1[n] + r_7_muMFF_SOcounts_1[J_7_SOcounts[n]] * Z_7_muMFF_SOcounts_1[n];
    }
    for (n in 1:N_SOcounts) {
      // add more terms to the linear predictor
      muFMM_SOcounts[n] += (bsp_muFMM_SOcounts[1]) * mo(simo_muFMM_SOcounts_1, Xmo_muFMM_SOcounts_1[n]) + r_3_muFMM_SOcounts_3[J_3_SOcounts[n]] * Z_3_muFMM_SOcounts_3[n] + r_8_muFMM_SOcounts_1[J_8_SOcounts[n]] * Z_8_muFMM_SOcounts_1[n] + r_9_muFMM_SOcounts_1[J_9_SOcounts[n]] * Z_9_muFMM_SOcounts_1[n] + r_10_muFMM_SOcounts_1[J_10_SOcounts[n]] * Z_10_muFMM_SOcounts_1[n];
    }
    for (n in 1:N_SOcounts) {
      // add more terms to the linear predictor
      muFFMM_SOcounts[n] += (bsp_muFFMM_SOcounts[1]) * mo(simo_muFFMM_SOcounts_1, Xmo_muFFMM_SOcounts_1[n]) + r_3_muFFMM_SOcounts_4[J_3_SOcounts[n]] * Z_3_muFFMM_SOcounts_4[n] + r_11_muFFMM_SOcounts_1[J_11_SOcounts[n]] * Z_11_muFFMM_SOcounts_1[n] + r_12_muFFMM_SOcounts_1[J_12_SOcounts[n]] * Z_12_muFFMM_SOcounts_1[n] + r_13_muFFMM_SOcounts_1[J_13_SOcounts[n]] * Z_13_muFFMM_SOcounts_1[n];
    }
    for (n in 1:N_IVSOint) {
      // add more terms to the linear predictor
      mu_IVSOint[n] += (bsp_IVSOint[1]) * mo(simo_IVSOint_1, Xmo_IVSOint_1[n]) + r_3_IVSOint_5[J_3_IVSOint[n]] * Z_3_IVSOint_5[n] + r_14_IVSOint_1[J_14_IVSOint[n]] * Z_14_IVSOint_1[n] + r_15_IVSOint_1[J_15_IVSOint[n]] * Z_15_IVSOint_1[n] + r_16_IVSOint_1[J_16_IVSOint[n]] * Z_16_IVSOint_1[n];
    }
    for (n in 1:N_SOcounts) {
      mu_SOcounts[n] = transpose([0, muMF_SOcounts[n], muMFF_SOcounts[n], muFMM_SOcounts[n], muFFMM_SOcounts[n]]);
    }
    for (n in 1:N_SOcounts) {
      target += multinomial_logit2_lpmf(Y_SOcounts[n] | mu_SOcounts[n]);
    }
    target += binomial_logit_lpmf(Y_IVSOint | trials_IVSOint, mu_IVSOint);
  }
  // priors including constants
  target += normal_lpdf(Intercept_muMF_SOcounts | 0,1);
  target += dirichlet_lpdf(simo_muMF_SOcounts_1 | con_simo_muMF_SOcounts_1);
  target += normal_lpdf(Intercept_muMFF_SOcounts | 0,1);
  target += dirichlet_lpdf(simo_muMFF_SOcounts_1 | con_simo_muMFF_SOcounts_1);
  target += normal_lpdf(Intercept_muFMM_SOcounts | 0,1);
  target += dirichlet_lpdf(simo_muFMM_SOcounts_1 | con_simo_muFMM_SOcounts_1);
  target += normal_lpdf(Intercept_muFFMM_SOcounts | 0,1);
  target += dirichlet_lpdf(simo_muFFMM_SOcounts_1 | con_simo_muFFMM_SOcounts_1);
  target += normal_lpdf(Intercept_IVSOint | 0,1);
  target += dirichlet_lpdf(simo_IVSOint_1 | con_simo_IVSOint_1);
  target += exponential_lpdf(sd_1 | 2);
  target += std_normal_lpdf(z_1[1]);
  target += exponential_lpdf(sd_2 | 2);
  target += std_normal_lpdf(z_2[1]);
  target += exponential_lpdf(sd_3 | 2);
  target += std_normal_lpdf(to_vector(z_3));
  target += lkj_corr_cholesky_lpdf(L_3 | 1);
  target += exponential_lpdf(sd_4 | 2);
  target += std_normal_lpdf(z_4[1]);
  target += exponential_lpdf(sd_5 | 2);
  target += std_normal_lpdf(z_5[1]);
  target += exponential_lpdf(sd_6 | 2);
  target += std_normal_lpdf(z_6[1]);
  target += exponential_lpdf(sd_7 | 2);
  target += std_normal_lpdf(z_7[1]);
  target += exponential_lpdf(sd_8 | 2);
  target += std_normal_lpdf(z_8[1]);
  target += exponential_lpdf(sd_9 | 2);
  target += std_normal_lpdf(z_9[1]);
  target += exponential_lpdf(sd_10 | 2);
  target += std_normal_lpdf(z_10[1]);
  target += exponential_lpdf(sd_11 | 2);
  target += std_normal_lpdf(z_11[1]);
  target += exponential_lpdf(sd_12 | 2);
  target += std_normal_lpdf(z_12[1]);
  target += exponential_lpdf(sd_13 | 2);
  target += std_normal_lpdf(z_13[1]);
  target += exponential_lpdf(sd_14 | 2);
  target += std_normal_lpdf(z_14[1]);
  target += exponential_lpdf(sd_15 | 2);
  target += std_normal_lpdf(z_15[1]);
  target += exponential_lpdf(sd_16 | 2);
  target += std_normal_lpdf(z_16[1]);
}
generated quantities {
  // actual population-level intercept
  real b_muMF_SOcounts_Intercept = Intercept_muMF_SOcounts - dot_product(means_X_muMF_SOcounts, b_muMF_SOcounts);
  // actual population-level intercept
  real b_muMFF_SOcounts_Intercept = Intercept_muMFF_SOcounts - dot_product(means_X_muMFF_SOcounts, b_muMFF_SOcounts);
  // actual population-level intercept
  real b_muFMM_SOcounts_Intercept = Intercept_muFMM_SOcounts - dot_product(means_X_muFMM_SOcounts, b_muFMM_SOcounts);
  // actual population-level intercept
  real b_muFFMM_SOcounts_Intercept = Intercept_muFFMM_SOcounts - dot_product(means_X_muFFMM_SOcounts, b_muFFMM_SOcounts);
  // actual population-level intercept
  real b_IVSOint_Intercept = Intercept_IVSOint - dot_product(means_X_IVSOint, b_IVSOint);
  // compute group-level correlations
  corr_matrix[M_3] Cor_3 = multiply_lower_tri_self_transpose(L_3);
  vector<lower=-1,upper=1>[NC_3] cor_3;
  // extract upper diagonal of correlation matrix
  for (k in 1:M_3) {
    for (j in 1:(k - 1)) {
      cor_3[choose(k - 1, 2) + j] = Cor_3[j, k];
    }
  }
}

