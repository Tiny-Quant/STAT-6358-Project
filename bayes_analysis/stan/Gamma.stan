data{
    int<lower=1> N; // Sample size. 
    int<lower=1> P; // Number of predictors. 
    matrix[N, P] X; // Design matrix.  
    vector<lower=0>[N] Y; // Outcome. 
    int<lower=1, upper=2> prior_beta; // 1 = normal, 2 = Laplace
    real<lower=0> sigma_beta; 
    real<lower=0> c; 
    real<lower=0> d; 
}

parameters {
    vector[P] beta;
    real<lower=0> alpha; 
}

transformed parameters {
    vector<lower=0>[N] mu_Y = exp(X * beta); // log-link. 
}

model{
    // Priors 
    if (prior_beta == 1){
        beta ~ normal(0, sigma_beta); 
    } else if (prior_beta == 2) {
        beta ~ double_exponential(0, sigma_beta); 
    }
    alpha ~ gamma(c, d); 

    // Likelihood 
    target += gamma_lpdf(Y | alpha, alpha / mu_Y); 
}

generated quantities {
    array[N] real<lower=0> Y_post; 
    Y_post = gamma_rng(alpha, alpha / mu_Y); 
}