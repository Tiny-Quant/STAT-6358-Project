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
    real<lower=0> sigma_2_Y; 
}

transformed parameters {
    vector[N] mu_Y = X * beta; 
    real<lower=0> sigma_Y = sqrt(sigma_2_Y); 
}

model { 
    // Priors 
    if (prior_beta == 1){
        beta ~ normal(0, sigma_beta); 
    } else if (prior_beta == 2) {
        beta ~ double_exponential(0, sigma_beta); 
    }
    sigma_2_Y ~ gamma(c, d); 

    // Likelihood 
    target += lognormal_lpdf(Y | mu_Y, sigma_Y); 
}

generated quantities {
    array[N] real<lower=0> Y_post; 
    Y_post = lognormal_rng(mu_Y, sigma_Y); 
}