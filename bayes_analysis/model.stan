data{
    int<lower=1> N; // Sample size. 
    int<lower=1> P; // Number of predictors. 
    matrix[N, P] X; // Design matrix.  
    vector<lower=0>[N] Y; // Outcome. 
    int<lower=1,upper=3> lik; // 1=log-Normal, 2=Gamma, 3=Weibull
    real<lower=0> sigma_beta; 
    real<lower=0> c; 
    real<lower=0> d; 
    real<lower=0> sigma_b; 
}

parameters {
    vector[P] beta;
    real<lower=0> a;
    real b; 
}

transformed parameters {
    vector[N] mu_Y = X * beta; // Link function. 
    vector<lower=0>[N] sigma_2_Y = a * pow(mu_Y, b); // Variance power law. 
    vector<lower=0>[N] sigma_Y = sqrt(sigma_2_Y); 
}

model {
    // Priors 
    beta ~ normal(0, sigma_beta); 
    a ~ gamma(c, d); 
    b ~ normal(0, sigma_b); 

    // Likelihood
    target += lognormal_lpdf(Y | mu_Y, sigma_Y); 
}

generated quantities {
    array[N] real<lower=0> Y_post; 
    Y_post = lognormal_rng(mu_Y, sigma_Y); 
}
