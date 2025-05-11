# Libraries
library(tidyverse)
library(kableExtra)
library(cmdstanr)

options(mc.cores = parallel::detectCores())

# Meta Data
response_vars <- c("count_sd", "p_value_sd", "effect_size_sd")

# Hyperparameters
n_chains <- 2
n_iter <- 500
burn_in <- 200

# Functions

define_model_grid <- function() {
    datasets <- list.files("data", pattern = "\\.csv", full.names = TRUE)

    models <- tibble(
        model_name = c("Normal", "logNormal", "Gamma"),
        stan_file = c(
            "stan/Normal.stan", "stan/logNormal.stan", "stan/Gamma.stan"
        )
    )

    priors <- tibble(
        prior_name = c("ridge", "lasso", "weak ridge", "weak lasso"),
        prior_beta = c(1, 2, 1, 2),
        prior_args = list(
            list(sigma_beta = 10, c = 0.1, d = 0.1),
            list(sigma_beta = 10, c = 0.1, d = 0.1),
            list(sigma_beta = 100, c = 0.01, d = 0.01),
            list(sigma_beta = 100, c = 0.01, d = 0.01)
        )
    )

    parameter_grid <- expand.grid(
        dataset = datasets,
        model_name = models$model_name,
        prior_name = priors$prior_name,
        stringsAsFactors = FALSE
    ) |>
        as_tibble() |>
        left_join(models, by = "model_name") |>
        left_join(priors, by = "prior_name")

    return(parameter_grid)
}

prepare_stan_data <- function(dataset, prior_beta, prior_args) {
    df <- read_csv(dataset)

    Y <- df |>
        select(any_of(response_vars)) |>
        pull()
    X <- df |>
        select(-any_of(response_vars)) |>
        as.matrix()

    stan_data_list <- list(
        N = nrow(df),
        P = ncol(X),
        Y = Y,
        X = X,
        prior_beta = prior_beta,
        sigma_beta = prior_args$sigma_beta,
        c = prior_args$c,
        d = prior_args$d
    )

    return(stan_data_list)
}

fit_stan_model <- function(model_grid_row, n_chains, n_iter, burn_in) {
    model <- cmdstan_model(model_grid_row$stan_file)

    stan_data_list <- prepare_stan_data(
        dataset = model_grid_row$dataset,
        prior_beta = model_grid_row$prior_beta,
        prior_args = model_grid_row$prior_args[[1]]
    )

    fit <- model$sample(
        data = stan_data_list,
        chains = n_chains, iter_sampling = n_iter,
        iter_warmup = burn_in,
        seed = 05042025
    )

    invisible(fit$summary())
    # invisible(fit$draws())

    attr(fit, "run_info") <- model_grid_row # Adds in meta data.

    return(fit)
}

beta_table <- function(stan_fit) {
    run_info <- attr(stan_fit, "run_info")

    # Extract dataset name.
    dataset_name <- case_when(
        str_detect(run_info$dataset, "count") ~ "count",
        str_detect(run_info$dataset, "p_value") ~ "p-value",
        str_detect(run_info$dataset, "effect_size") ~ "effect size",
        TRUE ~ "unknown"
    )

    caption <- paste(
        "Summary of", dataset_name, run_info$model_name, "model with",
        run_info$prior_name, "prior"
    )

    # Extract predictor names.
    pred_names <- read_csv(run_info$dataset) |>
        select(-any_of(response_vars)) |>
        colnames()

    beta_summary <- stan_fit$summary(variables = "beta") |>
        mutate(predictor = pred_names) |>
        mutate(avg_ess = (ess_bulk + ess_tail) / 2) |>
        mutate(avg_ess_percent = avg_ess / (n_iter * n_chains)) |>
        filter(sign(q5) == sign(q95)) |>
        filter(rhat <= 1.1) |>
        select(predictor, mean, sd, q5, q95, rhat, avg_ess_percent)

    table_1 <- beta_summary |>
        kbl(
            format = "latex", booktabs = T,
            longtable = T, linesep = "", align = "c",
            caption = caption, digits = 4
        )

    return(table_1)
}

post_pred_check_plot <- function(stan_fit) {
    run_info <- attr(stan_fit, "run_info")

    # Extract dataset name.
    dataset_name <- case_when(
        str_detect(run_info$dataset, "count") ~ "count",
        str_detect(run_info$dataset, "p_value") ~ "p-value",
        str_detect(run_info$dataset, "effect_size") ~ "effect size",
        TRUE ~ "unknown"
    )

    caption <- paste(
        "PPC for",
        dataset_name, run_info$model_name, "model with",
        run_info$prior_name, "prior"
    )

    post_draws <- stan_fit$draws(variables = "Y_post", format = "df") |>
        slice_sample(n = 2) |>
        pivot_longer(cols = everything(), values_to = "value") |>
        select(value) |>
        mutate(source = "posterior predictive") |>
        unname()

    obs <- read_csv(run_info$dataset) |>
        select(any_of(response_vars)) |>
        mutate(source = "observed") |>
        unname()

    ppc_plot_df <- rbind(post_draws, obs)
    colnames(ppc_plot_df) <- c("sample", "source")

    upper_bound <- ppc_plot_df |>
        filter(source == "observed") |>
        pull(sample) |>
        max()

    lower_bound <- ppc_plot_df |>
        filter(source == "observed") |>
        pull(sample) |>
        min()

    ppc_plot <- ppc_plot_df |>
        mutate(source = as.factor(source)) |>
        ggplot(aes(x = sample, col = source)) +
        geom_density(bounds = c(lower_bound, upper_bound)) +
        # geom_boxplot(
        #     outliers = FALSE,
        #     draw_quantiles = c(0.25, 0.5, 0.75),
        #     fill = "skyblue"
        # ) +
        labs(title = caption) +
        theme_bw()

    return(ppc_plot)
}

sampler_diag_test <- function(stan_fit) {
    return(stan_fit$diagnostic_summary())
}

loo_cv <- function(stan_fit) {
    # return(stan_fit$loo())
    return(stan_fit$metadata())
}
