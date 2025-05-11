# Created by use_targets().
# Follow the comments below to fill in this target script.
# Then follow the manual to check and run the pipeline:
#   https://books.ropensci.org/targets/walkthrough.html#inspect-the-pipeline

# Load packages required to define the pipeline:
library(targets)
library(tarchetypes) # Load other packages as needed.
cmdstanr::set_cmdstan_path(Sys.getenv("CMDSTAN"))

# Set target options:
tar_option_set(
    packages = c("tidyverse", "cmdstanr", "kableExtra"), # Packages that your targets need for their tasks.
    format = "qs", # Optionally set the default storage format. qs is fast.
    error = "null",
    #
    # Pipelines that take a long time to run may benefit from
    # optional distributed computing. To use this capability
    # in tar_make(), supply a {crew} controller
    # as discussed at https://books.ropensci.org/targets/crew.html.
    # Choose a controller that suits your needs. For example, the following
    # sets a controller that scales up to a maximum of two workers
    # which run as local R processes. Each worker launches when there is work
    # to do and exits if 60 seconds pass with no tasks to run.
    #
    controller = crew::crew_controller_local(workers = 36, seconds_idle = 60)
    #
    # Alternatively, if you want workers to run on a high-performance computing
    # cluster, select a controller from the {crew.cluster} package.
    # For the cloud, see plugin packages like {crew.aws.batch}.
    # The following example is a controller for Sun Grid Engine (SGE).
    #
    #   controller = crew.cluster::crew_controller_sge(
    #     # Number of workers that the pipeline can scale up to:
    #     workers = 10,
    #     # It is recommended to set an idle time so workers can shut themselves
    #     # down if they are not running tasks.
    #     seconds_idle = 120,
    #     # Many clusters install R as an environment module, and you can load it
    #     # with the script_lines argument. To select a specific verison of R,
    #     # you may need to include a version string, e.g. "module load R/4.3.2".
    #     # Check with your system administrator if you are unsure.
    #     script_lines = "module load R"
    #   )
    #
    # Set other options as needed.
)

# Run the R scripts in the R/ folder with your custom functions:
tar_source("./scripts/pipeline_functions.r")
# tar_source("other_functions.R") # Source other scripts as needed.

# Replace the target list below with your own:
list(
    tar_target(
        model_grid,
        define_model_grid()
    ),
    tar_target(
        stan_fit,
        fit_stan_model(
            model_grid,
            n_chains = n_chains, n_iter = n_iter, burn_in = burn_in
        ),
        pattern = map(model_grid),
        iteration = "list"
    ),
    tar_target(
        beta_tables,
        beta_table(stan_fit),
        pattern = map(stan_fit),
        iteration = "list"
    ),
    tar_target(
        ppc_plots,
        post_pred_check_plot(stan_fit),
        pattern = map(stan_fit),
        iteration = "list"
    ),
    # tar_target(
    #     diag_test,
    #     sampler_diag_test(stan_fit),
    #     pattern = map(stan_fit),
    #     iteration = "list"
    # ),
    # tar_target(
    #     loo_test,
    #     loo_cv(stan_fit),
    #     pattern = map(stan_fit),
    #     iteration = "list"
    # ),
    tar_quarto(
        report,
        path = "model_report.qmd",
        quiet = FALSE
    )
)
