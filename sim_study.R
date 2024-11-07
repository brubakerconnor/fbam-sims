# Simulation study ####
# Ensure the working directory is the fbam-sims directory
# perform many repetitions of the fbam optimization routine on independent
# realizations of a given model with given parameters
# arguments 7-10 are optional. if not specified, default parameter values will
# be used in running the genetic algorithm.
# [1] nsim              integer   number of simulation repetitions (e.g., 100)
# [2] ncores            integer   number of cores to use for parallelization
# [3] model_name        string    e.g., "model1", "model2a"
# [4] nrep              integer   number of replicate time series (per subpopulation)
# [5] len               integer   length of the time series epochs
# [6] results_dir       string    where to save the simulation results (.rds file)
args <- commandArgs(trailingOnly = T)
nsim <- as.integer(args[1])
ncores <- as.integer(args[2])
model_name <- args[3]
nrep <- as.integer(args[4])
len <- as.integer(args[5])
results_dir <- args[6]

# make sure results directory exists
if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)

# set seed for reproducibility and load fbam library
set.seed(451)
library(fbam)
source('sim_models.R')

# print parameter settings to log
cat("STUDY PARAMETERS:\n",
    "nsim: ", nsim, "\n",
    "ncores: ", ncores, "\n",
    "model_name: ", model_name, "\n",
    "nrep: ", nrep, "\n",
    "len: ", len, "\n",
    "results_dir: ", results_dir, "\n")

# path to results data file
data_fname <- file.path(results_dir,
  paste0(model_name, "_nrep=", nrep, "_len=", len, ".rds"))
cat("Results will be saved at", data_fname, "\n")

# run simulation nsim times and save results to disk at each iteration
SIM_START_TIME <- Sys.time()
output_data <- list()
nsuccess <- 0; nfail <- 0
while (nsuccess < nsim & nfail < nsim) {
  # attempt the simulation until the required number of successes are met
  # terminate once too many failures occur
  tryCatch({
    cat("\nReplicate starting at", format(Sys.time(), usetz = TRUE), "\n")
    cat("Number of successful runs:", nsuccess, "\n")
    cat("Number of failed runs:", nfail, "\n")
    cat("Generating data...\n")
    dat <- get(model_name)(nrep, len)

    cat("Running FBAM on generated data...\n")
    FBAM_START_TIME <- Sys.time()
    out <- fbam(dat$x, nbands = 2:6, nsubpop = 2:6, parallel = ncores)
    FBAM_RUNTIME <- as.numeric(Sys.time() - FBAM_START_TIME, units = "secs")
    cat("Completed in", FBAM_RUNTIME, " seconds.\n")

    nsuccess <- nsuccess + 1
    cat("Saving results to disk...\n")
    output_data[[nsuccess]] <- list(data = dat, fbam_out = out, time = FBAM_RUNTIME)
    save(output_data, file = data_fname)
  }, error = function(e) {
    cat("Failed with error\n\n", str(e), "\n")
    nfail <- nfail + 1
  })
}
cat("\n\nRUN COMPLETED", format(Sys.time(), usetz = TRUE), "\n")
cat("TOTAL RUNTIME", format(Sys.time() - SIM_START_TIME, usetz = TRUE), "\n\n")

