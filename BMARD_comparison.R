# Comparison of FBAM with a BMARD-based approach for frequency band estimation
# Ensure the working directory is the fbam-sims directory
# command line arguments
# [1] nsim              integer   number of simulation repetitions (e.g., 100)
# [2] ncores            integer   number of cores to use for parallelization
# [3] model_name        string    e.g., "model1", "model2a"
# [4] nrep              integer   number of replicate time series (per subpopulation)
# [5] len               integer   length of the time series epochs
# [6] nclust            integer   number of clusters to use in clustering on
#                                 cosine basis expansion
# [7] results_dir       string    where to save the simulation results (.rds file)
args <- commandArgs(trailingOnly = T)
nsim <- as.integer(args[1])
ncores <- as.integer(args[2])
model_name <- args[3]
nrep <- as.integer(args[4])
len <- as.integer(args[5])
nclust <- as.integer(args[6])
results_dir <- args[7]

# BMARD parameters
nsamples <- 500     # number of MCMC samples
nchains <- 3        # number of MCMC chains
sample_freq <- 1    # sample frequency

# seed for reproducibility and load required libraries
set.seed(541)
library(fda)
library(mclust)
library(fbam)
source('sim_models.R')

# confirm settings in the log
cat("BMARD STUDY PARAMETERS:\n",
    "nsim: ", nsim, "\n",
    "ncores: ", ncores, "\n",
    "model_name: ", model_name, "\n",
    "nrep: ", nrep, "\n",
    "len: ", len, "\n",
    "nclust: ", nclust, "\n",
    "nsamples: ", nsamples, "\n",
    "nchains: ", nchains, "\n",
    "results_dir: ", results_dir, "\n")

# source required functions
Rcpp::sourceCpp("BMARD/CPPcode/BMARD_V112020.cpp")
source("BMARD/Rcode_auxiliary/ExtractionBMARDmaincomponentsmodes.R")
source("BMARD_comparison_functions.R")

# path to the results data file
DATA_FILE_NAME <- file.path(results_dir, paste0(model_name, "_nrep=", nrep, "_len=", len, ".rda"))
cat("Results will be saved at", DATA_FILE_NAME, "\n")

# run simulation nsim times and save results to disk at each iteration
SIM_START_TIME <- Sys.time()
output_data <- list()
nsuccess <- 0; nfail <- 0
while (nsuccess < nsim & nfail < nsim) {
  tryCatch({
    cat("\nReplicate starting at", format(Sys.time(), usetz = TRUE), "\n")
    cat("Number of successful runs:", nsuccess, "\n")
    cat("Number of failed runs:", nfail, "\n")
    cat("Generating data...\n")
    dat <- get(model_name)(nrep, len)
    nfreq <- length(dat$mtfreq)

    cat("Running cluster -> BMARD -> band estimation procedure...\n")
    BMARD_START_TIME <- Sys.time()
    # k-means clustering of cosine basis expansion coefficients
    kmlabels <- cluster_cosine_basis(dat$mtfreq, dat$mtspec, nclust, nbasis = 15)

    # BMARD for peak detection
    ## compute periodogram for each replicate
    ## result is each periodogram replicated nchains times for next step
    periodograms <- unlist(apply(dat$x, MARGIN = 2, FUN = function(x) {
      len <- length(x); x <- (x - mean(x)) / sd(x)   # standardize time series
      per <- Mod(fft(x))^2; per <- 2 * per[2:floor(len / 2)] / len # periodogram
      out <- list(); for (n in 1:nchains) out[[n]] <- per; return(out)
    }), recursive = FALSE)

    ## MCMC sampling
    BDP <- parallel::mclapply(periodograms, function(x) {
      multirun(dat$mtfreq, x, nsamples, sample_freq)
    }, mc.cores = ncores)

    ## process BMARD samples for each replicate
    cat("Processing samples...\n")
    out <- data.frame()
    for (i in 1:ncol(dat$x)) {
      rep_start <- (nchains * (i - 1) + 1) # first chain of current replicate TS
      rep_end <- rep_start + (nchains - 1) # last chain of current replicate TS
      BDPout <- BDP[rep_start:rep_end] # all MCMC chains for the current replicate TS
      try({
        listmodes <- modes_extraction_gaussianmix(dat = BDPout, threshold_weight = .01,
                                                  chains = nchains, Nsamp = nsamples,
                                                  freq_rate = nrow(dat$x), quant = .9,
                                                  percsamp = .6, qlow = .025, qsup = .975)
      })
      peaks <- c(0, sort(listmodes$globalcenter[1,]), 0.5)
      endpoints <- peaks[1:(length(peaks) - 1)] + diff(peaks) / 2
      endpoints_index <- rep(0, length(endpoints))
      for(b in 1:length(endpoints)) {
        endpoints_index[b] <- which.min(abs(endpoints[b] - dat$mtfreq))
      }
      if (endpoints_index[1] == 1) endpoints_index[1] <- 2
      if (endpoints_index[length(endpoints)] == nfreq + 1) {
        endpoints_index[length(endpoints)] <- nfreq
      }

      nfreq <- nrow(dat$mtspec)
      clust_spec <- dat$mtspec[, kmlabels == kmlabels[i], drop = FALSE]
      collapsed <- fbam:::avg_summary(clust_spec, endpoints_index)
      expanded <- rep(collapsed, diff(c(1, endpoints_index, nfreq + 1)))
      clust_loss <- sum((clust_spec - expanded)^2) / (2 * (nfreq + 1))

      new.row <- data.frame(
        peaks = I(list(peaks)),
        endpoints = I(list(endpoints)),
        endpoints_index = I(list(endpoints_index)),
        clust_loss = clust_loss,
        label = kmlabels[i]
      )
      out <- rbind(out, new.row)
    }

    ## aggregate by cluster and compute average endpoints and bandwidths
    endpoints <- list()
    clust_loss <- rep(0, 3)
    out$n_endpoints <- sapply(out$endpoints, length)
    for (j in 1:3) {
      outj <- out[out$label == j, ]
      component_count <- table(outj$n_endpoints)
      mode <- which.max(component_count); mode <- as.integer(names(mode))
      outj <- outj[outj$n_endpoints == mode, ]
      min_loss <- outj[which.min(outj$clust_loss),]
      endpoints[j] <- min_loss$endpoints_index
      clust_loss[j] <- min_loss$clust_loss
    }
    BMARD_TIME <- as.numeric(Sys.time() - BMARD_START_TIME, units = "secs")
    cat("Finished in", BMARD_TIME, "seconds.\n")

    FBAM_START_TIME <- Sys.time()
    cat("Running FBAM on generated data...\n")
    fbam_out <- fbam(dat$x, nbands = 2:6, nsubpop = nclust, parallel = ncores)
    FBAM_TIME <- as.numeric(Sys.time() - FBAM_START_TIME, units = "secs")
    cat("FBAM finished in", FBAM_TIME, "seconds.\n")

    cat("Saving results to disk...\n")
    output_data[[nsuccess + 1]] <- list(
      data = dat,
      kmlabels = kmlabels,
      bmard_endpoints = endpoints,
      bmard_obj = sum(clust_loss),
      fbam_out = fbam_out,
      bmard_time = BMARD_TIME,
      fbam_time = FBAM_TIME
    )
    save(output_data, file = DATA_FILE_NAME)
    nsuccess <- nsuccess + 1
  }, error = function(e) {
    assign("nfail", nfail + 1, env=globalenv())
    cat('\n', str(e), '\n')
  })
}
cat("\n\nRun completed at", format(Sys.time(), usetz = TRUE), "\n")
cat("Total runtime:", format(Sys.time() - SIM_START_TIME, usetz = TRUE), "\n\n")

