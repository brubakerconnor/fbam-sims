# Model 1 ###################################################################
# piecewise smooth step function
# x         inputs to at which to evaluate function
# values    vector of length L
# breaks    vector of length L-1
# delta     double, spacing around breaks for cubic spline
piecewise_smooth_step <- function(x, values, breaks, delta = 0.025) {
  coefs <- matrix(nrow = length(breaks), ncol = 4)
  for (i in 1:length(breaks)) {
    left <- (breaks[i] - delta/2); right <- (breaks[i] + delta/2)
    A <- matrix(c(
      left^3, left^2, left, 1,
      right^3, right^2, right, 1,
      3 * left^2, 2 * left, 1, 0,
      3 * right^2, 2 * right, 1, 0
    ), nrow = 4, ncol = 4, byrow = TRUE)
    coefs[i, ] <- solve(A, c(values[i], values[i + 1], 0, 0))
  }
  y <- rep(values[length(values)], length(x))
  marks <- sort(c(0, breaks - (delta/2), breaks + (delta/2), 0.5))
  for (i in 1:length(breaks)) {
    ind <- 2 * (i - 1) + 1
    y[x >= marks[ind] & x < marks[ind + 1]] <- values[i]
    x_gap <- x[x >= marks[ind + 1] & x < marks[ind + 2]]
    y[x %in% x_gap] <- matrix(c(
      x_gap^3, x_gap^2, x_gap, rep(1, length(x_gap))
    ), ncol = 4) %*% coefs[i, ]
  }
  return(y)
}

# generate data from model 1
# nrep - number of replicates per cluster
# len - length of realizations
model1 <- function(nrep, len) {
  # locations of marginal breakpoints
  breaks <- matrix(c(
    0.1, 0.25, #lower
    0.2, 0.3, #middle
    0.25, 0.4 #upper
  ), nrow = 3, byrow = TRUE)
  delta1 <- 0.025  # spacing for cubic spline

  # marginal values of the constant segments
  values <- matrix(c(
    15, 7.5, 2, #lower
    25, 12.5, 4, #middle
    35, 17.5, 6 #upper
  ), nrow = 3, byrow = TRUE)
  delta2 <- 2  # spacing around average summary measures

  # underlying spectra
  labels <- rep(1:3, each = nrep)
  freq <- seq(0, 0.5, length = floor(len / 2))
  spec <- matrix(nrow = floor(len / 2), ncol = 3 * nrep)
  for (j in 1:3) {
    for (k in 1:nrep) {
      spec[, nrep * (j - 1) + k] <- piecewise_smooth_step(
        freq,
        values[j, ] + runif(ncol(values), min = -delta2, max = delta2),
        breaks[j, ],
        delta1
      )
    }
  }
  x <- apply(rbind(spec, spec[dim(spec)[1]:1, ]), 2, spec_sim)
  mtout <- fbam::sine_mt(x, ntapers = floor(sqrt(len / 2)))
  return(list(x = x, labels = labels, freq = freq, spec = spec,
              mtfreq = mtout$mtfreq, mtspec = mtout$mtspec))
}

# Model 2 ####
# ar2 spectrum
ar2_spec <- function(freq, phi1, phi2, sd) {
  sd^2 / (1 + phi1^2 + phi2^2 -
            2 * phi1 * (1 - phi2) * cos(2 * pi * freq) -
            2 * phi2 * cos(4 * pi * freq))
}

# nrep:         number of replicates per subpopulation
# len:          length of time series realizations
# peaks:        vector of 3 marginal peak locations
# bw:           vector of 3 marginal bandwidths
model2 <- function(nrep, len, peaks, bw, sd = 2) {
  labels <- rep(1:3, each = nrep)
  freq <- seq(0, 0.5, length = 250) # true frequencies
  spec <- matrix(nrow = 250, ncol = 3 * nrep) # true spectra
  x <- matrix(nrow = len, ncol = 3 * nrep) # time series data
  for(i in 1:nrep) {
    # realizations of parameters using peak/bw parameterization from
    # Granados-Garcia et al. (2022)
    peaks_rep <- peaks + runif(3, -0.02, 0.02)
    bw_rep <- bw + runif(3, -0.005, 0.005)
    phi1 <- 2 * cos(2 * pi * peaks_rep) * exp(-bw_rep)
    phi2 <- -exp(-2 * bw_rep)
    # underlying spectra
    spec[, i] <- ar2_spec(freq, phi1[1], phi2[1], sd)
    spec[, nrep + i] <- ar2_spec(freq, phi1[2], phi2[2], sd)
    spec[, 2*nrep + i] <- ar2_spec(freq, phi1[3], phi2[3], sd)
    # time series data
    x[, i] <- arima.sim(list(ar = c(phi1[1], phi2[1])), n = len, sd = sd)
    x[, nrep + i] <- arima.sim(list(ar = c(phi1[2], phi2[2])), n = len, sd = sd)
    x[, 2 * nrep + i] <- arima.sim(list(ar = c(phi1[3], phi2[3])), n = len, sd = sd)
  }
  mtout <- fbam::sine_mt(x)
  return(list(x = x, labels = labels, freq = freq, spec = spec,
              mtfreq = mtout$mtfreq, mtspec = mtout$mtspec))
}

# model2a - most overlap between peaks
model2a <- function(nrep, len) {
  peaks <- c(0.23, 0.25, 0.27)
  bw <- rep(0.15, 3)
  return(model2(nrep, len, peaks, bw))
}

# model2b - medium overlap between peaks
model2b <- function(nrep, len) {
  peaks <- c(0.21, 0.25, 0.29)
  bw <- c(0.155, 0.15, 0.155)
  return(model2(nrep, len, peaks, bw))
}

# model2c - least overlap between peaks
model2c <- function(nrep, len) {
  peaks <- c(0.19, 0.25, 0.31)
  bw <- c(0.165, 0.15, 0.165)
  return(model2(nrep, len, peaks, bw))
}

# Model 3 ###################################################################
# AR(1) theoretical spectrum
# freq - frequencies to evaluate on
# phi1 - double; autoregressive parameter
# sd - positive double; standard deviation of the innovation process
ar1_spec <- function(x, phi, sd) {
  sd^2 / (1 + phi^2 - 2 * phi * cos(2 * pi * x))
}

# model3 - three clusters of AR(1) processes that mimic gait data
model3 <- function(nrep, len) {
  # model parameters
  phi1_lower <- c(0.3, 0.3, 0.7); phi1_upper <- c(0.4, 0.4, 0.75)
  sd1_lower <- sqrt(c(0.5, 1.75, 1.25)); sd1_upper <- sqrt(c(1, 2.25, 1.5))

  # theoretical spectra/ts data
  labels <- rep(1:3, each = nrep)
  freq <- seq(0, 0.5, length = 250)
  spec <- matrix(nrow = 250, ncol = 3 * nrep)
  x <- matrix(nrow = len, ncol = 3 * nrep)

  for(i in 1:nrep) {
    # draw realizations of each of the parameters
    phi1_ <- runif(length(phi1_lower), phi1_lower, phi1_upper)
    sd_ <- runif(length(sd1_lower), sd1_lower, sd1_upper)

    # theoretical spectra
    spec[, i] <- ar1_spec(freq, phi1_[1], sd_[1])
    spec[, nrep + i] <- ar1_spec(freq, phi1_[2], sd_[2])
    spec[, 2*nrep + i] <- ar1_spec(freq, phi1_[3], sd_[3])

    # generate time series realization of length len
    x[, i] <- arima.sim(list(ar = phi1_[1]), n = len, sd = sd_[1])
    x[, nrep + i] <- arima.sim(list(ar = phi1_[2]), n = len, sd = sd_[2])
    x[, 2 * nrep + i] <- arima.sim(list(ar = phi1_[3]), n = len, sd = sd_[3])
  }

  ## multitaper spectral estimates with n_tapers = floor(sqrt(len))
  mtout <- fbam::sine_mt(x)

  return(list(x = x, labels = labels, freq = freq, spec = spec,
              mtfreq = mtout$mtfreq, mtspec = mtout$mtspec))
}

# Simulate from underlying spectrum ####
# simulate time series realization from a single theoretical spectrum
# See: Guo and Dai (2006) Multivariate time-dependent spectral analysis using
# Cholesky decomposition
# spec - vector of spectrum values
spec_sim <- function(spec) {
  nx <- length(spec)
  sd <- sqrt(1 / (2 * nx))
  z <- vector("complex", nx); y <- vector("complex")
  i <- complex(imaginary = 1)
  for (j in 1:nx) {
    if (j / nx == 0.5 || j / nx == 1) {
      z[j] <- rnorm(1, sd = sd)
    } else if (j < floor(nx / 2) + 1) {
      z[j] <- complex(real = rnorm(1, sd = sd), imaginary = rnorm(1, sd = sd))
    } else {
      z[j] <- Conj(z[nx - j])
    }
  }
  for (t in 1:nx) {
    y[t] <- sum(sqrt(spec) * exp(2 * pi * i * (1:nx) * t / nx) * z)
  }
  return(Re(y))
}

# Simulation study ####
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
# [7] pmutate           float     mutation probability (optional)
args <- commandArgs(trailingOnly = T)
nsim <- as.integer(args[1])
ncores <- as.integer(args[2])
model_name <- args[3]
nrep <- as.integer(args[4])
len <- as.integer(args[5])
results_dir <- args[6]
pmutate <- as.double(args[7])

# make sure results directory exists
if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)

# set seed for reproducibility and load fbam library
set.seed(451)
library(fbam)

# print parameter settings to log
cat("STUDY PARAMETERS:\n",
    "nsim: ", nsim, "\n",
    "ncores: ", ncores, "\n",
    "model_name: ", model_name, "\n",
    "nrep: ", nrep, "\n",
    "len: ", len, "\n",
    "results_dir: ", results_dir, "\n",
    "pmutate: ", ifelse(is.na(pmutate), "not provided", pmutate), "\n")

# path to results data file
data_fname <- file.path(results_dir,
  paste0(model_name, "_nrep=", nrep, "_len=", len,
         ifelse(is.na(pmutate), "", paste0("_pmutate=", pmutate)), ".rda"))
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
    out <- fbam(dat$x, nbands = 2:6, nsubpop = 2:6, parallel = ncores,
                pmutate = ifelse(is.na(pmutate), 0.1, pmutate))
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

