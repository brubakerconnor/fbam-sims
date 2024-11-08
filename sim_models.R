# Model 1 ###################################################################
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
# peaks:        vector of 3 marginal peak locations
# bw:           vector of 3 marginal bandwidths
model2 <- function(nrep, len, peaks, bw, sd = 2.25) {
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
  # bw <- c(0.155, 0.15, 0.155)
  bw <- c(0.15, 0.16, 0.2)
  return(model2(nrep, len, peaks, bw))
}

# model2c - least overlap between peaks
model2c <- function(nrep, len) {
  peaks <- c(0.19, 0.25, 0.31)
  # bw <- c(0.165, 0.15, 0.165)
  bw <- c(0.15, 0.165, 0.205)
  return(model2(nrep, len, peaks, bw))
}

model3 <- function(nrep, len) {
  freq <- seq(0, 0.5, length = 250)

  # ar1 processes
  x <- matrix(nrow = len, ncol = 3 * nrep)
  xspec <- matrix(nrow = 250, ncol = 3 * nrep)
  for (i in 1:nrep) {
    peaks_rep <- rep(0, 3)
    bw_rep <- rep(0.5, 3) + runif(3, -0.02, 0.02)
    phi1 <- 2 * cos(2 * pi * peaks_rep) * exp(-bw_rep)
    phi2 <- -exp(-2 * bw_rep)
    xspec[, i] <- ar2_spec(freq, phi1[1], phi2[1], 1)
    xspec[, nrep + i] <- ar2_spec(freq, phi1[2], phi2[2], 1)
    xspec[, 2 * nrep + i] <- ar2_spec(freq, phi1[3], phi2[3], 1)
    x[, i] <- arima.sim(list(ar = c(phi1[1], phi2[1])), n = len, sd = 1)
    x[, nrep + i] <- arima.sim(list(ar = c(phi1[2], phi2[2])), n = len, sd = 1)
    x[, 2 * nrep + i] <- arima.sim(list(ar = c(phi1[3], phi2[3])), n = len, sd = 1)
  }

  # ar2 processes
  y <- matrix(nrow = len, ncol = 3 * nrep)
  yspec <- matrix(nrow = 250, ncol = 3 * nrep)
  for (i in 1:nrep) {
    peaks_rep <- c(0.25, 0.3, 0.35) + runif(3, -0.015, 0.015)
    bw_rep <- c(0.15, 0.15, 0.15) #+ runif(3, -0.005, 0.005)
    phi1 <- 2 * cos(2 * pi * peaks_rep) * exp(-bw_rep)
    phi2 <- -exp(-2 * bw_rep)
    yspec[, i] <- ar2_spec(freq, phi1[1], phi2[1], 1)
    yspec[, nrep + i] <- ar2_spec(freq, phi1[2], phi2[2], 1)
    yspec[, 2 * nrep + i] <- ar2_spec(freq, phi1[3], phi2[3], 1)
    y[, i] <- arima.sim(list(ar = c(phi1[1], phi2[1])), n = len, sd = 2)
    y[, nrep + i] <- arima.sim(list(ar = c(phi1[2], phi2[2])), n = len, sd = 2)
    y[, 2 * nrep + i] <- arima.sim(list(ar = c(phi1[3], phi2[3])), n = len, sd = 2)
  }

  # add ar1 and ar2 processes
  labels <- rep(1:3, each = nrep)
  z <- x + y
  zspec <- xspec + yspec
  mtout <- fbam::sine_mt(z)
  return(list(x = z, labels = labels, freq = freq, spec = zspec,
              mtfreq = mtout$mtfreq, mtspec = mtout$mtspec))
}

# Utility Functions ####
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

# ar1 spectrum
ar1_spec <- function(x, phi, sd) {
  sd^2 / (1 + phi^2 - 2 * phi * cos(2 * pi * x))
}

# ar2 spectrum
ar2_spec <- function(x, phi1, phi2, sd) {
  sd^2 / (1 + phi1^2 + phi2^2 -
            2 * phi1 * (1 - phi2) * cos(2 * pi * x) -
            2 * phi2 * cos(4 * pi * x))
}

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
