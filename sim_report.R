library(ggplot2)
library(patchwork)
library(dplyr)
library(tidyr)
library(fbam)
source('sim_models.R')

args <- commandArgs(trailingOnly = TRUE)
results_dir <- args[1]
rda_files <- list.files(results_dir, pattern = "\\.rds$", full.names = TRUE)
report <- data.frame()

get_mode <- function(v) {
  uniq_vals <- unique(v)
  uniq_vals[which.max(tabulate(match(v, uniq_vals)))]
}

print("Generating .csv simulation report...")
for (fname in rda_files) {
  print(fname)
  load(fname)
  model_name <- stringr::str_extract(basename(fname), "^[^_]+")
  nrep <- as.integer(stringr::str_extract(basename(fname), "(?<=nrep=)\\d+"))
  len <- as.integer(stringr::str_extract(basename(fname), "(?<=len=)\\d+"))
  nsim <- length(output_data)

  ari <- nsubpop <- nbands <- obj <- runtime <- rep(0, nsim)
  for (i in 1:nsim) {
    ari[i] <- mclust::adjustedRandIndex(
      output_data[[i]]$data$labels,
      output_data[[i]]$fbam_out$selected_solution$labels
    )
    nsubpop[i] <- output_data[[i]]$fbam_out$selected_solution$params$nsubpop
    nbands[i] <- output_data[[i]]$fbam_out$selected_solution$params$nbands
    obj[i] <- output_data[[i]]$fbam_out$selected_solution$objective
    runtime[i] <- output_data[[i]]$time
  }
  new_row <- data.frame(
    model = model_name,
    nsim = nsim,
    nrep = nrep,
    len = len,
    mean_ari = mean(ari, na.rm = T),
    se_ari = sd(ari, na.rm = T),
    mean_nsubpop = mean(nsubpop, na.rm = T),
    se_nsubpop = sd(nsubpop, na.rm = T),
    mean_nbands = mean(nbands, na.rm = T),
    se_nbands = sd(nbands, na.rm = T),
    mean_obj = mean(obj),
    se_obj = sd(obj),
    mean_runtime = mean(runtime),
    se_runtime = sd(runtime)
  )
  report <- rbind(report, new_row)
  rm(output_data)
}
write.csv(report, paste0(results_dir, "/summary_results.csv"),
          row.names = FALSE)


# endpoint distribution plot
print("Plotting endpoint distributions by model...")
rds_files_dist <- list.files(path = results_dir,
                             pattern = "_nrep=30_len=1000\\.rds$",
                             full.names = TRUE)
models <- c("model1", "model2a", "model2b", "model2c", "model3")
model_names <- c("Model 1", "Model 2(a)", "Model 2(b)", "Model 2(c)", "Model 3")
nrep <- 20; len <- 200
dat <- lapply(1:length(models), function(i) {
  get(models[i])(nrep, len)
})

# top row: underlying spectra
spec_plots <- lapply(1:length(models), function(i) {
  spec <- data.frame(cbind(dat[[i]]$freq, dat[[i]]$spec))
  names(spec) <- c("freq", 1:ncol(dat[[i]]$x))
  spec <- spec %>%
    pivot_longer(cols = -freq, names_to = "rep", values_to = "spec") %>%
    mutate(rep = as.integer(rep)) %>% arrange(rep) %>%
    mutate(label = rep(dat[[i]]$labels, each = length(unique(freq))))
  ggplot(spec, aes(x = freq, y = spec, group = rep, color = as.factor(label),
                   linetype = as.factor(label))) +
    geom_line(linewidth = 0.5) +
    scale_linetype_manual(values = c("solid", "dashed", "dotted")) +
    scale_color_manual(values = c("black", "red", "blue")) +
    labs(x = "Frequency", y = "Power", title = model_names[i]) +
    theme_bw() +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5),
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())
})

# use the solution corresponding to nsubpop = 3 and the most commonly
# selected number of bands across the runs
param_grid <- expand.grid(nbands = 2:6, nsubpop = 2:6)
endpoint_plots <- lapply(rds_files_dist, function(fname) {
  print(fname)
  load(fname)
  model_name <- stringr::str_extract(basename(fname), "^[^_]+")
  nrep <- as.integer(stringr::str_extract(basename(fname), "(?<=nrep=)\\d+"))
  len <- as.integer(stringr::str_extract(basename(fname), "(?<=len=)\\d+"))
  nsim <- length(output_data)
  if (model_name == "model1") model_name_plot <- "Model 1"
  if (model_name == "model2a") model_name_plot <- "Model 2(a)"
  if (model_name == "model2b") model_name_plot <- "Model 2(b)"
  if (model_name == "model2c") model_name_plot <- "Model 2(c)"
  if (model_name == "model3") model_name_plot <- "Model 3"

  # determine most commonly selected value of nbands
  nbands <- unlist(lapply(output_data, function(x) {
    x$fbam_out$selected_solution$params$nbands
  }))
  nbands <- as.integer(names(which.max(table(nbands))))

  # endpoints by label
  all_endpoints <- lapply(1:nsim, function(i) {
    # all models contain 3 true subpopulations
    ind <- which(param_grid$nbands == nbands & param_grid$nsubpop == 3)

    # solve the linear sum assignment problem to re-label the estimated
    # labels so they are consistent across all replications
    X <- output_data[[i]]$data$labels
    Y <- output_data[[i]]$fbam_out$all_solutions[[ind]]$labels
    assignment <- clue::solve_LSAP(table(X, Y), maximum = TRUE)
    mapping <- cbind(1:length(assignment), as.integer(assignment))
    endpoints <- as.vector(t(output_data[[i]]$fbam_out$all_solutions[[ind]]$endpoints))
    label <- rep(1:length(assignment), each = nbands - 1)
    for (i in 1:length(label)) label[i] <- mapping[which(mapping[, 2] == label[i]), 1]
    return(data.frame(endpoints = endpoints, label = label))
  })
  all_endpoints <- do.call(rbind, all_endpoints)
  ggplot(data = all_endpoints,
         mapping = aes(x = endpoints, color = as.factor(label),
                       fill = as.factor(label))) +
    geom_density(alpha = 0.2, bw = 0.005) +
    # geom_rug() +
    scale_linetype_manual(values = c("solid", "dashed", "twodash")) +
    scale_colour_manual(values = c("black", "red", "blue")) +
    scale_fill_manual(values = c("black", "red", "blue")) +
    xlim(0, 0.5) +
    theme_bw() +
    labs(x = "Frequency", y = "Density") +
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = "none")
})
pdf(paste0(results_dir, 'endpoint_distributions.pdf'), width = 20, height = 4)
wrap_plots(c(spec_plots, endpoint_plots), nrow = 2, heights = c(1, 0.5))
dev.off()
