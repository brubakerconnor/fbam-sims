library(ggplot2)
library(patchwork)

args <- commandArgs(trailingOnly = TRUE)
results_dir <- args[1]
rds_files <- list.files(path = results_dir,
                        pattern = "_nrep=30_len=1000\\.rds$",
                        full.names = TRUE)

# use selected solution; no coloring of histograms according to label
histograms <- lapply(rds_files, function(fname) {
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

  all_endpoints <- data.frame(endpoints = unlist(lapply(1:nsim, function(i) {
    as.vector(output_data[[i]]$fbam_out$selected_solution$endpoints)
  })))
  ggplot(data = all_endpoints, aes(x = endpoints)) +
    geom_histogram(breaks = seq(0, 0.5, 0.005), fill = "grey80", color = "black") +
    # geom_density(fill = "grey80") +
    # geom_rug() +
    xlim(0, 0.5) +
    theme_bw() +
    labs(x = "Estimated Boundaries", y = "Density",
         title = model_name_plot) +
    theme(plot.title = element_text(hjust = 0.5))
})
pdf(paste0(results_dir, 'endpoint_histograms.pdf'), width = 10, height = 2)
wrap_plots(histograms, nrow = 1)
dev.off()

# use the solution corresponding to nsubpop = 3 and the most commonly
# selected number of bands across the runs
param_grid <- expand.grid(nbands = 2:6, nsubpop = 2:6)
histograms1 <- lapply(rds_files, function(fname) {
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

  nbands <- unlist(lapply(output_data, function(x) {
    x$fbam_out$selected_solution$params$nbands
  }))
  nbands <- as.integer(names(which.max(table(nbands))))

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
  ggplot(data = all_endpoints, aes(x = endpoints, color = as.factor(label),
                                   linetype = as.factor(label))) +
    geom_histogram(fill="white", alpha=0.3, position="identity",
                   breaks = seq(0, 0.5, 0.005)) +
    # geom_density(fill = "grey80") +
    # geom_rug() +
    scale_linetype_manual(values = c("solid", "dashed", "twodash")) +
    xlim(0, 0.5) +
    theme_bw() +
    labs(x = "Estimated Boundaries", y = "Density",
         title = model_name_plot) +
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = "none")
})
pdf(paste0(results_dir, 'endpoint_histograms1.pdf'), width = 10, height = 2)
wrap_plots(histograms1, nrow = 1)
dev.off()
