rm(list = ls())
library(tidyverse)
library(patchwork)
set.seed(6)
source("sim_models.R")
if (!dir.exists('figures/')) dir.create('figures/')
models <- c("model1", "model2a", "model2b", "model2c", "model3")
model_names <- c("Model 1", "Model 2(a)", "Model 2(b)", "Model 2(c)", "Model 3")
nrep <- 20; len <- 200
dat <- lapply(1:length(models), function(i) {
  get(models[i])(nrep, len)
})
colors <- c("black", "red", "blue")

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
    scale_color_manual(values = colors) +
    labs(x = "Frequency", y = "Power", title = model_names[i]) +
    theme_bw() +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5))
})
wrap_plots(spec_plots, nrow = 1)

# next three rows: time series plot from each of the subpopulations of each model
# each row corresponds to one subpopulation label
x_plots <- unlist(lapply(1:3, function(i) {
  lapply(1:length(models), function(j) {
    ind <- sample(which(dat[[j]]$labels == i), size = 1)
    x <- data.frame(time = 1:nrow(dat[[j]]$x), x = dat[[j]]$x[, ind])
    ggplot(x, aes(x = time, y = x)) +
      geom_line(color = colors[i]) +
      labs(x = "Time", y = "", title = paste0("Subpopulation ", i)) +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5))
  })
}), recursive = FALSE)

pdf('figures/sim_models.pdf', width = 20, height = 16)
wrap_plots(c(spec_plots, x_plots), nrow = 4)
dev.off()



