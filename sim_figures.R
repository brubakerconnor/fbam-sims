rm(list = ls())
set.seed(6)
source("sim_study.R")

if (!dir.exists('figures/')) dir.create('figures/')

models <- c("model1", "model2a", "model2b", "model2c", "model3")
model_names <- c("Model 1", "Model 2(a)", "Model 2(b)", "Model 2(c)", "Model 3")
dat <- list()
for (i in 1:length(models)) {
  dat[[i]] <- get(models[i])(20, 500)
}

colors <- c("grey70", "blue", "red")
ltys <- c(1, 2, 6)

pdf('figures/sim_models.pdf', width = 18, height = 16)
par(mfrow = c(4, 5))
for (i in 1:5) matplot(dat[[i]]$freq, dat[[i]]$spec,
                       col = colors[dat[[i]]$labels], type = 'l',
                       lty = ltys[dat[[i]]$labels],
                       xlab = 'Frequency', ylab = "Power",
                       main = model_names[i])
for (i in 1:3) {
  for (j in 1:5) {
    ts.plot(dat[[j]]$x[, which(dat[[j]]$labels == i)[1]],
            xlab = 'Time', ylab = '',
            ylim = c(min(dat[[j]]$x), max(dat[[j]]$x)),
            main = paste0(model_names[j], ', Subpopulation ', i))
  }
}
dev.off()


