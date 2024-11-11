args <- commandArgs(trailingOnly = TRUE)
results_dir <- args[1]
rda_files <- list.files(results_dir, pattern = "\\.rds$", full.names = TRUE)
report <- data.frame()

get_mode <- function(v) {
  uniq_vals <- unique(v)
  uniq_vals[which.max(tabulate(match(v, uniq_vals)))]
}

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
    mode_nsubpop = get_mode(nsubpop),
    percent_nsubpop = sum(nsubpop == get_mode(nsubpop)),
    mode_nbands = get_mode(nbands),
    percent_nbands = sum(nbands == get_mode(nbands)),
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
