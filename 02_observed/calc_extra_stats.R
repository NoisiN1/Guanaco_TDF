#!/usr/bin/env Rscript
.libPaths("~/R/libs")
suppressMessages({
  library(pegas)
  library(adegenet)
  library(hierfstat)
})
source("/home/apena/Splatche/ARP2pegas.R")

args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 3) {
  cat("Uso: Rscript calc_extra_stats.R <scenario> <start> <end>\n")
  quit(status=1)
}
scen  <- args[1]
start <- as.integer(args[2])
end   <- as.integer(args[3])

outdir <- paste0("/home/apena/Splatche/datasetC/outscenario", scen, "_C_PANEL953")
if (!dir.exists(outdir)) {
  cat("Error: El directorio", outdir, "no existe.\n")
  quit(status=1)
}

param <- read.table("/home/apena/Splatche/datasetC/param.txt",
                    col.names = c("simid","m","K1","K2"))
col_K <- if (scen == "1") "K1" else "K2"

# Función para calcular estadísticos extra
calc_extra <- function(d) {
  bs <- basic.stats(d)
  
  fst_locus <- bs$perloc$Fst
  fst_locus <- fst_locus[!is.na(fst_locus)]
  
  var_fst <- if (length(fst_locus) > 1) var(fst_locus) else NA
  fst_quant <- if (length(fst_locus) > 0) {
    quantile(fst_locus, probs = c(0.05, 0.5, 0.95), na.rm = TRUE)
  } else {
    rep(NA, 3)
  }
  
  Hs <- bs$overall["Hs"]
  
  he_locus <- bs$perloc$He
  he_locus <- he_locus[!is.na(he_locus)]
  pi_global <- if (length(he_locus) > 0) mean(he_locus) else NA
  
  priv_total <- NA
  
  return(c(Hs = as.numeric(Hs),
           Fst_var = var_fst,
           Fst_05 = fst_quant[1],
           Fst_50 = fst_quant[2],
           Fst_95 = fst_quant[3],
           Pi = pi_global,
           Priv = priv_total))
}

# Procesar rango
results <- data.frame(simid = integer(),
                      K = numeric(),
                      m = numeric(),
                      H = numeric(),
                      FST = numeric(),
                      Hs_extra = numeric(),
                      Fst_var = numeric(),
                      Fst_05 = numeric(),
                      Fst_50 = numeric(),
                      Fst_95 = numeric(),
                      Pi = numeric(),
                      Priv = numeric())

for (i in start:end) {
  f <- file.path(outdir, sprintf("out%s_%d.arp", scen, i))
  if (!file.exists(f)) next
  
  tryCatch({
    d <- ARP2pegas(f, otherpackage = "hierfstat")
    
    # Corregir nombres de filas si contienen NA
    if (any(is.na(rownames(d$genind@tab)))) {
      rownames(d$genind@tab) <- make.names(rownames(d$genind@tab), unique = TRUE)
    }
    
    bs <- basic.stats(d)
    H <- as.numeric(bs$overall["Hs"])
    FST <- as.numeric(bs$overall["Fst"])
    
    extra <- calc_extra(d)
    
    results <- rbind(results, data.frame(
      simid = i,
      K     = param[i, col_K],
      m     = param[i, "m"],
      H     = H,
      FST   = FST,
      Hs_extra = extra["Hs"],
      Fst_var = extra["Fst_var"],
      Fst_05 = extra["Fst_05"],
      Fst_50 = extra["Fst_50"],
      Fst_95 = extra["Fst_95"],
      Pi = extra["Pi"],
      Priv = extra["Priv"]
    ))
  }, error = function(e) {
    cat("ERROR en", basename(f), ":", conditionMessage(e), "\n", file = stderr())
  })
}

outfile <- sprintf("/home/apena/Splatche/datasetC/extra_stats_scen%s_chunk_%06d_%06d.txt",
                   scen, start, end)
write.table(results, outfile, row.names = FALSE, quote = FALSE, sep = "\t")
cat("Chunk guardado en:", outfile, "\n")
