#!/usr/bin/env Rscript
.libPaths("~/R/libs")
suppressMessages({
  library(pegas)
  library(adegenet)
  library(hierfstat)
})
source("/home/apena/Splatche/ARP2pegas.R")

args <- commandArgs(trailingOnly=TRUE)
scen  <- args[1]           # "1" o "2"
start <- as.integer(args[2])
end   <- as.integer(args[3])

outdir <- paste0("/home/apena/Splatche/datasetC/outscenario", scen, "_C")
param  <- read.table("/home/apena/Splatche/datasetC/param.txt",
                     col.names=c("simid","m","K1","K2"))

col_K <- if (scen=="1") "K1" else "K2"

results <- data.frame(simid=integer(), K=numeric(), m=numeric(),
                      H=numeric(), FST=numeric())

for (i in start:end) {
  f <- file.path(outdir, sprintf("out%s_%d.arp", scen, i))
  if (!file.exists(f)) next
  tryCatch({
    d  <- ARP2pegas(f, otherpackage="hierfstat")
    bs <- basic.stats(d)$overall
    results <- rbind(results, data.frame(
      simid = i,
      K     = param[i, col_K],
      m     = param[i, "m"],
      H     = as.numeric(bs["Hs"]),
      FST   = as.numeric(bs["Fst"])
    ))
  }, error=function(e){
    cat("ERROR", basename(f), ":", conditionMessage(e), "\n", file=stderr())
  })
}

# escribir chunk con nombre único
outfile <- sprintf("/home/apena/Splatche/datasetC/chunks_scen%s/chunk_%06d_%06d.txt",
                   scen, start, end)
write.table(results, outfile, row.names=FALSE, quote=FALSE, sep="\t")
