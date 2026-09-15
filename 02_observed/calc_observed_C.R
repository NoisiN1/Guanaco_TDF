#!/usr/bin/env Rscript
.libPaths("~/R/libs")
suppressMessages({
  library(adegenet)
  library(hierfstat)
})

# 1. Leer el .ped de PLINK
cat("Leyendo .ped...\n")
ped <- read.table("/home/apena/Splatche/datasetC/datasetC_bi_plink.ped",
                  header=FALSE, stringsAsFactors=FALSE)
cat("Individuos:", nrow(ped), "\n")
cat("Columnas totales:", ncol(ped), " (6 info + 2*SNPs alelos)\n")
cat("SNPs:", (ncol(ped)-6)/2, "\n")

# 2. Extraer nombres de individuos y prefijos de localidad
ind_names <- ped$V1
pops <- sub("_.*", "", ind_names)
cat("\nLocalidades detectadas:\n")
print(table(pops))

# 3. Reconstruir genotipos como pares "A/T" por locus
n_snps <- (ncol(ped) - 6) / 2
geno_matrix <- matrix(NA_character_, nrow=nrow(ped), ncol=n_snps)
for (j in 1:n_snps) {
  a1 <- ped[, 6 + 2*j - 1]
  a2 <- ped[, 6 + 2*j]
  # Missing en PLINK es "0", convertir a NA
  a1[a1=="0"] <- NA
  a2[a2=="0"] <- NA
  geno_matrix[, j] <- ifelse(is.na(a1) | is.na(a2), NA, paste(a1, a2, sep="/"))
}

# 4. Construir data.frame para adegenet
df_geno <- as.data.frame(geno_matrix, stringsAsFactors=FALSE)
colnames(df_geno) <- paste0("L", 1:n_snps)

# 5. Crear objeto genind
cat("\nConstruyendo genind...\n")
gi <- df2genind(df_geno, sep="/", ind.names=ind_names, pop=pops, ploidy=2, NA.char=NA)
cat("Genind creado. Loci:", nLoc(gi), ", Individuos:", nInd(gi), ", Pops:", nPop(gi), "\n")

# 6. Convertir a hierfstat y calcular basic.stats
hf <- genind2hierfstat(gi)
cat("\nCalculando basic.stats...\n")
bs <- basic.stats(hf)

cat("\n=== ESTADÍSTICOS OBSERVADOS datasetC (877 SNPs bialélicos) ===\n")
print(bs$overall)

# 7. Escribir observed_C.obs con formato compatible con el pipeline ABC
obs <- data.frame(H = as.numeric(bs$overall["Hs"]),
                  FST = as.numeric(bs$overall["Fst"]))
write.table(obs, "/home/apena/Splatche/datasetC/observed_C.obs",
            row.names=FALSE, quote=FALSE, sep="\t")
cat("\n>>> Guardado en observed_C.obs:\n")
cat("H (Hs) =", obs$H, "\n")
cat("FST    =", obs$FST, "\n")

# 8. También guardar el basic.stats completo para tenerlo de referencia
saveRDS(bs, "/home/apena/Splatche/datasetC/basic_stats_observed_C.rds")
cat("\nBasic stats completos guardados en basic_stats_observed_C.rds\n")
