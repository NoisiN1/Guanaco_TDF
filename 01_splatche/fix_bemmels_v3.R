#!/usr/bin/env Rscript
.libPaths("~/R/libs")
suppressMessages({
  library(terra)
  library(ggplot2)
  library(patchwork)
})

base_dir <- "/home/apena/Splatche/datasetC"
scen1_dir <- file.path(base_dir, "Splatche_scen1_deme_TDF_C/datasets")
scen2_dir <- file.path(base_dir, "Splatche_scen2_hum_TDF_C/datasets")

build_K_raster <- function(volcano_asc, veg2K_file, K_refuge_value=36) {
  r <- rast(volcano_asc)
  lines <- readLines(veg2K_file)
  lines <- lines[nchar(trimws(lines)) > 0]
  cats <- integer(); Ks <- integer()
  for (ln in lines) {
    parts <- strsplit(trimws(ln), "\\s+")[[1]]
    cat_i <- as.integer(parts[1])
    K_str <- parts[2]
    K_i <- if (K_str == "PARAM") K_refuge_value else as.integer(K_str)
    cats <- c(cats, cat_i); Ks <- c(Ks, K_i)
  }
  list(r_K = classify(r, cbind(cats, Ks)),
       r_cat = r)  # también retornamos las categorías originales
}

raster_to_df <- function(rK, rCat, phase_label) {
  df <- as.data.frame(rK, xy=TRUE, na.rm=FALSE)
  colnames(df) <- c("x","y","K")
  df$cat <- values(rCat)
  df$phase <- phase_label

  # Asignar clase discreta según K
  df$class <- cut(df$K,
    breaks=c(-Inf, 0, 100, 3000, 4000, Inf),
    labels=c("No habitat (K = 0)",
             "Refugium (K = 36)",
             "High ash (K = 2432)",
             "Low/moderate ash (K = 3892)",
             "Baseline habitat (K = 4865)"),
    right=TRUE)
  df
}

K_refuge <- 36
K_corridor <- 22150

phases_scen1 <- list(
  list(time=1, label="Pre-glacial"),
  list(time=2, label="Late Glacial"),
  list(time=3, label="Hudson H1 (~7,750 BP)"),
  list(time=4, label="Post-eruption recovery")
)
phases_scen2 <- list(
  list(time=3, label="Hudson H1 (baseline)"),
  list(time=5, label="Post-glacial corridor (early)"),
  list(time=6, label="Post-glacial corridor (late)")
)

df_s1 <- do.call(rbind, lapply(phases_scen1, function(p) {
  rl <- build_K_raster(file.path(scen1_dir, "Volcano_scen1.asc"),
                       file.path(scen1_dir, sprintf("veg2K_pop1_time_%d.txt", p$time)),
                       K_refuge)
  raster_to_df(rl$r_K, rl$r_cat, p$label)
}))
df_s1$phase <- factor(df_s1$phase, levels=sapply(phases_scen1, function(x) x$label))

scen2_volc <- list.files(scen2_dir, pattern="^Volcano.*\\.asc$", full.names=TRUE)[1]
df_s2 <- do.call(rbind, lapply(phases_scen2, function(p) {
  veg_file <- file.path(scen2_dir, sprintf("veg2K_pop1_time_%d.txt", p$time))
  if (!file.exists(veg_file)) return(NULL)
  rl <- build_K_raster(scen2_volc, veg_file, K_corridor)
  df <- raster_to_df(rl$r_K, rl$r_cat, p$label)
  # Añadir clase "Corridor" para el escenario 2 en fases 5 y 6
  df$class <- as.character(df$class)
  if (p$time %in% c(5, 6)) {
    corridor_mask <- df$cat == 10 | (df$K > 5000 & df$K < 50000)
    df$class[corridor_mask] <- sprintf("Human corridor (K = %d)", K_corridor)
  }
  df$class <- factor(df$class)
  df
}))
df_s2$phase <- factor(df_s2$phase, levels=sapply(phases_scen2, function(x) x$label))

# Ordenar niveles del factor para leyenda coherente
all_levels <- c("No habitat (K = 0)",
                "Refugium (K = 36)",
                "High ash (K = 2432)",
                "Low/moderate ash (K = 3892)",
                "Baseline habitat (K = 4865)",
                sprintf("Human corridor (K = %d)", K_corridor))

df_s1$class <- factor(df_s1$class, levels=all_levels)
df_s2$class <- factor(df_s2$class, levels=all_levels)

# Paleta discreta:
palette_map <- c(
  "No habitat (K = 0)"          = "#F5F5F5",  # gris muy claro
  "Refugium (K = 36)"           = "#D62828",  # rojo vibrante (destaca)
  "High ash (K = 2432)"         = "#F4A261",  # naranja
  "Low/moderate ash (K = 3892)" = "#E9C46A",  # amarillo
  "Baseline habitat (K = 4865)" = "#2A9D8F",  # verde teal
  "Human corridor (K = 22150)"  = "#264653"   # azul oscuro
)

# Coordenadas del refugio para anotar con flecha
refuge_x <- 36.5
refuge_y <- 5.5

make_panel <- function(df_phase, title, highlight_refuge=FALSE, highlight_corridor=FALSE) {
  p <- ggplot(df_phase, aes(x=x, y=y, fill=class)) +
    geom_raster() +
    scale_fill_manual(values=palette_map, name="Habitat class",
                      drop=FALSE, na.value="white") +
    coord_fixed(expand=FALSE) +
    labs(title=title, x=NULL, y=NULL) +
    theme_void(base_size=11) +
    theme(plot.title=element_text(hjust=0.5, size=11, face="bold"),
          legend.position="right",
          plot.margin=margin(2, 2, 2, 2))

  # Añadir contorno rojo alrededor de la celda del refugio en la fase H1
  if (highlight_refuge) {
    p <- p +
      annotate("rect", xmin=refuge_x-0.5, xmax=refuge_x+0.5,
               ymin=refuge_y-0.5, ymax=refuge_y+0.5,
               fill=NA, colour="black", linewidth=0.6) +
      annotate("segment", x=refuge_x+3, y=refuge_y-3,
               xend=refuge_x+0.7, yend=refuge_y-0.3,
               arrow=arrow(length=unit(0.15, "cm")),
               colour="black", linewidth=0.4) +
      annotate("text", x=refuge_x+3.2, y=refuge_y-3.3,
               label="Refugium", hjust=0, size=3, fontface="bold")
  }
  p
}

# Fila 1: Scen1 con highlight en H1 (fase 3)
plots_s1 <- list()
for (i in seq_along(phases_scen1)) {
  p <- phases_scen1[[i]]
  highlight <- (p$time == 3)
  plots_s1[[i]] <- make_panel(subset(df_s1, phase == p$label), p$label,
                              highlight_refuge=highlight)
}

# Fila 2: Scen2 con highlight en fase corridor (5 y 6)
plots_s2 <- list()
for (i in seq_along(phases_scen2)) {
  p <- phases_scen2[[i]]
  plots_s2[[i]] <- make_panel(subset(df_s2, phase == p$label), p$label,
                              highlight_corridor=(p$time %in% c(5,6)))
}
empty_panel <- ggplot() + theme_void()
plots_s2_padded <- c(plots_s2, list(empty_panel))

row1 <- wrap_plots(plots_s1, nrow=1, widths=c(1,1,1,1))
row2 <- wrap_plots(plots_s2_padded, nrow=1, widths=c(1,1,1,1))

combined <- (row1 / row2) +
  plot_layout(guides="collect", heights=c(1,1)) +
  plot_annotation(
                              ", Isla Grande de Tierra del Fuego and adjacent mainland")),
    tag_levels=list(c("A","","","","B","","","")),
    theme=theme(plot.title=element_text(size=15, face="bold"),
                plot.subtitle=element_text(size=11),
                legend.position="right",
                plot.tag=element_text(size=14, face="bold"))
  )

dir.create("figuras/mapas_bemmels", recursive=TRUE, showWarnings=FALSE)
ggsave("figuras/mapas_bemmels/Fig_Combined_landscapes.png", combined,
       width=16, height=8, dpi=300, bg="white")
ggsave("figuras/mapas_bemmels/Fig_Combined_landscapes.pdf", combined,
       width=16, height=8, bg="white")

cat("=== Figura regenerada con refugio destacado ===\n")
system("ls -la figuras/mapas_bemmels/Fig_Combined_landscapes.png")
