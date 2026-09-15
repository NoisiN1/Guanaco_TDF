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

sample_coords <- data.frame(
  loc = c("eceste","ecoeste","est","porv","prim","russ","sanseb","sg"),
  row = c(6, 4, 6, 8, 9, 7, 8, 11),
  col = c(33, 29, 29, 26, 28, 28, 30, 27)
)

build_K_raster <- function(volcano_asc, veg2K_file, K_refuge_value=36) {
  r <- rast(volcano_asc)
  lines <- readLines(veg2K_file)
  lines <- lines[nchar(trimws(lines)) > 0]
  cats  <- integer()
  Ks    <- integer()
  for (ln in lines) {
    parts <- strsplit(trimws(ln), "\\s+")[[1]]
    cat_i <- as.integer(parts[1])
    K_str <- parts[2]
    if (K_str == "PARAM") {
      K_i <- K_refuge_value
    } else {
      K_i <- as.integer(K_str)
    }
    cats <- c(cats, cat_i)
    Ks   <- c(Ks, K_i)
  }
  r_K <- classify(r, cbind(cats, Ks))
  return(r_K)
}

raster_to_df <- function(r, phase_label) {
  df <- as.data.frame(r, xy=TRUE, na.rm=FALSE)
  colnames(df) <- c("x","y","K")
  df$phase <- phase_label
  df$suitability <- df$K / 4865
  return(df)
}

K_refuge <- 36
K_corridor <- 22150

phases_scen1 <- list(
  list(time=1, label="Pre-glacial"),
  list(time=2, label="Late Glacial"),
  list(time=3, label="Hudson H1 (~7,750 BP)"),
  list(time=4, label="Post-eruption recovery")
)

df_all_s1 <- data.frame()
for (p in phases_scen1) {
  veg_file <- file.path(scen1_dir, sprintf("veg2K_pop1_time_%d.txt", p$time))
  r_K <- build_K_raster(file.path(scen1_dir, "Volcano_scen1.asc"),
                        veg_file, K_refuge_value=K_refuge)
  df <- raster_to_df(r_K, p$label)
  df_all_s1 <- rbind(df_all_s1, df)
}
df_all_s1$phase <- factor(df_all_s1$phase,
                          levels=sapply(phases_scen1, function(x) x$label))

scen2_volc <- list.files(scen2_dir, pattern="^Volcano.*\\.asc$", full.names=TRUE)[1]
cat("Scen2 volcano file:", scen2_volc, "\n")

phases_scen2 <- list(
  list(time=3, label="Hudson H1 (baseline)"),
  list(time=5, label="Post-glacial corridor (early)"),
  list(time=6, label="Post-glacial corridor (late)")
)

df_all_s2 <- data.frame()
for (p in phases_scen2) {
  veg_file <- file.path(scen2_dir, sprintf("veg2K_pop1_time_%d.txt", p$time))
  if (!file.exists(veg_file)) next
  r_K <- build_K_raster(scen2_volc, veg_file, K_refuge_value=K_corridor)
  df <- raster_to_df(r_K, p$label)
  df_all_s2 <- rbind(df_all_s2, df)
}
df_all_s2$phase <- factor(df_all_s2$phase,
                          levels=sapply(phases_scen2, function(x) x$label))

nrows_grid <- 38
sample_coords$x_plot <- sample_coords$col
sample_coords$y_plot <- sample_coords$row

bemmels_palette <- c("#F5DEB3", "#EDD8AC", "#F0D08C", "#E8B96C", "#D9B04E",
                     "#B8C64B", "#8CB94C", "#5AA83E", "#2E8B37", "#1A5A2E")

make_panel <- function(df_phase, title, show_samples=TRUE) {
  p <- ggplot(df_phase, aes(x=x, y=y, fill=suitability)) +
    geom_raster() +
    scale_fill_gradientn(colours=bemmels_palette,
                         limits=c(0, 1),
                         na.value="white",
                         name="Habitat\nsuitability") +
    coord_fixed(expand=FALSE) +
    labs(title=title, x=NULL, y=NULL) +
    theme_void(base_size=10) +
    theme(plot.title=element_text(hjust=0.5, size=10, face="bold"),
          legend.position="right")
  if (show_samples) {
    p <- p + geom_point(data=sample_coords,
                        aes(x=x_plot, y=y_plot),
                        inherit.aes=FALSE,
                        shape=21, fill="red", colour="black",
                        size=2.5, stroke=0.5)
  }
  return(p)
}

dir.create("figuras/mapas_bemmels", recursive=TRUE, showWarnings=FALSE)

plots_s1 <- lapply(phases_scen1, function(p) {
  df_phase <- subset(df_all_s1, phase == p$label)
  make_panel(df_phase, p$label)
})

fig_s1 <- wrap_plots(plots_s1, ncol=length(plots_s1), guides="collect") +
  plot_annotation(
    title="Scenario 1 - Demographic refugium",
    subtitle=sprintf("Habitat suitability across temporal phases (K_refuge = %d)", K_refuge),
    theme=theme(plot.title=element_text(size=14, face="bold"),
                plot.subtitle=element_text(size=11))
  ) & theme(legend.position="right")

ggsave("figuras/mapas_bemmels/Fig_Scen1_landscapes.png", fig_s1,
       width=14, height=4, dpi=300, bg="white")

plots_s2 <- lapply(phases_scen2, function(p) {
  df_phase <- subset(df_all_s2, phase == p$label)
  make_panel(df_phase, p$label)
})

fig_s2 <- wrap_plots(plots_s2, ncol=length(plots_s2), guides="collect") +
  plot_annotation(
    title="Scenario 2 - Human-mediated connectivity",
    subtitle=sprintf("Habitat suitability with corridor active (K_corridor = %d)", K_corridor),
    theme=theme(plot.title=element_text(size=14, face="bold"),
                plot.subtitle=element_text(size=11))
  ) & theme(legend.position="right")

ggsave("figuras/mapas_bemmels/Fig_Scen2_landscapes.png", fig_s2,
       width=12, height=4, dpi=300, bg="white")

combined <- (fig_s1 / fig_s2) +
  plot_annotation(
    title="Habitat suitability landscapes under competing demographic scenarios",
    subtitle="Lama guanicoe, Isla Grande de Tierra del Fuego - 8 sampling localities in red",
    theme=theme(plot.title=element_text(size=16, face="bold"),
                plot.subtitle=element_text(size=12))
  )

ggsave("figuras/mapas_bemmels/Fig_Combined_landscapes.png", combined,
       width=14, height=8, dpi=300, bg="white")

cat("\n=== FIGURAS GENERADAS ===\n")
system("ls -la figuras/mapas_bemmels/")
