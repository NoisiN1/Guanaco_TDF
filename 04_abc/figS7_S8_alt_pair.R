#!/usr/bin/env Rscript
# =====================================================================
# Figuras para el manuscrito - panel 1,773 SNPs, .par con MAF empirico
# Genera: Fig 6 (principal), Fig S6, S7, S8, S9 (suplementario)
# =====================================================================

suppressMessages(library(abc)); set.seed(1)
setwd(path.expand("~/Splatche/datasetC"))

s1 <- read.table("sim1_ext.txt", header=TRUE, check.names=FALSE)
s2 <- read.table("sim2_ext.txt", header=TRUE, check.names=FALSE)
ob <- read.table("observed_ext.txt", header=TRUE, check.names=FALSE)
pa <- read.table("param.txt", col.names=c("ID","m","K1","K2"))

fp  <- setdiff(grep("^FST_[a-z]", names(ob), value=TRUE), c("FST_glob","FST_IM"))
isl <- fp[!grepl("_sg$|^FST_sg", fp)]
cont<- fp[ grepl("_sg$|^FST_sg", fp)]

R <- list(
  "R1" = function(d) d$FST_IM/d$FST_glob,
  "R2" = function(d) rowMeans(d[,isl,drop=FALSE])/d$FST_IM,
  "R3" = function(d) rowMeans(d[,isl,drop=FALSE])/rowMeans(d[,cont,drop=FALSE]),
  "R4" = function(d) d$FST_glob/d$HT,
  "R5" = function(d) d$FST_IM/d$HT)
RLAB <- c(R1="F[ST](I-M)/F[ST](global)", R2="F[ST](island)/F[ST](I-M)",
          R3="F[ST](island)/F[ST](I-C)", R4="F[ST](global)/H[T]",
          R5="F[ST](I-M)/H[T]")

M1 <- sapply(R, function(f) f(s1)); M2 <- sapply(R, function(f) f(s2))
tg <- sapply(R, function(f) f(as.data.frame(ob)))
mod <- c(rep("Refugium", nrow(M1)), rep("HumanBridge", nrow(M2)))
SS  <- rbind(M1, M2)

CREF <- "#B03A2E"; CHB <- "#2C3E50"; COBS <- "#E67E22"

# ---------------------------------------------------------------- Fig 6
# Panel principal: distribuciones de los 2 estadisticos retenidos + GOF
png("FIG6_main.png", 2200, 1650, res=200)
layout(matrix(c(1,2,3,4), 2, 2, byrow=TRUE))
par(mar=c(4.5,4.5,3,1))

for (k in c("R1","R2")) {
  d1 <- density(M1[,k]); d2 <- density(M2[,k])
  xr <- range(c(d1$x, d2$x, tg[k]))
  plot(NA, xlim=xr, ylim=c(0, max(d1$y,d2$y)*1.12),
       xlab=parse(text=RLAB[k]), ylab="Density",
       main=bquote(bold(.(ifelse(k=="R1","(a)","(b)")))~.(gsub("\\[|\\]","",RLAB[k]))))
  polygon(d1, col=adjustcolor(CREF,.35), border=CREF, lwd=2)
  polygon(d2, col=adjustcolor(CHB,.35),  border=CHB,  lwd=2)
  abline(v=tg[k], col=COBS, lwd=3, lty=2)
  legend("topright", c("Refugium","Human Bridge","Observed"),
         fill=c(adjustcolor(CREF,.35), adjustcolor(CHB,.35), NA),
         border=c(CREF,CHB,NA), lty=c(NA,NA,2), lwd=c(NA,NA,3),
         col=c(NA,NA,COBS), bty="n", cex=.85)
}

# GOF de cada escenario
for (i in 1:2) {
  nm <- c("Refugium","HumanBridge")[i]
  M  <- if (i==1) M1[,c("R1","R2")] else M2[,c("R1","R2")]
  g  <- gfit(target=tg[c("R1","R2")], sumstat=M, statistic=median,
             nb.replicate=200, tol=0.05)
  sm <- summary(g)
  hist(g$dist.sim, breaks=40, col="grey85", border="grey60",
       xlim=range(c(g$dist.sim, sm$dist.obs))*1.05,
       xlab="Distance to model median", ylab="Frequency",
       main=bquote(bold(.(ifelse(i==1,"(c)","(d)")))~.(nm)~~ (italic(p)==.(sprintf("%.2f",sm$pvalue)))))
  abline(v=sm$dist.obs, col=COBS, lwd=3)
  text(sm$dist.obs, par("usr")[4]*.85, "observed", col=COBS, pos=if(i==1) 4 else 2, cex=.85)
}
dev.off(); cat("FIG6_main.png\n")

# ------------------------------------------------------------- Fig S6
# Los 5 cocientes: observado vs distribuciones simuladas
png("FIGS6_all_ratios.png", 2400, 1500, res=200)
par(mfrow=c(2,3), mar=c(4.3,4.3,3,1))
for (k in names(R)) {
  d1 <- density(M1[,k]); d2 <- density(M2[,k])
  xr <- range(c(quantile(M1[,k],c(.001,.999)), quantile(M2[,k],c(.001,.999)), tg[k]))
  plot(NA, xlim=xr, ylim=c(0,max(d1$y,d2$y)*1.1),
       xlab=parse(text=RLAB[k]), ylab="Density", main=k)
  polygon(d1, col=adjustcolor(CREF,.35), border=CREF, lwd=1.8)
  polygon(d2, col=adjustcolor(CHB,.35),  border=CHB,  lwd=1.8)
  abline(v=tg[k], col=COBS, lwd=2.5, lty=2)
}
plot.new()
legend("center", c("Refugium","Human Bridge","Observed"),
       fill=c(adjustcolor(CREF,.35), adjustcolor(CHB,.35), NA),
       border=c(CREF,CHB,NA), lty=c(NA,NA,2), lwd=c(NA,NA,3),
       col=c(NA,NA,COBS), bty="n", cex=1.2)
dev.off(); cat("FIGS6_all_ratios.png\n")

# ------------------------------------------------------------- Fig S7
# Estadisticos crudos: observado vs simulado
png("FIGS7_raw_stats.png", 2400, 1500, res=200)
B <- c("HT","FST_glob","FST_IM","H_isl","H_main","FIS_isl")
BL <- c("H[T]","F[ST]~(global)","F[ST]~(I-M)","H[isl]","H[main]","F[IS]~(isl)")
par(mfrow=c(2,3), mar=c(4.3,4.3,3,1))
for (i in seq_along(B)) {
  s <- B[i]; d1 <- density(s1[[s]]); d2 <- density(s2[[s]])
  xr <- range(c(d1$x,d2$x,ob[[s]]))
  plot(NA, xlim=xr, ylim=c(0,max(d1$y,d2$y)*1.1),
       xlab=parse(text=BL[i]), ylab="Density", main=parse(text=BL[i]))
  polygon(d1, col=adjustcolor(CREF,.35), border=CREF, lwd=1.8)
  polygon(d2, col=adjustcolor(CHB,.35),  border=CHB,  lwd=1.8)
  abline(v=ob[[s]], col=COBS, lwd=2.5, lty=2)
}
dev.off(); cat("FIGS7_raw_stats.png\n")

# ------------------------------------------------------------- Fig S8
# Prior vs posterior de m y K
p1 <- pa[match(s1$simid, pa$ID), c("m","K1")]; names(p1) <- c("m","K")
r  <- abc(target=tg[c("R1","R2")], param=p1, sumstat=M1[,c("R1","R2")],
          tol=0.05, method="rejection")
po <- r$unadj.values
png("FIGS8_posteriors.png", 2000, 900, res=200)
par(mfrow=c(1,2), mar=c(4.5,4.5,3,1))
for (q in c("m","K")) {
  pri <- if (q=="m") pa$m else pa$K1
  dpo <- density(po[,q]); dpr <- density(pri)
  plot(NA, xlim=range(pri), ylim=c(0,max(dpo$y,dpr$y)*1.1),
       xlab=if(q=="m") "Migration rate (m)" else expression(K[refuge]),
       ylab="Density", main=if(q=="m") "(a) Migration rate" else "(b) Refugial carrying capacity")
  polygon(dpr, col=adjustcolor("grey60",.3), border="grey50", lwd=2, lty=2)
  polygon(dpo, col=adjustcolor(CREF,.4), border=CREF, lwd=2.5)
  abline(v=median(po[,q]), col=CREF, lwd=2, lty=3)
  legend("topright", c("Prior","Posterior"),
         fill=c(adjustcolor("grey60",.3), adjustcolor(CREF,.4)),
         border=c("grey50",CREF), bty="n", cex=.9)
}
dev.off(); cat("FIGS8_posteriors.png\n")

# ------------------------------------------------------------- Fig S9
# Robustez: GOF y P(Refugium) para cada conjunto de estadisticos
rb <- tryCatch(read.table("ROBUST_resultados.txt", header=TRUE, sep="\t"),
               error=function(e) NULL)
if (!is.null(rb)) {
  rb <- rb[!is.na(rb$GOF_Ref),]
  png("FIGS9_robustness.png", 2000, 1100, res=200)
  par(mfrow=c(1,2), mar=c(9,4.5,3,1))
  bp <- barplot(rbind(rb$GOF_Ref, rb$GOF_HB), beside=TRUE,
                col=c(CREF,CHB), border=NA, ylim=c(0,1.05),
                ylab="Goodness-of-fit p-value", main="(a) Model adequacy",
                names.arg=rep("", nrow(rb)))
  abline(h=0.05, lty=2, col="grey40", lwd=2)
  text(colMeans(bp), -0.03, rb$conjunto, srt=45, adj=1, xpd=TRUE, cex=.65)
  legend("topright", c("Refugium","Human Bridge"), fill=c(CREF,CHB),
         border=NA, bty="n", cex=.85)
  bp2 <- barplot(rb$P_Refugium, col=CREF, border=NA, ylim=c(0,1.05),
                 ylab="P(Refugium)", main="(b) Posterior model probability",
                 names.arg=rep("", nrow(rb)))
  abline(h=0.5, lty=2, col="grey40", lwd=2)
  text(bp2, -0.03, rb$conjunto, srt=45, adj=1, xpd=TRUE, cex=.65)
  dev.off(); cat("FIGS9_robustness.png\n")
}

cat("\nFiguras generadas.\n")
