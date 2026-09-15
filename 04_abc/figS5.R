library(abc)
a <- read.table("sim1_multi_NEW.txt", header=TRUE)
b <- read.table("sim2_multi_NEW.txt", header=TRUE)
a <- a[complete.cases(a),]; b <- b[complete.cases(b),]
OBS <- 0.098301
XL  <- c(0.04, 0.15)                       # misma escala en ambos paneles
png("FigS5_GOF.png", 2200, 1000, res=250)
par(mfrow=c(1,2), mar=c(4.4,4.4,2.8,1))
for (nm in c("Refugium","Human bridge")) {
  d <- if (nm=="Refugium") a$FST_glob else b$FST_glob
  acc <- d[order(abs(d-OBS))][1:round(0.05*length(d))]
  br <- seq(min(d), max(d), length.out=80)
  hist(d, breaks=br, col="grey80", border=NA, main="", xlim=XL, ylim=c(0,1200),
       xlab=expression(F[ST]), ylab="Frequency")
  hist(acc, breaks=br, col="#2C6E9B", border=NA, add=TRUE)
  abline(v=OBS, lwd=2.4)
  set.seed(1)                              # semilla justo antes de gfit
  g <- gfit(target=OBS, sumstat=d, nb.replicate=1000, tol=0.01)
  p <- summary(g)$pvalue
  mtext(sprintf("%s   (p = %.3f)", nm, p), side=3, adj=0, line=0.8,
        font=2, cex=.9)
  cat(nm, "p =", round(p,3), "| max simulado =", round(max(d),4), "\n")
}
dev.off()
