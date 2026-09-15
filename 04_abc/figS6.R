library(abc); library(MASS)
a <- read.table("sim1_multi_NEW.txt", header=TRUE)
p <- read.table("param.txt")
a$K <- p[a$simid,3]; a$m <- p[a$simid,2]
a <- a[complete.cases(a),]
OBS <- 0.098301
tol <- 0.05
acc <- a[order(abs(a$FST_glob-OBS))[1:round(tol*nrow(a))],]
po  <- abc(target=OBS, param=a[,c("K","m")], sumstat=a$FST_glob,
           tol=0.01, method="loclinear")
ci  <- quantile(po$adj.values[,"m"], c(.025,.5,.975))
png("FigS6_JointPosterior.png", 1500, 1300, res=260)
par(mar=c(4.4,4.4,1.2,1.2))
plot(acc$K, acc$m, pch=16, cex=.45, col="#2C6E9B80",
     xlab="K (refugial deme)", ylab="m (migration rate)",
     xlim=c(0,200), ylim=c(0.025,0.10))
d <- kde2d(acc$K, acc$m, n=80)
contour(d, add=TRUE, drawlabels=FALSE, col="grey35", lwd=1.0,
        nlevels=5)
abline(h=ci[c(1,3)], lty=2, col="firebrick", lwd=1.4)
abline(h=ci[2], lty=1, col="firebrick", lwd=1.8)
legend("topright", c("accepted simulations","density contours",
       "m: median and 95% CI"), col=c("#2C6E9B","grey25","firebrick"),
       pch=c(16,NA,NA), lty=c(NA,1,1), lwd=c(NA,1.1,2), bty="n", cex=.72)
dev.off()
cat(sprintf("m: mediana %.4f  IC [%.4f, %.4f]\n", ci[2], ci[1], ci[3]))
