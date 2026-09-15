#!/usr/bin/env Rscript
# =====================================================================
# ANALISIS DEFINITIVO — Panel B, estadisticos H_isl + R1 = FST(I-M)/FST(global)
#
# Una sola corrida, semilla fija. Todos los gfit se calculan UNA VEZ
# y se reutilizan en tablas y figuras, de modo que los p-valores
# coinciden en todo el manuscrito.
#
# Salidas: FINAL_H_isl_R1_*.txt y FINAL_H_isl_R1_*.png
# =====================================================================

suppressMessages(library(abc))
SEED <- 20260811
set.seed(SEED)
setwd(path.expand("~/Splatche/datasetC"))

TOLS <- c(0.02, 0.05, 0.10, 0.20)
NBREP <- 200          # pseudo-datasets para gfit
NVAL  <- 200          # pseudo-observados para cv4postpr
S     <- c("H_isl","R1")
BASE  <- c("HT","FST_glob","FST_IM","H_isl","H_main","FIS_isl")

s1 <- read.table("sim1_ext_thinB.txt", header=TRUE, check.names=FALSE)
s2 <- read.table("sim2_ext_thinB.txt", header=TRUE, check.names=FALSE)
ob <- read.table("observed_ext_thinB.txt", header=TRUE, check.names=FALSE)
obd <- as.data.frame(ob)
pa <- read.table("param.txt", col.names=c("ID","m","K1","K2"))
p1 <- pa[match(s1$simid, pa$ID), c("m","K1")]; names(p1) <- c("m","K")
p2 <- pa[match(s2$simid, pa$ID), c("m","K2")]; names(p2) <- c("m","K")
s1$R1  <- s1$FST_IM / s1$FST_glob
s2$R1  <- s2$FST_IM / s2$FST_glob
obd$R1 <- obd$FST_IM / obd$FST_glob
k1 <- is.finite(s1$R1); k2 <- is.finite(s2$R1)
s1 <- s1[k1,]; p1 <- p1[k1,]
s2 <- s2[k2,]; p2 <- p2[k2,]

M1 <- as.matrix(s1[,S]); M2 <- as.matrix(s2[,S])
tg <- as.numeric(obd[1,S]); names(tg) <- S
mod <- c(rep("Refugium", nrow(M1)), rep("HumanBridge", nrow(M2)))
SS  <- rbind(M1, M2)

wr <- function(x,f){ write.table(x,f,sep="\t",row.names=FALSE,quote=FALSE)
                     cat("   ->",f,"\n") }

sink("FINAL_H_isl_R1_log.txt", split=TRUE)
cat("=====================================================\n")
cat("Panel B — H_isl + R1   ·  seed =", SEED, "\n")
cat("=====================================================\n")
cat("Replicates: Refugium", nrow(M1), " Human Bridge", nrow(M2), "\n")
cat(sprintf("Correlation between statistics: Refugium %.3f, Human Bridge %.3f\n",
            cor(M1[,1],M1[,2]), cor(M2[,1],M2[,2])))

# ------------------------------------------------ GOF: UNA sola vez
cat("\n=== TABLE 7 — Goodness of fit ===\n")
GOF <- list()
for (sc in c("Refugium","HumanBridge")) {
  M <- if (sc=="Refugium") M1 else M2
  set.seed(SEED)                       # misma semilla para cada escenario
  GOF[[sc]] <- gfit(target=tg, sumstat=M, statistic=median,
                    nb.replicate=NBREP, tol=0.05)
}
t7 <- do.call(rbind, lapply(names(GOF), function(sc){
  sm <- summary(GOF[[sc]])
  data.frame(Scenario=sc, p_value=round(sm$pvalue,4),
             dist_observed=round(sm$dist.obs,3),
             null_median=round(median(GOF[[sc]]$dist.sim),3),
             Fit=ifelse(sm$pvalue>0.05,"adequate","rejected"))
}))
print(t7, row.names=FALSE); wr(t7,"FINAL_H_isl_R1_Table7_GOF.txt")

# ------------------------------------------------ Table 5
cat("\n=== TABLE 5 — Summary statistics ===\n")
t5 <- do.call(rbind, lapply(BASE, function(x){
  o <- obd[[x]]; q1 <- quantile(s1[[x]],c(.025,.975)); q2 <- quantile(s2[[x]],c(.025,.975))
  data.frame(Statistic=x, Used=ifelse(x %in% S,"yes","no"), Observed=round(o,4),
             Refugium=sprintf("%.4f - %.4f",q1[1],q1[2]),
             In_Ref=ifelse(o>=q1[1]&o<=q1[2],"Yes","No"),
             HumanBridge=sprintf("%.4f - %.4f",q2[1],q2[2]),
             In_HB=ifelse(o>=q2[1]&o<=q2[2],"Yes","No"))
}))
print(t5, row.names=FALSE); wr(t5,"FINAL_H_isl_R1_Table5.txt")

# ------------------------------------------------ Table 6 — priors
cat("\n=== TABLE 6 — Priors ===\n")
t6 <- data.frame(
  Parameter = c("m","K_refuge","K_corridor"),
  Prior     = c("U(0.01, 0.10)","U(2, 200)","U(312, 44061)"),
  Scenario  = c("both","Refugium","Human Bridge"),
  Source    = c("Marin et al. 2013; this study",
                "Zubillaga et al. 2018, scaled to refugial cell",
                "Zubillaga et al. 2018, scaled to corridor cell"))
print(t6, row.names=FALSE); wr(t6,"FINAL_H_isl_R1_Table6_priors.txt")

# ------------------------------------------------ Table 8 — model choice
cat("\n=== TABLE 8 — Model choice ===\n")
set.seed(SEED)
t8 <- data.frame()
for (met in c("rejection","mnlogistic","neuralnet")) for (tl in TOLS) {
  r <- tryCatch(postpr(target=tg,index=mod,sumstat=SS,tol=tl,method=met),
                error=function(e) NULL)
  if (is.null(r)) next
  sm <- summary(r, print=FALSE)
  pp <- if (!is.null(sm$Prob)) sm$Prob else sm[[met]]$Prob
  if (is.null(pp)) pp <- sm$rejection$Prob
  pr <- as.numeric(pp["Refugium"]); pb <- as.numeric(pp["HumanBridge"])
  if (!length(pr) || is.na(pr)) next
  nHB <- round(pb * tl * nrow(SS))
  t8 <- rbind(t8, data.frame(Method=met, Tolerance=tl,
    P_Refugium=round(pr,4), P_HumanBridge=round(pb,4),
    BayesFactor=if (pb>0) sprintf("%.1f",pr/pb) else "not bounded",
    n_HB_accepted=nHB))
}
print(t8, row.names=FALSE); wr(t8,"FINAL_H_isl_R1_Table8_modelchoice.txt")
cat("\nNote: at tolerances <= 0.10 no Human Bridge replicate is accepted,\n")
cat("so the Bayes factor is not bounded. Only tol = 0.20 is interpretable.\n")

# ------------------------------------------------ Table 9 — parameters
cat("\n=== TABLE 9 — Parameter estimates ===\n")
set.seed(SEED)
t9 <- data.frame()
for (sc in c("Refugium","HumanBridge")) {
  M <- if (sc=="Refugium") M1 else M2
  pp <- if (sc=="Refugium") p1 else p2
  for (tl in TOLS) {
    r <- abc(target=tg,param=pp,sumstat=M,tol=tl,method="rejection")
    a <- r$unadj.values
    for (q in colnames(a)) {
      v <- a[,q]; d <- density(v); qq <- quantile(v,c(.025,.5,.975))
      t9 <- rbind(t9, data.frame(Scenario=sc,Parameter=q,Tolerance=tl,
        Mode=signif(d$x[which.max(d$y)],4), Median=signif(qq[2],4),
        CI95=sprintf("%.4g - %.4g",qq[1],qq[3])))
    }
  }
}
print(t9, row.names=FALSE); wr(t9,"FINAL_H_isl_R1_Table9_parameters.txt")

# ------------------------------------------------ Tables S
cat("\n=== TABLE S1 — Prediction error ===\n")
set.seed(SEED)
ts <- data.frame()
for (sc in c("Refugium","HumanBridge")) {
  M <- if (sc=="Refugium") M1 else M2
  pp <- if (sc=="Refugium") p1 else p2
  e <- summary(cv4abc(param=pp,sumstat=M,nval=100,tols=TOLS,method="rejection"))
  ts <- rbind(ts, data.frame(Scenario=sc,Tolerance=TOLS,
    err_m=round(e[,"m"],4), err_K=round(e[,"K"],4)))
}
print(ts, row.names=FALSE); wr(ts,"FINAL_H_isl_R1_TableS1_prederror.txt")
cat("\nPrediction error = MSE(estimate, true) / Var(true), leave-one-out\n")
cat("over 100 pseudo-observed datasets; equivalent to 1 - R^2.\n")

cat("\n=== TABLE S2 — Cross-validation of model choice ===\n")
set.seed(SEED)
cv <- cv4postpr(mod, SS, nval=NVAL, tol=0.05, method="rejection")
print(summary(cv))
capture.output(print(summary(cv)), file="FINAL_H_isl_R1_TableS2_crossval.txt")

# ------------------------------------------------ FIGURES
CREF <- "#B03A2E"; CHB <- "#2C3E50"; COBS <- "#E67E22"
labs <- list(expression(H[island]), expression(R[1]==F[ST]*(I-M)/F[ST]*(global)))

png("FINAL_H_isl_R1_Fig_main.png", 2400, 1600, res=200)
layout(matrix(1:4,2,2,byrow=TRUE)); par(mar=c(4.5,4.5,3,1))
for (j in 1:2) {
  d1 <- density(M1[,j]); d2 <- density(M2[,j])
  xr <- range(c(quantile(M1[,j],c(.002,.998)),quantile(M2[,j],c(.002,.998)),tg[j]))
  plot(NA,xlim=xr,ylim=c(0,max(d1$y,d2$y)*1.12),xlab=labs[[j]],ylab="Density",
       main=paste0("(",c("a","b")[j],") ",S[j]))
  polygon(d1,col=adjustcolor(CREF,.35),border=CREF,lwd=2)
  polygon(d2,col=adjustcolor(CHB,.35),border=CHB,lwd=2)
  abline(v=tg[j],col=COBS,lwd=3,lty=2)
  if (j==1) legend("topleft",c("Refugium","Human Bridge","Observed"),
    fill=c(adjustcolor(CREF,.35),adjustcolor(CHB,.35),NA),border=c(CREF,CHB,NA),
    lty=c(NA,NA,2),lwd=c(NA,NA,3),col=c(NA,NA,COBS),bty="n",cex=.85)
}
# GOF: reutiliza los objetos ya calculados -> p-valores identicos a Table 7
for (i in 1:2) {
  sc <- c("Refugium","HumanBridge")[i]; g <- GOF[[sc]]; sm <- summary(g)
  hist(g$dist.sim,breaks=40,col="grey85",border="grey60",
       xlim=c(0, max(c(g$dist.sim,sm$dist.obs))*1.15),
       xlab="Distance to model median",ylab="Frequency",
       main=sprintf("(%s) %s  (p = %.3f)",c("c","d")[i],sc,sm$pvalue))
  abline(v=sm$dist.obs,col=COBS,lwd=3)
}
dev.off(); cat("\n   -> FINAL_H_isl_R1_Fig_main.png\n")

png("FINAL_H_isl_R1_Fig_posteriors.png", 2000, 900, res=200)
par(mfrow=c(1,2),mar=c(4.5,4.5,3,1))
set.seed(SEED)
r <- abc(target=tg,param=p1,sumstat=M1,tol=0.05,method="rejection")
po <- r$unadj.values
for (q in c("m","K")) {
  pri <- if (q=="m") pa$m else pa$K1
  dpo <- density(po[,q]); dpr <- density(pri)
  plot(NA,xlim=range(pri),ylim=c(0,max(dpo$y,dpr$y)*1.1),
       xlab=if(q=="m") "Migration rate (m)" else expression(K[refuge]),
       ylab="Density",
       main=if(q=="m") "(a) Migration rate" else "(b) Refugial carrying capacity")
  polygon(dpr,col=adjustcolor("grey60",.3),border="grey50",lwd=2,lty=2)
  polygon(dpo,col=adjustcolor(CREF,.4),border=CREF,lwd=2.5)
  abline(v=median(po[,q]),col=CREF,lwd=2,lty=3)
  legend("topright",c("Prior","Posterior"),
         fill=c(adjustcolor("grey60",.3),adjustcolor(CREF,.4)),
         border=c("grey50",CREF),bty="n",cex=.9)
}
dev.off(); cat("   -> FINAL_H_isl_R1_Fig_posteriors.png\n")

cat("\n=====================================================\n")
cat("DONE — seed", SEED, "; all p-values in tables and figures\n")
cat("come from the same gfit objects.\n")
cat("=====================================================\n")
sink()
