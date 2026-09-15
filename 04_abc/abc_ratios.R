suppressMessages(library(abc)); set.seed(1)
setwd(path.expand("~/Splatche/datasetC"))
s1 <- read.table("sim1_ext.txt",header=TRUE,check.names=FALSE)
s2 <- read.table("sim2_ext.txt",header=TRUE,check.names=FALSE)
ob <- read.table("observed_ext.txt",header=TRUE,check.names=FALSE)
fst <- grep("^FST_[a-z]", names(ob), value=TRUE); fst <- setdiff(fst,c("FST_glob","FST_IM"))
isl <- fst[!grepl("_sg$|^FST_sg", fst)]
mk <- function(d) cbind(R1=d$FST_IM/d$FST_glob, R2=rowMeans(d[,isl])/d$FST_IM)
S1 <- mk(s1); S2 <- mk(s2)
tg <- c(R1=ob$FST_IM/ob$FST_glob, R2=mean(as.numeric(ob[1,isl]))/ob$FST_IM)
ok <- function(M) apply(M,1,function(r) all(is.finite(r)))
k1 <- ok(S1); k2 <- ok(S2); S1<-S1[k1,]; S2<-S2[k2,]
pa <- read.table("param.txt", col.names=c("ID","m","K1","K2"))
p1 <- pa[match(s1$simid[k1],pa$ID),c("m","K1")]; names(p1)<-c("m","K")
p2 <- pa[match(s2$simid[k2],pa$ID),c("m","K2")]; names(p2)<-c("m","K")
cat("replicas: s1=",nrow(S1)," s2=",nrow(S2),"\n\n=== GOF ===\n")
for (nm in c("Refugium","HumanBridge")) {
  S <- if(nm=="Refugium") S1 else S2
  r <- tryCatch(gfit(target=tg,sumstat=S,statistic=median,nb.replicate=200,tol=0.05),error=function(e)NULL)
  if(!is.null(r)) cat(sprintf("%-12s p=%.4f  dist=%.4f\n",nm,summary(r)$pvalue,summary(r)$dist.obs))
}
mod <- c(rep("Refugium",nrow(S1)),rep("HumanBridge",nrow(S2))); SS <- rbind(S1,S2)
cat("\n=== VALIDACION CRUZADA ===\n")
cv <- tryCatch(cv4postpr(mod,SS,nval=200,tol=0.05,method="mnlogistic"),error=function(e)NULL)
if(!is.null(cv)) print(summary(cv))
cat("\n=== SELECCION DE MODELO ===\n")
for (met in c("rejection","mnlogistic","neuralnet")) for (tl in c(0.02,0.05,0.10,0.20)) {
  r <- tryCatch(postpr(target=tg,index=mod,sumstat=SS,tol=tl,method=met),error=function(e)NULL)
  if(is.null(r)) next
  pp <- summary(r,print=FALSE)$Prob; pr<-as.numeric(pp["Refugium"]); pb<-as.numeric(pp["HumanBridge"])
  if(!length(pr)||is.na(pr)) next
  cat(sprintf(" %-11s tol=%.2f P_Ref=%.4f P_Bridge=%.4f BF=%s\n",met,tl,pr,pb,format(if(pb>0)pr/pb else Inf,digits=4)))
}
