#Explore the ability of the SEM to correlate with ChIP scores

results <- read.delim("~/SEM_project/src/ChIP_Plot/results/FOXA1/BASELINE/Enumerated_kmer_filtered.chiphit", header=TRUE)
colnames(results) <- c("chr", "start", "end", "kmer", "PWM", "SEM", "ChIP", "peak")
plot(results$SEM, results$ChIP)
std <- function(x) sd(x)/sqrt(length(x))

####### new threshold method
bins <- as.matrix(table(results$kmer))
bins<-cbind(bins, rep(NA,length(bins)), rep(NA,length(bins)), rep(NA,length(bins)), rep(NA,length(bins)))
colnames(bins) <- c("Count", "PWM", "SEM", "ChIP", "Err")
threshold.val <- 13.5446
for(i in 1:nrow(bins)) {
  resultSet <- results[which(results$kmer == rownames(bins)[i]),]
  if(nrow(resultSet) < 1) {
  } else {
    bins[i,"PWM"] <- mean(resultSet$PWM) #all the same
    bins[i,"SEM"] <- mean(resultSet$SEM) #all the same
    bins[i,"ChIP"] <- mean(resultSet$ChIP)
    bins[i,"Err"] <- std(resultSet$ChIP)
  }
}

bins.thresh <- bins[which(bins[,"Count"] > 50),]

#Figure for paper ---START
plot(bins.thresh[,"PWM"], bins.thresh[,"ChIP"], xlab="PWM kmers", ylab="Avg ChIP-seq signal",main="FOXA1 ChIP vs PWM")
arrows(bins.thresh[,"PWM"], bins.thresh[,"ChIP"]-bins.thresh[,"Err"], bins.thresh[,"PWM"], bins.thresh[,"ChIP"]+bins.thresh[,"Err"], length = 0.05, angle = 90, code=3, col = "slategray")
abline(v=threshold.val)
txt <- round(cor(bins.thresh[,"PWM"], bins.thresh[,"ChIP"])^2, 2)
txt <- bquote(r^2 == .(txt))
text(10, 3, txt, cex=1)
txt <- round(cor(bins.thresh[,"PWM"], bins.thresh[,"ChIP"], method="spearman"), 2)
txt <- bquote(rho == .(txt))
text(10, 2.9, txt, cex=1)

plot(2^bins.thresh[,"SEM"], bins.thresh[,"ChIP"], xlab="SEM kmers", ylab="Avg ChIP-seq signal",main="FOXA1 ChIP vs SEM")
arrows(2^bins.thresh[,"SEM"], bins.thresh[,"ChIP"]-bins.thresh[,"Err"], 2^bins.thresh[,"SEM"], bins.thresh[,"ChIP"]+bins.thresh[,"Err"], length = 0.05, angle = 90, code=3, col = "slategray")
abline(v=1)
txt <- round(cor(2^bins.thresh[,"SEM"], bins.thresh[,"ChIP"])^2, 2)
txt <- bquote(r^2 == .(txt))
text(.5, 3, txt, cex=1)
txt <- round(cor(2^bins.thresh[,"SEM"], bins.thresh[,"ChIP"], method="spearman"), 2)
txt <- bquote(rho == .(txt))
text(.5, 2.9, txt, cex=1)

#Figure for paper ---End

# 
# #correlation with SEM cutoff -  this improves R^2
# cor(2^bins.thresh[(which(bins.thresh[,"SEM"]>=-.5)),"SEM"], bins.thresh[(which(bins.thresh[,"SEM"]>=-.5)),"ChIP"])^2
# 
# #correlation with ChIP cutoff - this doesnt change much
# cor(2^bins.thresh[(which(bins.thresh[,"ChIP"]>=-1)),"SEM"], bins.thresh[(which(bins.thresh[,"ChIP"]>=-1)),"ChIP"])^2
# cor(2^bins.thresh[(which(bins.thresh[,"ChIP"]>=-1)),"PWM"], bins.thresh[(which(bins.thresh[,"ChIP"]>=-1)),"ChIP"])^2

rbPal <- colorRampPalette(c('red','blue'))
pwm.col <- rbPal(10)[as.numeric(cut(bins.thresh[,"PWM"], breaks=10))]
plot(2^bins.thresh[,"SEM"], bins.thresh[,"ChIP"], col=pwm.col, xlab="SEM kmers", ylab="Avg ChIP-seq signal",main="FOXA1 ChIP vs SEM", type="n")
arrows(2^bins.thresh[,"SEM"], bins.thresh[,"ChIP"]-bins.thresh[,"Err"], 2^bins.thresh[,"SEM"], bins.thresh[,"ChIP"]+bins.thresh[,"Err"], length = 0.05, angle = 90, code=3, col = "slategray")
points(2^bins.thresh[,"SEM"], bins.thresh[,"ChIP"], col=pwm.col)
#points by count size
#points(2^bins.thresh[,"SEM"], bins.thresh[,"ChIP"], col=pwm.col, cex=bins.thresh[,"Count"]/mean(bins.thresh[,"Count"]))
abline(v=1)
abline(v=0.8501, lty=2)
txt <- round(cor(2^bins.thresh[,"SEM"], bins.thresh[,"ChIP"])^2, 2)
txt <- bquote(r^2 == .(txt))
text(.5, 3, txt, cex=1)
txt <- round(cor(2^bins.thresh[,"SEM"], bins.thresh[,"ChIP"], method="spearman"), 2)
txt <- bquote(rho == .(txt))
text(.5, 2.9, txt, cex=1)
txt <- round(cor(2^bins.thresh[(which(2^bins.thresh[,"SEM"]>=0.8501)),"SEM"], bins.thresh[(which(2^bins.thresh[,"SEM"]>=0.8501)),"ChIP"])^2, 2)
txt <- bquote(thresh_r^2 == .(txt))
text(0.7, 3, txt, cex=1)
txt <- round(cor(2^bins.thresh[(which(2^bins.thresh[,"SEM"]>=0.8501)),"SEM"], bins.thresh[(which(2^bins.thresh[,"SEM"]>=0.8501)),"ChIP"], method="spearman"), 2)
txt <- bquote(thresh_rho == .(txt))
text(0.7, 2.9, txt, cex=1)

# 
# rbPal <- colorRampPalette(c('red','blue'))
# pwm.col <- rbPal(10)[as.numeric(cut(bins.thresh[,"PWM"], breaks=10))]
# plot(bins.thresh[,"SEM"], bins.thresh[,"ChIP"], xlab="SEM kmers", ylab="Avg ChIP-seq signal",main="FOXA1 ChIP vs SEM", type="n")
# arrows(bins.thresh[,"SEM"], bins.thresh[,"ChIP"]-bins.thresh[,"Err"], bins.thresh[,"SEM"], bins.thresh[,"ChIP"]+bins.thresh[,"Err"], length = 0.05, angle = 90, code=3, col = "slategray")
# abline(v=0,lty=1)
# abline(v=-0.234,lty=2)
# points(bins.thresh[,"SEM"], bins.thresh[,"ChIP"],col=pwm.col)
# txt <- round(cor(bins.thresh[,"SEM"], bins.thresh[,"ChIP"])^2, 2)
# txt <- bquote(r^2 == .(txt))
# text(-1, 3, txt, cex=1)
# txt <- round(cor(2^bins.thresh[,"SEM"], bins.thresh[,"ChIP"], method="spearman"), 2)
# txt <- bquote(rho == .(txt))
# text(-1, 2.9, txt, cex=1)
# txt <- round(cor(bins.thresh[(which(bins.thresh[,"SEM"]>=-.234)),"SEM"], bins.thresh[(which(bins.thresh[,"SEM"]>=-.234)),"ChIP"])^2, 2)
# txt <- bquote(thresh_r^2 == .(txt))
# text(-0.7, 3, txt, cex=1)
# txt <- round(cor(bins.thresh[(which(bins.thresh[,"SEM"]>=-.234)),"SEM"], bins.thresh[(which(bins.thresh[,"SEM"]>=-.234)),"ChIP"], method="spearman"), 2)
# txt <- bquote(thresh_rho == .(txt))
# text(-0.7, 2.9, txt, cex=1)

