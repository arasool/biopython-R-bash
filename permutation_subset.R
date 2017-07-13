library(dplyr)

#the purpose of this script is to calculate the difference of mean between two groups (treatment and control) and 
#apply permutations to difference of means. If the observed difference is higher than calculated 
#difference from permutations, then that group is taken into further consideration for the experiment. The group into 
#consideration with the most significant p-value is the one we are most interested in.

control1 <- NO.dmr1$value
treatment1 <- YES.dmr1$value
obsdiff1 <- mean(control1)-mean(treatment1)
obsdiff1

N <- 3995
avgdiff1 <- replicate(10000, {
  all <- sample(c(control1,treatment1))
  newcontrols1 <- all[1:N]
  newtreatments1 <- all[(N+1):(2*N)]
  return(mean(newcontrols1) - mean(newtreatments1))
})
hist(avgdiff1,
     xlim = c(-0.03, 0.03),
     col = "white",
     breaks = 100,
     main = "Diff.Mean | Sperm.v.Blood | pval= 0.0014",
     las = 2,
     xlab = "Mean Difference Distribution (Sperm.v.Blood)")
axis(1, at=seq(-0.03,0.03,by=0.01), labels=seq(-0.03,0.03,by=0.01), las=2 )
axis(2, at=seq(0,300,by=25), labels=seq(0,300,by=25), las=1 )
abline(v=obsdiff1, col = "red", lwd = 2)
#p-value, two tailed test
a <- avgdiff1[(abs(avgdiff1)) >= (abs(obsdiff1))]
length(a)
#13 values are >= obsdiff1
p1 <- (sum(abs(avgdiff1) >= abs(obsdiff1)) + 1) / (length(avgdiff1) + 1)



control2 <- NO.dmr2$value
treatment2 <- YES.dmr2$value
obsdiff2 <- mean(control2)-mean(treatment2)
obsdiff2

N <- 3995
avgdiff2 <- replicate(10000, {
  all <- sample(c(control2,treatment2))
  newcontrols2 <- all[1:N]
  newtreatments2 <- all[(N+1):(2*N)]
  return(mean(newcontrols2) - mean(newtreatments2))
})
hist(avgdiff2,
     xlim = c(-0.06, 0.06),
     col = "white", 
     breaks = 200,
     main = "Diff.Mean | Oocyte.v.Blood | pval= 0.00009",
     las = 2,
     xlab = "Diff.Mean Difference Distribution (Oocyte.v.Blood)")
axis(1, at=seq(-0.06,0.04,by=0.01), labels=seq(-0.06,0.04,by=0.01), las=2 )
axis(2, at=seq(0,2000,by=25), labels=seq(0,2000,by=25), las=1 )
abline(v=obsdiff2, col = "red", lwd = 2)
b <- avgdiff2[(abs(avgdiff2)) >= ((obsdiff2))]
length(b)
#0 values are >= obsdiff2
p2 <- (sum(abs(avgdiff2) >= abs(obsdiff2)) + 1) / (length(avgdiff2) + 1)



control3 <- NO.dmr3$value
treatment3 <- YES.dmr3$value
obsdiff3 <- mean(control3)-mean(treatment3)
obsdiff3

N <- 3995
avgdiff3 <- replicate(10000, {
  all <- sample(c(control3,treatment3))
  newcontrols3 <- all[1:N]
  newtreatments3 <- all[(N+1):(2*N)]
  return(mean(newcontrols3) - mean(newtreatments3))
})
hist(avgdiff3,
     xlim = c(-0.04, 0.12),
     col = "white",
     breaks = 150,
     las = 2,
     main = "Diff.Mean | Blastocyst.v.Blood | p-val = 0.00009",
     xlab = "Mean Difference Distribution (Blastocyst.v.Blood)")
axis(1, at=seq(-0.04, 0.12,by=0.01), labels=seq(-0.04, 0.12,by=0.01), las=2 )
axis(2, at=seq(0,300,by=50), labels=seq(0,300,by=50), las=1 )
abline(v=obsdiff3, col = "red", lwd = 2)
c <- avgdiff3[(abs(avgdiff3)) >= (abs(obsdiff3))]
length(c)
# 0 values are >= obsdiff3.
p3 <- (sum(abs(avgdiff3) >= abs(obsdiff3)) + 1) / (length(avgdiff3) + 1)

par(oma=c(3,3,0,0))
par(mfrow=c(1,3))
par(new=F)
par(mar=c(3,3,2,2))

mtext(text="Permutation Distribution - Difference in Mean of Methlation Fraction in Tissues",side=1,line=1,outer=TRUE)
mtext(text="Frequency",side=2,line=0,outer=TRUE)

##########################################################################################
#permutation test for median

control4 <- NO.dmr1$value
treatment4 <- YES.dmr1$value
obsdiff4 <- median(control4)-median(treatment4)
obsdiff4

N <- 3995
avgdiff4 <- replicate(10000, {
  all <- sample(c(control4,treatment4))
  newcontrols4 <- all[1:N]
  newtreatments4 <- all[(N+1):(2*N)]
  return(median(newcontrols4) - median(newtreatments4))
})
hist(avgdiff4,
     xlim = c(-0.007, 0.007),
     col = "white",
     breaks = 150,
     las = 2,
     main = "Diff.Median | Sperm.v.Blood | p-val = 0.0003",
     xlab = "Distribution of Difference in Median (Sperm.v.Blood)")
axis(1, at=seq(-0.007, 0.007,by=0.005), labels=seq(-0.007, 0.007,by=0.01), las=2 )
axis(2, at=seq(0,300,by=50), labels=seq(0,300,by=50), las=1 )
abline(v=obsdiff4, col = "blue", lwd = 2)
#legend("topleft", "actual.diff = 0.005", cex = 0.4)
d <- avgdiff4[(abs(avgdiff4)) >= ((obsdiff4))]
length(d)
# 2 values are >= obsdiff4.
p4 <- (sum(abs(avgdiff4) >= abs(obsdiff4)) + 1) / (length(avgdiff4) + 1)


control5 <- NO.dmr2$value
treatment5 <- YES.dmr2$value
obsdiff5 <- median(control5)-median(treatment5)
obsdiff5

N <- 3995
avgdiff5 <- replicate(10000, {
  all <- sample(c(control5,treatment5))
  newcontrols5 <- all[1:N]
  newtreatments5 <- all[(N+1):(2*N)]
  return(median(newcontrols5) - median(newtreatments5))
})
hist(avgdiff5,
     xlim = c(-0.01, 0.01),
     col = "white",
     breaks = 150,
     las = 2,
     main = "Diff.Median | Oocyte.v.Blood | p-val = 0.009",
     xlab = "Distribution of Difference in Median (Oocyte.v.Blood)")
axis(1, at=seq(-0.01, 0.01,by=0.005), labels=seq(-0.01, 0.01,by=0.01), las=2 )
axis(2, at=seq(0,300,by=50), labels=seq(0,300,by=50), las=1 )
abline(v=obsdiff5, col = "blue", lwd = 2)
e <- avgdiff5[(abs(avgdiff5)) >= ((obsdiff5))]
length(e)
# 100 values are >= obsdiff5.
p5 <- (sum(abs(avgdiff5) >= abs(obsdiff5)) + 1) / (length(avgdiff5) + 1)



control6 <- NO.dmr3$value
treatment6 <- YES.dmr3$value
obsdiff6 <- median(control6)-median(treatment6)
obsdiff6

N <- 3995
avgdiff6 <- replicate(10000, {
  all <- sample(c(control6,treatment6))
  newcontrols6 <- all[1:N]
  newtreatments6 <- all[(N+1):(2*N)]
  return(median(newcontrols6) - median(newtreatments6))
})
hist(avgdiff6,
     xlim = c(-0.1, 0.3),
     col = "white",
     breaks = 150,
     las = 2,
     main = "Diff.Median | Blastocyst.v.Blood | p-val = 0.00009",
     xlab = "Distribution of Difference in Median (Blastocyst.v.Blood)")
axis(1, at=seq(-0.1, 0.3,by=0.05), labels=seq(-0.1, 0.3,by=0.05), las=2 )
axis(2, at=seq(0,300,by=50), labels=seq(0,300,by=50), las=1 )
abline(v=obsdiff6, col = "blue", lwd = 2)
f <- avgdiff6[(abs(avgdiff6)) >= (abs(obsdiff6))]
length(f)
# 0 values are >= obsdiff6.
p6 <- (sum(abs(avgdiff6) >= abs(obsdiff6)) + 1) / (length(avgdiff6) + 1)

par(oma=c(3,3,0,0))
par(mfrow=c(1,3))
par(new=F)
par(mar=c(3,3,2,2))

mtext(text="Permutation Distribution - Difference in Median of Methylation Fraction in Tissues",side=1,line=1,outer=TRUE)
mtext(text="Frequency",side=2,line=0,outer=TRUE)