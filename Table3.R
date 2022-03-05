########
## function to compute the mcmc output
para_summary <- function(mcmc,a,b,print){
  y <- matrix(NA,ncol(mcmc),4)
  for (i in 1:ncol(mcmc)){
    y[i,1:3] <- quantile(mcmc[,i],c(0.5,0.025,0.975))
    y[i,4] <- sum(diff(mcmc[,i])!=0)/nrow(mcmc)
  }
  layout(matrix(1:(a*b),nrow=a,byrow=T))
  par(mar=c(2,4,1,1))
  if (print==1){
    for (i in 1:ncol(mcmc)){
      plot(mcmc[,i],type="l")
    }
  }
  return(y)
}

incmat <- matrix(NA,4,14)
infmat <- matrix(NA,6,22)

incmat[1,] <- c(0.091,0.16,0.19,0.18,0.15,0.1,0.061,0.032,0.015,0.0061,0.0022,7.2e-4,2.1e-4,5.3e-5)
incmat[2,] <- c(0.058,0.11,0.14,0.16,0.15,0.13,0.1,0.068,0.044,0.026,0.014,0.0072,0.0034,0.0015)
incmat[3,] <- c(0.043,0.079,0.11,0.12,0.13,0.12,0.11,0.088,0.07,0.052,0.037,0.025,0.016,0.0098)
incmat[4,] <- c(0.04,0.065,0.082,0.093,0.098,0.098,0.095,0.088,0.080,0.071,0.061,0.052,0.043,0.035)
incmat <- incmat/rowSums(incmat)
infmat[1,] <- c(rep(1.0,8),0.8,0.6,0.4,0.2,0.1,rep(0,9))
infmat[2,] <- c(rep(1.0,10),0.8,0.8,0.6,0.4,0.2,0.1,rep(0,6))
infmat[3,] <- c(rep(1.0,12),0.8,0.8,0.6,0.6,0.4,0.4,0.3,0.3,0.1,0.1)
infmat[4,] <- c(rep(1.0,6),0.9,0.8,0.7,0.6,0.4,0.2,0.1,rep(0,9))
infmat[5,] <- c(rep(1.0,7),0.9,0.8,0.8,0.7,0.7,0.6,0.4,0.2,0.1,rep(0,6))
infmat[6,] <- c(rep(1.0,8),0.9,0.9,0.8,0.8,0.7,0.7,0.6,0.6,0.4,0.4,0.3,0.3,0.1,0.1)

out1 <- matrix(NA,25,12)

loadid <- c("/Users/timtsang/SPH Dropbox/Tim Tsang/COVID_shandong/program_rcpp/2021_03_03/v4/model21/0.0277723572216928.Rdata",
            "/Users/timtsang/SPH Dropbox/Tim Tsang/COVID_shandong/program_rcpp/2021_03_03/v4/model23/0.0520897761452943.Rdata",
            "/Users/timtsang/SPH Dropbox/Tim Tsang/COVID_shandong/program_rcpp/2021_03_03/v4/model41/0.390996989794075.Rdata",
            "/Users/timtsang/SPH Dropbox/Tim Tsang/COVID_shandong/program_rcpp/2021_03_03/v4/model43/0.265370195731521.Rdata")

for (i in 1:4){
temp <- new.env()
load(file=loadid[i],temp)
zzz <- temp$z1
yyy <- temp$tt[[1]][5001:15000,]
out1[7,3*(i-1)+1:3] <- as.matrix(zzz[9,1:3])
out1[9,3*(i-1)+1:3] <- as.matrix(zzz[17,1:3])
out1[25,3*(i-1)+1:3] <- as.matrix(zzz[8,1:3])

out1[c(17,19),3*(i-1)+1:3] <- as.matrix(zzz[15:16,1:3])


## change the ref group to the oldest one
agem <- cbind(yyy[,10]-yyy[,12],-yyy[,12],yyy[,11]-yyy[,12])
tttt <- para_summary(agem,4,3,0)
out1[2,3*(i-1)+1:3] <- as.matrix(tttt[1,1:3])
out1[3,3*(i-1)+1:3] <- as.matrix(tttt[2,1:3])
out1[4,3*(i-1)+1:3] <- as.matrix(tttt[3,1:3])

out1[22,3*(i-1)+1:3] <- as.matrix(zzz[13,1:3])
out1[23,3*(i-1)+1:3] <- as.matrix(zzz[14,1:3])

for (j in 0:3){
tempvec <- log(yyy[,3+j]/yyy[,2])
out1[12+j,3*(i-1)+1:3] <- quantile(tempvec,c(0.5,0.025,0.975))
}

}



out1 <- exp(out1)
out1 <- round(out1,2)

tableout <- matrix(NA,nrow(out1),4)
for (j in 1:4){
tableout[,j] <- paste(out1[,3*(j-1)+1]," (",out1[,3*(j-1)+2],", ",out1[,3*(j-1)+3],")",sep="")  
}

tableout[grepl("NA",tableout)] <- " "
tableout[c(5,11,18,21),] <- "Ref" 

write.csv(tableout,"/Users/timtsang/Dropbox/COVID_shandong/summary/2021_03_03_table3.csv",row.names = F)


testing <- matrix(NA,5,2)
testing[,1] <- c(2,13,60,12,2)
testing[,2] <- c(10,19,68,4,1)

testing2 <- matrix(c(14,89-14,5,102-5),nrow=2)