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

out1 <- matrix(NA,14,12)

loadid <- c("/Users/timtsang/SPH Dropbox/Tim Tsang/COVID_shandong/program_rcpp/2021_03_03/v4_null/model21/0.22613280499354.Rdata",
            "/Users/timtsang/SPH Dropbox/Tim Tsang/COVID_shandong/program_rcpp/2021_03_03/v4_null/model23/0.500689913751557.Rdata",
            "/Users/timtsang/SPH Dropbox/Tim Tsang/COVID_shandong/program_rcpp/2021_03_03/v4_null/model41/0.876324594719335.Rdata",
            "/Users/timtsang/SPH Dropbox/Tim Tsang/COVID_shandong/program_rcpp/2021_03_03/v4_null/model43/0.0220032485667616.Rdata")

for (i in 1:4){
temp <- new.env()
load(file=loadid[i],temp)
zzz <- temp$z1
yyy <- temp$tt[[1]][5001:15000,]
out1[1,3*(i-1)+1:3] <- as.matrix(zzz[1,1:3])

for (j in 0:4){
out1[9+j,3*(i-1)+1:3] <- as.matrix(zzz[2+j,1:3])
tempvec <- yyy[,2+j]*exp(yyy[,8])
out1[3+j,3*(i-1)+1:3] <- quantile(tempvec,c(0.5,0.025,0.975))
}

out1[14,3*(i-1)+1:3] <- as.matrix(exp(zzz[8,1:3]))

}

out1[1,] <- out1[1,]*10^4
out1[2:13,] <- out1[2:13,]*10^2

out1 <- round(out1,2)

tableout <- matrix(NA,nrow(out1),4)
for (j in 1:4){
tableout[,j] <- paste(out1[,3*(j-1)+1]," (",out1[,3*(j-1)+2],", ",out1[,3*(j-1)+3],")",sep="")  
}

tableout[grepl("NA",tableout)] <- " "


write.csv(tableout,"/Users/timtsang/Dropbox/COVID_shandong/summary/2021_05_10_table2.csv",row.names = F)