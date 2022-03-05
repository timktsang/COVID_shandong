
rm(list = ls())

library(Rcpp)
library(RcppParallel)
library(mvtnorm)
print(defaultNumThreads())
#setThreadOptions(numThreads = 8)
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

########
## incubation and infectiousness

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

incvec <- incmat[2,]
infvec <- infmat[3,]

########
# set link
#setwd("/Users/timtsang/Dropbox/COVID_shandong/program_rcpp/2021_03_03/v4")
#setwd("Z:/COVID_shandong/program_rcpp/2021_03_03/v4_null/model23")
##
data <- read.csv("2021_01_11_data.csv",as.is = T)

#primarylist <- data[data[,3]==1,1]
#data <- data[data[,5]%in%primarylist|data[,6]%in%primarylist|data[,7]%in%primarylist|data[,3]==1,]

for (i in 1:ncol(data)){
data[,i] <- as.numeric(data[,i])
}
data$agegp <- as.numeric(cut(data$age,c(0,9,19,29,39,49,59,69,999)))-1

## 58 and 84 primary with no onset date
data$bound[data[,1]==4043] <- 55
data$bound[data[,1]==4044] <- 52
# 4171 wrong
data$contact_start_day[data[,1]==4171] <- NA
data$contact_end_day[data[,1]==6358] <- NA

## Shandong started quarantine on Feb 3. Before that, people without symptoms were quarantined at home.
## ref date (11/30/2019)
## feb3 = 65

data$upperbound_inf <- pmin(data$contact_end_day,data$quarantine_start_day,data$bound,data$onset_day-1,na.rm=T)
data$lowerbound_inf <- data$contact_start_day

# for inf time
# for primary, not other info, just assume at most 14 days incubation
data[data$primary==1,]$lowerbound_inf <- data[data$primary==1,]$upperbound_inf - 14
# for non-primary, conditional on infection must be later than primary
#data[is.na(data$lowerbound_inf)&data$pos==1&data$primary==0,]
for (i in 1:nrow(data)){
if (is.na(data$lowerbound_inf[i])&data$pos[i]==1&data$primary[i]==0){
temp <- data[data[,2]==data[i,2]&data[,3]==1,]
data$lowerbound_inf[i] <- min(temp$onset_day)-5
}
if (data$pos[i]==1){  
if (data$lowerbound_inf[i]>data$upperbound_inf[i]){
  data$lowerbound_inf[i] <-  data$upperbound_inf[i] - 14 
} 
}
}
data[which(data$lowerbound_inf>data$upperbound_inf),]

data$upperbound_exposure <- pmin(data$contact_end_day,data$quarantine_start_day,na.rm=T)
data$lowerbound_exposure <- data$contact_start_day

## lowerbound of exposure
data[is.na(data$lowerbound_exposure)&data$primary==0,]
for (i in 1:nrow(data)){
  if (is.na(data$lowerbound_exposure[i])&data$primary[i]==0){
    temp <- data[data[,2]==data[i,2]&data[,3]==1,]
    data$lowerbound_exposure[i] <- min(temp$onset_day,na.rm=T)-5
    }
}

data[is.na(data$upperbound_exposure)&data$primary==0,]
for (i in 1:nrow(data)){
  if (is.na(data$upperbound_exposure[i])&data$primary[i]==0){
    temp <- data[data[,2]==data[i,2],]
    data$upperbound_exposure[i] <- max(temp$upperbound_exposure,na.rm=T)
  }
}

testset <- which(data$upperbound_exposure<data$lowerbound_exposure)
length(testset)
## the contact end date may be wrong
data$upperbound_exposure[testset] <- data$quarantine_start_day[testset]

testset <- which(data$upperbound_exposure<data$lowerbound_exposure)
length(testset)
## missing upperbound again
testset <- which(is.na(data$upperbound_exposure)&data$primary==0)
data$upperbound_exposure[testset] <- 120

## for those inf=0, bound_inf = bound_exposure 
data[data$pos==0,c("lowerbound_inf","upperbound_inf")] <- data[data$pos==0,c("lowerbound_exposure","upperbound_exposure")]

# exposure is about contact
# use quarantine start day to determine the expsoure from community
for (i in 1:length(unique(data$clusterID))){
temp <- data[data$clusterID==unique(data$clusterID)[i],]  
data$quarantine_start_day[data$clusterID==unique(data$clusterID)[i]&is.na(data$quarantine_start_day)] <- max(data[data$clusterID==unique(data$clusterID)[i],c("quarantine_start_day","contact_end_day")],na.rm=T)
}
data$quarantine_start_day <- pmax(data$quarantine_start_day,data$upperbound_exposure,na.rm=T)

data$inf_time <- -1
## set the data set to following format
data1 <- data[,c("ID","primary","clusterID","contact_type","contact1","contact2","contact3","location","male","agegp"
                 ,"pos","onset_day","lowerbound_inf","upperbound_inf","inf_time","lowerbound_exposure","upperbound_exposure","quarantine_start_day","primary_severity","medppl")]

summary(data1[data1$primary==0,])
data1[is.na(data1$upperbound_exposure)&data1$primary==0,]


data1$inf_time[data1$upperbound_inf==data1$lowerbound_inf&data1$pos==1] <- data1$lowerbound_inf[data1$upperbound_inf==data1$lowerbound_inf&data1$pos==1]



data1[is.na(data1)] <- -1

data1 <- as.matrix(data1)


# a contact b means (a,b)=(b,a)=1
contactmatrix <- matrix(0,nrow(data1),nrow(data1))
for (i in 1:nrow(data1)){

    if (data1[i,2]==0){
    for (v in 1:3){
     contactmatrix[i,which(data1[,1]==data1[i,4+v])] <- 1
    }
    # get all the same household
    if (data1[i,4]==1){
    for (v in 1:3){  
      for (w in 1:3){
      contactmatrix[i,which(data1[,4]==1&data1[,4+w]==data1[i,4+v]&data1[,4+w]!=-1&data1[,4+v]!=-1)] <- 1
      }
      }  
    }
    }
}

contactmatrix2 <- matrix(-1,nrow(data1),max(rowSums(contactmatrix)))
for (i in 1:nrow(data1)){
vv <- which(contactmatrix[i,]==1)  
vv <- vv[vv!=i]
if (length(vv)>0){
contactmatrix2[i,1:length(vv)] <- vv
}
}

# para 1. from community
# para 2. from household
# para 3. from hosp 
# para 4. from work 
# para 5. from flight 
# para 6. other contacts
# para 7. asym vs sym
# para 8. pre-sym vs sym
# para 9. sex effect
# para 10-12. age effect 1-3
# para 13-14. location effect
# para 15-16. primary case severity
# para 17. med person or not

## edit cpp to test it


# para2: 1-6 sex dist for 3 location, 7-30: age groups for 3 location, 31-32: medppl
para <- c(0.00001,rep(0.1,5),rep(0,11),2)
para2 <- c(rep(0.5,6),rep(1/8,24),rep(0.5,2))

#para <- read.csv("para1.csv")[,1]
#para2 <- read.csv("para2.csv")[,1]
para[1:6] <- 0.03
# normalized para2
para2[1:2] <- para2[1:2]/sum(para2[1:2])
para2[3:4] <- para2[3:4]/sum(para2[3:4])
para2[5:6] <- para2[5:6]/sum(para2[5:6])
para2[7:14] <- para2[7:14]/sum(para2[7:14])
para2[15:22] <- para2[15:22]/sum(para2[15:22])
para2[23:30] <- para2[23:30]/sum(para2[23:30])
para2[31:32] <- para2[31:32]/sum(para2[31:32])
para[c(7,9:17)] <- 0
sourceCpp("shandong.cpp")
sim <- sim_data(data1,contactmatrix2,matrix(0,nrow(data1),1),para,para2,infvec,incvec,c(0,0))


testdata <- sim[[1]]
summary(testdata)
table(testdata[,11]==1,testdata[,2])
table(testdata[testdata[,10]==-1,8])
sim[[4]][testdata[,10]==-1,] 

#testing <- loglik2(testdata,contactmatrix2,para,para2,infvec,incvec)



fisher.test(table(testdata[,9],testdata[,11]==1))
para[9]

hazard <- sim[[2]]
table(testdata[testdata[,11]==1,12])
table(testdata[testdata[,11]==1,15])

sigma <- abs(para+0.1)/5
move <- rep(1,length(para))
move[c(7)] <- 0
aaaaa1 <- Sys.time()
tt <- mcmc(data1,contactmatrix2,para,para2,infvec,incvec,15000,move,sigma)
aaaaa2 <- Sys.time()
print(aaaaa2-aaaaa1)


testing <- loglik2(tt[[4]],contactmatrix2,as.matrix(tt[[6]][1,]),para,para2,infvec,incvec)




id <- runif(1,0,1)
inc <- 5000+1:10000
z1 <- para_summary((tt[[1]][inc,]),4,4,1)
print(z1)
z2 <- para_summary((tt[[2]][inc,]),4,4,1)
print(z2)

write.csv(z1,paste("mcmc_summary1_",id,".csv"),row.names = F)
write.csv(z2,paste("mcmc_summary2_",id,".csv"),row.names = F)



recenter <- function(data1){
  tempdata <- data1
  cutid <- tempdata[tempdata[,2]==1&tempdata[,12]==-1,3]
  tempdata <- tempdata[!tempdata[,3]%in%cutid,]
  tempdata[tempdata[,12]==-1,12] <- NA
  for (i in 1:length(unique(tempdata[,3]))){
    tempdata[tempdata[,3]==unique(tempdata[,3])[i],12] <- tempdata[tempdata[,3]==unique(tempdata[,3])[i],12] - min(tempdata[tempdata[,3]==unique(tempdata[,3])[i]&tempdata[,2]==1,12],na.rm=T)
  }
  return(tempdata)
}

primarylist <- data[data[,3]==1,1]
data1pri <- data1[data[,5]%in%primarylist|data[,6]%in%primarylist|data[,7]%in%primarylist|data[,3]==1,]

validate <- matrix(0,10000,40)
validate2 <- matrix(0,10000,40)
for (i in 1:10000){
  curpara1 <- tt[[1]][5000+i,]
  curpara2 <- tt[[2]][5000+i,]
  rm <- as.matrix(tt[[6]][5000+i,])
  sim2 <- sim_data(data1,contactmatrix2,rm,curpara1,curpara2,infvec,incvec,c(1,1))
  sim22 <- recenter(sim2[[1]])
  sim22pri <- recenter(sim2[[1]][data[,5]%in%primarylist|data[,6]%in%primarylist|data[,7]%in%primarylist|data[,3]==1,])
  
  for (k in -10:29){
    validate[i,k+11] <- sum(sim22[sim22[,2]==0,12]==k,na.rm = T)  
    validate2[i,k+11] <- sum(sim22pri[sim22pri[,2]==0,12]==k,na.rm = T)  
  }

  print(i)
}

data1re <- recenter(data1)
data1prire <- recenter(data1pri)
obs <- matrix(0,5,40)
obs2 <- matrix(0,5,40)
for (i in 1:40){
  obs[1,i] <- sum(data1re[data1re[,2]==0,12]==i-11,na.rm = T)  
  obs[2,i] <- mean(validate[,i])
  obs[3,i] <- quantile(validate[,i],0.025)
  obs[4,i] <- quantile(validate[,i],0.975)
  obs[5,i] <- 1*(  obs[1,i]>=  obs[3,i]&&obs[1,i]<=  obs[4,i])
  
  obs2[1,i] <- sum(data1prire[data1prire[,2]==0,12]==i-11,na.rm = T)  
  obs2[2,i] <- mean(validate2[,i])
  obs2[3,i] <- quantile(validate2[,i],0.025)
  obs2[4,i] <- quantile(validate2[,i],0.975)
  obs2[5,i] <- 1*(  obs2[1,i]>=  obs2[3,i]&&obs2[1,i]<=  obs2[4,i])
}

write.csv(obs,paste("pred_",id,".csv"),row.names = F)
write.csv(obs2,paste("pred2_",id,".csv"),row.names = F)

save.image(file=paste(id,".Rdata",sep=""))

colMeans(tt[[6]]) -> kk
kk <- kk[kk!=0]



data$nsec <- 0
for (i in 1:length(unique(data$clusterID))){
data$nsec[data$clusterID==unique(data$clusterID)[i]]  <- sum(data$pos[data$clusterID==unique(data$clusterID)[i]&data$primary==0])
}

nsec <- data$nsec[data$primary==1]
