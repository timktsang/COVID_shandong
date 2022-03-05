## close contacts distribution

data1 <- read.csv("/Users/timtsang/Dropbox/COVID_shandong/data/2021_01_11_data_all.csv")

temp <- data1[data1$primary==1,]
temp <- temp[!duplicated(temp$clusterID),]
quantile(temp $cluster_size,c(0,0.25,0.5,0.75,1))

data11 <- data1[data1$pos==1,]

data11[,c("contact_all","contact_type1","contact_type2","contact_type3","contact_type4","contact_type5")] <- NA

for (i in 1:nrow(data11)){
data11$contact_all[i] <- sum(data1[,c("contact1","contact2","contact3")]==data11$ID[i],na.rm=T)
data11$contact_type1[i] <- sum(data1[,c("contact1","contact2","contact3")]==data11$ID[i]&data1$contact_type==1,na.rm=T)
data11$contact_type2[i] <- sum(data1[,c("contact1","contact2","contact3")]==data11$ID[i]&data1$contact_type==2,na.rm=T)
data11$contact_type3[i] <- sum(data1[,c("contact1","contact2","contact3")]==data11$ID[i]&data1$contact_type==3,na.rm=T)
data11$contact_type4[i] <- sum(data1[,c("contact1","contact2","contact3")]==data11$ID[i]&data1$contact_type==4,na.rm=T)
data11$contact_type5[i] <- sum(data1[,c("contact1","contact2","contact3")]==data11$ID[i]&data1$contact_type==5,na.rm=T)
}

mean(data11$contact_type1)
mean(data11$contact_type2)
mean(data11$contact_type3)
mean(data11$contact_type4)
mean(data11$contact_type5)

AIC_matrix <- matrix(NA,6,3)
para_m <- matrix(NA,6,3)
summarystat <- matrix(NA,6,6)
## confirm nb dist is the best
library(fitdistrplus)
library(loo)
out_t <- matrix(NA,6,13)



for (i in 1:6){
  a <- fitdist(data11[,21+i],"geom")
  out_t[i,1:3] <- a$estimate + c(0,-1.96,1.96)*a$sd
  out_t[i,4] <- a$aic
  out_t[i,5] <- a$bic
  a <- fitdist(data11[,21+i],"nbinom")
  out_t[i,6:8] <- a$estimate[1] + c(0,-1.96,1.96)*a$sd[1]
  out_t[i,9:11] <- a$estimate[2] + c(0,-1.96,1.96)*a$sd[2]
  out_t[i,12] <- a$aic
  out_t[i,13] <- a$bic
}
out_t <- pmax(round(out_t,3),0.001)

out_t2 <- matrix(NA,6,7)
out_t2[,1] <- paste0(out_t[,1]," (",out_t[,2],", ",out_t[,3],")")
out_t2[,2:3] <- out_t[,4:5]
out_t2[,4] <- paste0(out_t[,6]," (",out_t[,7],", ",out_t[,8],")")
out_t2[,5] <- paste0(out_t[,9]," (",out_t[,10],", ",out_t[,11],")")
out_t2[,6:7] <- out_t[,12:13]

write.csv(out_t2,"/Users/timtsang/SPH Dropbox/Tim Tsang/China COVID/COVID_shandong/summary/epidemics_R1/nb_fit.csv",row.names = F)

for (i in 1:6){
  a <- fitdist(data11[,21+i],"pois")
  AIC_matrix[i,1] <- summary(a)$aic
  a <- fitdist(data11[,21+i],"nbinom")
  AIC_matrix[i,2] <- summary(a)$aic
  para_m[i,1:2] <- a$estimate
  a <- fitdist(data11[,21+i],"geom")
  AIC_matrix[i,3] <- summary(a)$aic
  para_m[i,3] <- a$estimate
  get <- data11[,21+i]
  summarystat[i,] <- c(mean(get),quantile(get,c(0,0.25,0.5,0.75,1)))
}

library(gridExtra)
library(grid)
library(ggplot2)

nb <- data.frame(matrix(NA,151,7))
nb[,1] <- 0:150
for (i in 1:6){
nb[,1+i] <- dnbinom(0:150,size=para_m[i,1],mu=para_m[i,2])*191
}

geo <- data.frame(matrix(NA,151,7))
geo[,1] <- 0:150
for (i in 1:6){
  geo[,1+i] <- dgeom(0:150,para_m[i,3])*191
}

p1 <- ggplot(data11,aes(x=contact_all))
p1 <- p1 + theme(
  panel.background = element_blank(),
  axis.line = element_line(colour = "black"),
  legend.position = c(.95, .95),
  legend.justification = c("right", "top"),
  legend.box.just = "right",
  legend.margin = margin(6, 6, 6, 6)
)

p1 <- p1 + geom_histogram(binwidth=1) 
p1 <- p1 + geom_line(data=geo,aes(x=X1,y=X2,colour="Geometric (AIC=1557)")) 
p1 <- p1 + geom_line(data=nb,aes(x=X1,y=X2,colour="Negative binomial (AIC=1430), k=0.34"))  +labs(colour=" ",y= "Count", x = "All Settings") 
p1 <- p1 + scale_color_manual(values=c("purple", "red"))
#p1 <- p1  + annotate(geom="text", x=70, y=30, label="AIC")
#p1 <- p1  + annotate(geom="text", x=70, y=20, label="Geometric = ")
#p1 <- p1  + annotate(geom="text", x=70, y=10, label="Negative binomial = ")


p2 <- ggplot(data11,aes(x=contact_type1))
p2 <- p2 + theme(
  panel.background = element_blank(),
  axis.line = element_line(colour = "black"),
  legend.position = c(.95, .95),
  legend.justification = c("right", "top"),
  legend.box.just = "right",
  legend.margin = margin(6, 6, 6, 6)
)
p2 <- p2 + geom_histogram(binwidth=1) 
p2 <- p2 + geom_line(data=geo[1:32,],aes(x=X1,y=X3,colour="Geometric (AIC=950)")) 
p2 <- p2 + geom_line(data=nb[1:32,],aes(x=X1,y=X3,colour="Negative binomial (AIC=885), k=0.34"))  +labs(colour=" ",y= "Count", x = "Household")
p2 <- p2 + scale_color_manual(values=c("purple", "red"))
#p2 <- p2 + geom_line(data=geo[1:32,],aes(x=X1,y=X3),colour="purple") 
#p2 <- p2 + geom_line(data=nb[1:32,],aes(x=X1,y=X3),colour="red")  +labs(y= "Count", x = "Contacts in households") 
 

p3 <- ggplot(data11,aes(x=contact_type2))
p3 <- p3 + theme(
  panel.background = element_blank(),
  axis.line = element_line(colour = "black"),
  legend.position = c(.95, .95),
  legend.justification = c("right", "top"),
  legend.box.just = "right",
  legend.margin = margin(6, 6, 6, 6)
)
p3 <- p3 + geom_histogram(binwidth=1) 
p3 <- p3 + geom_line(data=geo[1:92,],aes(x=X1,y=X4,colour="Geometric (AIC=849)")) 
p3 <- p3 + geom_line(data=nb[1:92,],aes(x=X1,y=X4,colour="Negative binomial (AIC=523), k=0.08"))  +labs(colour=" ",y= "Count", x = "Healthcare Facility")
p3 <- p3 + scale_color_manual(values=c("purple", "red"))

#p3 <- p3 + geom_line(data=geo[1:92,],aes(x=X1,y=X4),colour="purple") 
#p3 <- p3 + geom_line(data=nb[1:92,],aes(x=X1,y=X4),colour="red")  +labs(y= "Count", x = "Contacts in medical-related facilities") 

p4 <- ggplot(data11,aes(x=contact_type3))
p4 <- p4 + theme(
  panel.background = element_blank(),
  axis.line = element_line(colour = "black"),
  legend.position = c(.95, .95),
  legend.justification = c("right", "top"),
  legend.box.just = "right",
  legend.margin = margin(6, 6, 6, 6)
)
p4 <- p4 + geom_histogram(binwidth=1) 
p4 <- p4 + geom_line(data=geo[1:54,],aes(x=X1,y=X5,colour="Geometric (AIC=632)")) 
p4 <- p4 + geom_line(data=nb[1:54,],aes(x=X1,y=X5,colour="Negative binomial (AIC=358), k=0.05"))  +labs(colour=" ",y= "Count", x = "Workplaces")
p4 <- p4 + scale_color_manual(values=c("purple", "red"))

#p4 <- p4 + geom_line(data=geo[1:54,],aes(x=X1,y=X5),colour="purple") 
#p4 <- p4 + geom_line(data=nb[1:54,],aes(x=X1,y=X5),colour="red")  +labs(y= "Count", x = "Contacts in workplaces") 

p5 <- ggplot(data11,aes(x=contact_type4))
p5 <- p5 + theme(
  panel.background = element_blank(),
  axis.line = element_line(colour = "black"),
  legend.position = c(.95, .95),
  legend.justification = c("right", "top"),
  legend.box.just = "right",
  legend.margin = margin(6, 6, 6, 6)
)
p5 <- p5 + geom_histogram(binwidth=1) 
p5 <- p5 + geom_line(data=geo[1:74,],aes(x=X1,y=X6,colour="Geometric (AIC=862)")) 
p5 <- p5 + geom_line(data=nb[1:74,],aes(x=X1,y=X6,colour="Negative binomial (AIC=255), k=0.01"))  +labs(colour=" ",y= "Count", x = "Air transportation")
p5 <- p5 + scale_color_manual(values=c("purple", "red"))

#p5 <- p5 + geom_line(data=geo[1:74,],aes(x=X1,y=X6),colour="purple") 
#p5 <- p5 + geom_line(data=nb[1:74,],aes(x=X1,y=X6),colour="red")  +labs(y= "Count", x = "Contacts in air transportation") 

p6 <- ggplot(data11,aes(x=contact_type5))
p6 <- p6 + theme(
  panel.background = element_blank(),
  axis.line = element_line(colour = "black"),
  legend.position = c(.95, .95),
  legend.justification = c("right", "top"),
  legend.box.just = "right",
  legend.margin = margin(6, 6, 6, 6)
)
p6 <- p6 + geom_histogram(binwidth=1) 
p6 <- p6 + geom_line(data=geo[1:141,],aes(x=X1,y=X6,colour="Geometric (AIC=1236)")) 
p6 <- p6 + geom_line(data=nb[1:141,],aes(x=X1,y=X6,colour="Negative binomial (AIC=944), k=0.16"))  +labs(colour=" ",y= "Count", x = "Other Settings")
p6 <- p6 + scale_color_manual(values=c("purple", "red"))

#p6 <- p6 + geom_line(data=geo[1:141,],aes(x=X1,y=X7),colour="purple") 
#p6 <- p6 + geom_line(data=nb[1:141,],aes(x=X1,y=X7),colour="red")  +labs(y= "Count", x = "Other types of contacts") 


pdf("/Users/timtsang/SPH Dropbox/Tim Tsang/China COVID/COVID_shandong/summary/epidemics_R1/figureS1.pdf",width=11, height=7)
grid.arrange(p1,p2,p3,p4,p5,p6,nrow=2)
dev.off()
