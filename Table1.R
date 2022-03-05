## it is to compute the different version of data-based secondary attack rate

## primary, secondary and close contact

data1 <- read.csv("/Users/timtsang/SPH Dropbox/Tim Tsang/China COVID/COVID_shandong/data/2021_01_11_data_all.csv")

tableout <- matrix(NA,33,5)

#tableout[1,] <- c(" ","Primary","Secondary","Close contact","Data-based SAR")
tableout[,1] <- c("N","Age","Age groups"," 0-19"," 20-39"," 40-59"," 60+","Sex"," Female"," Male",
                  "Location"," Jinan"," Jining"," Qingdao","Contact type"," households"," hospital"," work"," flight"," Others",
                  "Severity of primary case"," Asymptomatic or mild"," Moderate"," Severe or critical",
                  "Severity"," Asymptomatic"," Mild"," Average"," Severe"," Critial","Medical-related personnel"," Yes"," No")

cond <- list(data1$primary==1,data1$primary==0&data1$pos==1,data1$primary==0)

out <- matrix(NA,33,12)


for (i in 1:3){
temp <- data1[cond[[i]],]  

out[1,3*i-2] <- nrow(temp)
out[2,3*i-2] <- median(temp$age,na.rm=T)
out[2,3*i-1] <- quantile(temp$age,0.25,na.rm = T)
out[2,3*i] <- quantile(temp$age,0.75,na.rm = T)
out[3,3*i-2] <- sum(temp$age<=19,na.rm=T)
out[4,3*i-2] <- sum(temp$age>19&temp$age<=39,na.rm=T)
out[5,3*i-2] <- sum(temp$age>39&temp$age<=59,na.rm=T)
out[6,3*i-2] <- sum(temp$age>59,na.rm=T)

out[8,3*i-2] <- sum(temp$male==0,na.rm=T)
out[9,3*i-2] <- sum(temp$male==1,na.rm=T)

out[11,3*i-2] <- sum(temp$location==1)
out[12,3*i-2] <- sum(temp$location==2)
out[13,3*i-2] <- sum(temp$location==3)

out[15,3*i-2] <- sum(temp$contact_type==1)
out[16,3*i-2] <- sum(temp$contact_type==2)
out[17,3*i-2] <- sum(temp$contact_type==3)
out[18,3*i-2] <- sum(temp$contact_type==4)
out[19,3*i-2] <- sum(temp$contact_type==5)

out[22,3*i-2] <- sum(temp$primary_severity==1)
out[23,3*i-2] <- sum(temp$primary_severity==2)
out[24,3*i-2] <- sum(temp$primary_severity==3)

out[26,3*i-2] <- sum(temp$severity==1)
out[27,3*i-2] <- sum(temp$severity==2)
out[28,3*i-2] <- sum(temp$severity==3)
out[29,3*i-2] <- sum(temp$severity==4)
out[30,3*i-2] <- sum(temp$severity==5)

out[32,3*i-2] <- sum(temp$medppl==1,na.rm = T)
out[33,3*i-2] <- sum(temp$medppl==0,na.rm=T)
}


out[15:19,1] <- NA

# add index contact

for (i in 1:5){
  primaryid <- data1$ID[data1$primary==1]
  temp <- data1[data1$contact_type==i,5:7]
  out[14+i,1] <-  sum(primaryid%in%temp[,1] | primaryid%in%temp[,2] | primaryid%in%temp[,3])
}


for (i in 1:3){
out[3:6,3*i-1]  <- sum(out[3:6,3*i-2] )
out[8:9,3*i-1]  <- sum(out[8:9,3*i-2] ) 
out[11:13,3*i-1]  <- sum(out[11:13,3*i-2] )
out[15:19,3*i-1]  <- sum(out[15:19,3*i-2] )
out[22:24,3*i-1]  <- sum(out[22:24,3*i-2] )
out[26:30,3*i-1]  <- sum(out[26:30,3*i-2] )
out[32:33,3*i-1]  <- sum(out[32:33,3*i-2] ) 
if (i == 1){
  out[15:19,2] <- 97
}
out[-2,3*i] <- round(100*(out[-2,3*i-2]/out[-2,3*i-1]))
}

for (i in 1:nrow(out)){
if (!is.na(out[i,7])&i!=2){  
out[i,10] <- out[i,4]/out[i,7]
temp2 <- binom.test(out[i,4],out[i,7])
out[i,11:12] <- temp2$conf.int
}
}
out[-2,10:12] <- round(out[-2,10:12]*100,2)



tableout[,2] <- paste(out[,1],"/",out[,2]," (",out[,3],"%)",sep="")
tableout[,3] <- paste(out[,4],"/",out[,5]," (",out[,6],"%)",sep="")
tableout[,4] <- paste(out[,7],"/",out[,8]," (",out[,9],"%)",sep="")
tableout[,5] <- paste(out[,10],"% (",out[,11],"-",out[,12],")",sep="")
tableout[grepl("NA",tableout)] <- " "
tableout[1,2:4] <- out[1,c(1,4,7)]
tableout[2,2] <- paste(out[2,1]," (",out[2,2],", ",out[2,3],")",sep="")
tableout[2,3] <- paste(out[2,4]," (",out[2,5],", ",out[2,6],")",sep="")
tableout[2,4] <- paste(out[2,7]," (",out[2,8],", ",out[2,9],")",sep="")

write.csv(tableout,"/Users/timtsang/SPH Dropbox/Tim Tsang/China COVID/COVID_shandong/summary/epidemics_R1/2022_02_21_table1_v2.csv",row.names = F)

# overall SAR
binom.test(102,3158)

# severity
fisher.test(matrix(c(2,10,89-2,102-10),nrow=2))


# age
fisher.test(matrix(c(29,32,1261-29,854-32),nrow=2))

fisher.test(matrix(c(29,23,1261-29,308-23),nrow=2))

fisher.test(matrix(c(29,18,1261-29,315-18),nrow=2))

# sex
fisher.test(matrix(c(42,60,1543-42,1407-60),nrow=2))

# location
fisher.test(matrix(c(28,34+40,540-28,1237+1381-40-34),nrow=2))

fisher.test(matrix(c(28,40,540-28,1381-40),nrow=2))


# contact type
fisher.test(matrix(c(68,14+8,674-68,504+244-8-14),nrow=2))


# contact type
fisher.test(matrix(c(68,14,674-68,504-14),nrow=2))

fisher.test(matrix(c(68,8,674-68,244-8),nrow=2))

# severity of primary cases
fisher.test(matrix(c(32,62,678-32,2199-62),nrow=2))

fisher.test(matrix(c(32,8,678-32,281-8),nrow=2))

# being a medical personal
fisher.test(matrix(c(5,69,333-5,2299-69),nrow=2))

