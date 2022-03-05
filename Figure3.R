####
data <- read.csv("/Users/timtsang/SPH Dropbox/Tim Tsang/COVID_shandong/summary/plotdata.csv")

## plot function

plotfunction <- function(table,pred,title,panel_index,legend){
  
  par(mar=c(0,0,0,0))
  
  plot(0,0,xlab="",ylab="", main="", axes=F, xlim=c(0,7), ylim=c(-1,15),type="n")
  
  axis(1,at=1:7,labels=0:6*2,pos = 1.5)
  axis(2,at=1:7*2,labels=1:7*2-2,pos = 0.7)
  
  lines(c(1,7),c(2,14),lty=2,lwd=0.8)
  lines(c(1,7),c(2,14)+1,lty=3,lwd=0.8)
  lines(c(1,7)+0.5,c(2,14),lty=3,lwd=0.8)
  
  mtext("Observed number of secondary cases",side = 1,line = -2.5,cex = 0.9,font = 1.1 , at=4)
  mtext("Model-predicted number of secondary cases",side = 2,line = -1.6,cex = 0.9,font = 1.1, at = 8)
  
  colbox <- c("red","blue","orange","black","purple")
  if (pred==1) {
    for (i in 1:5) {
      points((table[table$index_type==i,2]/2)+1,table[table$index_type==i,3]+2,col=colbox[i],pch=16)}
  }
  if (pred==2){
    for (i in 1:5) {
      points((table[table$index_type==i,2]/2)+1,table[table$index_type==i,4]+2,col=colbox[i],pch=16)}
  }
  if (pred==3){
    for (i in 1:5) {
      points((table[table$index_type==i,2]/2)+1,table[table$index_type==i,5]+2,col=colbox[i],pch=16)
    }
  }
  if (pred==4){
    for (i in 1:5) {
      points((table[table$index_type==i,2]/2)+1,table[table$index_type==i,6]+2,col=colbox[i],pch=16)
    }
  }
  
  if (legend==1){
    legend(1,14.3,c("Household", "Healthcare Facility", "Workplace", "Air Transportation","Other Settings"),pch=16,bty="n",col=colbox)
  }
  mtext(title,side = 3,line = -1.2,font = 1.7)
  text(0.2,15.2,labels = panel_index,font = 1.7,cex = 1.6)
  
}

###plot
pdf("/Users/timtsang/SPH Dropbox/Tim Tsang/COVID_shandong/2021_09_11_scatter plot.pdf",width=9.5, height=9)

layout(matrix(1:4, nrow=2,byrow=T))

plotfunction(data,pred = 1,'Incubation period: 5 days, Infectious period: 13 days','A',1)
plotfunction(data,pred = 2,'Incubation period: 5 days, Infectious period: 21 days','B',0)
plotfunction(data,pred = 3,'Incubation period: 7 days, Infectious period: 13 days','C',0)
plotfunction(data,pred = 4,'Incubation period: 7 days, Infectious period: 21 days','D',0)

dev.off()

