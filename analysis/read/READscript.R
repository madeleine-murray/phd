########################READ###############################

args = commandArgs(trailingOnly=TRUE)

mydata <- read.table("READ_output_ordered", sep="\t", header = T)

means_2AlleleDifference_percentage<-aggregate(P0~PairIndividuals, mydata, mean)


if (length(args)<1){
	print("No standardization method specified!")
} else if (args[1] =='median'){
	normalization_val = median(means_2AlleleDifference_percentage$P0, na.rm = TRUE)
} else if (args[1] =='mean'){
	normalization_val = mean(means_2AlleleDifference_percentage$P0, na.rm = TRUE)
} else if (args[1] =='max'){
	normalization_val = max(means_2AlleleDifference_percentage$P0, na.rm = TRUE)
} else if (args[1] == 'value'){
	if (length(args)<2){
	print("Value standardization selected, but no valid standardization value specified!")}
	normalization_val = as.numeric(args[2])
}else{
	print("No valid standardization method specified!")}


mydata[,"Normalized2AlleleDifference"] <- mydata[,"P0"]/normalization_val

means_2AlleleDifferenceNormalized<-aggregate(Normalized2AlleleDifference~PairIndividuals, mydata, mean)
serr_2AlleleDifferenceNormalized<-aggregate(Normalized2AlleleDifference~PairIndividuals, mydata, function(x){return(sd(x)/sqrt(length(x)))})
unnorm_serr_2AlleleDifferenceNormalized<-aggregate(P0~PairIndividuals, mydata, function(x){return(sd(x)/sqrt(length(x)))})
unnorm_2AlleleDifferenceNormalized<-aggregate(P0~PairIndividuals, mydata, mean)
means_2AlleleDifferenceNormalized$StandardError<-serr_2AlleleDifferenceNormalized$Normalized2AlleleDifference
means_2AlleleDifferenceNormalized$NonNormalizedP0<-unnorm_2AlleleDifferenceNormalized$P0
means_2AlleleDifferenceNormalized$NonNormalizedStandardError<-unnorm_serr_2AlleleDifferenceNormalized$P0




if ((length(means_2AlleleDifference_percentage$P0) %% 2)==1){
	index=which(means_2AlleleDifference_percentage$P0==normalization_val)[1]
	seb2_nn=means_2AlleleDifferenceNormalized$NonNormalizedStandardError[index]
}else{
	indices=which(abs(means_2AlleleDifference_percentage$P0-normalization_val)== min(abs(means_2AlleleDifference_percentage$P0-normalization_val),na.rm=TRUE))
	sorted_nn=means_2AlleleDifferenceNormalized$NonNormalizedStandardError[order(means_2AlleleDifference_percentage$P0)]
	se1_nn=sorted_nn[length(means_2AlleleDifference_percentage$P0)/2]
	se2_nn=sorted_nn[1+length(means_2AlleleDifference_percentage$P0)/2]
	seb2_nn=mean(c(se1_nn,se2_nn))
}




write.table(means_2AlleleDifferenceNormalized,"meansP0_AncientDNA_normalized",row.names=F,quote=FALSE)

tab2=means_2AlleleDifferenceNormalized[order(means_2AlleleDifferenceNormalized$NonNormalizedP0),]

error.bars <- function(y,z){
	n=length(y)
	x=1:n
	for (i in 1:n){
		arrows(x[i],y[i]-z[i],x[i],y[i]+z[i],code=3,angle=90,length=0)
	}
	points(x,y,pch=16)
}

add.alpha <- function(col, alpha=1){
  apply(sapply(col, col2rgb)/255, 2, 
                     function(x) 
                       rgb(x[1], x[2], x[3], alpha=alpha))  
}




pdf('READ_results_plot.pdf',height=6,width=8+0.1*nrow(tab2))

layout(mat=t(1:2),widths=c(3,1.5))
par(mar=c(10, 4, 2, 2))

plot(1,type='n',xlim=c(0,nrow(tab2)+1),ylim=c(0.8*min(means_2AlleleDifferenceNormalized$NonNormalizedP0),1.1*max(means_2AlleleDifferenceNormalized$NonNormalizedP0)),axes=F,xlab='',ylab='Average pairwise P0 (\u00B1 2SE)')
axis(2)

ci_factor=1.96

polygon(c(-100,nrow(tab2)+100,nrow(tab2)+100,-100),c((normalization_val-ci_factor*seb2_nn)*0.90625,(normalization_val-ci_factor*seb2_nn)*0.90625,(normalization_val+ci_factor*seb2_nn)*0.90625,(normalization_val+ci_factor*seb2_nn)*0.90625),col=add.alpha('gray70',0.5),border='NA')

polygon(c(-100,nrow(tab2)+100,nrow(tab2)+100,-100),c((normalization_val-ci_factor*seb2_nn)*0.8125,(normalization_val-ci_factor*seb2_nn)*0.8125,(normalization_val+ci_factor*seb2_nn)*0.8125,(normalization_val+ci_factor*seb2_nn)*0.8125),col=add.alpha('gray70',0.5),border='NA')

polygon(c(-100,nrow(tab2)+100,nrow(tab2)+100,-100),c((normalization_val-ci_factor*seb2_nn)*0.625,(normalization_val-ci_factor*seb2_nn)*0.625,(normalization_val+ci_factor*seb2_nn)*0.625,(normalization_val+ci_factor*seb2_nn)*0.625),col=add.alpha('gray70',0.5),border='NA')

error.bars(tab2$NonNormalizedP0,tab2$NonNormalizedStandardError*2)

abline(h=normalization_val)
abline(h=normalization_val*0.90625,lty=2)
text(2,normalization_val*0.89,labels='2nd degree',col='blue')
abline(h=normalization_val*0.8125,lty=2)
text(2,normalization_val*0.8,labels='1st degree',col='blue')
abline(h=normalization_val*0.625,lty=2)
text(2,normalization_val*0.61,labels='Identical/Twins',col='blue')

par(las=2)
axis(1,at=1:nrow(tab2),labels=tab2$PairIndividuals,cex.axis=0.7)

par(las=0)
hist(tab2$NonNormalizedP0,xlab='Average pairwise P0',main='',breaks=30,col='blue')
abline(v=normalization_val)
abline(v=normalization_val*0.90625,lty=2)
text(normalization_val*0.885,0,labels='2nd degree',srt = 90,pos=4,col='gray40')
abline(v=normalization_val*0.8125,lty=2)
text(normalization_val*0.795,0,labels='1st degree',srt = 90,pos=4,col='gray40')
abline(v=normalization_val*0.625,lty=2)
text(normalization_val*0.605,0,labels='Identical/Twins',srt = 90,pos=4,col='gray40')

dev.off()
