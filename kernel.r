library(distr)
library(distrEx)
library(mixtools)
library(MASS)
library(mclust)

gaussianMetric <- function(dataFile1, dataFile2, index){
	tab1 = read.csv(dataFile1)
	tab2= read.csv(dataFile2)
	col1 = as.vector(tab1[[index]])
	col2 = as.vector(tab2[[index]])

	mclustModel = Mclust(c(col1,col2),G=2)
        model = normalmixEM(c(col1,col2),k=2,maxit=3000,sigma=mclustModel$parameters$variance$sigmasq,mu=mclustModel$parameters$mean) 
      
        distribution1 = Norm(mean=model$mu[1],sd=model$sigma[1]) 
	distribution2 = Norm(mean=model$mu[2],sd=model$sigma[2]) 

	delta = abs(model$mu[1]/model$sigma[1]-model$mu[2]/model$sigma[2])
	png("histogram.png")
	hist(col1,col=rgb(0,0,1,0.5), breaks=15,main=expression(paste(H[n], " distribution, ", "bins = 180")), xlab="bins")
	hist(col2,col=rgb(1,0,0,0.5), breaks=15, add=T)
	png("mixtools.png")
	plot(model, which=2,breaks=20)
	print(model$lambda)

	return(c(HellingerDist(distribution1,distribution2),KolmogorovDist(distribution1,distribution2),delta))
}

confusionMatrix <- function(dataFile1, dataFile2, index){
	tp = 0
	fp = 0
	fn = 0
	tn = 0
	
        tab1 = read.csv(dataFile1)
        tab2= read.csv(dataFile2)
        col1 = as.vector(tab1[[index]])
        col2 = as.vector(tab2[[index]])
	
	model = normalmixEM(c(col1,col2),k=2)
	posterior = model$posterior
		
	for(i in 1:(length(col1))){
		if(posterior[i,1]>=posterior[i,2]){
			if(mean(col1)< mean(col2))
				tp = tp + 1
			else
				fp = fp + 1 
		}
		else{
			if(mean(col1)< mean(col2))
				fp = fp + 1
			else
				tp = tp +1
		}
	}
	for(i in length(col1):nrow(posterior)){
                if(posterior[i,1]>=posterior[i,2]){
                        if(mean(col1)> mean(col2))
                                tn = tn + 1
                        else
                                fn = fn + 1
                }
                else{
                        if(mean(col1)> mean(col2))
                                fn = fn + 1
                        else
                                tn = tn +1
                }
        } 
        hist(col1,col=rgb(1,0,0,0.5), breaks=30,probability=TRUE)
	hist(col2,col=rgb(0,0,1,0.5), breaks=30, add=T, probability=TRUE)
        min = min(c(col1,col2))
        max = max(c(col1,col2))
        xseq<-seq(min,max,.01)
 	lines(xseq,dnorm(xseq,model$mu[1],model$sigma[1]))
	lines(xseq,dnorm(xseq,model$mu[2],model$sigma[2]))
	plot(model,which=2,breaks=35)
	png("histogram.pdf")

	return(matrix(c(tp,fp,fn,tn),nrow=2,ncol=2,byrow = TRUE))	
}
#confusionMatrix("output/ellipticals.csv","output/spirals.csv","C1")



v = gaussianMetric("output/r1.csv", "output/r2.csv", "sS3")
write(v,'routput.txt')
