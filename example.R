### example of the use of the package

## load the package

load(simugo)

## 1st example

## create a network
set.seed(32)
mockgeneset=makeset(nbgo=30,nbgenes=100,nbgopergenes=5)


## create an experiment where longer genes are more active

distribution="normal" # parameters needed to simulate the bias in the measure of the activity
sdBiasGo=0.5
sigma=5
meanBias=10
advantage=0.5
goforce=0.9
nbgoselected=0.5

experiment <- simuExp(mockgeneset,meanBias,sdBiasGo,sigma, distribution,advantage,goforce,nbgoselected)
experiment # this object contains the data for a simulated experieents, with a list of active go and a list of active genes and their simulated covariate


## test for enrichment
result=dotest(experiment) # for each gene set return the adjusted pvalue corresponding to different methods and the mean value of the bias of its genes

#png("example1.png")
plot(result$pval_classic,col=result$active+1,ylab="pvalue",xlab="gene set index") # shows the pvalue color coded by the true activity of the gene set - the pvalue are corrected for mutiple testing with FDR
abline(0.05,0,lty=2) # pvalue threshold at 0.05
text(5,0.07,"pvalue = 0.05 threshold")
#dev.off()
# in this example, only the truly active go are significantly active, but many active go are not picked by the method as significantly active


## 2nd example
set.seed(32)
mockgeneset=makeset(nbgo=30,nbgenes=1000,nbgopergenes=5)

## create an experiment where longer genes are more active

distribution="normal" # parameters needed to simulate the bias in the measure of the activity
sdBiasGo=0.5
sigma=5
meanBias=10
advantage=0.5
goforce=0.9
nbgoselected=0.2
experiment <- simuExp(mockgeneset,meanBias,sdBiasGo,sigma, distribution,advantage,goforce,nbgoselected)
experiment # this object contains the data for a simulated experieents, with a list of active go and a list of active genes and their simulated covariate

## test for enrichment
result=dotest(experiment) # for each gene set return the adjusted pvalue corresponding to different methods and the mean value of the bias of its genes

#png("example2.png")
plot(result$pval_classic,col=result$active+1,ylab="pval",xlab="index of gene set") # shows the pvalue color coded by the true activity of the gene set
abline(0.05,0,lty=2) # pvalue threshold at 0.05
#dev.off()
# in this second, we find similar result as the first one but because we increase the number of genes per gene set, we obtain more true positive

## 3rd example
set.seed(33)
mockgeneset=yeastnetwork()
## in this example we use the real yeast go

distribution="normal" # parameters needed to simulate the bias in the measure of the activity
sdBiasGo=0.5
sigma=5
meanBias=10
advantage=0.5
goforce=0.9
nbgoselected=0.2
experiment <- simuExp(mockgeneset,meanBias,sdBiasGo,sigma, distribution,advantage,goforce,nbgoselected)
#experiment # this object contains the data for a simulated experieents, with a list of active go and a list of active genes and their simulated covariate

## test for enrichment
result=dotest(experiment) # for each gene set return the adjusted pvalue corresponding to different methods and the mean value of the bias of its genes
plot(result$pval_classic,col=result$active+1) # shows the pvalue color coded by the true activity of the gene set
abline(0.05,0,lty=2) # pvalue threshold at 0.05


png("example3.png")
hist(result$pval_classic[which(!result$active)],col=rgb(1,0,0,0.4),xlab="pvalue",ylab="number of go",main="",breaks=seq(0,1,by=0.1))
hist(result$pval_classic[which(result$active==1)],col=rgb(0,0,1,0.4),add=T,breaks=seq(0,1,by=0.1))
legend('top',c("active go","not active go"),col=c("blue","red"),box.col="white",fill=c("blue","red"))
dev.off()


