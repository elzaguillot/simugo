### example of the use of the package


## create a network
set.seed(32)
mockgeneset=makeset(nbgo=30,nbgenes=100,nbgopergenes=5)

## create an experiment where longer genes are more active

distribution="normal" # parameters needed to simulate the bias in the measure of the activity
sdBias=10
meanBias=10
advantage=0.5
goforce=0.9
nbgoselected=0.5

experiment <- simuExp(mockgeneset,meanBias,sdBias, distribution,popGoExp,advantage,goforce,nbgoselected)
experiment # this object contains the data for a simulated experieents, with a list of active go and a list of active genes and their simulated covariate


## test for enrichment 
dotest(experiment) # for each go return a pvalue corresponding to different methods and the mean value of the bias of its genes
