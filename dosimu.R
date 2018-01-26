## functions to do GO simulation in R
## from theoretical distributions and empirical distribution of go size en genelength

dosimu<-function(k, meanGLLT, sdGLLT, nbsimu, nbgenes, nbgo, goID, geneID, model, distribution)
{
    library(MASS)
    parameters = rep(0, 5)
    geneID = paste("GENE", 1:nbgenes, sep = ":")
    goID = paste("GO", 1:nbgo, sep = ":")
    simugenes = data.frame("gene_id" = geneID, "genelength" = rep(0, nbgenes ))
    simugenes$go_id = sample(goID, nbgenes, replace = T)
    goterms = unique(simugenes$go_id)
    rownames(simugenes) = simugenes$gene_id
    simugenes$selected = rep(0, nbgenes)
    ## Random choice of paramters
    ## proportion of go selected vs non selected
    nbgoselected = runif(1, 0, 0.2)
    ## linear factor between probability of being selected and gene length
    advantage = runif(1, 0, 1)
    goforce = runif(1, 0, 1)
    goselected = rbinom(length(goterms), 1, nbgoselected)
    ## save for later
    parameters[1] = nbgoselected
    parameters[2] = advantage
    parameters[3] = goforce
    parameters[4] = distribution
    parameters[5] = model
    ## create a fake distribution of mean of each go    
    gogldistribution = rnorm(length(goterms), meanGLLT, sdGLLT/10)
    ##    gogldistribution = rgamma(length(goterms), meanGLLT/5, scale = 5)
    ##        print("YO")
    ##   gogldistribution = rcauchy(length(goterms), meanGLLT, sdGLLT/10)
    ##    print(distribution)
    ## simulated go datasetm,  ID ,  SELECTED,  GENE LENGTH
    simugo = data.frame("go_id" = goterms, "selected" = goselected, "glteo" = gogldistribution, "glemp" = 0, "nbgenes" = 0)
    rownames(simugo) = goterms
    ##
    ## Simulate a pvalue for each gene
    for(i in 1:nbgenes)
    {
        ## find the go
        thisgo = simugenes$go_id[i][1]
        tgo = which(simugo$go_id == thisgo)
        ## simulate gene length of the gene depending on GO mean GL
        if (distribution == "normal")
            simugenes$genelength[i] = max(0.1, rnorm(1, simugo[tgo, 3], sd = sdGLLT))
        if (distribution == "gamma")
        {
            simugenes$genelength[i] = max(0.1, rgamma(1, simugo[tgo, 3]/5, scale = 5))
                                        #            gogldistribution = rgamma(length(goterms), meanGLLT/5, scale = 5)
        }
        if (distribution == "cauchy")
        {
            simugenes$genelength[i] = max(0.1, rcauchy(1, simugo[tgo, 3], scale = 0.1))
                                        #            print(paste(simugenes$genelength[i], simugo[tgo, 3]))
        }
        simugo[thisgo, ]$glemp = simugo[thisgo, ]$glemp+simugenes$genelength[i]
        simugo[thisgo, ]$nbgenes = simugo[thisgo, ]$nbgenes+1
    }
    simugo$glemp = simugo$glemp*1.0/simugo$nbgenes+1
    print(sum(simugenes$genelength<0))
    divisor1<-(max(simugenes$genelength)*advantage)/(1-goforce)
    divisor2 = (max(exp(simugenes$genelength))*advantage)/(1-goforce)/2    
    for(i in 1:nbgenes)
    {
        ## find the go
        thisgo = simugenes$go_id[i][1]
        tgo = which(simugo$go_id == thisgo)
        
        ## is that go selected ?
        sel = simugo$selected[tgo]        
        ##
        if (sel == 1) ## if yes            
        {
            if (floor(model/2.0) == 0)
            {
                simugenes$selected[i] = rbinom(1, 1, prob = min(goforce+simugenes$genelength[i]/divisor1*advantage, 1))

            }
            else
            {
                simugenes$selected[i] = rbinom(1, 1, prob = advantage)#min(goforce+exp(simugenes$genelength[i])/divisor2*advantage, 1))
            }
        }
        else
        {
            if((model%%2) == 0)
                simugenes$selected[i] = rbinom(1, 1, prob = min(0+simugenes$genelength[i]/divisor1*advantage/2, 1))
            if((model%%2) == 1)
                simugenes$selected[i] = rbinom(1, 1, prob = 0.1)
        }        
    }    
    print("simu done")
    print(cor(simugenes[, c(2, 4)]))
    print(parameters)
    selectset = which(simugenes$go_id %in% simugo$go_id[(simugo$selected == 1)])
    print(cor(simugenes[selectset, c(2, 4)]))
    ## save for each go terms: ID,  pval1,  pval2,  meanGL    
    ##    gofinal = data.frame("ID" = goterms, "pval1" = rep(1, length(goterms)), "pval2" = rep(1, length(goterms)), "pva;l3" = rep(1, length(goterms)))
    ##    print(table((gofinal[, 3]<0.05)-simugo[, 2]))
    save(simugo, simugenes, parameters, file = paste("simusept/simu", paste(parameters[c(4, 5)], collapse = "-"), nbgo, nbgenes, format(Sys.time(),  "%H:%M:%S"), ".Rdata", sep = "-"))
    return(c(cor(simugenes[selectset, c(2, 4)])[1, 2], parameters))
}



doempsimunew<-function(k, simugenes, simugo, model, subsamplesize, foldername)
{    ## pb if two genes in two different gos
    if(subsamplesize>dim(simugenes)[1])
    {
        print(paste("ERROR the subsample", subsamplesize, " is bigger than the size of the dataset", dim(simugenes)[1]))
        return(0)
    }
    ##    simugenes = simugenes[sample(dim(simugenes)[1], subsamplesize), ]
    ##    simugenes$gene_id = simugenes$ensembl_gene_id
    simugo = simugo %>% filter(go_id %in% unique(simugenes$go_id))
    nbgenes = length(simugenes[, 1])
    parameters = rep(0, 5)
    geneID = simugenes$ensembl_gene_id
    goID = simugo$go_id
    goterms = unique(simugenes$go_id)
    ##    rownames(simugenes) = simugenes$ensembl_gene_id
    simugenes$selected = rep(0, nbgenes)
    ## Random choice of paramters
    ## proportion of go selected vs non selected
    nbgoselected = runif(1, 0, 0.5)
    ## linear factor between probability of being selected and gene length
    advantage = runif(1, 0, 1)
    goforce = runif(1, 0, 1)
    nbgo = dim(simugo)[1]
    goselected = rbinom(nbgo, 1, nbgoselected)
    simugo$selected = goselected
    ## save for later
    parameters[1] = nbgoselected
    parameters[2] = advantage
    parameters[3] = goforce
    ##    parameters[4] = distribution
    parameters[5] = model
    ##
    ## Simulate a pvalue for each gene
    divisor1<-(max(simugenes$genelength)*advantage)/(1-goforce)
    divisor2 = (max(exp(simugenes$genelength))*advantage)/(1-goforce)/2    
                                        #    for(i in 1:nbgenes)
    for(g in 1:nbgo)
    {
        ## find the go
        thisgo = simugo$go_id[g]	
        tgenes = which(simugenes$go_id == thisgo)
        ## is that go selected ?
        sel = simugo$selected[g]        
        ##
        if(sel == 1) ## if yes            
        {
            simugenes$selected[tgenes] = rbinom(length(tgenes), 1, prob = min(goforce+simugenes$genelength[tgenes]/divisor1*advantage, 1))
        }
        else
        {
            simugenes$selected[tgenes] = rbinom(length(tgenes), 1, prob = min(0+simugenes$genelength[tgenes]/divisor1*advantage/2, 1))
        }
    }
    uniquegenes = unique(simugenes$ensembl_gene_id)
    for(t in length(uniquegenes))
    {
        thesegenes = which(simugenes$gene_id == uniquegenes[t])
        if(length(thesegenes)>1)
        {
            simugenes$selected[thesegenes] = sample(simugenes$selected[thesegenes], 1)
        }
    }
                                        #    simugenes <- simugenes %>% distinct(gene_id) #remove the doublons
    print("simu done")
    selectset = which(simugenes$go_id %in% simugo$go_id[(simugo$selected == 1)])
    ## save for each go terms: ID,  pval1,  pval2,  meanGL    
    ##    gofinal = data.frame("ID" = goterms, "pval1" = rep(1, length(goterms)), "pval2" = rep(1, length(goterms)), "pva;l3" = rep(1, length(goterms)))
    save(simugo, simugenes, parameters, file = paste(foldername, "/simu", paste(parameters[c(4, 5)], collapse = "-"), nbgo, nbgenes, format(Sys.time(),  "%H:%M:%S"), ".Rdata", sep = ""))
    return(parameters)
    ##    return(parameters)
}

