

#### Compute Equilibrium and their stability - for all scenarios for different b values ####

  rm(list=ls())
  library(prodlim)
  numSims = 1000;
  source('parameterizeEffort.R')
  source('Functions_equilibrium.R')

#initialize variables
  delta = 100*r
  bMax = 10^5 #largest demand curve slope tested
  numb = 9 #number of slopes tested
  bSeq = 10^seq(1, log10(bMax), length.out = numb) #vector of slopes tested (log10 scale)
  p0v = array(NA, dim=c(numSims, numb) ) #note p0 depends b, and price demand data point, so must be generated for each scenario, using the formula below, this line of code initializes the value
    
  N.eq.SQ = matrix( NA, numSims, numb)
  N.numRoots.SQ = matrix( NA, numSims, numb)
  beta.SQ = matrix( NA, numSims, numb)
  
  N.eq.M = matrix( NA, numSims, numb)
  N.numRoots.M = matrix( NA, numSims, numb)
  beta.M = matrix( NA, numSims, numb)
  
  N.eq.FE = matrix( NA, numSims, numb)
  N.numRoots.FE = matrix( NA, numSims, numb)
  beta.FE = matrix( NA, numSims, numb)
  
  roots.FE = list()
  stab.FE = list()
  
  roots.M = list()
  stab.M = list()
  
  roots.SQ = list()
  stab.SQ = list()
  
#compute equilibria and their stability 
  for(j in 1:numb){
    
    p0v[,j]  = currentPriceV[1:numSims] + bSeq[j]*currentDemandV[1:numSims];  
    roots.FE[[j]] = list()
    stab.FE[[j]] = list()
    roots.M[[j]] = list()
    stab.M[[j]] = list()
    roots.SQ[[j]] = list()
    stab.SQ[[j]] = list()
    
    for(i in 1:numSims){
      
      # set parameters
      parsFE = c(rv[i], kv[i], qv[i], cpv[i], s0v[i], bSeq[j], s1v[i], muv[i], 
               pcv[i], ccv[i], p0v[i,j], area, 1,  1, stockDrip, delta) #zv[i]=1
      parsM = c(rv[i], kv[i], qv[i], cpv[i], s0v[i], bSeq[j], 0, muv[i], 
                pcv[i], ccv[i], p0v[i,j], area, 1,  1, stockDrip, delta) #zv[i]=1
      parsSQ = c(rv[i], kv[i], qv[i], cpv[i], s0v[i], bSeq[j], 0, muv[i], 
                pcv[i], ccv[i], p0v[i,j], area, 1,  1, 0, delta) #zv[i]=1
      tempFE = calculate.equilibria(parsFE, N0)
      tempM = calculate.equilibria(parsM, N0)
      tempSQ = calculate.equilibria(parsSQ, N0)
      
      N.eq.FE[i,j] = tempFE[[1]]
      N.numRoots.FE[i,j] = tempFE[[2]]
      beta.FE[i,j] = tempFE[[3]]
      
      N.eq.M[i,j] = tempM[[1]]
      N.numRoots.M[i,j] = tempM[[2]]
      beta.M[i,j] = tempM[[3]]
      
      N.eq.SQ[i,j] = tempSQ[[1]]
      N.numRoots.SQ[i,j] = tempSQ[[2]]
      beta.SQ[i,j] = tempSQ[[3]]
      
      roots.FE[[j]][[i]] = tempFE[[4]]
      stab.FE[[j]][[i]] = tempFE[[5]]
      roots.M[[j]][[i]] = tempM[[4]]
      stab.M[[j]][[i]] = tempM[[5]]
      roots.SQ[[j]][[i]] = tempSQ[[4]]
      stab.SQ[[j]][[i]] = tempSQ[[5]]
      
    }
  }
  
#### For plotting: Store 'Proportion of sucessful scenarios' variables for each management strategy ####  
  par(mfrow=c(3,3));for(j in 1:numb){ hist(N.eq.SQ[,j]/kv)}
  par(mfrow=c(3,3));for(j in 1:numb){ hist(N.eq.M[,j]/kv)}
  par(mfrow=c(3,3));for(j in 1:numb){ hist(N.eq.FE[,j]/kv)}
  
  prop.SQ = apply(N.eq.SQ > N0, 2, sum)/numSims
  prop.M = apply(N.eq.M > N0, 2, sum)/numSims
  prop.FE = apply(N.eq.FE > N0, 2, sum)/numSims
  

       
          
          
          