
#### Compute Equilibrium  for Demand Reduction and punishment probs  ####

rm(list=ls())
numSims = 1;
source('parameterizeEffort.R')
source('Functions_equilibrium.R')
source('ODEelephantFunctions_Eff.R')

#other parameters (not in parameterizeEffort)
  delta = 100*r
  np = 500 #number of parameter values tested
  pmax = .02
  
  
  ## compute Equilibria in the baseline scenario ####
  

  y0 = c(N0, E0/10)
  pars = c( r, k, q0, cp, s0, b, 0, mu, pc, cc, p0, area, eta, 1, 0, delta)
  parsFE = c( r, k, q0, cp, s0, b, s1, mu, pc, cc, p0, area, eta, 1, stockDrip, delta)
  parsNFE = c( r, k, q0, cp, s0, b, 0, mu, pc, cc, p0, area, eta, 1, stockDrip, delta)
  out =  calculate.equilibria(pars, N0)[[1]]
  outFE = calculate.equilibria(parsFE, N0)[[1]]
  outNFE = calculate.equilibria(parsNFE, N0)[[1]]  
  
  
  
#Punishment probability and Demand shift [1D sensitivity - for combination see other file]  
  punish.prob = seq(0, pmax, length.out = np) #vector of slopes tested (log10 scale)
  pcSeq = punish.prob/s0
  
  shift = seq(0,p0/b, length.out= np) 

#initialize variables for stored output
  roots.punish=list()
  stab.punish=list()
  N.eq.punish=rep(NULL,np)
  N.numRoots.punish=rep(NULL,np)

  roots.demand=list()
  stab.demand=list()
  N.eq.demand=rep(NULL,np)
  N.numRoots.demand=rep(NULL,np)
  

#compute equilibria and their stability 
  for(j in 1:np){
      
      # Equilibrium as function of punishment
        pars.punish = c(r, k, q0, cp, s0, b, s1, mu, 
                   pcSeq[j], cc, p0, area, 1,  1, 0, delta) #zv[i]=1
        
        temp.punish = calculate.equilibria(pars.punish, N0)
        
        N.eq.punish[j] = temp.punish[[1]]
        N.numRoots.punish[j] = temp.punish[[2]]
        roots.punish[[j]] = temp.punish[[4]]
        stab.punish[[j]] = temp.punish[[5]]
        
      # Equilibrium as function of demand shift
        pars.demand = c(r, k, q0, cp, s0, b, s1, mu, 
                        pc, cc, p0-shift[j]*b, area, 1,  1, 0, delta) #zv[i]=1
        
        temp.demand = calculate.equilibria(pars.demand, N0)
        
        N.eq.demand[j] = temp.demand[[1]]
        N.numRoots.demand[j] = temp.demand[[2]]
        roots.demand[[j]] = temp.demand[[4]]
        stab.demand[[j]] = temp.demand[[5]]  
  
  }

plot(punish.prob, N.eq.punish, type = 'l', ylim = c(0,k))
plot(shift, N.eq.demand, type = 'l', ylim = c(0,k))



