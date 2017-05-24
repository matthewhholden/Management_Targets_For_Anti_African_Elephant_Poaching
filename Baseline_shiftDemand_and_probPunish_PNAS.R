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
punish.prob = seq(0, pmax, length.out = np) #vector of slopes tested (log10 scale)
pcSeq = punish.prob/s0
shift = seq(0,p0/b, length.out= np) 

#initialize variables for stored output

N.eq = array(NA, dim = c(np, np) )

#compute equilibria and their stability 
for(i in 1:np){
  for(j in 1:np){  
    # Equilibrium as function of punishment
    pars = c(r, k, q0, cp, s0, b, s1, mu, 
                    pcSeq[j], cc, p0-shift[i]*b, area, 1,  1, 0, delta) #zv[i]=1
    temp = calculate.equilibria(pars, N0)
    N.eq[i,j] = temp[[1]]
  
  }
}

image(N.eq)


