#### Compute Equilibrium  for Demand Reduction and punishment probs  ####

rm(list=ls())
numSims = 10000
source('parameterizeEffort.R')
source('Functions_equilibrium.R')
source('ODEelephantFunctions_Eff.R')

#other parameters (not in parameterizeEffort)
delta = 100*r
np = 200 #number of parameter values tested
pmax = .09
punish.prob = seq(0, pmax, length.out = np) #vector of slopes tested (log10 scale)
pcSeq = punish.prob/s0
shift = seq(0,1, length.out= np) 
sim.time = 5.3e-4*np^2*numSims

#warn the user how long the simulation is expected to take on 
#a 4-core desktop computer 
print( paste('predicted simulation time is', sim.time,
             ' secs, aka', sim.time/60, 
             ' mins, aka', sim.time/(60^2),' hours, on my desktop, yours may vary') )
#initialize variables for stored output

N.eq = array(NA, dim = c(np, np, numSims) )

tic = Sys.time()
  #compute equilibria and their stability 
  for(i in 1:np){
    print(i/np)
    for(j in 1:np){ 
      for(w in 1:numSims){
        # Equilibrium as function of punishment
        pars = c(rv[w], kv[w], qv[w], cpv[w], s0, bv[w], s1v[w], muv[w], 
                 punish.prob[j]/s0, ccv[w], p0v[w]*(1-shift[i]), area, 1,  1, 0, delta) #zv[w][i]=1
        temp = calculate.equilibria(pars, N0)
        N.eq[i,j, w] = temp[[1]]
      }
    }
  }
print( Sys.time() - tic )

prop = apply(N.eq > N0, c(1,2), sum)/numSims

save.image('Punish_and_shift_sensitivity.Rdata')

#sumarize a sensitivity analysis
# y = N.eq[10,10,]/(kv+1e-8)
# summary( betareg(y ~ ccv+qv+cpv+p0v+bv+kv+rv,link='logit' ))
