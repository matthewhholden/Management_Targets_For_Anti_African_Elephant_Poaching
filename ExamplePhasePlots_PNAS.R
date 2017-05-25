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


### Identify choice indeces for plotting ####
  ## This code is commented out after identifying typical examples

  #### display the unique types of equilibria situations that can occur in these models ####
  # 1 = stable, 2 = unstable saddle,
  # 3= unstable node or spiral, 0=no more equilibria

  stabCheckFE = array(0, dim = c(numSims*numb, 5) )
  stabCheckM = array(0, dim = c(numSims*numb, 5) )
  stabCheckSQ = array(0, dim = c(numSims*numb, 5) )

  stabCheckFE2 = array(0, dim = c(numSims,numb, 5) )
  stabCheckM2 = array(0, dim = c(numSims,numb, 5) )
  stabCheckSQ2 = array(0, dim = c(numSims,numb, 5) )

  for(i in 1:numSims){
    for(j in 1:numb){
      stabCheckFE[(i-1)*numb + j, 1:N.numRoots.FE[i,j] ] = stab.FE[[j]][[i]]
      stabCheckM[(i-1)*numb + j, 1:N.numRoots.M[i,j] ] = stab.M[[j]][[i]]
      stabCheckSQ[(i-1)*numb + j, 1:N.numRoots.SQ[i,j] ] = stab.SQ[[j]][[i]]

      stabCheckFE2[i, j, 1:N.numRoots.FE[i,j] ] = stab.FE[[j]][[i]]
      stabCheckM2[i, j, 1:N.numRoots.M[i,j]  ] = stab.M[[j]][[i]]
      stabCheckSQ2[i, j, 1:N.numRoots.SQ[i,j]  ] = stab.SQ[[j]][[i]]
    }
  }


  #For both status quo and legal trade dynamics only 3 qualitatively different types of
  #phase portraits are possible (1) 3 equilibria where the middle one is stable, row 1
  #(2) (0,0) and (k,0) equilibria where carrying capacity is stable, row 3 and
  #(3) 5 equilibria where the middle equilibrium is an unstable saddle, initial elephant
  #densities below this equilibria approach the smaller second equilibrium and above this
  #equilibrium approach the larger fourth equilibrium, both N=0 and N=k are also unstable
  #saddles in this case. See Auger et al (2014) for a detailed mathematical analysis of
  #potential dynamics from this model.
  stab.types.M = unique(stabCheckM, MARGIN = 1)
  stab.types.SQ = unique(stabCheckSQ, MARGIN = 1)

  #when funding frome sales increases enforcement this leads to the posability of 2
  #additional types of dyanmics including (1) row 5: 4 equilibria where only carrying capacity is stable
  #and (2) row 4, 5 equilbria, where the smallest nonzero one is an unstable source
  #rather than a saddle
  stab.types.FE = unique(stabCheckFE, MARGIN = 1)

  #find all the indices where different stability types occur
  # 1 = 3 equilibria (0 = saddle, stable , k = saddle)
  # 2 = 5 equilibria (0 = saddle, stable , saddle, stable, k = saddle)
  # 3 = 2 equilibria (0 = saddle, k = stable)
  # 4 = 5 equilibria (0 = saddle, unstable, saddle, k = stable)
  # 5 = 4 equilibria (0 = saddle, unstable, saddle, stable, k = saddle)


  categories.SQ = apply( stabCheckSQ2, c(1,2), function(x) row.match(x, stab.types.SQ) )
  categories.M = apply( stabCheckM2, c(1,2), function(x) row.match(x, stab.types.M) )
  categories.FE = apply( stabCheckFE2, c(1,2), function(x) row.match(x, stab.types.FE) )

  num.cat.b.SQ = apply(categories.SQ, 1, function(x) (length(unlist(unique(x)))) )
  num.cat.b.M = apply(categories.M, 1, function(x) (length(unlist(unique(x)))) )
  num.cat.b.FE = apply(categories.FE, 1, function(x) (length(unlist(unique(x)))) )



  # pick an example scenario i, and plot null clines and phase diagram to verify the model's correctness
  ind1FE = which( categories.FE == 1, arr.ind=TRUE );
  ind2FE = which( categories.FE == 2, arr.ind=TRUE );
  ind3FE = which( categories.FE == 3, arr.ind=TRUE );
  ind4FE = which( categories.FE == 4, arr.ind=TRUE );
  ind5FE = which( categories.FE == 5, arr.ind=TRUE );

  ind1M = which( categories.M == 1, arr.ind=TRUE );
  ind2M = which( categories.M == 2, arr.ind=TRUE );
  ind3M = which( categories.M == 3, arr.ind=TRUE );

  ind1SQ = which( categories.SQ == 1, arr.ind=TRUE );
  ind2SQ = which( categories.SQ == 2, arr.ind=TRUE );
  ind3SQ = which( categories.SQ == 3, arr.ind=TRUE );

#uncomment this code to plot the nullclines of all studies for all scenarios. There 
#are only 5 possibible combinations of roots and there stability. You can inspect
#These graphs to pick the cases that illustrate the type of dynamic behavior.
  
  # Plot.all.types(l=dim(ind1M)[1], ind1M, ty = 1, blim = c(6,7), strat = 'M')
  # dev.off(dev.list()["RStudioGD"])
  # Plot.all.types(l=dim(ind2M)[1], ind2M, ty = 2, blim = c(6,7), 'M')
  # Plot.all.types(l=dim(ind3M)[1], ind3M, ty = 3, blim = c(6,7), 'M')
  # 
  # Plot.all.types(dim(ind1SQ)[1], ind1SQ, 1, blim = c(6,7), 'SQ')
  # Plot.all.types(dim(ind2SQ)[1], ind2SQ, 2, blim = c(6,7), 'SQ')
  # Plot.all.types(dim(ind3SQ)[1], ind3SQ, 3, blim = c(6,7), 'SQ')
  # 
  # Plot.all.types(dim(ind1FE)[1], ind1FE, 1, blim = c(6,7), strat = 'FE')
  # Plot.all.types(dim(ind2FE)[1], ind2FE, 2, blim = c(6,7), 'FE')
  # Plot.all.types(dim(ind3FE)[1], ind3FE, 3, blim = c(6,7), 'FE')
  # Plot.all.types(dim(ind4FE)[1], ind4FE, 4, blim = c(6,7), 'FE')
  # Plot.all.types(dim(ind5FE)[1], ind5FE, 5, blim = c(1,9), 'FE')

## Plot the Phase planes of typical examples (by simulating the ODES ) ####

#indeces of typical examples - choose them from visually inspecting the
#roots from the indeces provided above (by plotting the nullclines - see
#commented out code above)
  
  ind.typical.FE = matrix( c(72, 7,
                         962, 7,
                         495,7,
                         772,7,
                         313, 9), 5, 2, byrow=TRUE)
  
  ind.typical.M = matrix( c(313, 7,
                             515, 7,
                             818,7,
                             327,4), 4, 2, byrow=TRUE)

#initialize variables
  init1 = c(N0,E0)
  init2 = c(.03,.0015)
  t.end = 40000
  lt = t.end*2+1
  times = seq(0, t.end, length.out = lt)
  out1FE = array(NA, dim = c(5,lt,3))
  out2FE = array(NA, dim = c(5,lt,3))
  out1M = array(NA, dim = c(4,lt,3))
  out2M = array(NA, dim = c(4,lt,3))

## Solve for typical ODE trajectories in the cases where legal ivory trade funds enforcement
  for(counter in 1:5){
    print(counter)
    i = ind.typical.FE[counter,1] 
    j = ind.typical.FE[counter,2]
    pars = c(rv[i], kv[i], qv[i], cpv[i], s0v[i], bSeq[j], s1v[i], muv[i], 
             pcv[i], ccv[i], p0v[i,j], area, 1,  1, stockDrip, delta)
    out1FE[counter,,] = lsoda(init1, times = times, 
                 func = ElephantPoacherODE.eff, 
                 parms = pars, hmax=.002)
    out2FE[counter,,] = lsoda(init2, times = times, 
                 func = ElephantPoacherODE.eff, 
                 parms = pars, hmax=.0005)
  }

## Solve for typical ODE trajectories in the cases where legal ivory DOESN'T fund enforcement
  par(mfrow=c(3,2))
  for(counter in 1:4){
    i = ind.typical.M[counter,1] 
    j = ind.typical.M[counter,2]
    pars = c(rv[i], kv[i], qv[i], cpv[i], s0v[i], bSeq[j], 0, muv[i], 
             pcv[i], ccv[i], p0v[i,j], area, 1,  1, stockDrip, delta)
    out1M[counter,,] = lsoda(init1, times = times, 
                 func = ElephantPoacherODE.eff, 
                 parms = pars, hmax=.01)
    out2M[counter,,] = lsoda(init2, times = times, 
                 func = ElephantPoacherODE.eff, 
                 parms = pars, hmax=.01)
  }

save.image('phasePlotData.Rdata')



## plot the trajectories in the phase plane
col.size = 3.4252
pdf('PhasePlots.pdf', width=2*col.size, height=2*col.size)
  col.traj ='gray25'
  lw.traj = 2.5
  cexa = 2
  
  titlesM = c('(A)','(B)','(C)','(D)')
  titlesFE = c('(E)','(F)','(G)','(H)','(I)')
  par(mfrow=c(3,3) , oma=c(4,4,0,0), mar = c(3,3,1.5,1) )
  options(scipen=999)
  step.a = round(.6*lt)
  step.a1 = round(.2*lt)
  ## M phase plots (no funding enforcement with legal sales)
    strat='M'
    #par(mfrow=c(3,2))
    for(counter in 1:4){
      i = ind.typical.M[counter,1] 
      j = ind.typical.M[counter,2]
      pars = c(rv[i], kv[i], qv[i], cpv[i], s0v[i], bSeq[j], 0, muv[i], 
               pcv[i], ccv[i], p0v[i,j], area, 1,  1, stockDrip, delta)
      fixedPoints = roots.M[[j]][[i]]
      stability = stab.M[[j]][[i]]
      PhasePlot(pars, strat, fixedPoints, stability, title = titlesM[counter])
      lines(out1M[counter,,2], out1M[counter,,3], col=col.traj, lwd=lw.traj, lty=2)
      arrow.step = ifelse(counter>1,step.a,step.a1)
      arrows(out1M[counter,lt-arrow.step,2], out1M[counter,lt-arrow.step,3],
            out1M[counter,lt,2], out1M[counter,lt,3], 
            lwd = lw.traj, col = col.traj, length = .15, xpd=TRUE)
    }
  
  ## FE phase plots (funding enforcement with legal sales)
    strat='FE'
    for(counter in 1:5){
      i = ind.typical.FE[counter,1] 
      j = ind.typical.FE[counter,2]
      pars = c(rv[i], kv[i], qv[i], cpv[i], s0v[i], bSeq[j], s1v[i], muv[i], 
               pcv[i], ccv[i], p0v[i,j], area, 1,  1, stockDrip, delta)
      fixedPoints = roots.FE[[j]][[i]]
      stability = stab.FE[[j]][[i]]
      PhasePlot(pars, strat, fixedPoints, stability, title = titlesFE[counter])
      lines(out1FE[counter,,2], out1FE[counter,,3], col=col.traj, lwd=lw.traj, lty=2)
      arrows(out1FE[counter,lt-step.a,2], out1FE[counter,lt-step.a,3],
             out1FE[counter,lt,2], out1FE[counter,lt,3], 
             lwd = lw.traj, col = col.traj, length = .15, xpd=TRUE)
    }
    
    mtext(expression( paste('Elephants / km'^'2') ), 
          outer=TRUE, side=1, cex = cexa, padj=.6)
    mtext(expression( paste('Poachers / km'^'2') ), 
          outer=TRUE, side=2, cex = cexa, padj=-.2)
  
dev.off()