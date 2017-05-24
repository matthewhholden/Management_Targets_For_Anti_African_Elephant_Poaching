
# Functions ####

dEdt.nulcline = function(pars, N){
  r = pars[1]; k = pars[2]; q0 = pars[3]; 
  cp = pars[4]; s0 = pars[5]; b = pars[6];
  s1 = pars[7]; mu = pars[8]; pc = pars[9];
  cc = pars[10]; p0 = pars[11]; area = pars[12];
  eta = pars[13]; z = pars[14]; 
  m = pars[15]; delta = pars[16];
  
  E = (-cp - b*m*mu*N*q0 + mu*N*p0*q0 - cc*pc*s0 + b*cc*m^2*pc*s1 - 
         cc*m*p0*pc*s1)/(area*b*N*q0*(mu*N*q0 - cc*m*pc*s1))
  
  return(E)
}

stability = function(eigvals){
  #imput eigenvalues and determine stability 
  #1: stable (nonde or spiral), 2: saddle 3: unstable 4:indeterminate form
  if( sum( Re(eigvals) < 0 ) == 2 ){
    stab = 1;
  } else if ( sum ( Re(eigvals) > 0 ) == 2){
    stab = 3;
  } else if ( eigvals[1]*eigvals[2]<0 ) {
    stab = 2;
  } else { stab = NA }
}

dNdt.nulcline = function(pars, N){
  r = pars[1]; k = pars[2]; q0 = pars[3]; 
  cp = pars[4]; s0 = pars[5]; b = pars[6];
  s1 = pars[7]; mu = pars[8]; pc = pars[9];
  cc = pars[10]; p0 = pars[11]; area = pars[12];
  eta = pars[13]; z = pars[14]; 
  m = pars[15]; delta = pars[16];
  
  E = (r/q0)*(1-N/k)
  
  return(E)
}

Jacobian = function(pars, x){
  r = pars[1]; k = pars[2]; q0 = pars[3]; 
  cp = pars[4]; s0 = pars[5]; b = pars[6];
  s1 = pars[7]; mu = pars[8]; pc = pars[9];
  cc = pars[10]; p0 = pars[11]; area = pars[12];
  eta = pars[13]; z = pars[14]; 
  m = pars[15]; delta = pars[16];
  n = x[1]; e = x[2];
  
  dfdn = (-e)*q0 - (n*r)/k + (1 - n/k)*r 
  dfde = (-n)*q0
  dgdn = e*((-area)*b*e*mu*n*q0^2 + mu*q0*(p0 - b*(m + area*e*n*q0)) + 
              area*b*cc*e*m*pc*q0*s1)
  dgde = -cp + mu*n*q0*(p0 - b*(m + area*e*n*q0)) + 
    e*((-area)*b*mu*n^2*q0^2 + area*b*cc*m*n*pc*q0*s1) - 
    cc*pc*(s0 + m*(p0 - b*(m + area*e*n*q0))*s1)
  
  J = matrix(c(dfdn, dfde, dgdn, dgde), 2, 2, byrow = TRUE)
  
  return(J)
  
}


dEdtCoefsN = function(pars){
  #coefficients for the polynomial (in N) dE/dt
  
  r = pars[1]; k = pars[2]; q0 = pars[3]; 
  cp = pars[4]; s0 = pars[5]; b = pars[6];
  s1 = pars[7]; mu = pars[8]; pc = pars[9];
  cc = pars[10]; p0 = pars[11]; area = pars[12];
  eta = pars[13]; z = pars[14]; 
  m = pars[15]; delta = pars[16];
  
  coef0 = -cp - cc*pc*s0 + b*cc*m^2*pc*s1 - cc*m*p0*pc*s1
  coef1 = (-b*m*mu*q0 + mu*p0*q0 + area*b*cc*m*pc*r*s1)
  coef2 = (-area*b*mu*q0*r - (area*b*cc*m*pc*r*s1)/k)
  coef3 = (area*b*mu*q0*r)/k
  return( c(coef0, coef1, coef2, coef3) )	
}

calculate.equilibria = function(pars, N0){
  r = pars[1]; k = pars[2]; q0 = pars[3]; 
  cp = pars[4]; s0 = pars[5]; b = pars[6];
  s1 = pars[7]; mu = pars[8]; pc = pars[9];
  cc = pars[10]; p0 = pars[11]; area = pars[12];
  eta = pars[13]; z = pars[14]; 
  m = pars[15]; delta = pars[16];
  
  
  beta = Beta(pars)
  root = polyroot( dEdtCoefsN(pars) )
  N.eq = NA;
  
  N.nonTrivRoots =  sort( Re( root[ 
    which( abs(Im(root)) < 1e-6 &
             Re(root) > 0 & Re(root) < k ) ] ))
  E.nonTrivRoots = dNdt.nulcline(pars, N.nonTrivRoots)
  N.roots = c(0, N.nonTrivRoots, k);
  E.roots = c(0, E.nonTrivRoots, 0);
  
  rl = length(N.roots)
  
  
  #compute stability of all roots
    N.stab = rep(NA, rl)
    for(i in 1:rl){
      N.stab[i] = stability(eigen( Jacobian(pars, c(N.roots[i], E.roots[i])) )$values )
    }
  
    
   # print(rl); print(N.stab); #print( (N.stab==1) == c(0,1,0,1,0) ) ;
    
    
  # Determine the stability of the roots
  if( (rl == 2 || rl == 3) && N.stab[2] == 1){
    N.eq = N.roots[2]
  } else if (rl == 5 && sum( (N.stab==1) == c(0,1,0,1,0)) ==5 ){
    N.eq = ifelse(N0 > N.roots[3], N.roots[4], N.roots[2])
  } else if (rl == 5 && sum( (N.stab==1) == c(0,0,0,1,0)) ==5 ){
    N.eq = N.roots[4]
  } else if (rl == 4 && N.stab[4] == 1 && 
             N.stab[3]>1 && N.stab[2]>1 && N.stab[1]>1){
    N.eq = N.roots[4]
  } 
  
  return(list(N.eq, rl, beta, N.roots, N.stab))
}



Beta = function(pars){
  
  r = pars[1]; k = pars[2]; q0 = pars[3]; 
  cp = pars[4]; s0 = pars[5]; b = pars[6];
  s1 = pars[7]; mu = pars[8]; pc = pars[9];
  cc = pars[10]; p0 = pars[11]; area = pars[12];
  eta = pars[13]; z = pars[14]; 
  m = pars[15]; delta = pars[16];
  
  beta = cp + cc*s0 + (p0 - b*m)*(m*cc*s1 - mu*q0*k)
  beta2 = r + cp + cc*pc*s0 + (p0 - b*m)*(m*cc*pc*s1 - mu*q0*k)
  return(beta2)
}



#################### Function to plot the nullclines for several examples
Plot.all.types =  function(l, ind, ty=NA, blim = c(1,numb), strat='SQ'){
  par(mfrow=c(4,4))
  sDrip = stockDrip
  s.new = rep(0, numSims)
  
  if(strat == 'FE'){ s.new = s1v} else if ( strat == 'SQ') {sDrip = 0}
  
  for(counter in 1:l){
    if(dim(ind)[1]>0){
      i = ind[counter,1] 
      j = ind[counter,2]
      
      if(j>=blim[1] && j<=blim[2]){
        
        if( strat == 'FE'){
          ass = s.new[i]*ccv[i]*pcv[i]*sDrip/(muv[i]*qv[i])
          buffer = .00002
          if (ass > 0 && ty == 4){
            N = c(seq(buffer, ass-buffer, by = min(.001,ass)), seq(ass+buffer, kv[i], by = .001))
          } else { N = c(seq(0, kv[i], by = .001)) } 
        } else {
          N = c(seq(0, kv[i], by = .0002))
        }
        
        #nullcline plots for legal trade funding enforcement
        pars = c(rv[i], kv[i], qv[i], cpv[i], s0v[i], bSeq[j], s.new[i], muv[i], 
                 pcv[i], ccv[i], p0v[i,j], area, 1,  1, sDrip, delta)
        E.E = dEdt.nulcline(pars, N)
        E.N = dNdt.nulcline(pars, N)
        plot(N, E.E,   
             xlim = c(0,max(1.05*N)),  ylim = c(0,max(2*E.N)), 
             type = 'l', main = paste(i,', ',j) )
        lines(N, E.N, col='blue')
        
        #equilibria
        if(strat =='FE'){  
          points( roots.FE[[j]][[i]], c(0, dNdt.nulcline(pars, roots.FE[[j]][[i]][-1])), 
                  lwd=2, col='green' )
          abline(v=ass, col = 'red', lty = 2)
        } else if(strat == 'M'){
          points( roots.M[[j]][[i]], c(0, dNdt.nulcline(pars, roots.M[[j]][[i]][-1])), 
                  lwd=2, col='green' )
        } else {
          points( roots.SQ[[j]][[i]], c(0, dNdt.nulcline(pars, roots.SQ[[j]][[i]][-1])), 
                  lwd=2, col='green' )
        }
      }
    } else { print(paste('no solutions of type ', ty, ' for strategy', strat)) }
  }
}

#plot typical phase plot ####
  PhasePlot = function(pars, strat, fp, stab, 
                       title='', lab.axis = c('',''), 
                       cexp=3, ca=1.2, lw=3){
    
    r = pars[1]; k = pars[2]; q0 = pars[3]; 
    cp = pars[4]; s0 = pars[5]; b = pars[6];
    s1 = pars[7]; mu = pars[8]; pc = pars[9];
    cc = pars[10]; p0 = pars[11]; area = pars[12];
    eta = pars[13]; z = pars[14]; 
    m = pars[15]; delta = pars[16];
    
    mult.k = 1.1;
    col.null.e = 'mediumseagreen'
    col.null.n = 'lightpink3'
    col.fp = 'black'
    N = seq(0, mult.k*k, by = .0001)
    symbols = 18*(stab == 1) + 1
    ylimit = c(0,.0022)
    xlimit = c(0,mult.k*k)
    
    if(strat=='FE'){
      ass = s1*cc*pc*m/(mu*q0)
      if(ass>0){ 
        buffer = 0.0001; 
        N1 = c(seq(buffer, ass-buffer, by = min(buffer, ass)))
        if(ass<mult.k*k) N2 = seq(ass+buffer,  mult.k*k, by = buffer)#, seq(ass+buffer, k, by = buffer))
      } 
    }
    
    if(strat=='FE' && ass>0 ){
      
      E.E1 = dEdt.nulcline(pars, N1)
      E.N1 = dNdt.nulcline(pars, N1)

      plot(N1, E.E1,   
           xlim = xlimit,  ylim = ylimit, 
           type = 'l', main = title, 
           xlab = lab.axis[1], ylab=lab.axis[2],
           xaxs="i", yaxs="i", lwd=lw, col=col.null.e)
      lines(N1, E.N1, col=col.null.n , lwd=lw)
      lines(xlimit, c(0,0), col=col.null.e, lwd=lw, xpd=TRUE)
      lines(c(0,0), ylimit, col=col.null.n, lwd=lw, xpd=TRUE)
      if( ass<mult.k*k ){
        E.E2 = dEdt.nulcline(pars, N2)
        E.N2 = dNdt.nulcline(pars, N2)
        lines(N2, E.E2,   
             type = 'l', main = paste(i,', ',j),
             xaxs="i", yaxs="i", lwd=lw, col=col.null.e)
        lines(N2, E.N2, col=col.null.n , lwd=lw)
      }
      points( fp, c(0, dNdt.nulcline(pars, fp[-1])), 
              cex = cexp, lwd=lw , pch = 16, col = 'white', xpd=TRUE)
      points( fp, c(0, dNdt.nulcline(pars, fp[-1])), 
              cex = cexp, lwd=lw , pch = symbols, col = col.fp, xpd=TRUE)
    } else {
      #print('got here')
      E.E = dEdt.nulcline(pars, N)
      E.N = dNdt.nulcline(pars, N)
      ylimit = c(0, max(E.N,E.E)+.0001)
      plot(N, E.E,   
           xlim = xlimit,  ylim = ylimit, 
           type = 'l', main = title,
           xlab = lab.axis[1], ylab=lab.axis[2],
           xaxs="i", yaxs="i", lwd=lw, col=col.null.e)
      lines(N, E.N, col=col.null.n , lwd=lw)
      lines(xlimit, c(0,0), col=col.null.e, lwd=lw, xpd=TRUE)
      lines(c(0,0), ylimit, col=col.null.n, lwd=lw, xpd=TRUE)
      points( fp, c(0, dNdt.nulcline(pars, fp[-1])), 
              cex = cexp, lwd=lw, pch = 16, col = 'white', xpd=TRUE )
      points( fp, c(0, dNdt.nulcline(pars, fp[-1])), 
              cex = cexp, lwd=lw, pch = symbols, col = col.fp, xpd=TRUE )
      #if( strat=='FE' && ass>0 ){abline(v=ass, col = 'grey', lty = 2, lwd=lw)}
    }
    
  }
