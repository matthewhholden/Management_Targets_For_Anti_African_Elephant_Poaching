#################################################################
### ODE Functions for elephant poaching supply and demand project
#################################################################

##ODE function used to simulate the poacher and elephant population dynamics - for use with lsoda and runsteady
	ElephantPoacherODE.eff = function(t, x, parms){
		#assign parameters and variables
			pars = parms
			N=x[1]; E=x[2];
			r = pars[1]; k = pars[2]; q0 = pars[3]; 
			cp = pars[4]; s0 = pars[5]; b = pars[6];
			s1 = pars[7]; mu = pars[8]; pc = pars[9];
			cc = pars[10]; p0 = pars[11]; area = pars[12];
			eta = pars[13]; z = pars[14]; 
			m = pars[15]; delta = pars[16];

		#price of an elephant as a function of demand	
			p = p0 - b*(q0*E*N*area + m);
			if(p<0){p=0};	
			s = s0 + s1*m*p
		
		#ODEs for elephant (N) and poacher effort (E) dynamics	
			dNdt = r * N * ( 1 - (N/k)^z ) - q0 * E * N;
			dEdt = delta * E * ( mu*p*q0*N - cp - cc*pc*s);
			
		return( list( c(dNdt, dEdt) ) );
			
	}

	
	##Same ODE as above but takes a combined parameter pcs = pc*s, and ignores the individual
	##values for s0, s1, and pc. This is just a quick and dirty way of running the ODE without 
	##changing the parameter set in the for loop in the script
	ElephantPoacherODE.Punish.eff = function(t, x, parms){
	  #assign parameters and variables
	  pars = parms
	  N=x[1]; E=x[2];
	  r = pars[1]; k = pars[2]; q0 = pars[3]; 
	  cp = pars[4]; s0 = pars[5]; b = pars[6];
	  s1 = pars[7]; mu = pars[8]; pc = pars[9];
	  cc = pars[10]; p0 = pars[11]; area = pars[12];
	  eta = pars[13]; z = pars[14]; 
	  m = pars[15]; delta = pars[16];pcs = pars[17];
	  
	  #price of an elephant as a function of demand	
	  p = p0 - b*(q0*E*N*area+m);
	  if(p<0){p=0};	
	  
	  #ODEs for elephant (N) and poacher effort (E) dynamics	
	  dNdt = r * N * ( 1 - (N/k)^z ) - q0 * E * N;
	  dEdt = delta * E * ( mu*p*q0*N - cp - cc*pcs );
	  return( list( c(dNdt, dEdt) ) );
	}
	
	
	
	#################################################
	####   Functions not used in the Manuscript - but used to validate models
	#################################################
	
	##Dynamic model with an expoential supply and demand curve
	ElephantPoacherODEexp = function(t, x, parms){
	  
	  #assign parameters and variables
	  pars = parms
	  N=x[1]; E=x[2];
	  r = pars[1]; k = pars[2]; q0 = pars[3]; 
	  cp = pars[4]; s0 = pars[5]; b = pars[6];
	  s1 = pars[7]; mu = pars[8]; pc = pars[9];
	  cc = pars[10]; p0 = pars[11]; area = pars[12];
	  eta = pars[13]; z = pars[14]; 
	  m = pars[15]; delta = pars[16];
	  
	  #price of an elephant as a function of demand	
	  p = p0*exp(- b*(q0*E*N*area + m));
	  if(p<0){p=0};	
	  s = s0 + s1*m*p
	  
	  #ODEs for elephant (N) and poacher effort (E) dynamics	
	  dNdt = r * N * ( 1 - (N/k)^z ) - q0 * E * N;
	  dEdt = delta * E * ( mu*p*(1 - pc*s)*q0*N - cp - cc*pc*s*q0*N );
	  
	  return( list( c(dNdt, dEdt) ) );
	  
	}
	