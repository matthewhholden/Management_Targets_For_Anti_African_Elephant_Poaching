######################################################################################
### Sets up baseline parameterization and global sensitivity analysis parameterization
######################################################################################

## Load functions
library(lhs)
library(deSolve)
library(rootSolve)
source('ODEelephantFunctions_Eff.R')


#Elephant demography parameters
  r = 0.067/365                     #growth rate (in days) - Milner-Gulland 1992, divided by days in year
  rLim = c(.5*r, 1.5*r)             #range of r for sensitivity analysis
  k = 0.198                         #carrying capacity - van Kooten 2008
  kLim = c(.15,.3)                  #range of k for sensitivity analysis
  area  = 3335827;                  #area of elephant habitat in km^2, from CITES 2014 
  areaLim = c(.7*area, area)        #small decrease in area tested for (since some clearing might have occured since 2014)

#probability of getting caught in Laungwa Valley, Zambia
  nKills = 3.54            #kills per expedition
  Te = 7                   #time spent on an expedition
  pL1 = 6/152              #proportion of gangs caught (method 2 in annex 3) Milner Gulland and Leader-Williams 1992
  pL2 = sum(25,22,14,19)/sum(395,412,432,453) #annex table 3.1 Milner Gulland, sum of arrests over gangs
  probCaughtL = mean(c(pL1, pL2))/Te #divide by number of days a gang poaches on average, Te and take mean
  patrolsLKm2Year = c(2, .2) #patrol denisty per km^2 per day, roughly approximated from histograms on pg 1077 of Leader-Williams et al 1990
  areaOfPatrolledLocs = c(400, 212+4630+3800) #from the corresponding table in Leader-Williams et al 1990
  s0L = sum(patrolsLKm2Year*areaOfPatrolledLocs/365)/sum(areaOfPatrolledLocs)
  probCaughtLperS0 = probCaughtL/s0L

#probability of getting caught in Serengeti, Tanzania
  probCaughtS = 0.0007               #knapp 2010, daily prob getting caught while poaching (mostly for bushmeat)
  areaS = 14763;                     #area of Serengeti National Park km^2
  patrolsS = 8.1;                    #num of patrols/day in serengeti park (average in 2004 & 05), Hilborn 2006, 
  s0S = patrolsS/areaS               #number of anti-poaching patrols per day per km^2
  probCaughtSperS0 = probCaughtS/s0S	

#probability of getting caught in Selous, Tanzania
  selousMaxS0 = (8*4)/44800          #the target which has not been met is 8 patols(personal comm.)


#General probability of getting caught parameters using the values above
  probCaught = mean(c(probCaughtLperS0, probCaughtSperS0)) #baseline
  probCaughtLim = c(.5*probCaughtSperS0, 1.5*probCaughtLperS0) #range for sensitivity analysis
  probConvicted = 0.86 #If caught, prob of getting convicted, from Leader-Williams 1990 pg 1073 from Laungwa Valley
  probConvictLim = c(.5,.9)
  pc = probCaught*probConvicted

#Cost parameters
  annualIncome = 762                         #PPP adjusted annual SubSaharan Africa income in USD
  cp = annualIncome/365;                     #Lost wages in USD/day
  cpLim = c(1, 5)                            #Cost of poaching in lost wages in USD/day, ranges for sensitivity analysis
  lifeExpectancy = 58;                       #life expectancy in years for average African male
  fineIvoryMax = cp* 20 *365                 #Kenya's largest prison term (corresponds to a life in prison punishment)
  fineIvoryMin = cp * 365                    #lost wages from a 12 month prison sentence 
  FineIvoryLim = c(fineIvoryMin, fineIvoryMax)
  cc = cp * 10 * 365                         #lost wages from a 10 year prison sentence is the baseline

#catch-ability of elephants (q) - using Selous data, then Serengetti data, then Laungwa Valley data
  qSelous = (44800/420)*(0.00025-log(0.186)/(365*(2013-2006)))         #q in Selous if all poachers were allways poaching
  qSelousNP = (44800/(420*.12))*(0.00025-log(0.186)/(365*(2013-2006))) #q in Selous if we factor in non-poaching time
  gamma = .25/160;                               #proportion of serengeti arrests going after elephant (~160k animals killed per year, ~0.25k elephants)
  qHilborn = 0.049;                              #Hilborn 2006, personal comm. to get fitted parameter value
  qSer = qHilborn*probCaughtSperS0/(365*gamma);  #qHilborn rescaled to our model
  NumExpLV = 432;                                #poaching expedition time length, Milner-Gulland 1992
  ExpLengthLV = 7;                               #poaching expedition time length, Milner-Gulland 1992
  numShootersLV = 2;                             #number of shooters per poaching expedition, Milner-Gulland 1992
  areaLV = 21492;                                #Leader Williams Mills 1990
  numElephantsLV = 25323;                        #number of Elephants in the park, Milner-Gulland 1992
  qLV = nKills*areaLV/(numShootersLV*numElephantsLV*ExpLengthLV);
  q0 = mean(c(qSelousNP, qLV));                  #baseline rate of killing elephants per poacher per day (average of selous and L Valley)
  qSD = sd(c(qSer, qSelous, qSelousNP, qLV))     #standard deviation of the four estimates for q
  qLim = c(qSelous, qSelousNP+2*qSD)               #range for sensitivtiy, contains all values of q from selous, serengeti, L valley 

#price and demand parameters
  mkg = 8.85		               #mass of ivory in kg per elephant
  pkg = 2100; 			           #market price of ivory [USD] per kg (during a demand of 16k per year)
  mu =.062;                    #proportion of money that goes to poachers .062, from Milner Gulland 1993, 
  muLim =c(.024,.10)           #range includes .032 from popsci interview of kenyan poacher 58/1800, https://www.thedodo.com/interview-with-an-elephant-poa-390317914.html more likely to lie down than up
  p = pkg * mkg                #current price of ivory per poached elephant
  currentPriceLim = c(500, 3000)*mkg    #range for current price of ivory
  currentDemandLim = c(10000, 50000)/365#demand between 6000 and 50000 elephants per year at current price         

#stockpile and legal ivory sale parameters
  n = 10;						 	           #number of years ivory is sold for
  stock = 800000;              	 #stockpile in kg
  stockDrip = stock/(365*mkg*n)  #amount of ivory sold per day which would deplete the stockpile in n years
  m = stockDrip;          		   #rate of elephants provided to the legal market [ 0 or (stock/mkg)/365/n]
  mLim = c(.25*m,2*m);
  muc = .5;  
  mucLim = c(0.1,0.8)
  patrolCost = 150;              #cost of anti-poaching patrols in USD per day#proportion of revenue that funds elephant sales
  patrolCostLim = c(100, 300);   #
  s1 = muc/(patrolCost*area);    #number of poaching patrols per km^2 added per USD from legal ivory sales
  s1Lim = c(.1*s1, 2*s1);        #a range that corresponds to a patrol cost of 50 - 500 USD 
  s1Lim2 = mucLim/(patrolCostLim*area)
  
#Supply and demand	
  bVanKooten = 0.0005;           #decrease in price (for which a seller would actually poach) per additional kg of ivory demanded
  b = (bVanKooten*mkg^2*365);    #rate of decrease in demand (elephants/day) per dollar increase in price per elephant
  QD = 20000/365;                #quantity demanded per day at price of p (16186 from old paper), 20000 in 2013
  p0 = p + b*QD;                 #ivory demanded [in # elephants] when free
  bLim = c(.1*b, 10*b)

#patrol Density
  selousMinS0 = (8*1)/44800   #estimate for Selous park min patrol density as low as 1 patrol/day for 8 regions of park 
  s0 = mean(c(s0L, s0S))      #mean of patrol densities in L Valley, and Serengeti
  s0Lim = c(1e-4, s0L)        #this range includes the patrol density in Serengeti, the two Selous densities above and the L Valley density	

### parameters to simulate dynamics (not needed to calculate beta)				
#Initial conditions	
  N0 = .13;                     #initial elephant density
  E0 = (420*.12)/44000;         #initial density of poachers in seleous (.12 is the proportion of days poachers poach, 420 in num poachers)
  delta = 100*r;                 #measure of how quickly poachers adjust behaviour to market equilibrium 
  deltaLim = c(20*r, 200*r)	
  eta = N0/k;
  z = 7;                        #Fowler 1981, for exponential model
  zLim = c(1,10);
  

###############################################################################
##### Latin Hypercube Sample (random parameter values used in global sensitivity analysis)
###############################################################################

set.seed(1)
numParms = 19;
lhs = randomLHS(n = numSims, k = numParms + 1)

# note the parameter order in the matrix below is:
# 1. r, 2. k, 3. q, 4. cp, 5. s0, 
# 6. delta, 7. currentPrice, 8. currentDemand 
# 9. b, 10. m, 11. s1, 
# 12. area 13. mu
# 14. probCaught per unit s0 15. probConvict 
# 16. Ivoryfine 17. tresspassfine 18. probCaughtKill

parLims = matrix( c(rLim, kLim, qLim, cpLim, s0Lim,deltaLim, currentPriceLim, 
                    currentDemandLim, bLim, mLim, s1Lim, areaLim, muLim, 
                    probCaughtLim, probConvictLim, FineIvoryLim, zLim, patrolCostLim, mucLim),       
                  ncol = 2, nrow = numParms, byrow = TRUE
)

limDif = parLims[,2] - parLims[,1];		

parsLHS = array(0, dim=c(numSims, numParms));
for(i in 1:numSims){
  parsLHS[i,] = parLims[,1] + lhs[i,-(numParms+1)] * limDif; #last column removed of LHS sample
}

rv = parsLHS[,1]
kv = parsLHS[,2]
qv = parsLHS[,3]
cpv = parsLHS[,4]
s0v = parsLHS[,5]
dv = parsLHS[,6]
currentPriceV = parsLHS[,7]
currentDemandV = parsLHS[,8]
bv = parsLHS[,9]
mv = parsLHS[,10]
s1v = parsLHS[,11]
areav = parsLHS[,12]
muv = parsLHS[,13]
pCaughtv = parsLHS[,14]
pConv = parsLHS[,15]
ccv = parsLHS[,16]
zv = parsLHS[,17]
patCostv = parsLHS[,18]
mucv = parsLHS[,19]

#composite parameters
  s1v = mucv/(patCostv*area)
  pcv = pCaughtv*pConv;
  p0v = currentPriceV + bv*currentDemandV;