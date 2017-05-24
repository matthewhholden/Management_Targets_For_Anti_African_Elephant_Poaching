# Management_Targets_For_Anti_African_Poaching

The files allow one to reproduce the work in manuscript
"Management Targets to Save the African Elephant From Poaching"

Note all files should be stored in one directory, and the user
should set the working directory to this location

parameterizeEffort.R 

  Sets the baseline parameter values and generates the values 
  for the sensitivity analysis as described in the supplement


ODEelephantFunctions_Eff.R

  Contains the ordinary differential equation functions which 
  can then be run with function `lsoda`. These are
  used to simulate the dynamics, and used to verify the
  dynamics behaves as described in the equilibrium analysis


Functions_equilibrium.R

  Contains functions to compute the equilibria of the model,
  taking adavantage of the fact that the equilibria are
  the solutions of a 3rd or 5th order polynomial, using the
  function polyroot. There are also functions to compute the 
  stability of the equilibria and clasify which equilibrium
  will be approached from a given initial condition


Baseline_Analysis_for_bar_plot_and_1D_analysis.R

  Runs the analysis to produce the data in Fig 1A and 1B. To plot
  the data run the portion of the script `PNAS_figures_1abc_2ab.R`
  for the corresponding figure. Note it also produces data on just
  punishment probability as well (instead of demand shift)
  
  
Baseline_shiftDemand_and_probPunish_PNAS.R

  Runs the analysis to produce the data in Fig 1C. To plot
  the data run the portion of the script `PNAS_figures_1abc_2ab.R`
  for the corresponding figure.


Sensitivity_of_equilibrium_to_b2_PNAS.R

  Runs the analysis to produce the data in Fig 2A. To plot
  the data run the portion of the script `PNAS_figures_1abc_2ab.R`
  for the corresponding figure


Sensitivity_shiftDemand_and_probPunish_PNAS.R
  Runs the analysis to produce the data in Fig 2B. To plot
  the data run the portion of the script `PNAS_figures_1abc_2ab.R`
  for the corresponding figure

ExamplePhasePlots_PNAS.R

Generates the 9 typical phase plots in Fig S1 in the supplement
