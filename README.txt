Use of the files in this directory:
----------------------------------------
All macros are intended to be run compiled (i.e. with the ".C+" option).

"time_spent_nucleus_glauber_vs_energy.C"
The original macro used to generate the cross sections from the Arleo model.
NO LONGER USED (methods are now in calculate_sigma_breakup.C, see below), but 
included here because it generated the file:
"time_spent_glauber_vs_energy_out.txt" used to get the best fit to the worlds breakup 
cross section data that is used extensively below.

"fit_sigma_vs_tau_vs_energy.C"
Reads the data file: 
"time_spent_glauber_vs_energy_out.txt"
written by 
"time_spent_nucleus_glauber.C"
The file contains the results of a calculation of chisq from a comparison of calculated breakup cross sections with the worlds shadowing corrected 
breakup cross section data (McGlinchey et. al. PRC ). Results are for a systematic scan of breakup model fit parameters (sigma1, r0, vcc). 
Finds the minimum chisq value and prints out the best fit parameters, as well as chisq values 
around the optimum for each parameter, for error estimation.

"fit_parabola.C"
Takes the results from "fit_sigma_vs_tau_vs_energy.C" 
Fits a parabola to the chisq results as each parameter is varied around its optimum, to get sigma for that parameter.
Sigma is the parameter value where the chisq increases by 1.0.

"calculate_sigma_breakup.C"
A set of methods used to calculate the breakup cross section from the model, given
the collision energy, rapidity of the Jpsi, pT and rT values.
The calculation is done by calling the method:
"return_sigma_breakup()"
(See later for more)
This macro is "#included" in "calculate_breakup_modifications.C", to provide access to the methods. 

"calculate_breakup_modification.C"
contains "#include calculate_sigma_breakup.C"
The main method "calculate_breakup_modification(double sigma1 = 7.2, double r0 = 0.16, double vcc = 1.0, int process = 0)"
calculates the modification in p+Al (#define PAL) or p+Au (#define PAU) collisions using the supplied model parameters, 
or defaults to the  best fit parameters.
Uses the appropriate rT distributions from the PHENIX Glauber (from Jamie Nagle).
Calls the "return_sigma_breakup(Ebeam, mstate, y[irap], pt[ipt], mtarget, rT[irt])" method in "calculate_sigma_breakup.C" to 
get the breakup cross section as a function of rT and Jpsi pT. 
Then calculates the breakup modification vs pT averaged over the rT distribution for each centrality bin.
Also calculates the modification averaged over pT.

"estimate_breakup_modification_errors.C"
contains "#include calculate_breakup_modification.C"
Used to estimate the errors for the breakup modification from the model parameter uncertainties obtained in "fit_parabola.C"
s Randomly varies the model fit parameters assuming a gaussian distribution with the width of the fit parameter uncertainties.
Calls "calculate_breakup_modification()" to get the modifications for each model parameter combination.
"calculate_breakup_modification()" writes the modifications to a file using the process number to name the file.
Repeats this for 1000 parameter combinations.

i"runit_bu_estimate_errors"
Script to run "estimate_breakup_modification_errors.C".
Usage is
runit_bu_estimate_errors <process number>
Called by "submit_condor_estimate"

"submit_condor_estimate"
Calls "runit_bu_estimate_errors" 1000 times.

"process_bu_mod_error_estimate_output.C" 
Reads in the condor output files written by 
"estimate_breakup_modification_errors.C"
and makes a histogram for every rapidity and centrality bin (only pT integrated results have been made so far).

Update on 7/19/2021:
---------------------------
Added the method "get_sigma_psi2s()" to "calculate_sigma_breakup.C". 
If "jpsi = true" is set in "calculate_breakup_modification.C" 
-- The mass is set to 3.4 GeV (average of three states)
-- Passing that mass to "return_sigma_breakup()" causes "get_sigma()" to be called and that returns the psi(1S) sigma
If "jpsi = false" is set
.-- The mass is et to 3.7 GeV (the psi(2S) mass)
-- Passing that mass to "return_sigma_breakup()" causes "get_sigma_psi2s()" to be called and that returns the ps(2S) sigma.
The psi(2S) is assumed to have the same pT and rapidity distribution as the psi(1S). 
-- Should scale pT distribution with pT/M? see arXiv:2006.15446 . Look into it.
-- What about dN/dy?

I also added a plot of the arm integrated modification to "process_bu_mod_error_estimate.C".
