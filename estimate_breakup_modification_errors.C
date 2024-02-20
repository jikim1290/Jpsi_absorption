// macro uses calculate breakup modification with breakup
// model parameters varied within the parameter uncertainties.
// it is intended to be called by condor.

// Uncertainties:
// from the macro used in the paper: /phenix/hhj/frawley/old_time_spent_nucleus/fit_sigma_vs_tau_vs_energy.C
// Uncertainties: find by varying each parameter up and down and finding where the chisq increased by 1
// sigma1 
//    sigma1  6.8  17.78
//    sigma1  7.0  15.91
//    sigma1  7.2  15.10
//    sigma1  7.4  15.34
//    sigma1  7.6  16.64
// from fit_parabola.C, mean = 7.25 and sigma = 0.278
// r0: 
//    r0 0.14 ch2 19.86
//    r0 0.16 ch2 15.10
//    r0 0.18 ch2  21.03
// from fit_parabola.C, mean = 0.159 and sigma = 0.0087 
// vcc:
//    vcc 0.95 ch2 16.22
//    vcc 1.00 ch2 15.10
//    vcc 1.05 ch2 15.64
//    vcc 1.10 ch2 17.70
// from fit_parabola.C, mean = 1.009 and sigma = 0.056

// Procedure for getting uncertainty in modification for each collision system:
//   Randomly sample each uncertainty distribution for 10 runs
//   Get RMS of modification at each of 4 rapidities

#include <TF1.h>
#include <TRandom3.h>

#include "calculate_breakup_modification.C"

void estimate_breakup_modification_errors(int process = 0)
{
  TRandom3 *ran = new TRandom3(0);  // want a new seed every time

  // errors come from fit_parabola.C
  double sigma1_mean = 7.2;   double sigma1_err = 0.258;
  double r0_mean = 0.16; double r0_err = 0.0087;
  double vcc_mean = 1.0; double vcc_err = 0.056;

  cout << "Input parameters:" << endl;
  cout << " sigma1 " << sigma1 << " error " << sigma1_err << endl;
  cout << " r0 " << r0 << " error " << r0_err << endl;
  cout << " vcc " << vcc << " error " << vcc_err << endl << endl;

  double sigma1_trial = ran->Gaus(sigma1_mean, sigma1_err);
  cout << " calculated sigma_1 = " <<  sigma1_trial << endl;

  double r0_trial = ran->Gaus(r0_mean, r0_err);
  cout << " calculated r0 = " <<  r0_trial << endl;

  double vcc_trial = ran->Gaus(vcc_mean, vcc_err);
  cout << " calculated vcc = " <<  vcc_trial << endl;


  calculate_breakup_modification(sigma1_trial, r0_trial, vcc_trial, process);

}
