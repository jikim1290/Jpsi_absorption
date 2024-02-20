#include <fstream>
#include <iostream>
#include <TMath.h>

void fit_sigma_vs_tau_vs_energy()
{

  bool fix_vcc = false;    // set this to true to fix vcc to the value determined from fit_psip_psi_ratio.C, false to fit inclusive data only
  // This value was arrived at by iterating between this macro and fit_psip_psi_ratio.C
  double vcc_fixed_value = 1.45;

  // Reads calculations for a range of values of sigma1, r_0 and v_cc and written to this file by "time_spent_nucleus_glauber_vs_energy.C"
  ifstream fin("time_spent_glauber_vs_energy_out.txt");

  // Performs a chisq fit and writes the best values to arrays ready to paste into "time_spent_nucleus.C" for plotting

  static const int NDATA = 19;

  // Put all of the data into a single set of arrays
  // These are the data values that determine the fit results, so they have been carefully checked on 2/23/13 against the sources

  double sigma_exp[NDATA] = {6.7, 5.3, 4.8, 4.6, 2.6, 1.7, 0.9,                  // PHENIX - checked 2/23/13 dAu_chisq_contours/plot_sigma_comparison.C 
			     4.35,                                               // HERA-B - checked Lourenco 2/23/13 Tab. 6 EKS98
			     4.67, 5.39, 4.98, 4.36, 3.95, 3.13, 2.76,           // E866  - checked Lourenco 2/23/13, Tab. 5 EKS98 and Fig. 7 W/Be
			     7.25, 7.24,                                         // NA50-450, NA50-400 - checked Lourenco 2/23/13 Tab. 6 EKS98
			     7.82,                                               // NA3 - checked Lourenco 2/23/13 Tab. 6 EKS98
			     9.3 };                                              // NA60 - checked arXiv:0907.0504 Sec. 2 EKS98

  double sigma_exp_up[NDATA] = {1.3, 1.1, 1.1, 1.1, 1.2, 1.1, 1.0,               // PHENIX - checked 2/23/13 dAu_chisq_contours/plot_sigma_comparison.C 
				1.37,                                            // HERA-B - checked Lourenco 2/23/13 Tab. 6 EKS98
				0.92, 0.82, 0.78, 0.76, 0.75, 0.75, 0.45,        // - checked Lourenco 2/23/13, Tab. 5 EKS98 and Fig. 7 W/Be 
				1.03, 0.73,                                      // NA50-450, NA50-400 - checked Lourenco 2/23/13 Tab. 6 EKS98
				0.90,                                            // NA3 - checked Lourenco 2/23/13 Tab. 6 EKS98
				0.99 };                                           // NA60 - checked arXiv:0907.0504 Sec. 2 EKS98 (all errors combined)

  double sigma_exp_dn[NDATA] = {0.8, 0.8, 0.8, 0.9, 1.2, 1.2, 1.3,               // PHENIX - checked 2/23/13 dAu_chisq_contours/plot_sigma_comparison.C 
				1.03,                                            // HERA-B - checked Lourenco 2/23/13 Tab. 6 EKS98
				0.85, 0.76, 0.73, 0.71, 0.71, 0.71, 0.45,        // - checked Lourenco 2/23/13, Tab. 5 EKS98 and Fig. 7 W/Be 
				1.03, 0.73,                                      // NA50-450, NA50-400 - checked Lourenco 2/23/13 Tab. 6 EKS98
				0.84,                                            // NA3 - checked Lourenco 2/23/13 Tab. 6 EKS98
				0.99 };                                           // NA60 - checked arXiv:0907.0504 Sec. 2 EKS98 (all errors combined)


  int NR0 = 5;
  int NCC = 20;
  int NSG = 30;

  int idata;
  double tau[NDATA];
  double sigma[NDATA];

  double chisq_min = 100000;
  double sigma1_min;
  double r0_min;
  double vcc_min;
  double sigma_min[NDATA];

  double sigma1,r0,vcc;

  /*
  // these are for looking at combinations other than the minimum
  // increasing vcc makes the psi'/psi ratio more suppressed
  // then decreasing r0 partly compensates in the inclusive fit
  //==============================================================
  double sigma1_test = 7.2;
  double r0_test = 0.16;
  double vcc_test = 1.0;
  double chisq_min_test = 100000;
  double sigma1_min_test;
  double r0_min_test;
  double vcc_min_test;
  double sigma_min_test[NDATA];
  */
  
  int icomb = 0;
  do 
    {
      icomb++;

      // starting a combination of fit parameters
      fin >> sigma1 >> r0 >> vcc;

      //cout << "Combination " << icomb << " Testing sigma1 " << sigma1 << " r0 " << r0 << " vcc " << vcc << "   ";

      double chisq = 0.0;

      for(int idata=0;idata<NDATA;idata++)
	{
	  
	  fin >> idata >> tau[idata] >> sigma[idata];

	  // calculate the contribution to chisq at this idata
	  // asymmetric error bars for most data points, calculate up and down contributions separately
	  double chisq_temp = pow( (sigma_exp[idata] - sigma[idata]) / sigma_exp_up[idata],2);    	    
	  chisq_temp += pow( (sigma_exp[idata] - sigma[idata]) / sigma_exp_dn[idata],2);
	  // average the up and down chisq values

	  chisq += chisq_temp/2.0;

	  //cout << "idata " << idata << " tau " << tau[idata] << " sigma " << sigma[idata] << " chisq " << chisq << endl;
	}
	//cout << " chisq " << chisq << endl;

      if(sigma1 == 7.2 && r0 == 0.16 && chisq < (15.01 + 3.0) )
      //if(sigma1 == 7.2 && vcc == 1.0 && r0 > 0.1 && r0 < 0.22 )
      //if(vcc == 1.0 && r0 == 0.16 && chisq < (15.01 + 3.0) )
	cout << "sigma1 " << sigma1 << " r0 " << r0 << " vcc " << vcc << " ch2 " << chisq << endl;


      /*
      if(r0 == r0_test && vcc == vcc_test && sigma1 == sigma1_test)
	{
	  chisq_min_test = chisq;
	  sigma1_min_test = sigma1;
	  r0_min_test = r0;
	  vcc_min_test = vcc;
	  for(int idata=0;idata<NDATA;idata++)
	    {
	      sigma_min_test[idata] = sigma[idata];
	    }
	}
      */

	
      if(chisq < chisq_min)
	{
	  if(fix_vcc && vcc != vcc_fixed_value)
	    continue;

	  chisq_min = chisq;
	  sigma1_min = sigma1;
	  r0_min = r0;
	  vcc_min = vcc;
	  for(int idata=0;idata<NDATA;idata++)
	    {
	      sigma_min[idata] = sigma[idata];
	    }
	}
    

    } while(fin.good());

  cout << " chisq_min " << chisq_min
       << " sigma1_min " << sigma1_min
       << " r0_min " << r0_min
       << " vcc_min " << vcc_min
       << endl;
  
  for(int idata=0;idata<NDATA;idata++)
    {
      cout << " idata " << idata << " sigma " << sigma_min[idata] << endl;
    }
  
  double NDF = (double) NDATA - 3.0;
  cout << " chisq/NDF = " << chisq_min/NDF << endl;

  // sort the best fit values to get increasing tau values

  int gotone;
  do
    {
      gotone = 1;
      for(int idata = 0;idata<NDATA-1;idata++)
	{
	  if(tau[idata] < tau[idata+1])
	    {
	      double tau_temp = tau[idata];
	      double sigma_min_temp = sigma_min[idata];
	      //double sigma_min_test_temp = sigma_min_test[idata];
	      
	      tau[idata] = tau[idata+1];
	      sigma_min[idata]=sigma_min[idata+1];
	      //sigma_min_test[idata]=sigma_min_test[idata+1];
	      tau[idata+1] = tau_temp;
	      sigma_min[idata+1]=sigma_min_temp;
	      //sigma_min_test[idata+1]=sigma_min_test_temp;

	      gotone = 0;
	    }
	}
    } while(gotone == 0);

  
  // Capture the best fit as arrays
  cout << endl << "Best fit is sigma1 = " << sigma1_min << " r0 = " << r0_min << " vcc = " << vcc_min << " chisq_min " << chisq_min << endl;
  cout << "double sigma_fitted1[" << NDATA << "] = {";
    for (int idata = 0;idata<NDATA-1;idata++)
      {
	cout << sigma_min[idata] << ", ";
      } 
    cout << sigma_min[NDATA-1] << "};" << endl;

  cout << "double tau_fitted1[" << NDATA << "] = {";
    for (int idata = 0;idata<NDATA-1;idata++)
      {
	cout << tau[idata] << ", ";
      } 
    cout << tau[NDATA-1] << "};" << endl;

    /*
  // Capture the test case as arrays
    cout << endl << "Test case for sigma1 = " << sigma1_test << " r0 = " << r0_test << " vcc = " << vcc_test << " chisq_min " << chisq_min_test << endl;
    cout << "double sigma_fitted1[" << NDATA << "] = {";
    for (int idata = 0;idata<NDATA-1;idata++)
      {
	cout << sigma_min_test[idata] << ", ";
      } 
    cout << sigma_min_test[NDATA-1] << "};" << endl;

  cout << "double tau_fitted1[" << NDATA << "] = {";
    for (int idata = 0;idata<NDATA-1;idata++)
      {
	cout << tau[idata] << ", ";
      } 
    cout << tau[NDATA-1] << "};" << endl;
    */

}
