#include <TCanvas.h>
#include <TH1.h>
#include <TH1F.h>
#include <TF1.h>
#include <TFile.h>
#include <TGraph.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TLegend.h>
#include <TLine.h>
#include <TLatex.h>
#include <TPolyLine.h>
#include <iostream>
#include <TMath.h>
#include <fstream.h>

using namespace std;

// Value of the nuclear density at r = 0, from "get_WS_rho_0()";
double rho_0 = 0.169;

// from the glauber_dau.inp file that I used                                                                          
double ws_radius = 6.38;
double ws_diff = 0.54;

double get_WS_rho_0();
double get_tau(double yrel, double MJpsi, double pt, double L);
double get_ccbar_sqrts(double yrel, double MJpsi, double pt);
double get_z_wt(double rT, double z, double zstep);
double get_path_thickness(double b, double zstart);
double get_thickness_integral(double b);
double get_sigma(double ccbar_roots, double tau, double sigma1=8.0, double r0=0.15, double vcc=1.85);



void time_spent_nucleus_glauber_vs_energy()
{
  // Calculation of sigma vs tau using Arleo model of color neutral expanding psi meson
  // Calculation is done for a range of parameters sigma1, r_0 and v_cc
  // The results are written to a text file "time_spent_glauber_vs_energy_out.txt"
  // The best parameter set is chosen from a chisq fit to data in "fit_sigma_vs_tau_vs_energy.C"
  // The fit is plotted vs data in "time_spent_in_nucleus.C"

  // This calculation is for PHENIX data at backward/mid rapidity
  // and lower energy data at mid rapidity

  // Assumptions:

  // 1) Assumes that the rT distribution at all energies can be approximated by the rT distribution of Au
  // from the dAu Glauber at 200 GeV. I have checked the mean tau values obtained from "time_spent_in_nucleus.C"
  // (which uses the average thickness scaled by A^(1/3)) and the mean tau values found by this macro 
  // "time_spent_nucleus_glauber_vs_energy.C" (which averages the tau values by sampling the PHENIX WS 
  // to get the z weighted contribution at each rT, then using the rT distribution from the Glauber to 
  // average over rT). The agreement is better than 2.5% at all energies and rapidities used here. OK.

  // 2) Assumes that the tau value can be calculated assuming that the psi has the mean pT at each point (rT,z).
  // I have checked this by comparing at 200 GeV, y=0, the mean tau and sigma values obtained in that way with those obtained 
  // by sampling the pT distribution from 200 GeV p+p. There is an option in "calculate_sigma_vs_roots.C" to 
  // sample the 200 GeV pT distribution instead of just using the mean pT. At y=0 at 200 GeV this makes a difference
  // of only 1% in tau and much less than that in the average sigma. Note that the mean pT value at 200 GeV is 
  // obtained from the pT distribution used in the check, so the check is purely of the assumption that the mean pT 
  // can be used to get the average tau and sigma accurately.  

  // 3) The mean pT values are obtained from (see note below for more details):
  // PHENIX d+Au and p+p pT distributions at 200 GeV for the 200 GeV Au target case
  // A HERA-B paper (Eur.Phys.J. C60 (2009) 525) containing a function describing the HERA-B data at 920 GeV for the 920 GeV W target case
  // The same HERA-B paper provides <pT^2> values for heavy masses from E866(800) and NA50(450)
  // The same HERA-B paper shows that the <pT^2> is linear in roots^2 for 450, 800 and 920 GeV
  // I derived the formula for the straight line from the 450 and 920 GeV <pT^2> values for heavy targets
  // I took the ratio of <pT>/sqrt(<pT^2 >) for the two cases 200 GeV and 920 GeV, and found them very similar
  // I took that ratio at 920 GeV and assumed it works for all energies below 920 GeV
  // The resulting <pT> values are given below in "mean_pT"

  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  ofstream fout("time_spent_glauber_vs_energy_out.txt");

  char *hname[5]={"hbAu0","hbAu1","hbAu2","hbAu3","hsum"};
  TH1F *hrT[5];
  
  TFile *fhist = new TFile("rT_distribution_dAu.root");
  if(!fhist)
    {
      cout << "Failed to open rT distribution file " << endl;
      return;
    }

  for(int icent=0;icent<5;icent++)
    {
      cout << " icent " << icent << endl;
      cout << "Getting histogram " << hname[icent] << endl;

      hrT[icent] = (TH1F *) fhist->Get(hname[icent]);
      if(!hrT[icent])
	{
	  cout << "Did not get histogram " << hname[icent] << endl;
	  return;
	}

    }

  // Assumed average mass over all charmonium states by Arleo
  double MJpsi = 3.4;
  double Mproton = 0.938;

  static const int NDATA = 19;  // number of energies at which to calculate sigma
  // list of beam energies and mass values for different energies. Counting down from RHIC.
  // RHIC(Au), HERA-B(W), E866(W), NA50-450(W), NA50-400(Pb), NA3(Pt), NA60(Pb)
  double Ebeam[NDATA] = {100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0,   // PHENIX
			 920.0,                                             // HERA-B
			 800.0, 800.0, 800.0, 800.0, 800.0, 800.0, 800.0,          // E866
			 450.0, 400.0,                                      // NA50
			 200.0,                                             // NA3
			 158.0};                                            // NA60
  double ycm[NDATA] = {-2.075,-1.825,-1.575,-1.325, -0.3, 0.0, 0.3, 
		       0.0,
		       -0.390, -0.115, 0.161, 0.433, 0.680, 0.896, 1.091,
		       0.0, 0.0, 
		       0.0, 
		       0.3};

  //=================
  // <pT> estimates:
  //=================
  // I know <pT^2> and <pT> for 
  // PHENIX from PPG125
  //       <pT^2> = 5.10 for Au
  //       from p+p case: <pT> = 0.872 * <pT^2>. So <pT> = 5.10 * 0.872 = 1.90 for Au        
  // From I. Abt et al. (the HERA-B collaboration), Eur.Phys.J. C60 (2009) 525-542 (arXiv:0812.0734) which has a parameterization for p-W 
  //       dN/dpT ~ pT(1 + 1/(beta-2) * pT^2/<pT^2>)^(-beta)
  //       with beta = 8.13, <pT^2> = 2.432
  //       I made a TF1 and checked, and this formula gives <pT> = 1.36 and <pT^2> = 2.433 for HERA-B W data
  // I know <pT^2> for NA50 (Pb) E866(Au) and HERA-B(W) because HERA-B has a plot in their paper showing it is linear with roots^2
  // I derived the formula from:
  //   Take HERA-B <pT^2> for W as 2.432 at S = 41.6^2 = 1731
  //   Take NA50 <pT^2> for Pb as 1.96 at S = 29.086^2 = 846
  //   Then we get for the slope and offset from HERAB and NA50:
  //     slope = (2.432 - 1.96)/(1731 - 846) = 0.472 / 885 = 0.00053333  
  //     offset = 2.432 - 1731*0.00053333 = 1.509 
  // Now, to get <pT> from <pT^2>:
  // for HERA-B, <pT> / sqrt(<pT^2>) = 1.36 / sqrt(2.432) = 0.872
  // for PHENIX (p+p), <pT> / sqrt(<pT^2>) = 1.76 / sqrt(4.38) = 0.841
  // very similar, probably OK to use HERA-B ratio at all lower energies
  // So for the lower energy data take <pT> = 0.872 * sqrt(<pT^2>)
  // This gives for <pT^2> and >pT> :
  //---------------------------------
  // dataset   roots       <pT^2>     <pT>
  // PHENIX    200 GeV     5.10       1.90  (pT scaled from measured <pT^2> using ratio of <pT> / <pT^2> from p+p)
  // HERA-B    41.6 GeV    2.432      1.36  both calculated from functional form in Eur.Phys.J. C60 (2009) 525 (arXiv:0812.0734)
  // E866      38.8 GeV    2.310      1.32  <pT^2> calculated from formula above, then <pT> = <pT^2> * 0.872
  // NA50-450  29.08 GeV   1.960      1.22  <pT^2> calculated from formula above, then <pT> = <pT^2> * 0.872
  // NA50-400  27.43 GeV   1.910      1.20  <pT^2> calculated from formula above, then <pT> = <pT^2> * 0.872
  // NA3       19.42 GeV   1.710      1.14  <pT^2> calculated from formula above, then <pT> = <pT^2> * 0.872
  // NA60      17.3 GeV    1.669      1.12  <pT^2> calculated from formula above, then <pT> = <pT^2> * 0.872
  //=======================

  double mean_pt[NDATA] = {1.90, 1.90, 1.90, 1.90, 1.90, 1.90, 1.90, 
			   1.36, 
			   1.32, 1.32, 1.32, 1.32, 1.32, 1.32, 1.32,
			   1.22, 1.20,
			   1.14, 
			   1.12};

  double ybeam[NDATA];
  for(int idata=0;idata<NDATA;idata++)
    {
      double pz = sqrt(pow(Ebeam[idata],2) - pow(Mproton,2));
      ybeam[idata] = 0.5 * log((Ebeam[idata]+pz)/(Ebeam[idata]-pz));
      cout << " idata " << idata
	   << " Ebeam " << Ebeam[idata]
	   << " pz " << pz
	   << " ybeam " << ybeam[idata]
	   << endl;
    }

  double roots[NDATA];



    
  // Specify number of steps
  static const int NL = 40;
  static const int NRT=40;
  static const int NPT=20;

  double rT[NRT];
  double thick[NRT];

  double rT_max = 10.0;
  double rT_step = rT_max / (double) (NRT);
 
  for(int irt=0;irt<NRT;irt++)
    {
      rT[irt] = rT_step/2.0 + (double) (irt) * rT_step;
      //cout << " rT " << rT[irt] << endl;
    }

  double sigma_avge[NDATA][NRT];
  double sigma_wt[NDATA][NRT];
  double tau_avge[NDATA][NRT];
  double L_avge[NDATA][NRT];
  double ccbar_roots_avge[NDATA][NRT];
    
  for(int idata=0;idata<NDATA;idata++)
    for(int irt=0;irt<NRT;irt++)
      {
	L_avge[idata][irt]=0.0;
	tau_avge[idata][irt]=0.0;
	sigma_avge[idata][irt]=0.0;
	sigma_wt[idata][irt]=0.0;
      }

  // Store these for later

  double ccbar_roots[NDATA];
  double yrel[NDATA];

  for(int idata=0;idata<NDATA;idata++)
    {
      // collider or fixed target
      if(idata < 7)
	yrel[idata] = ycm[idata] + ybeam[idata];
      else
	yrel[idata] = ycm[idata] + ybeam[idata]/2.0;

      ccbar_roots[idata] = get_ccbar_sqrts(yrel[idata],MJpsi,mean_pt[idata]);
    }
  
  // Use the rT distribution calculated from my d+Au glauber in "/Users/frawley/work/rdau_from_rpau/"
  // to get the rT weights
  
  double rT_weight[NRT];
  
  for(int irt=0;irt<NRT;irt++)
    {
      float binlo = rT[irt]-rT_step/2.0;
      float binhi = rT[irt]+rT_step/2.0;
      
      int ilo = hrT[4]->FindBin(binlo);
      int ihi = hrT[4]->FindBin(binhi);
      
      rT_weight[irt] = (double) hrT[4]->Integral(ilo,ihi);
    }
  
  double LrT_avge[NRT];
  double zwt_ccbar[NRT][NL];
  double path[NRT][NL];
  
  double zmax = 10.0;    // distribution is zero beyond 10 fm radius
  double zstep = zmax*2.0 / ((double) (NL+1) );
  
  double z[NL];
  for(int iL=0;iL<NL;iL++)
    {
      z[iL] = -zmax + zstep * (double) (iL+1);
    }
  
  for(int irt=0;irt<NRT;irt++)
    {
      // We want to get the average distance before a produced ccbar pair hits a nucleon
      // From: M.C. Abreu et al., PLB 410 (1997) 337
      // L(b) = 1/rho_0 * Integral{rho * dz}
      // where rho_0 is taken by them as 0.17 /fm^3 (note, I get 0.1691 / fm^3 and use that to normalize rho)

      // To get the average nuclear length at rT, use the integral over all z
      LrT_avge[irt] = get_thickness_integral(rT[irt]) / rho_0;

      for(int iL=0;iL<NL-1;iL++)
	{
	  // Get the weight for this z value from the WS z distribution at rT
	  zwt_ccbar[irt][iL] = get_z_wt(rT[irt],z[iL],zstep);
	  
	  // Get the average value of L seen by a ccbar pair produced at z	  
	  path[irt][iL] = get_path_thickness(rT[irt],z[iL]) / rho_0;
	}

    }

  //===========================
  // OK, here we go
  //===========================
  


  int NSG = 30;
  double sigstep = 0.2;
  double siglow = 5.0;

  // Reduced the step size to 0.05 from 0.2 when adding a fit to the psi'/psi ratio on 3/16/13
  int NVCC = 20;
  double vcclow = 0.8;
  double vccstep = 0.05;

  int NR0 = 10;
  double r0low = 0.00;
  double r0step = 0.02;

  // Loop over the trial r0 parameters
  for(int ir=0;ir<NR0;ir++)
    {
      double r0 = r0low + r0step * (double) ir;
  
      cout << "Trial value of r0 = " << r0 << endl;

      // Loop over the trial vcc parameters
      for(int iv=0;iv<NVCC;iv++)
	{
	  double vcc = vcclow + vccstep * (double) iv;

	  cout << "  r0 = " << r0 << " trial value of vcc = " << vcc << endl;

	  // Loop over the trial sigma1 parameters
	  for(int isig=0;isig<NSG;isig++)
	    {
	      // specify the model parameters
	      double sigma1 = siglow + sigstep * (double) isig;

	      cout << "    r0 = " << r0 << " vcc = " << vcc << " trial value of sigma1 = " << sigma1 << endl;

	      fout << sigma1 << " "
		   << r0 << " " 
		   << vcc << endl;

	      for(int idata=0;idata<NDATA;idata++)
		{
		  //cout << "Begin y = " << y[idata] << endl;
	  
		  // Now loop over the r_T values and calculate the length, and then the sigma
	  
		  for(int irt=0;irt<NRT;irt++)
		    {
		      // Since sigma is not linear in tau, using the average path length is not good enough
		      // So we subdivide the z distribution at rT into NL steps
		      // calculate the production probability at each step
		      // calculate the remaining path length at each step
		      // calculate tau from the remaining path length
		      // Calculate sigma from tau
		      // Average sigma over all of the steps in z at this rT
	      
		      // These will be the averages over the line
		      double tau_line = 0.0;
		      double sigma_line = 0.0;
		      double LrT_line = 0.0;
		      double wt_line = 0.0;
	      
		      // Going only to NL-1 keeps the path length from going to zero
		      for(int iL=0;iL<NL-1;iL++)
			{

			  // Weight path by the ccbar production probability and add it to the average over all iL values for this irt value
			  // This gives us the average over the whole line of z values at rT

			  // We use mean pT
    
			  double tau = get_tau(yrel[idata], MJpsi, mean_pt[idata], path[irt][iL]); 
			  tau_line += tau*zwt_ccbar[irt][iL];
			  // Use the Arleo form for the tau dependence of the psi radius
			  double sigma = get_sigma(ccbar_roots[idata], tau, sigma1, r0, vcc);
			  sigma_line += sigma*zwt_ccbar[irt][iL];
			  LrT_line += path[irt][iL]*zwt_ccbar[irt][iL];
			  wt_line += zwt_ccbar[irt][iL];
			}
		      // averages for this rT

		      tau_line = tau_line/wt_line;
		      sigma_line = sigma_line/wt_line;
		      LrT_line = LrT_line/wt_line;	      
		      L_avge[idata][irt] = LrT_line;

		      tau_avge[idata][irt] = tau_line;
		      sigma_avge[idata][irt] = sigma_line;
		      sigma_wt[idata][irt] = wt_line;
		    }
		}
  
	      // Get the average sigma over all rT at each rapidity
  
	      double L[NDATA];
	      double L_avge_tot[NDATA];
	      double tau[NDATA];
	      double sigma[NDATA];
	      for(int idata=0;idata<NDATA;idata++)
		{
		  // integrate the modification over the rT weights for MB
	    
		  L[idata]=0.0;
		  tau[idata]=0.0;
		  sigma[idata]=0.0;
		  double wt = 0.0;
		  L_avge_tot[idata]=0.0;
 
		  for(int irt=0;irt<NRT;irt++)
		    {
		      L_avge_tot[idata] += LrT_avge[irt]*rT_weight[irt];
		      L[idata] += L_avge[idata][irt]*rT_weight[irt];
		      tau[idata] += tau_avge[idata][irt]*rT_weight[irt];
		      sigma[idata] += sigma_avge[idata][irt]*rT_weight[irt];
		      wt += rT_weight[irt];
		    }
		  sigma[idata] = sigma[idata] / wt;      
		  tau[idata] = tau[idata] / wt;      
		  L[idata] = L[idata] / wt;      
		  L_avge_tot[idata] = L_avge_tot[idata] / wt;

		  fout << idata << " "
		       << tau[idata] << " "
		       << sigma[idata] << endl;

		  cout << "idata " << idata << " sigma " << sigma[idata] << " tau " << tau[idata] << " L " << L[idata]  << endl; 

		}
	      cout << endl << "Average path length for MB Au is 1/2 of the total length = " << L_avge_tot[0]/2.0 << endl;

	    }  // end loop over sigma1
	}  // end loop over vcc
    } // end loop over r0
}

double get_WS_rho_0()
{
  // The WS distribution  
  TF1 *ws = new TF1("f1","1.0/(1+exp((x-[0])/[1]))",0.,20.);

  ws->SetParameters(ws_radius,ws_diff);

  // calculate the volume integral of the Woods Saxon distribution
  // The volume element is:
  // dV = dr * d(theta) * {Sin(theta) * d(phi)} = r^2 Sin(theta) dr d(theta) d(phi)
  // Then:
  // Integral = Int{f(r) r^2 dr}Int{sin(theta) d(theta)}Int{d(phi)}
  // Int{Sin(theta)d(theta)} = [-cos(theta)](0,pi) = 1.0 -(-1) = 2.0
  // Int{d(phi)} = [phi](0,2pi) = 2 pi
  // Int{r^2/(1+exp((r-b)/a)) dr} - do it numerically!

  double rstep = 0.1;
  double ws_integral = 0.0;

  for(int ir=0;ir<120;ir++)
    {
      double r = rstep * (double) ir;
      ws_integral += pow(r,2) / (1.0+exp((r-ws_radius)/ws_diff)) * rstep;
    }

  //  Now multiply by the theta and phi terms, which equal (2 * 2 * pi)

  ws_integral = 2.0 * 2.0 * 3.14159 * ws_integral;

  // we need rho_0 properly normalized so that it makes the integral equal to 197
  double rho0 = 197.0 / ws_integral;

  // The product should be 197 if all is done correctly!
  //cout << " rho0 " << rho0 << " WS integral = " << ws_integral << " product " << rho * ws_integral << endl;

  return rho0;
}

double get_thickness_integral(double b)
{
  // The WS distribution  
                                                                                              
  TF1 *ws = new TF1("f1","1.0/(1+exp((x-[0])/[1]))",0.,20.);

  // from the glauber_dau.inp file that I used  
                                                                        
  ws->SetParameters(ws_radius,ws_diff);

  // Now calculate the longitudinal density at the impact parameter b 
                                                  
  // Go from -15 fm to +15 fm  
                                                                                         
  int nz = 60;
  double dz = 0.25;
  double lint = 0.0;
  for(int iz=-nz; iz<nz;iz++)
    {
      double z = (double) iz * dz;
      double r = sqrt(pow(z,2)+pow(b,2));
      double rho = ws->Eval(r);

      lint += rho*dz;
    }

  // Number of nucleons / fm^3 in Au (mass 197) is roughly  197 / {(4/3) pi R_WS^3} = 0.181 / fm^3
  // Note that PLB 410 (1997) 337 says the standard is 0.17 / fm^3
  // Ramona says rho_0 = 0.1693 /fm^3, so the same

  // \rho_0 = A/(\int d^3r \rho_A(r)) = A/(\int d2r_T dz \rho_A(r_T,z))

  // Normalize to nucleons/ fm^2 by multiplying by rho_0
  // rho_0 = get_WS_rho_0();
  lint = lint * rho_0;

  return lint;
}

double get_path_thickness(double b, double zstart)
{
  // The WS distribution  
                                                                                              
  TF1 *ws = new TF1("f1","1.0/(1+exp((x-[0])/[1]))",0.,20.);

  ws->SetParameters(ws_radius,ws_diff);

  // Now calculate the longitudinal density at the impact parameter b 
                                                  
  // Go from zstart fm to +15 fm  
                                                                                         
  int nz = 60;
  double dz = 0.25;
  double lint = 0.0;
  for(int iz=-nz; iz<nz;iz++)
    {
      double z = (double) iz * dz;

      if(z < zstart)
	continue;

      double r = sqrt(pow(z,2)+pow(b,2));
      double rho = ws->Eval(r);

      lint += rho*dz;
    }

  // Number of nucleons / fm^3 in Au (mass 197) is roughly  197 / {(4/3) pi R_WS^3} = 0.181 / fm^3
  // Note that PLB 410 (1997) 337 says the standard is 0.17 / fm^3
  // Ramona says rho_0 = 0.1693 /fm^3, so the same:
  // \rho_0 = A/(\int d^3r \rho_A(r)) = A/(\int d2r_T dz \rho_A(r_T,z))

  // Normalize to nucleons/ fm^2 by multiplying by rho_0
  // double rho_0 = get_WS_rho_0();
  lint = lint * rho_0;

  return lint;
}

double get_z_wt(double rT, double z, double zstep)
{
  double zmax = 10.0;

  TF1 *zdist = new TF1("zdist","1.0/(1+exp((sqrt(pow(x,2)+pow([0],2))-[1])/[2]))",-zmax,zmax);
  zdist->SetParameters(rT, ws_radius, ws_diff);

  double zlo = z - zstep/2.0;
  double zhi = z + zstep/2.0;

  //double rho_0 = get_WS_rho_0();
  double zwt = rho_0 * zdist->Integral(zlo,zhi);

  /* 
 cout << " rT " << rT 
       << " z " << z
       << " zstep " << zstep
       << " zlo " << zlo
       << " zhi " << zhi
       << " zwt " << zwt
       << endl;
  */

  return zwt;


}

double get_ccbar_sqrts(double yrel, double MJpsi, double pt)
{
  // write the four vectors p1 and p2 in terms of E and the three vectors p1' and p2'
  // p1+p2 is Lorentz invariant, its magnitude is s
  // in any frame:
  // (p1+p2) =  (E1+E2,p1'+p2') 
  // s = (E1+E2)^2 - (p1+p2)^2
  // In the CM, by definition p1'+p2' = 0, so s = (E1+E2)^2 
  // In the lab, p2' = 0, so:
  // s = (E1+m2)^2 -(p1'+0)^2


  // we need Ebeam and pbeam for the ccbar
  double Mtarget = 0.938;
  double Ecm = sqrt(pow(pt,2) + pow(MJpsi,2));
  double Ebeam = Ecm * sqrt(1+pow(sinh(yrel),2));
  double pbeam = sqrt(pow(Ebeam,2) - pow(MJpsi,2));
  double s = pow(Ebeam + Mtarget,2) - pow(pbeam,2);
  double sqrts = sqrt(s);

  return sqrts;
}

double get_tau(double yrel, double MJpsi, double pt, double L)
{
  // E, m and p are in GeV, GeV/c^2 and GeV/c repectively
  // Generally, if we use those units for E, p, m and v, we consider c=1 and ignore it
  // p = m * v 
  // E = m * c^2 
  // So p/E = v/c^2 = v/c * 1/c is in units of (GeV/c) / GeV = 1/c 
  // ----- ie. p/E is v in units of 1/c, ie. p/E = beta

  // We need the J/psi Lorentz gamma factor in the target frame, which is p/E in the target frame:
  // p_z = m_T sinh(y)
  // m_T = sqrt(m_0^2+p_T^2)
  double Ecm = sqrt(pow(pt,2) + pow(MJpsi,2));
  double pz = Ecm * sinh(yrel);

  double p = sqrt(pow(pt,2)+pow(pz,2));
  // E = m_T cosh(y) = m_T sqrt(1+sinh^2(y)))
  double E = Ecm * sqrt(1+pow(sinh(yrel),2));
  // p/E is v in units of 1/c, ie. p/E = beta
  double beta = p/E;
  
  // The numerator is pz, the denominator is energy
  double beta_z = pz/E;

  // We need the J/psi Lorentz gamma factor in the target frame calculated with beta, not beta_z
  double gamma = 1.0/sqrt(1-pow(beta,2));

  // proper time spent in nucleus by J/psi:
  // This is (v in units of 1/c) times (L in units of fm), so it has units of fm/c = time
  // If beta = v/c = 0.99 and L = 4.3 fm, then tau = 4.26 fm/c before correcting for gamma
  // Have to correct apparent distance traveled in nucleus by J/psi for the gamma of J/psi in nuclear rest frame to get proper time
  // gamma would be 7.09 if beta = 0.99, so in our example tau becomes 0.6 fm/c due to time dilation
  double tau = beta_z *L/gamma;
  
  return tau;
}

double get_sigma(double ccbar_roots, double tau, double sigma1, double r0, double vcc)
{
  // There are three states - J/psi, psi' and chi_C
  // These are the numbers used by Arleo
  double rpsi[3] = {0.43,0.87,0.67};

  // These are from PPG104
  double fdfrac[3] = {0.58,0.1,0.32};

  double sigma = 0.0;
  for(int istate=0;istate<3;istate++)
    {
      // This is the model of Arleo PRC61 where rccbar is linear in tau
      double rccbar = r0 + vcc*tau;

      if(rccbar > rpsi[istate])
	rccbar = rpsi[istate];

      // From Arleo PRC61 
      // The reference radius is the J/psi one, so that the cross section is larger for the fully formed excited states
      double sigma_state = sigma1* pow(ccbar_roots/10.0,0.4) * pow(rccbar/rpsi[0],2);
      sigma += sigma_state*fdfrac[istate];

      /*
      cout << " sigma1 " << sigma1 
	   << " ccbar_roots " << ccbar_roots
	   << " istate " << istate
	   << " rpsi " << rpsi[istate]
	   << " rccbar " << rccbar
	   << " fdfrac " << fdfrac[istate]
	   << " sigma_state " << sigma_state 
	   << " sigma " << sigma 
	   << endl;
      */
    }
  //cout << "           average over states is sigma = " << sigma << endl;
      
  return sigma;      
}
