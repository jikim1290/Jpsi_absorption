#include <TCanvas.h>
#include <TH1.h>
#include <TH1F.h>
#include <TF1.h>
#include <TFile.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TLegend.h>
#include <TLine.h>
#include <TLatex.h>
#include <TPolyLine.h>
#include <iostream>
#include <fstream>


using namespace std;

double get_thickness_integral(double b);
//double get_modification(int npow, double A,double thick,double rg);
double get_WS_rho_0();

// Note: ws_radius, ws_diff, rho_0 and mtarget are defined in the file included below 

// we will need to call functions from this macro
#include "calculate_sigma_breakup.C"

// only one of these at a a time!
#define PAU
//#define PAL
//#define HEAU


// 1: System with AuAu
// 2: System with PbPb
// 4: Upsilon


// set parameters here so macro can be called by condor with different parameters, for error estimation
void calculate_breakup_modification(double sigma1 = 7.2, double r0 = 0.16, double vcc = 1.0, int process = 0) 
{
 cout << "ChangeOptions in breakup_mod: " << ChangeOptions << endl;
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  // choose the state (jpsi = true for J/psi,  jpsi = false for psi(2S)
  bool jpsi = true;

  // Set the breakup parameters in "calculate_sigma_breakup")
  set_breakup_parameters(sigma1, r0, vcc);
  if( ChangeOptions&4 ){
	set_breakup_parameters(sigma1, 0, 0.28);
  }


  cout << " Breakup model parameters set to sigma 1 = " << sigma1 << " r0 " << r0 << " vcc " << vcc << endl;

  // initially define target parameter, will be set to correct value later
  int targ_index = 0;  // 0 for Au, 1 for Al
  
#ifdef PAU
  // has to be 0 for Au target
  targ_index = 0;
#endif
  
#ifdef PAL
  // has to be 1 for Al target
  targ_index = 1;
#endif

#ifdef HEAU
  targ_index = 0;
#endif
  
  static const int NRT=100;
  static const int NPT=40;
  double ptstep = 0.5;
  
  // set up the rapidity 
  // currently assumed to be the same for psi(1S) and psi(2S)
  static const int NRAP = 9;
  double y[NRAP] = {-2.075, -1.825 , -1.575, -1.325, 0, 1.325, 1.575, 1.825, 2.075};   
  double BdNdy[NRAP] = {0.325, 0.515, 0.68, 0.90, 1, 0.66, 0.53, 0.41, 0.30};
  
  static const int NARMS = 2;
  double ylow[NARMS]={-2.2, 1.2};
  double yhigh[NARMS]={-1.2, 2.2};
  
  // the default in calculate_sigma_breakup.C is a Au target
  double ws_radius_list[2] = {6.38, 3.34};
  double ws_diff_list[2] = {0.54, 0.58};
  double rho_0_list[2] = {0.169, 0.133};
  // Au or Al target - projectile does not matter for this
  double mtarget_list[2] = {197, 27};
  
  double use_ws_radius = ws_radius_list[targ_index];
  double use_ws_diff = ws_diff_list[targ_index];
  double use_rho_0 = rho_0_list[targ_index];
  double use_mtarget = mtarget_list[targ_index];
  
  set_WS_parameters(use_ws_radius, use_ws_diff, use_rho_0, use_mtarget);
  double check_rho_0 = get_WS_rho_0();
  cout << "WS parameters set to mtarget " << mtarget << " radius " << ws_radius << " diffusion " << ws_diff << " rho_0 = " << rho_0 << " check rho_0 " << check_rho_0 << endl;
  
  // The rT distribution should be obtained from a Glauber model of the specific collision system, with the trigger efficiency handled properly
  
#ifdef PAU

  cout << "Using measured pAu rT distributions" << endl;
  
  // get the pAu rT distributions from a file made by Jamie
  //===================================
  /*
    for p+Au with 9 centrality bins (the standard 9 we have used, i.e. 0-5, 5-10, 10-20, 20-30, 30-40, 40-50, 50-60, 60-70, 70-84). 
  */
  TFile* fhist;
  fhist = new TFile("rT_distributions_pau.root");
  if( ChangeOptions&1 ){
	fhist = new TFile("rT_distributions_auau.root");
  } else if( ChangeOptions&2 ){
	fhist = new TFile("rT_distributions_pbpb.root");
  }
  if(!fhist)
    {
      cout << "Failed to open rT distribution file " << endl;
      return;
    }

  static const int NCENT_in = 9;
  char hname[NCENT_in][100]={"himpact_cent0","himpact_cent1","himpact_cent2","himpact_cent3","himpact_cent4", "himpact_cent5", "himpact_cent6", "himpact_cent7", "himpact_cent8"};
  TH1D *hrT_in[NCENT_in];
  for(int icent=0;icent<NCENT_in;icent++)
    {
      cout << "Getting pAu histogram " << hname[icent] << endl;
      
      hrT_in[icent] = (TH1D *) fhist->Get(hname[icent]);
      if(!hrT_in[icent])
	{
	  cout << "Did not get histogram " << hname[icent] << endl;
	  return;
	}
      cout << "  icent " << icent << " hrT_in integral " << hrT_in[icent]->Integral() << endl;
    }

  // combine glauber centralities to get the experimental ones
  //========================================
  static const int NCENT = 6;  // includes MB and 0-20%
  TH1D *hrT[NCENT];
  hrT[0] = (TH1D*) hrT_in[0]->Clone();   // 0-5
  hrT[1] =  (TH1D*) hrT_in[1]->Clone();   // 5-10
  hrT[2] =  (TH1D*) hrT_in[2]->Clone();   // 10-20
  hrT[3] =  (TH1D*) hrT_in[3]->Clone();   // 20-30
  hrT[3]->Add( (TH1D*) hrT_in[4]->Clone());  // 30-40
  hrT[4] =  (TH1D*) hrT_in[5]->Clone();   // 40-50
  hrT[4]->Add( (TH1D*) hrT_in[6]->Clone());  // 50-60
  hrT[5] =  (TH1D*) hrT_in[8]->Clone();   // is this  actually 60-84, based on counts?
  hrT[5]->Add( (TH1D*) hrT_in[7]->Clone());  // 60-70

  double ncoll[NCENT] = {9.7, 8.4, 7.4, 6.1, 4.4, 2.6};
  int col[NCENT] = {kRed, kGreen, kBlue, kMagenta, kViolet, kRed };
  int centlow[NCENT] = {0, 5, 10, 20, 40, 60 };
  int centhigh[NCENT] = {5,10,20,40,60,84 };

#endif

#ifdef PAL

  cout << "Using measured pAl rT distributions" << endl;

  // get the pAl rT distributions from a file made by Jamie
  //===================================
  /*
    for p+Al with 9 centrality bins (the standard 9 we have used, i.e. 0-5, 5-10, 10-20, 20-30, 30-40, 40-50, 50-60, 60-70, 70-84). 
  */

  TFile *fhist = new TFile("rT_distributions_pal.root");
  if(!fhist)
    {
      cout << "Failed to open rT distribution file " << endl;
      return;
    }

  static const int NCENT_in = 9;
  char hname[NCENT_in][100]={"himpact_cent0","himpact_cent1","himpact_cent2","himpact_cent3","himpact_cent4", "himpact_cent5", "himpact_cent6", "himpact_cent7", "himpact_cent8"};
  TH1D *hrT_in[NCENT_in];
  for(int icent=0;icent<NCENT_in;icent++)
    {
      cout << "Getting pAl histogram " << hname[icent] << endl;

      hrT_in[icent] = (TH1D *) fhist->Get(hname[icent]);
      if(!hrT_in[icent])
	{
	  cout << "Did not get histogram " << hname[icent] << endl;
	  return;
	}
      cout << "  icent " << icent << " hrT_in integral " << hrT_in[icent]->Integral() << endl;
    }

  // combine glauber centralities to get the experimental ones
  //========================================
  static const int NCENT = 6;  // includes MB and 40-60% and 60-72%
  TH1D *hrT[NCENT];
  // 0-20%
  hrT[0] = (TH1D*) hrT_in[0]->Clone();   // 0-5
  hrT[0] ->Add( (TH1D*) hrT_in[1]->Clone());   // 5-10
  hrT[0]->Add((TH1D*) hrT_in[2]->Clone());   // 10-20
  // 20-40%
  hrT[1] =  (TH1D*) hrT_in[3]->Clone();   // 20-30
  hrT[1]->Add( (TH1D*) hrT_in[4]->Clone());  // 30-40
  // 40-72%
  hrT[2] =  (TH1D*) hrT_in[5]->Clone();   // 40-50
  hrT[2]->Add( (TH1D*) hrT_in[6]->Clone());  // 50-60
  hrT[2]->Add( (TH1D*) hrT_in[7]->Clone());  // 60-70
  //hrT[2] =  (TH1D*) hrT_in[8]->Clone();  // 70-84  // the counts here seem too large - Jamie said to ignore this bin
  // MB
  hrT[3] = (TH1D*) hrT[0]->Clone();
  hrT[3]->Add(hrT[1]);
  hrT[3]->Add(hrT[2]);
  // 40-60%
  hrT[4] = (TH1D*) hrT_in[5]->Clone();
  hrT[4]->Add( (TH1D*) hrT_in[6]->Clone() );
  // 60-72%
  hrT[5] = (TH1D*) hrT_in[7]->Clone();

  double ncoll[NCENT] = {3.35, 2.3, 1.7, 2.1};
  int col[NCENT] = {kRed, kGreen, kBlue, kBlack};
  int centlow[NCENT] = {0, 20, 40, 0};
  int centhigh[NCENT] = {20,40,72,100};

#endif

#ifdef HEAU

  cout << "Using PHENIX Glauber HeAu rT distributions" << endl;
  
  // get the HeAu rT distributions from a file made by Jamie
  //===================================
  /*
    for He+Au with 9 centrality bins (the standard 9 we have used, i.e. 0-5, 5-10, 10-20, 20-30, 30-40, 40-50, 50-60, 60-70, 70-84). 
  */
  TFile *fhist = new TFile("rT_distributions_3heau.root");
  if(!fhist)
    {
      cout << "Failed to open rT distribution file " << endl;
      return;
    }

  static const int NCENT_in = 9;
  char hname[NCENT_in][100]={"himpact_cent0","himpact_cent1","himpact_cent2","himpact_cent3","himpact_cent4", "himpact_cent5", "himpact_cent6", "himpact_cent7", "himpact_cent8"};
  TH1D *hrT_in[NCENT_in];
  for(int icent=0;icent<NCENT_in;icent++)
    {
      cout << "Getting 3He+Au histogram " << hname[icent] << endl;
      
      hrT_in[icent] = (TH1D *) fhist->Get(hname[icent]);
      if(!hrT_in[icent])
	{
	  cout << "Did not get histogram " << hname[icent] << endl;
	  return;
	}
      cout << "  icent " << icent << " hrT_in integral " << hrT_in[icent]->Integral() << endl;
    }

  // combine glauber centralities to get the experimental ones
  //========================================
  static const int NCENT = 6;  // includes MB and 40-60% and 60-84%
  TH1D *hrT[NCENT];
  hrT[0] = (TH1D*) hrT_in[0]->Clone();   // 0-5
  hrT[0] ->Add( (TH1D*) hrT_in[1]->Clone());   // 5-10
  hrT[0]->Add((TH1D*) hrT_in[2]->Clone());   // 10-20
  hrT[1] =  (TH1D*) hrT_in[3]->Clone();   // 20-30
  hrT[1]->Add( (TH1D*) hrT_in[4]->Clone());  // 30-40
  hrT[2] =  (TH1D*) hrT_in[5]->Clone();   // 40-50
  hrT[2]->Add( (TH1D*) hrT_in[6]->Clone());  // 50-60
  hrT[2]->Add( (TH1D*) hrT_in[7]->Clone());  // 60-70
  //hrT[2] =  (TH1D*) hrT_in[8]->Clone();  // 70-84  // the counts here seem too large - Jamie said to ignore this bin
  // MB
  hrT[3] = (TH1D*) hrT[0]->Clone();
  hrT[3]->Add(hrT[1]);
  hrT[3]->Add(hrT[2]);
  // 40-60%
  hrT[4] = (TH1D*) hrT_in[5]->Clone();
  hrT[4]->Add( (TH1D*) hrT_in[6]->Clone() );
  // 60-88%
  hrT[5] = (TH1D*) hrT_in[7]->Clone();
  
  double ncoll[NCENT] = {};
  int col[NCENT] = {kRed, kGreen, kBlue, kBlack};
  int centlow[NCENT] = {0, 20, 40, 0};
  int centhigh[NCENT] = {20,40,88,100};
  
#endif

#ifndef PAU
#ifndef PAL
#ifndef HEAU

  cout << "Using default rT distribution from WS parameters ws_radius " << ws_radius << " ws_diffusion " << ws_diff  
       << " rho_0 " << rho_0<< endl;

  // make something up for now based on the Glauber parameters.
  // assume a linear dependence on rT with a sharp cut at ws_radius

  static const int NCENT = 7;

  TH1D *hrT[NCENT];
  for(int icent = 0; icent < NCENT; ++icent)
    {
      char hname[500];
      sprintf(hname, "hrt%i", icent);
      hrT[icent] = new TH1D(hname, hname, NRT, 0, ws_radius*1.2);
    }

  double rT_step = ws_radius*1.2 / (double) NRT;
  for(int irT=0; irT < NRT; ++irT)
    {
      double rt = 0.0 + (double) irT * rT_step;

      if(rt > ws_radius) continue;

      int ncounts = (int) (rt * 200.0);
      hrT[NCENT-1]->SetBinContent(irT+1, ncounts);
      for(int icent = 0; icent < NCENT-1; ++icent)
      {
	hrT[icent]->SetBinContent(irT+1, ncounts);
      }

    }
  double ncoll[NCENT] = {3.35, 2.3, 1.7, 2.1};
  int col[NCENT] = {kRed, kGreen, kBlue, kBlack};
  int centlow[NCENT] = {0, 20, 40, 0};
  int centhigh[NCENT] = {20,40,72,100};

#endif
#endif
#endif

  double rT[NRT];
  double thick[NRT];
  double rT_weight[NCENT][NRT];
  for(int icent=0;icent<NCENT;icent++)
    {  
      cout << "  icent " << icent << " hrT integral " << hrT[icent]->Integral() << endl;
      for(int irt=0;irt<NRT;irt++)
	{
	  rT_weight[icent][irt] = (double) hrT[icent]->GetBinContent(irt+1);
	  rT[irt] = hrT[icent]->GetBinCenter(irt+1);
	  //cout << "irt " << irt << " rT "<< rT[irt] << endl;
	}
    }
  for(int irt=0;irt<NRT;irt++)
    {
      // this is essentially TA
      thick[irt] = get_thickness_integral(rT[irt]);
      //cout << " irt " << irt << " rT " << rT[irt] << " thick " << thick[irt] << endl;
    }

  // As a check, calculate the average rT within each centrality bin using the rT distributions for the centrality bins
  //   -- this agrees perfectly with the results of "calculate_sigma_breakup.C"
  double rT_avge[NCENT];
  double thick_avge[NCENT]; // average T_A weighted only by r_T
  double thick_avge_ncoll[NCENT];  // Ncoll (really TA) weighted average T_A
  for(int icent=0;icent<NCENT;icent++)
    {
      rT_avge[icent]=0.0;
      thick_avge[icent]=0.0;
      thick_avge_ncoll[icent]=0.0;
      double wt=0.0;
      double wt_ncoll=0.0;

      for(int irt=0;irt<NRT;irt++)
	{
	  rT_avge[icent] += rT[irt]*rT_weight[icent][irt];
	  thick_avge[icent] += thick[irt]*rT_weight[icent][irt];
	  thick_avge_ncoll[icent] += thick[irt]*rT_weight[icent][irt] * thick[irt];  // add weight for T_A scaling of hard process
	  wt += rT_weight[icent][irt];
	  wt_ncoll += rT_weight[icent][irt] * thick[irt];
	  //cout << " icent " << icent << " rT " << rT[irt] << " thick " << thick[irt] << " wt " << wt << " wt_ncoll " << wt_ncoll << " rT_weight " << rT_weight[icent][irt] << " rT_avge " << rT_avge[icent] << " thick_avge " << thick_avge[icent] << " thick_avge_ncoll " << thick_avge_ncoll[icent] << endl;
	}

      rT_avge[icent] = rT_avge[icent]/wt;
      thick_avge[icent] = thick_avge[icent]/wt;
      thick_avge_ncoll[icent] = thick_avge_ncoll[icent]/wt_ncoll;
      cout << " icent " << icent << " rT_avge " << rT_avge[icent] 
	   << " thick_avge " << thick_avge[icent] 
	   << " thick_avge_ncoll " << thick_avge_ncoll[icent]
	   << endl;
    }


  // breakup cross section at backward rapidity only
  // calculated from time spent in nucleus as per 2013 paper
  // so it is a function of pTand rapidity

  cout << "Make pT TF1's " << endl;

  TF1* fpt[2];

  double fpt_par[2][3];
  //define the fit parameters for each rapidity bin
  //Based on fits to J/psi Invariant Yields * Acceptance from run6+run8 p+p
  //backward/forward rapidity p+p, from Darren
  // currently assumed to be the same for psi(2S)
  //   -- data? maybe should scale with pT/M? see arXiv:2006.15446 

  fpt_par[0][0] = 5.02367e-09; 
  fpt_par[0][1] = 2.94113;
  fpt_par[0][2] = -4.27;
  fpt[0] = new TF1("fpt","[0]*x*(1+[1]*x^2)^([2])",0,20);  
  fpt[0]->SetParameter(0,fpt_par[0][0]);
  fpt[0]->SetParameter(1,fpt_par[0][1]);
  fpt[0]->SetParameter(2,fpt_par[0][2]);

  fpt_par[1][0] = 6.58295e-09;
  fpt_par[1][1] = 2.89819;
  fpt_par[1][2] = -4.19;
  fpt[1] = new TF1("fpt","[0]*x*(1+[1]*x^2)^([2])",0,20);  
  fpt[1]->SetParameter(0,fpt_par[1][0]);
  fpt[1]->SetParameter(1,fpt_par[1][1]);
  fpt[1]->SetParameter(2,fpt_par[1][2]);
  
  // We need to calculate the modification vs pT for each rT value
  // later we will get the value for each centrality using the appropriate 
  // rT weighting
  //===========================================

  cout << "Calculate pT values " << endl;

  double ptmin = 0.1;
  double pt[NPT];      
  double ta_keep[NRT];
  for(int ipt=0;ipt<NPT;ipt++)
    {
      pt[ipt] = ptmin + (double) ipt * ptstep;
    }

  // Calculate the effect of the absorption cross section
  //     - uses methods in "calculate_sigma_breakup.C" - included above

  cout << "Calculate sigbr " << endl;

  // Jpsi at RHIC
  double Ebeam = 100;
  if( ChangeOptions&2 ){
	Ebeam = 2510;
  }
  double mstate;
  if( !ChangeOptions&4 ){
	mstate = 3.4;
  } else{
	mstate = 10.0;
  }

  double sigbr[NRAP][NPT][NRT];
  for(int irap=0;irap < NRAP;++irap)
    {
      for(int ipt=0;ipt<NPT;++ipt)
	{
	  for(int irt=0;irt < NRT; ++irt)
	  {
	    sigbr[irap][ipt][irt] = return_sigma_breakup(Ebeam, mstate, y[irap], pt[ipt], mtarget, rT[irt]);
	  }
	}      
    }

  cout << "Calculate sigmod " << endl;

  double sigmod[NRAP][NPT][NRT];
  for(int irap=0;irap < NRAP;irap++)
    {
      for(int ipt = 0; ipt<NPT; ++ipt)
	{	
	  for(int irt=0;irt<NRT;irt++)
	    {
	      // Note that sigbr is in mb, multiplying it by 0.10 converts it to fm^2
	      // thick is T_A at rT[irt]. It has units of 1/fm^2. So the product of thick*(sigbr*0.1) is breakup probability
	      // Take the average thickness as being 1/2 the total thickness 
	      //   - i.e. the average absorption corresponds to a starting point of z = 0

	      sigmod[irap][ipt][irt] = exp(-sigbr[irap][ipt][irt]*0.10*thick[irt]*0.5);  
		  if( ChangeOptions&1 || ChangeOptions&2 )
			sigmod[irap][ipt][irt] = exp(-sigbr[irap][ipt][irt]*0.10*thick[irt]*1.);
//	      if(irap > 3)
//		sigmod[irap][ipt][irt] = 1.0;   // sigbr not applied at forward rapidity
	    }
	}
    }

  /*
  // as a check, get the average sigbr vs pT in each centrality bin
  //   -- this agrees perfectly with the results of "calculate_sigma_breakup.C"
  double sigbr_avge[NRAP][NCENT][NPT] = {0};
  for(int irap = 0; irap < NRAP; ++irap)
    {
      for(int icent = 0; icent <NCENT; ++icent)
	{
	  for(int ipt = 0; ipt < NPT; ++ipt)
	    {
	      double wt = 0.0;
	      for(int irt = 0; irt < NRT; ++irt)
		{
		  sigbr_avge[irap][icent][ipt] += sigbr[irap][ipt][irt] * rT_weight[icent][irt];
		  wt += rT_weight[icent][irt]; 
		}
	      sigbr_avge[irap][icent][ipt] /= wt;
	      //cout << " irap " << irap << " icent " << icent << " ipt " << ipt << " sigbr_avge " << sigbr_avge[irap][icent][ipt] << endl;
	    }
	}
    }
  */

  cout << "Average over centrality" << endl;

  // Get the average effect of the modification and absorption cross section within each centrality and  pT bin
  //=========================================================================
  double cent_pt_sigmod[NRAP][NCENT][NPT];

  for(int irap=0;irap < NRAP;irap++)
    {
      for(int icent=0;icent<NCENT;icent++)
	{
	  for(int ipt = 0; ipt < NPT; ++ipt)
	    {

	      // integrate the modification over the rT weights for this centrality
	      cent_pt_sigmod[irap][icent][ipt]=0.0;
	      double wt = 0.0;
	      
	      for(int irt=0;irt<NRT;irt++)
		{
		  // The modifications we have are correct at each value of rT
		  // Because we want to average over all Jpsi in this centrality bin, we have to consider that ncoll is different at different irt so the Jpsi yield weight should be different
		  // Number of collisions:
		  //    The rT_weight at irt to supply the number of projectile-target collisions in that irt bin - i.e. events with 1 OR MORE N-N collisions
		  // We have to weight the mod at each rT by the T_A at that rT. So we need for the shadowing:
		  //    The mod at irt
		  //    The T_A at irt to weight the mod for that bin by jpsi ncoll scaling
		  // Multiplying these two gives something proportional to the modified Jpsi yield

		  cent_pt_sigmod[irap][icent][ipt] += rT_weight[icent][irt]*sigmod[irap][ipt][irt] * thick[irt];
		  wt += rT_weight[icent][irt] * thick[irt];
		}
	      cent_pt_sigmod[irap][icent][ipt] = cent_pt_sigmod[irap][icent][ipt]/wt;
	    }
	}
    }      

  cout << "Average over rapidity in arms" << endl;

  // Now average all of the modifications over rapidity within the two arms
  //==================================================

  double cent_pt_sigmod_arm[NARMS][NCENT][NPT] = {0};
  for(int iarm = 0; iarm < NARMS; ++iarm)
    {
      for(int icent=0;icent<NCENT;++icent)
	{
	  for(int ipt = 0; ipt < NPT; ++ipt)
	    {
	      double wt = 0.0;
	      for(int irap = 0+iarm*4; irap < 4+iarm*4; ++irap)
		{
		  cent_pt_sigmod_arm[iarm][icent][ipt] += cent_pt_sigmod[irap][icent][ipt] * BdNdy[irap];
		  wt += BdNdy[irap];		    
		}
	      cent_pt_sigmod_arm[iarm][icent][ipt] /= wt;
	    }
	}
    }
  
  // output the MB sigmod
  cout << "double sigmod_MB_arm0[NPT] = {";
  for(int ipt = 0; ipt<NPT; ++ipt)
    {
      cout << cent_pt_sigmod_arm[0][NCENT-1][ipt];
      if(ipt < NPT-1)
	cout << ", ";
      else
	cout << "};" << endl;
    }
  cout << "double pt_bu[NPT] = {";
  for(int ipt = 0; ipt<NPT; ++ipt)
    {
      cout << pt[ipt];
      if(ipt < NPT-1)
	cout << ", ";
      else
	cout << "};" << endl;
    }  

  // now get the pT integrated modification vs centrality
  //====================================
  
  double cent_sigmod[NCENT][NRAP] = {0};
  
  for(int irap=0;irap < NRAP; ++irap)
    {
      for(int icent=0;icent<NCENT; ++icent)
	{
	  
	  double wt = 0;
	  for(int ipt = 0; ipt < NPT; ++ipt)
	    {
	      double ptwt;
	      if(irap < 4)
		ptwt  = fpt[0]->Eval(pt[ipt]);
	      else
		ptwt = fpt[1]->Eval(pt[ipt]);

	      wt += ptwt;
	      cent_sigmod[icent][irap] += cent_pt_sigmod[irap][icent][ipt] * ptwt;
	    }
	  cent_sigmod[icent][irap] = cent_sigmod[icent][irap] / wt;
	}
    }  

      for(int icent=0;icent<NCENT; ++icent)
	{
	  for(int irap=0;irap < NRAP; ++irap)
	    {
	      cout << "Centrality: " << icent << " irap " << irap << " cent_sigmod " << cent_sigmod[icent][irap] << endl;	      
	    }
	}

  
  // Average these over the rapidity within the two arms
  //===================================== 
  double cent_sigmod_arm[NARMS][NCENT] = {0};

  for(int iarm = 0; iarm < NARMS; ++iarm)
    {
      for(int icent=0;icent<NCENT;++icent)
	{
	  double wt = 0.0;
	  for(int irap = 0+iarm*4; irap < iarm*4 + 4; ++irap)
	    {
	      cent_sigmod_arm[iarm][icent] += cent_sigmod[icent][irap] * BdNdy[irap];
	      wt += BdNdy[irap];
	    }
	  cent_sigmod_arm[iarm][icent] /= wt;

	  cout << " iarm " << iarm << " icent " << icent << "  cent_sigmod_arm " << cent_sigmod_arm[iarm][icent] << endl;
	}
    }

  // output modification to a file for all centralities and backward rapidities
  char fname[500];
  sprintf(fname,"estimate_output/BU_mod_%i",process);
  ofstream fout(fname);
  for(int irap=0;irap<4;++irap)
    {
      for(int icent=0; icent<NCENT; ++icent)
	{
	  fout <<  y[irap] << "  " << icent << "  " << sigma1 << "  " << r0 << "  " << vcc << "  " << cent_sigmod[icent][irap] << "  " << BdNdy[irap] << endl;      
	}
    }
  fout.close();

  bool plot_results = true;

  if(plot_results)
    {
      // Plot the results vs pT for each arm and centrality bin
      //=====================================
      
      TCanvas *c1 = new TCanvas("c1","c1",20,20,1600,800);
      //c1->SetLeftMargin(0.15);
      c1->Divide(2,1);
      c1->cd(1);
      gPad->SetLeftMargin(0.13);
      gPad->SetBottomMargin(0.12);
      gPad->SetTopMargin(0.05);
      gPad->SetRightMargin(0.05);
      c1->cd(2);
      gPad->SetLeftMargin(0.13);
      gPad->SetBottomMargin(0.12);
      gPad->SetTopMargin(0.05);
      gPad->SetRightMargin(0.05);
      
      TH1F *htemplate = new TH1F("htemplate","htemplate",1000,0.0,7.0);
      htemplate->GetYaxis()->SetTitleSize(0.06);
      htemplate->GetXaxis()->SetTitleSize(0.06);
      htemplate->GetYaxis()->SetLabelSize(0.06);
      htemplate->GetYaxis()->SetNdivisions(502);
      htemplate->GetXaxis()->SetLabelSize(0.06);
      htemplate->GetYaxis()->CenterTitle();
      htemplate->GetYaxis()->SetTitleOffset(0.9);
      htemplate->GetXaxis()->SetTitleOffset(0.9);
      htemplate->SetXTitle("p_{T} (GeV/c)");
      htemplate->SetYTitle("R_{pAu}");
      htemplate->SetMaximum(1.5);
      htemplate->SetMinimum(0.0);
      
      TLine *unit = new TLine(0.0, 1.0, 7.0, 1.0);
      unit->SetLineStyle(2);
      
      TGraph *gsigmod[NARMS][NCENT];
      
      for(int iarm=0;iarm < NARMS;iarm++)
	{
	  c1->cd(iarm + 1);
	  htemplate->DrawCopy();
	  
	  for(int icent = 0; icent < NCENT; ++icent)   
	    {
	      
	      gsigmod[iarm][icent] = new TGraph(NPT, pt, cent_pt_sigmod_arm[iarm][icent]);
	      gsigmod[iarm][icent]->SetLineStyle(2);
	      gsigmod[iarm][icent]->SetLineWidth(2.0);
	      gsigmod[iarm][icent]->SetLineColor(kRed);
	      if(icent == NCENT-1) gsigmod[iarm][icent]->SetLineColor(kBlue); 
	      gsigmod[iarm][icent]->Draw("same");
	      
	    }
	  char name[100];
	  sprintf(name,"y=%.1f to %.1f",ylow[iarm],yhigh[iarm]);
	  TLatex *raplab = new TLatex(0.35,0.77,name);
	  raplab->SetNDC();
	  raplab->SetTextSize(0.08);
	  raplab->Draw();
	  
	  unit->Draw();
	}
      
      
      
      // Plot the results vs Ncoll
      //=================
      
      TCanvas *c2 = new TCanvas("c2","c2",20,20,1600,800);
      c2->Divide(2,1);
      c2->SetLeftMargin(0.15);
      c2->cd(1);
      gPad->SetLeftMargin(0.13);
      gPad->SetBottomMargin(0.12);
      gPad->SetTopMargin(0.05);
      gPad->SetRightMargin(0.05);
      c2->cd(2);
      gPad->SetLeftMargin(0.13);
      gPad->SetBottomMargin(0.12);
      gPad->SetTopMargin(0.05);
      gPad->SetRightMargin(0.05);
      
      TH1F *htemplate2 = new TH1F("htemplate2","htemplate2",1000,0.0,10.0);
      htemplate2->GetYaxis()->SetTitleSize(0.06);
      htemplate2->GetXaxis()->SetTitleSize(0.06);
      htemplate2->GetYaxis()->SetLabelSize(0.06);
      htemplate2->GetYaxis()->SetNdivisions(202);
      htemplate2->GetXaxis()->SetLabelSize(0.06);
      htemplate2->GetYaxis()->CenterTitle();
      htemplate2->GetYaxis()->SetTitleOffset(0.9);
      htemplate2->GetXaxis()->SetTitleOffset(0.9);
      htemplate2->GetXaxis()->SetTitle("N_{coll}");
      htemplate2->SetYTitle("R_{pAu}");
      htemplate2->SetMaximum(2.0);
      htemplate2->SetMinimum(0.0);
      
      TLine *unit2 = new TLine(0.0, 1.0, 10.0, 1.0);
      unit2->SetLineStyle(2);
      
      TGraph *gsigmodcent[NARMS];
      
      c2->cd(1);
      htemplate2->Draw();
      unit2->Draw();
      c2->cd(2);
      htemplate2->Draw();
      unit2->Draw();
      
      for(int iarm=0;iarm < NARMS;iarm++)
	{
	  c2->cd(iarm + 1);
	  
	  gsigmodcent[iarm] = new TGraph(NCENT-1, ncoll, cent_sigmod_arm[iarm]);
	  gsigmodcent[iarm]->SetLineColor(kRed);
	  gsigmodcent[iarm]->Draw("same");
	  
	  char name[100];
	  sprintf(name,"y=%.1f to %.1f",ylow[iarm],yhigh[iarm]);
	  TLatex *raplab = new TLatex(0.35,0.77,name);
	  raplab->SetNDC();
	  raplab->SetTextSize(0.08);
	  raplab->Draw();
	}
      
      // plot the results vs  rapidity in each centrality bin, integrated over pT
      //===============================================
      
      TCanvas *crap = new TCanvas("crap","crap", 5,5,800,800);
      gPad->SetLeftMargin(0.15);
      
      TH1D *hrap = new TH1D("hrap","hrap",100,-2.5, 2.5);
      hrap->GetYaxis()->SetTitle("Modification (absorption only)");
      hrap->GetYaxis()->SetTitleOffset(1.4);
      hrap->GetXaxis()->SetTitle("rapidity");
      hrap->Draw();
      
      TLegend * lrap = new TLegend(0.6, 0.6, 0.9, 0.9,"Centrality");
      
      TGraph *gr_rap[NCENT];
      for(int icent=0; icent<NCENT; ++icent)
	{
	  gr_rap[icent] = new TGraph(NRAP, y, cent_sigmod[icent]);      
	  gr_rap[icent]->SetLineColor(col[icent]);
	  gr_rap[icent]->SetMarkerStyle(20);
	  gr_rap[icent]->SetMarkerColor(col[icent]);
	  gr_rap[icent]->Draw("same p");
	  
	  char cname[500];
	  sprintf(cname,"%i - %i",centlow[icent], centhigh[icent]);
	  lrap->AddEntry(gr_rap[icent],cname,"p");
	}
      lrap->Draw();
      
      // pT dependence N
      
      TH1D *h = new TH1D("h","h",100, 0, 7.0);
      h->SetMinimum(0.0);
      h->SetMaximum(2.0);
      
      TCanvas *cdatN = new TCanvas("cdatN","cdatN", 5,5,1600,800);
      cdatN->Divide(3,3);
      for(int icent=0; icent < NCENT; ++icent)
	{
	  cdatN->cd(icent+1);
	  h->Draw();
	  
	  gsigmod[1][icent]->Draw();
	}
      
      // pT dependence S
      
      TCanvas *cdatS = new TCanvas("cdatS","cdatS", 5,5,1600,800);
      cdatS->Divide(3,3);
      for(int icent=0; icent < NCENT; ++icent)
	{
	  cdatS->cd(icent+1);
	  h->Draw();
	  
	  gsigmod[0][icent]->Draw();
	}
      
      
      TCanvas *crt = new TCanvas("crt","crt", 5,5,1200,800);
      hrT[NCENT-1] ->Draw();
      for(int icent=0; icent<NCENT-1; ++icent)
	{
	  hrT[icent]->SetLineColor(col[icent]);
	  hrT[icent]->Draw("same");
	}
      
    }
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

  // ws->Eval gives the relative prob of location/unit volume (relative because we have taken rho_0 = 1 here)
  // rho_0 normalizes the WS distribution so that the integral is the atomic number, 197 in this case
  // So multiply ws by rho_0 to get absolute prob
  // See get_WS_rho_0() for how to get rho_0
  // lint * rho_0 is TA at rT = b and has units of 1/fm^2
  lint = lint * rho_0;

  return lint;
}

/*
double get_modification(int npow,double A,double thick,double rg)
{
  // used only for cases where the mod is a power of the thickness
  // A is a normalization factor that depends on the power
  double mod = 1.0 - ((1-rg)/A) * pow(thick,npow);

  cout << "     get_modification: npow " << npow << " A " << A 
       << " thick " << thick << " rg " << rg << " mod " << mod 
       << endl;

  if(mod < 0.0)
    mod=0.0;

  return mod;

}
*/

double get_WS_rho_0()
{
  // The WS distribution  
  TF1 *ws = new TF1("f1","1.0/(1+exp((x-[0])/[1]))",0.,20.);

  ws->SetParameters(ws_radius,ws_diff);

  // calculate the volume integral of the Woods Saxon distribution
  // The volume element is:
  //      dV = dr * d(theta) * {Sin(theta) * d(phi)} = r^2 Sin(theta) dr d(theta) d(phi)
  // Then:
  // Integral = Int{f(r) r^2 dr}Int{sin(theta) d(theta)}Int{d(phi)}
  // where f(r) = 1/(1+exp( (r-ws_radius)/ws_diff))
  //     Int{Sin(theta)d(theta)} = [-cos(theta)](0,pi) = 1.0 -(-1) = 2.0
  //     Int{d(phi)} = [phi](0,2pi) = 2 pi
  //     Int{r^2/(1+exp((r-b)/a)) dr} - do it numerically!

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
  double rho0 = mtarget / ws_integral;

  // The product should be 197 if all is done correctly!
  cout << " mtarget " << mtarget << " rho0 " << rho0 << " WS integral = " << ws_integral << " product " << rho0 * ws_integral << endl;

  return rho0;
}


