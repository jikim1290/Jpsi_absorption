#include <iostream>
#include <TF1.h>
#include <TH1F.h>
#include <TFile.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLegend.h>

// This code calculates the J/psi breakup cross section 
// due to time spent in the target nucleus from the model
// and global fit parameters described in 
// Phys. Rev. C87 (2013) 05490 

double get_tau(double yrel, double MJpsi, double pt, double L);
double get_sigma(double ccbar_roots, double tau);
double get_sigma_psi2s(double ccbar_roots, double tau);
double get_ccbar_sqrts(double yrel, double MJpsi, double pt);
double return_sigma_breakup(double Ebeam, double mstate, double ystate, double pt, double mtarget, double L);
double get_L(double Targmass);
double get_z_wt(double rT, double z, double zstep);
double get_path_thickness(double b, double zstart);
void set_WS_parameters(double ws_radius, double ws_diffusion, double rho_0, double mtarget);
void set_breakup_parameters(double sigma1_in, double r0_in, double vcc_in);

// Default values are for Au target, can be overridden when using 
double rho_0 = 0.169;
double mtarget = 197;
// from the glauber_dau.inp file that I used                                                                          
double ws_radius = 6.38;
double ws_diff = 0.54;

// this one gives a chisq of 15..006, was used in our paper
double sigma1 = 7.2; double r0 = 0.16; double vcc = 1.0;
// this one gives a chisq of 15.076
//double sigma1 = 7.6; double r0 = 0.14; double vcc = 1.05;


int ChangeOptions = 1;
// 1: System with AuAu
// 2: System with PbPb
// 4: Upsilon


void calculate_sigma_breakup()
{
  // This main program exercises the function "return_sigma_breakup" 
  // to calculate the breakup cross section for p+Au collisions at a specified rapidity

  cout << "ChangeOptions: " << ChangeOptions << endl;

  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(true);

  double Ebeam = 100;  // GeV
  double mtarget = 197;
  //double ystate = -2.075;  // specified rapidity
  double ystate = 0.0;  // specified rapidity
  //double ystate = 2.075;  // specified rapidity
  double pt = 1.90;  // mean pT for d+Au J/psi
  double mstate = 3.4;

  /*
  double Ebeam = 2500;  // GeV
  double mtarget = 208;
  double ystate = -4.0;  // specified rapidity 
  double pt = 2.4;  // guess for LHC?
  */

  if( ChangeOptions&2 ){
    Ebeam = 2510;
    mtarget = 208;
    ystate = 0.0;  
    pt = 2.4;
    if( ChangeOptions&4 ){
		pt = 5.91;
		mstate = 10.0;
	}
  } else{
	if( ChangeOptions&4 ){
		pt = 3.07; 
		mstate = 10.0;
	}
  }


  
  // Get the rT distributions from the input file
  //==============================
  static const int NCENT_in = 9;
  char hname[NCENT_in][100]={"himpact_cent0","himpact_cent1","himpact_cent2","himpact_cent3","himpact_cent4", "himpact_cent5", "himpact_cent6", "himpact_cent7", "himpact_cent8"};
  TH1D *hrT_in[NCENT_in];
  
  TFile *fhist;

  fhist = new TFile("rT_distributions_pau.root");
  if( ChangeOptions&1 ){
	fhist = new TFile("rT_distributions_auau.root");
  }
  else if( ChangeOptions&2 ){
    fhist = new TFile("rT_distributions_pbpb.root");
  }
  if(!fhist)
    {
      cout << "Failed to open rT distribution file " << endl;
      return;
    }

  for(int icent=0;icent<NCENT_in;icent++)
    {
      cout << " icent " << icent << endl;
      cout << "Getting histogram " << hname[icent] << endl;

      hrT_in[icent] = (TH1D *) fhist->Get(hname[icent]);
      if(!hrT_in[icent])
	{
	  cout << "Did not get histogram " << hname[icent] << endl;
	  return;
	}
      cout << "  icent " << icent << " hrT_in integral " << hrT_in[icent]->Integral() << endl;
    }

  // combine centralities to get the experimental ones
  //===================================
  static const int NCENT = 6;
  TH1D *hrT[NCENT];
  hrT[0] = (TH1D*) hrT_in[0]->Clone();   // 0-5
  hrT[1] =  (TH1D*) hrT_in[1]->Clone();   // 5-10
  hrT[2] =  (TH1D*) hrT_in[2]->Clone();   // 10-20
  hrT[3] =  (TH1D*) hrT_in[3]->Clone();   // 20-30
  hrT[3]->Add( (TH1D*) hrT_in[4]->Clone());  // 30-40
  hrT[4] =  (TH1D*) hrT_in[5]->Clone();   // 40-50
  hrT[4]->Add( (TH1D*) hrT_in[6]->Clone());  // 50-60
  hrT[5] =  (TH1D*) hrT_in[8]->Clone();   // 60-70
  /*
  hrT[5] =  (TH1D*) hrT_in[7]->Clone();   // 60-70
  hrT[5]->Add( (TH1D*) hrT_in[8]->Clone());  // 70-84
  */

  static const int NRT = 100; //hard-coded, the number of hrt bins
  double rT[NRT];
 double rT_weight[NCENT][NRT];
  for(int icent=0;icent<NCENT;icent++)
    {  
      cout << "  icent " << icent << " hrT integral " << hrT[icent]->Integral() << endl;
      for(int irt=0;irt<NRT;irt++)
	{

	  rT_weight[icent][irt] = (double) hrT[icent]->GetBinContent(irt+1);
	  rT[irt] = hrT[icent]->GetBinCenter(irt+1);
	}
    }

  // as a check, get rT averages in the centrality bins
  double rT_avge[NCENT] = {0};
  for(int icent = 0; icent < NCENT; ++icent)
    {
      double wt = 0.0;
      for(int irt = 0; irt < NRT; ++irt)
	{
	  rT_avge[icent] += rT[irt] * rT_weight[icent][irt];	  
	  wt += rT_weight[icent][irt];
	}
      rT_avge[icent] /= wt;
      cout << " icent " << icent << " rT_avge " << rT_avge[icent] << endl;
    }

  static const int NPT = 40;
  double pt_array[NPT];
  double ptmin = 0.5;
  double ptstep = 0.5;
  for(int ipt = 0; ipt < NPT; ++ipt)
    {
      pt_array[ipt] = ptmin + (double) ipt * ptstep;
    }

  double sigma[NCENT][NPT];
  TGraph *grpt[NCENT] ;
  for(int icent = 0; icent < NCENT; ++icent)
    {
      for(int ipt = 0; ipt < NPT; ++ipt)
	{
	  double sigma_avge = 0.0;
	  double wt = 0.0;
	  for(int irt=0;irt<NRT;irt++)
	    {
	      double sigbreakup = return_sigma_breakup(Ebeam, mstate, ystate, pt_array[ipt], mtarget, rT[irt]);
	      //cout << " irt " << irt << " rT " << rT[irt] << " rT_weight " << rT_weight[irt] << " sigma  " << sigbreakup << endl;
	      
	      sigma_avge += sigbreakup * rT_weight[icent][irt];      
	      wt += rT_weight[icent][irt];
	    }
	  sigma[icent][ipt] = sigma_avge / wt;
	  cout << " icent " << icent << " pT " << pt_array[ipt] << " sigma_avge = " << sigma_avge / wt << endl;
	}
    }

  //int col[NCENT] = {kRed, kBlue, kGreen, kViolet, kYellow};
  int col[NCENT] = {kBlack, kCyan, kRed, kBlue, kViolet, kGreen};
  //char centlabel[NCENT][100] = {"0-20", "20-40","40-60","60-84","0-100"};
  char centlabel[NCENT][100] = {"0-5", "5-10","10-20","20-40","40-60","60-84"};
  TLegend *leg = new TLegend(0.65, 0.6, 0.89, 0.89,"Centrality","NDC");

  TCanvas *cpt = new TCanvas("cpt","cpt",5,5,800,800);
  TH1F *h = new TH1F("h","y = -1.7",100, 0, 20);
  h->GetYaxis()->SetTitle("#sigma_{abs} (mb)");
  h->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  h->SetMinimum(0);
  h->SetMaximum(25);
  h->Draw();
  for(int icent = 0; icent < NCENT; ++icent)
    {
      // make a TGraph for this centrality
      grpt[icent] = new TGraph(NPT, pt_array, sigma[icent]);
      grpt[icent]->SetLineColor(col[icent]);
      grpt[icent]->SetLineWidth(2.0);
      grpt[icent]->Draw("same");      
      leg->AddEntry(grpt[icent], centlabel[icent]);
    }
  leg->Draw();

  // Plot the rT distributions

  TCanvas *crt = new TCanvas("crt","rT distributions", 20, 20, 800, 800);
  gPad->SetLogy(true);
  for(int icent=0;icent<NCENT;++icent)
    {
      double scale = hrT[icent]->Integral();
      //hrT[icent]->Scale(20.0/scale);
      hrT[icent]->SetLineColor(col[icent]);
      hrT[icent]->SetLineWidth(2);
      //hrT[icent]->Rebin(2);
      if(icent==0)
	{
	  hrT[icent]->SetMaximum(400);
	  hrT[icent]->SetMinimum(0.5);
	  hrT[icent]->GetXaxis()->SetRangeUser(0, 18);
	  hrT[icent]->GetYaxis()->SetTitle("Yield");
	  hrT[icent]->SetTitle("r_{T} distributions");
	  hrT[icent]->GetXaxis()->SetTitle("r_{T} (fm)");
	  hrT[icent]->Draw("hist");
	}
	else
	  hrT[icent]->Draw("same hist");
    }
  leg->Draw();


}

double return_sigma_breakup(double Ebeam, double Mstate, double ystate, double pt, double mtarget, double rT)
{
  // ybeam is the full rapidity of one of the colliding beams, ystate is the rapidity of the ccbar pair
  double Mbeam = 0.938;
  double pz = sqrt(pow(Ebeam,2) - pow(Mbeam,2));
  double ybeam = 0.5 * log((Ebeam+pz)/(Ebeam-pz));
  double yrel = ybeam + ystate;
  yrel =0.0;
  double ccbar_roots = get_ccbar_sqrts(yrel,Mstate,pt);

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

  double zmax = 10.0;    // distribution is zero beyond 10 fm radius
  double NL = 40;
  double zstep = zmax*2.0 / ((double) (NL+1) );
  
  // Going only to NL-1 keeps the path length from going to zero
  for(int iL=0;iL<NL-1;iL++)
    {
      double z = -zmax + zstep * (double) (iL+1);

      // Get the weight for this z value from the WS z distribution at rT
      double zwt_ccbar = get_z_wt(rT, z, zstep);
      
      // Get the average value of L seen by a ccbar pair produced at z	  
      // get_path_thickness returns TA, which is rho_0 * {integral(rho.dz) along rT}
      // dividing out rho_0 leaves  {integral(rho.dz) along rT}
      // which is a probability weighted mean path length (check: for rT = 0 and z = 0, path comes out at (2.156/2) / 0.17 = 6.34 fm, the WS radius)
      double path = get_path_thickness(rT, z) / rho_0;
       
      // Weight path by the ccbar production probability and add it to the average over all iL values for this irt value
      // This gives us the average over the whole line of z values at rT
      double tau = get_tau(yrel, Mstate, pt, path); 
      //cout << " rT " << rT << " z " << " rho_0 " << rho_0 << " path " << path << " tau " << tau << endl;

      tau_line += tau*zwt_ccbar;
      // Use the Arleo form for the tau dependence of the psi radius

      double sigma;
//      if(Mstate < 3.5)   // J/psi
	  sigma = get_sigma(ccbar_roots, tau);
  //    else                           // psi(2S)
//	sigma = get_sigma_psi2s(ccbar_roots, tau);
	  // TODO: correct for feed-down corrections


      sigma_line += sigma*zwt_ccbar;
      LrT_line += path*zwt_ccbar;
      wt_line += zwt_ccbar;
    }
  // averages for this rT
  
  double tau_avge = tau_line / wt_line;
  double sigma_avge = sigma_line / wt_line;

  //cout << "rT " << rT << " ccbar_roots " << ccbar_roots << " tau_avge " << tau_avge << " sigma_avge " << sigma_avge << endl;

  return sigma_avge;

}

void set_breakup_parameters(double sigma1_in, double r0_in, double vcc_in)
{
  sigma1 = sigma1_in;
  r0 = r0_in;
  vcc = vcc_in;
}

void set_WS_parameters(double ws_radius_in, double ws_diffusion_in, double rho_0_in, double mtarget_in)
{
  ws_radius = ws_radius_in;
  ws_diff = ws_diffusion_in;
  rho_0 = rho_0_in;
  mtarget = mtarget_in;
}


double get_tau(double yrel, double MJpsi, double pt, double L)
{
  // E, m and p are in GeV, GeV/c^2 and GeV/c repectively
  // Generally, if we use those units for E, p, m and v, we consider c=1 and ignore it
  // p = m * v and E = m * c^2 
  // So p/E = v/c^2 = v/c * 1/c is in units of (GeV/c) / GeV = 1/c 
  // ----- ie. p/E is v in units of 1/c, ie. p/E = beta

  // We need the J/psi Lorentz gamma factor in the target frame, which is p/E in the target frame:
  // p_z = m_T sinh(y)  and m_T = sqrt(m_0^2+p_T^2)
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
  double tau = L / (gamma * beta_z);
//  double tau = L * beta_z / gamma;
  //  if(tau > 0.0002) 
  //cout << " tau " << tau << " beta_z " << beta_z << " L " << L << " gamma " << gamma << endl;
  return tau;
}

double get_sigma(double ccbar_roots, double tau)
{
  // best fit values to global data from PRC 83 (2013) 054910
 
  // There are three states - J/psi, psi' and chi_C
  // These are the numbers used by Arleo
  double rpsi[3] = {0.43,0.87,0.67};
  // These are from PPG104
  double fdfrac[3] = {0.58,0.1,0.32};

  if( ChangeOptions&4 ){
	rpsi[0] = 0.14; rpsi[1] = 0.28; rpsi[2] = 0.39;
	fdfrac[0] = 0.67; fdfrac[1] = 0.30; fdfrac[2] = 0.03;
  }

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
	   << " r0 " << r0
	   << " vcc " << vcc
	   << " tau " << tau 
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

double get_sigma_psi2s(double ccbar_roots, double tau)
{
  // best fit values to global data from PRC 83 (2013) 054910
 
  // The psi(2S)  is the highest mass charmonium bound state, no feeddown
  // These are the numbers used by Arleo
  double rpsi1s = 0.43;
  double rpsi2s = 0.87;

  // This is the model of Arleo PRC61 where rccbar is linear in tau
  double rccbar = r0 + vcc*tau;

  if(rccbar > rpsi2s)
    rccbar = rpsi2s;

  // From Arleo PRC61 
  // The reference radius is the J/psi one, so that the cross section is larger for the fully formed excited states
  double sigma = sigma1* pow(ccbar_roots/10.0,0.4) * pow(rccbar/rpsi1s,2);

      /*
      cout << " sigma1 " << sigma1 
	   << " r0 " << r0
	   << " vcc " << vcc
	   << " tau " << tau 
	   << " ccbar_roots " << ccbar_roots
	   << " rpsi1s " << rpsi1s
	   << " rccbar " << rccbar
	   << " sigma " << sigma 
	   << endl;
      */

  return sigma;      
}

double get_ccbar_sqrts(double yrel, double Mstate, double pt)
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
  double Ecm = sqrt(pow(pt,2) + pow(Mstate,2));
  double Ebeam = Ecm * sqrt(1+pow(sinh(yrel),2));
  double pbeam = sqrt(pow(Ebeam,2) - pow(Mstate,2));
  double s = pow(Ebeam + Mtarget,2) - pow(pbeam,2);
  double sqrts = sqrt(s);

  return sqrts;
}

double get_L(double Targmass)
{
  // This is the average length over all impact parameters, not the radius
  // However it should scale as the radius

  // Use A^/13 scaling assuming W -> L = 3.95
  //double C = 3.95 / pow(184,1.0/3.0);

  // Use A^/13 scaling assuming Au -> L = 4.36
  // That value of L comes from my Glauber calculation (see time_spent_nucleus_glauber.C)
  double C = 4.36 / pow(197,1.0/3.0);

  double L = C * pow(Targmass,1.0/3.0);

  return L;

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

double get_path_thickness(double b, double zstart)
{
  // The WS distribution  
  // We drop the rho_0 normalization for now and ad it back in at the end                                                                                              
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
  // I get 0.169 from numerical volume integral of WS distribution - see get_WS_rho_0() in calculate_pAu_CNM_modification.C

  // Normalize to nucleons/ fm^2 by multiplying by rho_0
  // double rho_0 = get_WS_rho_0();
  lint = lint * rho_0;

  return lint;
}


/* 
double return_sigma_breakup(double Ebeam, double Mstate, double ystate, double pt, double mtarget, double L)
{
  // ybeam is the full rapidity of one of the colliding beams, ystate is the rapidity of the ccbar pair
  double Mbeam = 0.938;
  double pz = sqrt(pow(Ebeam,2) - pow(Mbeam,2));
  double ybeam = 0.5 * log((Ebeam+pz)/(Ebeam-pz));

  double yrel = ybeam + ystate;

  double tau = get_tau(yrel, Mstate, pt, L);
  double ccbar_roots = get_ccbar_sqrts(yrel,Mstate,pt);
  double sigma = get_sigma(ccbar_roots, tau);
  
  cout << " ybeam " << ybeam
       << " rapidity " << ystate
       << " relative rapidity " << yrel
       << " tau " << tau
       << " ccbar_sqrts " << ccbar_roots
       << " sigma_breakup " << sigma
       << endl;

  return sigma;

}
*/
