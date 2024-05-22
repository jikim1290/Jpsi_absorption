double get_sigma(double ccbar_roots, double tau, int feedcor=3, bool isUps=false, int UpsState=0, bool useDR=false, double r0=0.16, double vcc=1, double sigma1=7.2)
{
 double rpsi[3] = {0.43,0.87,0.67};
 double fdfrac[3] = {0.58,0.1,0.32};
 if( isUps ){
	rpsi[0] = 0.14; rpsi[1] = 0.28; rpsi[2] = 0.39;
	fdfrac[0] = 0.67; fdfrac[1] = 0.30; fdfrac[2] = 0.03;
	if( UpsState==0 ){
		r0 = 0.0;
	} else if( UpsState==1 ){
		r0 = 0.093090909;
	} else if( UpsState==2 ){
		r0 = 0.054161983;
	}
	sigma1 = 4.2;
	vcc = (0.14 - r0) / 0.5;

	cout << "r0 and vcc: " << r0 << ", " << vcc << endl;
 }
 double sigma_state = 0;
 double rccbar;

 for(int i=0;i<feedcor;i++){
	rccbar = r0 + vcc*tau;
	if( rccbar > rpsi[i] ) rccbar = rpsi[i];
	if( feedcor>1 ){
		if( !useDR ){
			sigma_state += sigma1* pow(ccbar_roots/10.0,0.4) * pow(rccbar/rpsi[0],2) * fdfrac[i];
		} else if( useDR ){
			sigma_state += sigma1* pow(ccbar_roots/10.0,0.4) * pow(rccbar/rpsi[i],2) * fdfrac[i];
		}
	} else if( feedcor==1 ){
		sigma_state = sigma1* pow(ccbar_roots/10.0,0.4) * pow(rccbar/rpsi[0],2);
	}
 }
 return sigma_state;
}

double get_tau(double yrel, double MJpsi, double pt=1.9, double L=4.36)
{
  double Ecm = sqrt(pow(pt,2) + pow(MJpsi,2));
  double pz = Ecm * sinh(yrel);

  double p = sqrt(pow(pt,2)+pow(pz,2));
  double E = Ecm * sqrt(1+pow(sinh(yrel),2));
  double beta = p/E;
 
  double beta_z = pz/E;

  double gamma = 1.0/sqrt(1-pow(beta,2));

  double tau = L / (gamma * beta_z);
//  if( ChangeOptions&32 ) tau = L * beta_z / gamma;
  return tau;
}

double get_ccbar_sqrts(double yrel, double Mstate, double pt)
{
  double Mtarget = 0.938;
  double Ecm = sqrt(pow(pt,2) + pow(Mstate,2));
  double Ebeam = Ecm * sqrt(1+pow(sinh(yrel),2));
  double pbeam = sqrt(pow(Ebeam,2) - pow(Mstate,2));
  double s = pow(Ebeam + Mtarget,2) - pow(pbeam,2);
  double sqrts = sqrt(s); 

  return sqrts;
} 


TGraph* drawTauSigInd(int feedcor=3, bool isUps=false, int UpsState=0, bool useDR=false){
 double ystate = 0.0;
 double Ebeam = 100;
 double Mbeam = 0.938;
 double Mstate = 3.4;
 if( isUps ) Mstate = 9.46;

 double rapmin = 5;
 const int nrap = 200;

 double pz = sqrt(pow(Ebeam,2) - pow(Mbeam,2));
 double ybeam = 0.5 * log((Ebeam+pz)/(Ebeam-pz));

 double yrel;
 double ccbar_roots;
 double tau;
 double sigma;

 TGraph* gTauSig = new TGraph();

 for(int i=0;i<nrap;i++){
	ystate = -rapmin + 2.*rapmin*(double)i/nrap;
	yrel = ybeam + ystate;

	ccbar_roots = get_ccbar_sqrts(yrel,Mstate, 1.9);
	if( isUps ) ccbar_roots = get_ccbar_sqrts(yrel,Mstate, 3.07);

	tau = get_tau(yrel, Mstate);
    if( isUps ) tau = get_tau(yrel, Mstate, 3.07);

	sigma = get_sigma(ccbar_roots, tau, feedcor, isUps, UpsState, useDR);
	gTauSig->SetPoint( i, tau, sigma );
 }
 
 return gTauSig;
}

void drawTauSig(){
 TGraph* gjpsi = (TGraph*)drawTauSigInd(3);
 TGraph* gjpsi_nfd = (TGraph*)drawTauSigInd(1);
 TGraph* gjpsi_dr = (TGraph*)drawTauSigInd(3, false, 0, true);

 TGraph* gUps0 = (TGraph*)drawTauSigInd(3, true, 0);
 TGraph* gUps1 = (TGraph*)drawTauSigInd(3, true, 1);
 TGraph* gUps2 = (TGraph*)drawTauSigInd(3, true, 2);

 gjpsi->SetLineColor(1);
 gjpsi_nfd->SetLineColor(2);
 gjpsi_dr->SetLineColor(3);
 gjpsi->SetLineWidth(4);
 gjpsi_nfd->SetLineWidth(4);
 gjpsi_dr->SetLineWidth(4);


 gUps0->SetLineColor(1);
 gUps1->SetLineColor(2);
 gUps2->SetLineColor(3);

 gUps0->SetLineWidth(3);
 gUps1->SetLineWidth(3);
 gUps2->SetLineWidth(3);

 gUps0->SetLineStyle(4);
 gUps1->SetLineStyle(4);
 gUps2->SetLineStyle(4);

 
 TCanvas* c = new TCanvas("c","c",800,700);
 gPad->SetLeftMargin(0.15);
 gPad->SetBottomMargin(0.15);
 gPad->SetRightMargin(0.03);
 gPad->SetTopMargin(0.03); 
 gPad->SetTicks(); 
 gPad->SetLogx();
 
 TLegend* leg = new TLegend(0.199, 0.592, 0.697, 0.965);
 leg->SetFillStyle(0);
 leg->SetTextFont(43);
 leg->SetTextSize(22);
 leg->SetLineWidth(0); 
 leg->SetNColumns(2);
 gjpsi->SetTitle(";#tau (fm/#it{c});#sigma_{abs} (mb)"); 
 
 gjpsi->SetMinimum(0);
 gjpsi->SetMaximum(15);

 gjpsi->Draw("AC");
 gjpsi_nfd->Draw("C");
 gjpsi_dr->Draw("C");

 gUps0->Draw("C");
 gUps1->Draw("C");
 gUps2->Draw("C");

 TLine* ly0 = new TLine(0.0357022,0,0.0357022,15);
 ly0->SetLineWidth(2);
 ly0->SetLineStyle(2);
 ly0->Draw("same");

 leg->AddEntry( (TObject*)0, "J/#psi", "");
 leg->AddEntry( (TObject*)0, "#varUpsilon(1S)", "");
 leg->AddEntry( gjpsi, "default", "l");
 leg->AddEntry( gUps0, "r_{0} = 0 fm", "l");
 leg->AddEntry( gjpsi_nfd, "only J/#psi", "l");
 leg->AddEntry( gUps1, "r_{0} = 0.093 fm", "l");
 leg->AddEntry( gjpsi_dr, "modified radius", "l");
 leg->AddEntry( gUps2, "r_{0} = 0.054 fm", "l");

 leg->Draw();

 c->SaveAs("figs/tau_abs.pdf");

}
