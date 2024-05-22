int GetSerialColors(int colorbins = 1){
    const int NRGBs = 5;
    double stops[NRGBs] = {0.00, 0.34, 0.61, 0.84, 1.00};
    double red[NRGBs] = {0.00, 0.00, 0.87, 0.9 * 1.00, 0.51};
    double green[NRGBs] = {0.00, 0.81, 0.9 * 1.00, 0.20, 0.00};
    double blue[NRGBs] = {0.51, 0.9 * 1.00, 0.12, 0.00, 0.00};
    int FIf = TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue,
                                                 colorbins);
    return FIf;
}

void compare(){
 TFile* f200GeV = new TFile("data/NMF_absonly_out_0option.root","read");
 TFile* f5020GeV = new TFile("data/NMF_absonly_out_8option.root","read");

 const int ncent = 6;
 const int narm = 3;

 int RainbowColor[ncent];
 int sysColorPallet = GetSerialColors(ncent);
 for(int i=0;i<ncent;i++){
    RainbowColor[i] = sysColorPallet + ncent - i - 1;
 }

 const int cent_min[ncent] = {
	0,  5, 10, 20, 40, 60 };
 const int cent_max[ncent] = {
	5, 10, 20, 40, 60, 84 };

 TGraph* grapAu[ncent];
 TGraph* gptAu[narm][ncent];
 TGraph* grapPb[ncent];
 TGraph* gptPb[narm][ncent];


 for(int j=0;j<ncent;j++){
	cout << Form("abs_%dto%d_%doption",cent_min[j],cent_max[j],0) << endl;
	grapAu[j] = (TGraph*)f200GeV->Get(Form("abs_%dto%d_%doption",cent_min[j],cent_max[j],0));
	grapPb[j] = (TGraph*)f5020GeV->Get(Form("abs_%dto%d_%doption",cent_min[j],cent_max[j],8));

	grapAu[j]->SetMarkerColorAlpha( RainbowColor[j], 0.6 );
	grapPb[j]->SetMarkerColorAlpha( RainbowColor[j], 0.6 );

	grapAu[j]->GetXaxis()->SetTitleFont(43);
	grapAu[j]->GetYaxis()->SetTitleFont(43);
	grapAu[j]->GetXaxis()->SetLabelFont(43);
	grapAu[j]->GetYaxis()->SetLabelFont(43);

	grapAu[j]->GetXaxis()->SetTitleSize(32);
	grapAu[j]->GetYaxis()->SetTitleSize(32);
	grapAu[j]->GetXaxis()->SetLabelSize(28);
	grapAu[j]->GetYaxis()->SetLabelSize(28);

	grapAu[j]->SetTitle(";rapidity;NMF (absorption only)");
	grapAu[j]->SetMaximum(1.3);
	grapAu[j]->SetMinimum(0.0);

	grapPb[j]->SetMarkerSize(2);
	grapAu[j]->SetMarkerSize(2);

    grapAu[j]->SetMarkerStyle(28);
	grapPb[j]->SetMarkerStyle(29);


	for(int i=0;i<narm;i++){
		gptAu[i][j] = (TGraph*)f200GeV->Get(Form("NMF_pt_%drap_%d_%dcent_%doption",i,cent_min[j],cent_max[j],0));
		gptPb[i][j] = (TGraph*)f5020GeV->Get(Form("NMF_pt_%drap_%d_%dcent_%doption",i,cent_min[j],cent_max[j],8));

		gptAu[i][j]->SetMarkerColor( RainbowColor[j] );
		gptPb[i][j]->SetMarkerColor( RainbowColor[j] );

		gptAu[i][j]->SetMarkerColorAlpha( RainbowColor[j], 0.6 );
		gptPb[i][j]->SetMarkerColorAlpha( RainbowColor[j], 0.6 );

	    gptAu[i][j]->GetXaxis()->SetTitleFont(43);
	    gptAu[i][j]->GetYaxis()->SetTitleFont(43);
	    gptAu[i][j]->GetXaxis()->SetLabelFont(43);
	    gptAu[i][j]->GetYaxis()->SetLabelFont(43);

	    gptAu[i][j]->GetXaxis()->SetTitleSize(32);
	    gptAu[i][j]->GetYaxis()->SetTitleSize(32);
	    gptAu[i][j]->GetXaxis()->SetLabelSize(28);
	    gptAu[i][j]->GetYaxis()->SetLabelSize(28);

	    gptAu[i][j]->SetTitle(";#it{p}_{T} (GeV/#it{c});NMF (absorption only)");
	    gptAu[i][j]->SetMaximum(1.3);
	    gptAu[i][j]->SetMinimum(0.0);

	    gptPb[i][j]->SetMarkerSize(2);
	    gptAu[i][j]->SetMarkerSize(2);

	    gptAu[i][j]->SetMarkerStyle(28);
	    gptPb[i][j]->SetMarkerStyle(29);

	}
 }

 TCanvas* c = new TCanvas("c","c",800,700);
 gPad->SetLeftMargin(0.15);
 gPad->SetBottomMargin(0.15);
 gPad->SetRightMargin(0.03);
 gPad->SetTopMargin(0.03);
 gPad->SetTicks();

 TLegend* leg = new TLegend(0.597, 0.195, 0.997, 0.568);
 leg->SetFillStyle(0);
 leg->SetTextFont(43);
 leg->SetTextSize(28);
 leg->SetLineWidth(0); 
 grapAu[0]->Draw("AP");
 leg->AddEntry( grapAu[0], "200 GeV", "p");
 leg->AddEntry( grapPb[0], "5020 GeV", "p");

 for(int j=0;j<ncent;j++){
	grapAu[j]->Draw("P");
	grapPb[j]->Draw("P");

	leg->AddEntry( grapPb[j], Form("%d#font[122]{-}%d%%",cent_min[j],cent_max[j]), "p");	
 }
 leg->Draw();
 c->SaveAs("figs/nmf_rap.pdf");

 gptAu[1][0]->Draw("AP");
 for(int j=0;j<ncent;j++){
    gptAu[1][j]->Draw("P");
    gptPb[1][j]->Draw("P");
 }

 leg->Draw();
 c->SaveAs("figs/nmf_pt.pdf");

}
