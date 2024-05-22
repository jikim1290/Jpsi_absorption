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

void compareVar(){

 const int ncent = 6;
 double cent[ncent] = {
	355.646, 291.16, 233.436,
	138.747, 62.5276, 21.2094
 };
 int centmin[ncent] = {
	0,  5, 10, 20, 40, 60 };
 int centmax[ncent] = {
	5, 10, 20, 40, 60, 84 };

 const int nfiles = 3;

 double r0[nfiles] = {
	0, 0.054161983, 0.093090909 };

 int lc[3] = { 1, 30, 46};

 TCanvas* c = new TCanvas("c","c",800,700);
 gPad->SetLeftMargin(0.15);
 gPad->SetBottomMargin(0.15);
 gPad->SetRightMargin(0.03);
 gPad->SetTopMargin(0.03);
 gPad->SetTicks();

 TLegend* leg = new TLegend(0.2, 0.195, 0.6, 0.568);
 leg->SetFillStyle(0);
 leg->SetTextFont(43);
 leg->SetTextSize(28);
 leg->SetLineWidth(0);

 leg->AddEntry( (TObject*)0, "Au#font[122]{-}Au, 200 GeV, #varUpsilon(1S)", "");


 int RainbowColor[ncent];
 int sysColorPallet = GetSerialColors(ncent);
 for(int i=0;i<ncent;i++){
    RainbowColor[i] = sysColorPallet + ncent - i - 1;
 }
 TFile* fin;
 TGraph* gNMF_cent[nfiles];
 TF1* ffit[nfiles];
 TGraph* gAdd;
 for(int i=0;i<nfiles;i++){
	fin = new TFile(Form("data/NMF_absonly_out_r%d_5option.root",i),"read");
	gNMF_cent[i] = new TGraph();
	ffit[i] = new TF1("f1","pol3",0,400);
	ffit[i]->SetLineColor(i+1);
	for(int j=0;j<ncent;j++){
		gAdd = (TGraph*)fin->Get(Form("abs_%dto%d_5option",centmin[j],centmax[j]));
		gNMF_cent[i]->SetPoint( j, cent[j], gAdd->GetY()[4] ); 
	}
	gNMF_cent[i]->SetMarkerStyle(20);
	gNMF_cent[i]->SetMarkerColorAlpha(i+1, 0.6);
	gNMF_cent[i]->SetMarkerSize(2);

	gNMF_cent[i]->Fit( ffit[i] );
	leg->AddEntry( gNMF_cent[i], Form("r_{0} = %.3lf fm",r0[i]), "p");
 }

 gNMF_cent[0]->SetMinimum(0);
 gNMF_cent[0]->SetMaximum(1.3);
 gNMF_cent[0]->SetTitle(";#it{N}_{part};Nuclear modification factor (absorption only)");

 gNMF_cent[0]->Draw("AP");
 gNMF_cent[1]->Draw("P");
 gNMF_cent[2]->Draw("P");
 leg->Draw();

 c->SaveAs("figs/NMF_abs_only.pdf");


 TFile* fdata_np = new TFile("/Users/junleekim/work/SHINCHON/Final/data/ExpNMF_npart_Out.root","read");
 TGraphErrors* gSTAR_1S_stat = (TGraphErrors*)fdata_np->Get("gSTAR_1S_stat");
 TGraphErrors* gSTAR_1S_syst = (TGraphErrors*)fdata_np->Get("gSTAR_1S_syst");
 TGraphErrors* gSTAR_2S_stat = (TGraphErrors*)fdata_np->Get("gSTAR_2S_stat");
 TGraphErrors* gSTAR_2S_syst = (TGraphErrors*)fdata_np->Get("gSTAR_2S_syst");


 gSTAR_1S_syst->SetFillColorAlpha(lc[0],0.3);
 gSTAR_1S_syst->SetLineWidth(0);
 gSTAR_1S_stat->SetMarkerStyle(29);
 gSTAR_1S_stat->SetMarkerColor(lc[0]);
 gSTAR_1S_stat->SetLineColor(lc[0]);
 gSTAR_1S_stat->SetMarkerSize(2);

 gSTAR_2S_syst->SetFillColorAlpha(lc[1],0.3);
 gSTAR_2S_syst->SetLineWidth(0);
 gSTAR_2S_stat->SetMarkerStyle(29);
 gSTAR_2S_stat->SetMarkerColor(lc[1]);
 gSTAR_2S_stat->SetLineColor(lc[1]);
 gSTAR_2S_stat->SetMarkerSize(2);

 gSTAR_1S_syst->SetMaximum(1.3);
 gSTAR_1S_syst->SetMinimum(0);
 gSTAR_1S_syst->SetTitle(";#it{N}_{part};Nuclear modification factor");

 gSTAR_1S_syst->Draw("A2");
 gSTAR_1S_stat->Draw("P");

 TFile* fSHINCHON = new TFile("/Users/junleekim/work/SHINCHON/Final/data/FDCOR_out_AuAu_trento_p1.root","read");
 TGraphErrors* g_tp1_npart = (TGraphErrors*)fSHINCHON->Get("g_RAA_multdep_fdcor_1S_trento_p1");
 TGraphErrors* g_tp1_npart_addabs[nfiles];
 for(int f=0;f<nfiles;f++){
	g_tp1_npart_addabs[f] = new TGraphErrors();
	for(int i=0;i<g_tp1_npart->GetN();i++){
		g_tp1_npart_addabs[f]->SetPoint( i, g_tp1_npart->GetX()[i],
			g_tp1_npart->GetY()[i] * ffit[f]->Eval( g_tp1_npart->GetX()[i] ) );
	}
 } 
 g_tp1_npart->SetLineColor(30);
 g_tp1_npart->SetLineWidth(3);
 g_tp1_npart->Draw("C");
 for(int f=0;f<nfiles;f++){
	g_tp1_npart_addabs[f]->SetLineColor(f+1);
	g_tp1_npart_addabs[f]->SetLineWidth(3);
	g_tp1_npart_addabs[f]->Draw("C");
 }

 TLegend* leg1 = new TLegend(0.2, 0.595, 0.9, 0.968);
 leg1->SetFillStyle(0);
 leg1->SetTextFont(43);
 leg1->SetTextSize(22);
 leg1->SetLineWidth(0);
 leg1->SetNColumns(3);

 leg1->SetHeader("Au#font[122]{-}Au, 200 GeV, #varUpsilon(1S)");
 leg1->AddEntry( gSTAR_1S_stat, "STAR", "p");
 leg1->AddEntry( g_tp1_npart, "default SHINCHON", "l");
 leg1->AddEntry( (TObject*)0, "", "");
 
 for(int f=0;f<nfiles;f++){
	leg1->AddEntry( g_tp1_npart_addabs[f], Form("r_{0} = %.3lf fm",r0[f]), "l");
 }
 leg1->Draw();

 c->SaveAs("figs/NMF_wabs.pdf");
 
 
}
