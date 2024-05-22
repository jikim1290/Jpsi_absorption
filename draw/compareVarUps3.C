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

void compareVarUps3(){

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

 leg->AddEntry( (TObject*)0, "Au#font[122]{-}Au, 200 GeV", "");


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
	if( i==0 ) fin = new TFile(Form("data/NMF_absonly_out_17option_4mb.root"),"read");
	if( i!=0 ) fin = new TFile(Form("data/NMF_absonly_out_r2_%doption.root",9),"read");

	gNMF_cent[i] = new TGraph();
	ffit[i] = new TF1("f1","pol3",0,400);
	ffit[i]->SetLineColor(i+1);
	for(int j=0;j<ncent;j++){
		if( i==0 ) gAdd = (TGraph*)fin->Get(Form("abs_%dto%d_%doption",centmin[j],centmax[j],17));
		if( i!=0 ) gAdd = (TGraph*)fin->Get(Form("abs_%dto%d_%doption",centmin[j],centmax[j],9));
		gNMF_cent[i]->SetPoint( j, cent[j], gAdd->GetY()[4] ); 
	}
	gNMF_cent[i]->SetMarkerStyle(20);
	gNMF_cent[i]->SetMarkerColorAlpha(i+1, 0.6);
	gNMF_cent[i]->SetMarkerSize(2);

	gNMF_cent[i]->Fit( ffit[i] );
	if( i!=2 ) leg->AddEntry( gNMF_cent[i], Form("#varUpsilon(%dS)",i+1), "p");
 }

 gNMF_cent[0]->SetMinimum(0);
 gNMF_cent[0]->SetMaximum(1.3);
 gNMF_cent[0]->SetTitle(";#it{N}_{part};Nuclear modification factor (absorption only)");

 gNMF_cent[0]->Draw("AP");
 gNMF_cent[1]->Draw("P");
 leg->Draw();

// c->SaveAs("figs/NMF_abs_only_ups2.pdf");

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
 gSTAR_2S_syst->Draw("2");
 gSTAR_2S_stat->Draw("P");

 TFile* fSHINCHON = new TFile("/Users/junleekim/work/SHINCHON/Final/data/FDCOR_out_AuAu_trento_p1.root","read");
 TFile* fSHINCHON1 = new TFile("/Users/junleekim/work/SHINCHON/Final/data/FDCOR_out_AuAu_trento_p1_tform600.root","read");

 TGraphErrors* g_tp1_npart[nfiles];
 TGraphErrors* g_tp1_npart_addabs[nfiles];

 for(int f=0;f<nfiles;f++){
//	if( f==0 )
	g_tp1_npart[f] = (TGraphErrors*)fSHINCHON->Get(Form("g_RAA_multdep_fdcor_%dS_trento_p1",f+1));

//	if( f==1 )
//	g_tp1_npart[f] = (TGraphErrors*)fSHINCHON1->Get(Form("g_RAA_multdep_fdcor_%dS_trento_p1_tform600",f+1));
	if( f==2 )
	g_tp1_npart[f] = (TGraphErrors*)fSHINCHON1->Get(Form("g_RAA_multdep_fdcor_%dS_trento_p1_tform600",f));

	g_tp1_npart_addabs[f] = new TGraphErrors();
	if( f==0 ){
		for(int i=0;i<g_tp1_npart[f]->GetN();i++){
			g_tp1_npart_addabs[f]->SetPoint( i, g_tp1_npart[f]->GetX()[i],
				g_tp1_npart[f]->GetY()[i] * ffit[f]->Eval( g_tp1_npart[f]->GetX()[i] ) );
		}
	}
 }
 
 int color1[3] = {15, 8, 40};
 int color2[3] = {1, 30, 40};

 for(int f=0;f<nfiles;f++){
	if(f==0)
	g_tp1_npart[f]->SetLineStyle(6);
	g_tp1_npart[f]->SetLineColor(color1[f]);
	g_tp1_npart[f]->SetLineWidth(5);
	g_tp1_npart[f]->Draw("Cex0");

	if(f==0){
		g_tp1_npart_addabs[f]->SetLineColor(color2[f]);
		g_tp1_npart_addabs[f]->SetLineWidth(5);
		g_tp1_npart_addabs[f]->Draw("Cex0");
	}
 }

 TLegend* leg1 = new TLegend(0.166, 0.730, 0.974, 0.964);
 leg1->SetFillStyle(0);
 leg1->SetTextFont(43);
 leg1->SetTextSize(22);
 leg1->SetLineWidth(0);
 leg1->SetNColumns(2);

 leg1->SetHeader("Au#font[122]{-}Au, 200 GeV");
 leg1->AddEntry( gSTAR_1S_stat, "STAR, #varUpsilon(1S)", "p");
 leg1->AddEntry( gSTAR_2S_stat, "STAR, #varUpsilon(2S)", "p");

/*
 for(int f=0;f<nfiles;f++){
	if( f!=2 ){
		leg1->AddEntry( g_tp1_npart[f], Form("#varUpsilon(%dS)",f+1), "l");
		if( f==0 ) leg1->AddEntry( g_tp1_npart_addabs[f], Form("#varUpsilon(%dS)",f+1), "l");
	} else if( f==2 ){
        leg1->AddEntry( g_tp1_npart[f], Form("#varUpsilon(%dS), #tau_{form} = 1.2 fm/#it{c}",2), "l");
        if( f==0 ) leg1->AddEntry( g_tp1_npart_addabs[f], Form("#varUpsilon(%dS), #tau_{form} = 1.2 fm/#it{c}",2), "l");
	}
 }
*/
 leg1->AddEntry( g_tp1_npart[0], "#varUpsilon(1S), default", "l");
 leg1->AddEntry( g_tp1_npart_addabs[0], "#varUpsilon(1S), absorption (#sigma_{abs} = 4 mb)", "l");
 leg1->AddEntry( g_tp1_npart[1], "#varUpsilon(2S), #tau_{form} = 1.0 fm/#it{c}", "l");
 leg1->AddEntry( g_tp1_npart[2], "#varUpsilon(2S), #tau_{form} = 1.2 fm/#it{c}", "l");


 leg1->Draw();

 c->SaveAs("figs/NMF_wabs_Ups3.pdf");
}
