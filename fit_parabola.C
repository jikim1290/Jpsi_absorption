#include <TF1.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TH1D.h>

void fit_parabola()
{

//    vcc 0.95 ch2 16.22
//    vcc 1.00 ch2 15.10
//    vcc 1.05 ch2 15.64
//    vcc 1.10 ch2 17.70
// r0: 
//    r0 0.14 ch2 19.86
//    r0 0.16 ch2 15.10
//    r0 0.18 ch2  21.03
// sigma1 
//    sigma1  6.8  17.78
//    sigma1  7.0  15.91
//    sigma1  7.2  15.10
//    sigma1  7.4  15.34
//    sigma1  7.6  16.64

  TF1 *f = new TF1("f","[0]*(x-[1])*(x-[1])");

  int n[3] = {5, 3, 4};
  double mean[3] = {7.2, 0.16, 1.0};
  double p0[3] = {13.18, 13363.0, 268.0};

  double val[3][5] = {
    6.8, 7.0, 7.2, 7.4, 7.6,
    0.14, 0.16, 0.18, 0, 0,
    0.95, 1.0, 1.05, 1.10, 0
  };

  double err[3][5] = {
    17.78, 15.91, 15.10, 15.34, 16.64,
    19.86, 15.10, 21.03, 0, 0,
    16.22, 15.10, 15.64, 17.70, 0
  };

 
  for(int i=0;i<3;++i)
    {
      double min = 15.1;
      for(int j = 0; j<n[i];++j)
	{
	  err[i][j] -= min;
	  cout << " i " << i << " j " << j << " val " << val[i][j] << " err " << err[i][j] << endl;
	}
    }

  TGraph *gr[3];
  for(int i=0;i<3;++i)
    {
      gr[i] = new TGraph(n[i], val[i], err[i]);
    }

  TCanvas *c = new TCanvas("c","c",5,5,1200,800);
  c->Divide(3,1);

  double xmin[3] = {6.8, 0.13, 0.85};
  double xmax[3] = {7.8, 0.19, 01.15};

  for(int i=0;i<3; ++i)
    {
      c->cd(i+1);

      TH1D *h = new TH1D("h","h",100,xmin[i],xmax[i]);
      h->SetMaximum(6);
      h->DrawCopy();


      gr[i]->SetMarkerStyle(20);
      gr[i]->SetMarkerColor(kBlue);

      gr[i]->Draw("p");

      f->SetParameter(0,p0[i]);
      //f->FixParameter(1,mean[i]);
      f->SetParameter(1, mean[i]);

      gr[i]->Fit("f");
      cout << "results: for " << i << " are " << f->GetParameter(0) << "   " << f->GetParameter(1) << endl;

      // calculate the err = 1 points here
      double x = sqrt(1.0/f->GetParameter(0) * 1.0);
      cout << " x = " << x << endl;
    }

}
