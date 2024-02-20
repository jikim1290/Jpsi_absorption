#include <iostream>
#include <fstream>
#include <TCanvas.h>
#include <TH1D.h>

void process_bu_mod_error_estimate_output()
{
  static const int NCENT = 8;

  // seven centralities max
  // 4 rapidities in south arm
  TH1D *h[8][4];
  for(int icent=0; icent<NCENT; ++icent)
    {
      for(int irap=0;irap<4;irap++)
	{
	  char hname[500];
	  sprintf(hname,"h%i%i",icent,irap);
	  h[icent][irap] = new TH1D("hname", "", 100, 0.4, 1.0);
	}
    }

  TH1D *harm[8];
  for(int icent=0;icent<NCENT;icent++)
    {
      char hname[500];
      sprintf(hname,"harm%i",icent);
      harm[icent] = new TH1D("hname", "", 100, 0.4, 1.0);
    }
  
  for (int i=0;i<1000;++i)
    {
      char fname[500];
      sprintf(fname,"estimate_output/BU_mod_%i",i);
      cout << "Open file " << fname << endl;

      ifstream fin(fname);

      if(!fin.is_open())
	{
      	  cout << " failed to open file " << fname << endl;
	  continue;
	}
      
      double y, s1, r0, vcc, mod, BdNdy;
      int cent;

      double arm_mod[8] = {0,0,0,0,0,0,0,0};
      double arm_wt[8] = {0,0,0,0,0,0,0,0};
      for(int j=0;j<4;++j)
	{
	  for(int icent=0;icent<NCENT;++icent)
	    {
	      fin >> y >> cent >> s1 >> r0 >> vcc >> mod >> BdNdy;
	      cout << " y " << y << " icent " << icent << " cent " << cent << " s1 " << s1 << " r0 " << r0 << " vcc " << vcc << " mod " << mod << " BdNdy " << BdNdy << endl;  
	      h[icent][j]->Fill(mod);
	      arm_mod[icent] += mod*BdNdy;
	      arm_wt[icent] += BdNdy;
	    }
	}      

      for(int icent=0;icent<NCENT;++icent)
	harm[icent]->Fill(arm_mod[icent] / arm_wt[icent]);
    }

  TCanvas *c[NCENT];
  for(int icent=0;icent<NCENT; ++icent)
    {
      char cname[500];
      sprintf(cname,"c%i",icent);
      c[icent] = new TCanvas(cname,cname, 5,5,1600,1200);
      c[icent]->Divide(2,2);

      for(int i=0; i<4; ++i)
	{
	  c[icent]->cd(i+1);
	  h[icent][i]->Draw();
	}
    }

  TCanvas *carm = new TCanvas("carm","carm", 1600, 1200);
  harm[6]->Draw();  

}
