void draw_exclusions()
{
 
 //----------------------------------------
 
 string model="mhmodm"; string model_label="m_{h}^{mod-}";
 //string model="hmssm"; string model_label="hMSSM";
 
 string limitname="run2tb"; string decay="tb";
 
 string lumi="13.2 fb^{-1}";
 
 double minmass=200.0; double maxmass=1000.0;
 double mintanb=0.5; double maxtanb=60; //useful for tb
 
 int do_observed=1;
 
 //----------------------------------------
 
 cout<<"Drawing exclusions '"<<limitname<<"' for "<<model<<" decay to "<<decay<<endl;
 
 int obs_fillstyle=3004;
 
 string filename="../output/exclusion_"+limitname+"_"+model+"_"+decay+".root";
 TFile* file=TFile::Open(filename.c_str());
 if(!file)
  cout<<"### ERROR: File "<<filename<<" not found"<<endl;
 if(file)
 {
  cout<<"file "<<filename<<" opened"<<endl;
  
  TGraph* g_exp_low =(TGraph*)file->Get("expected_low");  g_exp_low->SetName("g_exp_low");
  TGraph* g_exp_high=(TGraph*)file->Get("expected_high"); g_exp_high->SetName("g_exp_high");
  
  TGraph* g_1up_low =(TGraph*)file->Get("sigma1up_low");  g_1up_low->SetName("g_1up_low");
  TGraph* g_1up_high=(TGraph*)file->Get("sigma1up_high"); g_1up_high->SetName("g_1up_high");
  
  TGraph* g_1dn_low =(TGraph*)file->Get("sigma1dn_low");  g_1dn_low->SetName("g_1dn_low");
  TGraph* g_1dn_high=(TGraph*)file->Get("sigma1dn_high"); g_1dn_high->SetName("g_1dn_high");
  
  TGraph* g_2up_low =(TGraph*)file->Get("sigma2up_low");  g_2up_low->SetName("g_2up_low");
  TGraph* g_2up_high=(TGraph*)file->Get("sigma2up_high"); g_2up_high->SetName("g_2up_high");
  
  TGraph* g_2dn_low =(TGraph*)file->Get("sigma2dn_low");  g_2dn_low->SetName("g_2dn_low");
  TGraph* g_2dn_high=(TGraph*)file->Get("sigma2dn_high"); g_2dn_high->SetName("g_2dn_high");
  
  //add dummy points to the obs to make that a complete area:
  double x,y;
  int n;
  
  TGraph* g_obs_low;
  TGraph* g_obs_high;
  TGraph* g_low_area;
  TGraph* g_high_area;
  
  if(do_observed)
  {
   g_obs_low =(TGraph*)file->Get("observed_low");  g_obs_low->SetName("g_obs_low");
   g_obs_high=(TGraph*)file->Get("observed_high"); g_obs_high->SetName("g_obs_high");
   g_obs_low->SetLineWidth(3);
   g_obs_high->SetLineWidth(3);
   
   g_low_area=(TGraph*) g_obs_low->Clone("g_low_area");
   g_high_area=(TGraph*) g_obs_high->Clone("g_high_area");
   n=g_obs_low->GetN();
  }
  else
  {
   g_low_area=(TGraph*) g_exp_low->Clone("g_low_area");
   g_high_area=(TGraph*) g_exp_high->Clone("g_high_area");
   n=g_exp_low->GetN();
  }
  
   g_low_area->SetPoint(n,maxmass,-1.0);
   g_low_area->SetPoint(n+1,minmass,-1.0);
   g_low_area->GetPoint(0,x,y);
   g_low_area->SetPoint(n+2,minmass,y);
   g_high_area->SetPoint(n,maxmass,1001);
   g_high_area->SetPoint(n+1,minmass,1001);
   g_high_area->GetPoint(0,x,y);
   g_high_area->SetPoint(n+2,minmass,y);
   g_low_area->SetFillColor(1);
   g_low_area->SetFillStyle(obs_fillstyle);
   g_low_area->SetLineWidth(0);
   g_high_area->SetFillColor(1);
   g_high_area->SetFillStyle(obs_fillstyle);
   g_high_area->SetLineWidth(0);
  
   
  TCanvas* can=new TCanvas("can","can",0,0,600,600);

  g_exp_low->Draw("al");
  g_exp_low->GetXaxis()->SetTitle("m_{H^{+}} [GeV]");
  g_exp_low->GetYaxis()->SetTitle("tan#beta");
  g_exp_low->GetXaxis()->SetRangeUser(minmass,maxmass);
  g_exp_low->GetYaxis()->SetRangeUser(mintanb,maxtanb);
  g_exp_low->SetLineWidth(3);  g_exp_low->SetLineStyle(2);
  g_exp_high->SetLineWidth(3); g_exp_high->SetLineStyle(2);
  
  //Bands:
  //1sigma
  
  TGraph* g_1sigma_low=new TGraph(); g_1sigma_low->SetName("g_1sigma_low");
  g_1sigma_low->SetFillColor(3);
  for(int i=0;i<g_1up_low->GetN();i++)
  {
   double x,y;
   g_1up_low->GetPoint(i,x,y);
   g_1sigma_low->SetPoint(i,x,y);
  }
  int a=0;
  for(int i=g_1dn_low->GetN()-1;i>=0;i--)
  {
   double x,y;
   g_1dn_low->GetPoint(i,x,y);
   g_1sigma_low->SetPoint(a+g_1up_low->GetN(),x,y);
   a++;
  }
   
  TGraph* g_1sigma_high=new TGraph(); g_1sigma_high->SetName("g_1sigma_high");
  g_1sigma_high->SetFillColor(3);
  for(int i=0;i<g_1up_high->GetN();i++)
  {
   double x,y;
   g_1up_high->GetPoint(i,x,y);
   g_1sigma_high->SetPoint(i,x,y);
  }
  a=0;
  for(int i=g_1up_high->GetN()-1;i>=0;i--)
  {
   double x,y;
   g_1dn_high->GetPoint(i,x,y);
   g_1sigma_high->SetPoint(a+g_1up_high->GetN(),x,y);
   a++;
  }
  
  //2sigma
  
  TGraph* g_2sigma_low=new TGraph(); g_2sigma_low->SetName("g_2sigma_low");
  g_2sigma_low->SetFillColor(5);
  for(int i=0;i<g_2dn_low->GetN();i++)
  {
   double x,y;
   g_2dn_low->GetPoint(i,x,y);
   g_2sigma_low->SetPoint(i,x,y);
  }
  a=0;
  for(int i=g_2dn_low->GetN()-1;i>=0;i--)
  {
   double x,y;
   g_2up_low->GetPoint(i,x,y);
   g_2sigma_low->SetPoint(a+g_2dn_low->GetN(),x,y);
   a++;
  }
   
  TGraph* g_2sigma_high=new TGraph(); g_2sigma_high->SetName("g_2sigma_high");
  g_2sigma_high->SetFillColor(5);
  for(int i=0;i<g_2dn_high->GetN();i++)
  {
   double x,y;
   g_2dn_high->GetPoint(i,x,y);
   g_2sigma_high->SetPoint(i,x,y);
  }
  a=0;
  for(int i=g_2dn_high->GetN()-1;i>=0;i--)
  {
   double x,y;
   g_2up_high->GetPoint(i,x,y);
   g_2sigma_high->SetPoint(a+g_2dn_high->GetN(),x,y);
   a++;
  }
  
  g_2sigma_low->Draw("fsame");  g_2sigma_high->Draw("fsame");
  g_1sigma_low->Draw("fsame");  g_1sigma_high->Draw("fsame");

  g_exp_low->Draw("lsame");     g_exp_high->Draw("lsame");
  
  g_low_area->Draw("fsame"); 
  g_high_area->Draw("fsame");

  if(do_observed)
  {
   g_obs_low->Draw("lsame");
   g_obs_high->Draw("lsame");
  }
  
  TLegend* leg=new TLegend(0.55,0.65,0.93,0.93);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetHeader(Form("H^{+}#rightarrow%s, %s, %s",decay.c_str(),model_label.c_str(),lumi.c_str()));
  TH1D* h_obs=new TH1D("h_obs","h_obs",1,0,1);
  h_obs->SetLineWidth(3);
  h_obs->SetFillColor(1);
  if(!do_observed) h_obs->SetLineStyle(2);
  h_obs->SetFillStyle(obs_fillstyle);
  if(do_observed)
  {
   leg->AddEntry(h_obs,"#bf{Observed exclusion}","f");
   leg->AddEntry(g_exp_low,"#bf{Expected exclusion}","l");
  }
  if(!do_observed) leg->AddEntry(h_obs,"#bf{Expected exclusion}","f");
  leg->AddEntry(g_1sigma_low,"#bf{#pm 1sigma}","f");
  leg->AddEntry(g_2sigma_low,"#bf{#pm 2sigma}","f");
  leg->Draw();
  
  can->RedrawAxis();
  
  ATLAS_LABEL_internal(0.2,0.9,0.05);
  
  can->Print(Form("exclusion_%s_%s.pdf",decay.c_str(),model.c_str()));
  cout<<"Plot safed as "<<Form("exclusion_%s_%s.pdf",decay.c_str(),model.c_str())<<endl;
  
 }
 
}


