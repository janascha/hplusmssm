#include <iterator>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TLegend.h"
#include "TFile.h"
#include <math.h>
#include <stdlib.h>

#include "hplus.h"

using namespace std;

void hplus::run()
{
 
 int choice;
 while(1)
 {
  cout<<endl;
  cout<<"What do you wish to do?"<<endl;
  cout<<" (1) Simple Output"<<endl;
  cout<<" (2) Plots (A vs. B)"<<endl;
  cout<<" (3) Exclusions"<<endl;
  cout<<" (0) EXIT"<<endl;
  cout<<"Enter number: "; cin>>choice;
  if(choice==0) break;
  if(choice==1) simple_output();
  if(choice==2) a_vs_b();
  if(choice==3) int a=exclusions();
 }
 
}

void hplus::init()
{
 
 m_debug=0;	
  
 ifstream file("setup.txt");
 string line;
 
 while(getline(file, line))
 {
  
  stringstream linestream(line);
  string word1,word2,word3;
  linestream >> word1 >> word2 >> word3;
  
  if(word1=="lhc_energy")
   energy=word2;
  
  if(word1=="model" && word2=="1")
   model_list.push_back(word3);

  if(word1=="BR" && word2=="1")
   BR_list.push_back(word3);

 } // while file
 
 cout<<" -> LHC energy "<<energy<<" TeV"<<endl;
 for(int m=0;m<model_list.size();m++)
  cout<<" -> use model "<<model_list[m]<<endl;
 for(int m=0;m<BR_list.size();m++)
  cout<<" -> use BR "<<BR_list[m]<<endl;
  
}

void hplus::init_xsec()
{
 
 string filename="../WG3_inputs/xsec_"+energy+"tev.txt";
 
 ifstream file(filename.c_str());
 string line;
 
 while(getline(file, line))
 {
  
  stringstream linestream(line);
  double xsec,down,up;
  double mH,tb;
  linestream >> mH >> tb >> xsec >> down >> up;
  
  xsec_val.push_back(xsec);
  xsec_up.push_back(xsec);
  xsec_dn.push_back(xsec);
  xsec_mH.push_back(mH);
  xsec_tanb.push_back(tb);
  
 } // while file
 
 if(xsec_val.size()==0)
  cout<<"### Error: Initialisation of raw xsec failed ###"<<endl;
 
 if(energy=="8")
 {
  for(int m=200;m<=1000;m+=20)
   xsec_input_mass.push_back(m);

 for(int t=1;t<=9;t++)
  xsec_input_tanb.push_back(0.1*t);
 for(int t=1;t<=20;t++)
  xsec_input_tanb.push_back(t);
 for(int t=25;t<=60;t+=5)
  xsec_input_tanb.push_back(t);

 }

 if(energy=="13")
 {
 for(int m=200;m<=500;m+=20)
  xsec_input_mass.push_back(m);
  
 for(int m=550;m<=1000;m+=50)
  xsec_input_mass.push_back(m);

 for(int m=1100;m<=2000;m+=100)
  xsec_input_mass.push_back(m);

 for(int t=1;t<=9;t++)
  xsec_input_tanb.push_back(0.1*t);
 for(int t=1;t<=60;t+=1)
  xsec_input_tanb.push_back(t);

 }
 
 
} //init


void hplus::init_model()
{
 
 for(int m=0;m<model_list.size();m++)
 {
 
  model thismodel(model_list[m],BR_list);
  thismodel.init_BR();
  thismodel.init_deltab();
  thismodel.init_topBR();
  mymodel.push_back(thismodel); 
 
 }
 
}


void hplus::simple_output()
{
 
 double input_tanb,input_mass;
 
 cout<<endl;
 cout<<"Please enter:"<<endl;
 
 cout<<" mH+ [GeV]\t";
 cin>>input_mass;
 
 input_tanb=-1;
 while(input_tanb<0.5 || input_tanb>60)
 {
  cout<<" tanb (0.5-60)\t";
  cin>>input_tanb;
 }
 
 cout<<endl;
 
 double xsec_raw=get_xsec(input_mass,input_tanb);
 cout<<"xsec 2HDM [pb]\t"<<xsec_raw<<endl; 
 
 for(int m=0;m<mymodel.size();m++)
 {
  cout<<endl;
 	cout<<"Model "<<mymodel[m].get_name()<<":"<<endl;
 	
 	//Width:
 	cout<<" -> width [GeV]\t"<<mymodel[m].get_BR("width",input_mass,input_tanb)<<endl;
  
 	//Width:
 	cout<<" -> mA [GeV]\t"<<mymodel[m].get_BR("mA",input_mass,input_tanb)<<endl;
  
  //xsec:
  double xsec_mssm=get_xsec_mssm(input_mass,input_tanb,m);
  cout<<" -> xsec [pb]\t"<<xsec_mssm<<endl;
  
  //BR:
 	for(int b=0;b<BR_list.size();b++)
 	 cout<<" -> BR "<<BR_list[b]<<"\t"<<mymodel[m].get_BR(BR_list[b],input_mass,input_tanb)<<endl;
  
  //top->H+b BR:
  cout<<" -> BR(t->H+b)\t"<<mymodel[m].get_topBR(input_mass,input_tanb)<<endl;
  
 }
 
}

double hplus::get_xsec_mssm(double mass, double tanb, int model_index)
{
 
 double xsec_mssm=-1;
 
 int dodeltab=mymodel[model_index].m_dodeltab;
 
 if(dodeltab)
 {
 	
  double deltab=mymodel[model_index].get_deltab(tanb);
  
  double tbeff = tanb/sqrt(1.0 + deltab);
  double xsec_tbeff=get_xsec(mass,tbeff);
    
  if(xsec_tbeff>0) 
   xsec_mssm=xsec_tbeff*(1.0/(1.0 + deltab));
  
  if(m_debug) cout<<"tbeff "<<tbeff<<" deltab "<<deltab<<" xsec_tbeff "<<xsec_tbeff<<" xsec_mssm "<<xsec_mssm<<endl;

 }
 
 if(!dodeltab)
 {
 	xsec_mssm=get_xsec(mass,tanb);
 }
 
 return xsec_mssm;
 
}


double hplus::get_xsec(double mass, double tanb)
{
 
 if(m_debug) cout<<"get_xsec()"<<endl;
 
 double xsec=-1.0;
 
 if(mass>=200 && (energy=="13" && mass<=2000) || (energy=="8" && mass<=1000) )
 {
 
  double m1,m2,t1,t2,t3,xsec11,xsec12,xsec21,xsec22,xsec13,xsec23;
  xsec11=xsec12=xsec21=xsec22=xsec13=xsec23=-1;
  
  for(int m=0;m<xsec_input_mass.size();m++)
  {
   int found_m=0;
   if(xsec_input_mass[m]>=mass && m!=0)
   {
    m1=m-1;
    m2=m;
    found_m=1;
   }
   if(xsec_input_mass[m]>=mass && m==0)
   {
    m1=m;
    m2=m+1;
    found_m=1;
   }
   if(found_m)
    m=xsec_input_mass.size();
  }
  
  if(m_debug) cout<<"m1 "<<m1<<" m2 "<<m2<<endl;
  
  
  //3-point quadratic interpolation:
  if(tanb<=60.0)
  {
   for(int t=0;t<xsec_input_tanb.size();t++)
   {
    int found_t=0;
    if(xsec_input_tanb[t]>=tanb && t!=0 && t<xsec_input_tanb.size())
    {
     t1=t-1;
     t2=t;
     t3=t+1;
     found_t=1;
    }
    if(xsec_input_tanb[t]>=tanb && t==xsec_input_tanb.size()-1)
    {
     t1=t-2;
     t2=t-1;
     t3=t;
     found_t=1;
    }
    if(xsec_input_tanb[t]>=tanb && t==0)
    {
     t1=t;
     t2=t+1;
     t3=t+2;
     found_t=1;
    }
    if(found_t)
     t=xsec_input_tanb.size();
   }
   
   if(m_debug) cout<<"t1 "<<t1<<" t2 "<<t2<<" t3 "<<t3<<endl; 
   
   xsec11=get_xsec_raw_exact(xsec_input_mass[m1],xsec_input_tanb[t1]);
   xsec12=get_xsec_raw_exact(xsec_input_mass[m1],xsec_input_tanb[t2]);
   xsec21=get_xsec_raw_exact(xsec_input_mass[m2],xsec_input_tanb[t1]);
   xsec22=get_xsec_raw_exact(xsec_input_mass[m2],xsec_input_tanb[t2]);
   xsec13=get_xsec_raw_exact(xsec_input_mass[m1],xsec_input_tanb[t3]);
   xsec23=get_xsec_raw_exact(xsec_input_mass[m2],xsec_input_tanb[t3]);
   
   if(m_debug)
   {
    cout<<"xsec11 "<<xsec11<<endl;
    cout<<"xsec12 "<<xsec12<<endl;
    cout<<"xsec21 "<<xsec21<<endl;
    cout<<"xsec22 "<<xsec22<<endl;
    cout<<"xsec13 "<<xsec13<<endl;
    cout<<"xsec23 "<<xsec23<<endl;
   }
  
   if(xsec11<0 || xsec12<0 || xsec21<0 || xsec22<0 || xsec13<0 || xsec23<0)
    xsec=-1;
   else
   { 
    xsec=interpolate_quad(mass,tanb,xsec_input_mass[m1],xsec_input_tanb[t1],xsec_input_mass[m2],
                        xsec_input_tanb[t2],xsec11,xsec12,xsec21,xsec22,xsec_input_tanb[t3],xsec13,xsec23);
   }
  }  //tanb<=60
  
  
  if(tanb>60.0 && tanb<64.0) //this is intended for the mhmod- problem for tanb=60
  {
   if(m_debug) cout<<"!!! Trying to recover xsec for tanb=60"<<endl;
   
   //extrapolate linearily from 59-60
   double xsec59_m1=get_xsec_raw_exact(xsec_input_mass[m1],59.0);
   double xsec59_m2=get_xsec_raw_exact(xsec_input_mass[m2],59.0);
   double xsec59=linear(xsec_input_mass[m1],xsec_input_mass[m2],xsec59_m1,xsec59_m2,mass);
   double xsec60_m1=get_xsec_raw_exact(xsec_input_mass[m1],60.0);
   double xsec60_m2=get_xsec_raw_exact(xsec_input_mass[m2],60.0);
   double xsec60=linear(xsec_input_mass[m1],xsec_input_mass[m2],xsec60_m1,xsec60_m2,mass);
   xsec=linear(59.0,60.0,xsec59,xsec60,tanb);
   //double hplus::linear(double x1,double x2,double y1,double y2,double x)
   
   if(m_debug) cout<<"xsec_input_mass[m1] "<<xsec_input_mass[m1]<<" xsec_input_mass[m2] "<<xsec_input_mass[m2]<<" mass "<<mass<<endl;
   if(m_debug) cout<<"xsec59_m1 "<<xsec59_m1<<" xsec59_m2 "<<xsec59_m2<<" xsec59 "<<xsec59<<endl;
   if(m_debug) cout<<"xsec60_m1 "<<xsec60_m1<<" xsec60_m2 "<<xsec60_m2<<" xsec60 "<<xsec60<<endl;
   if(m_debug) cout<<"-> xsec "<<xsec<<endl;
  }  //tanb>60
  
 } //mass is valid
 
 xsec*=2.0;
 
 return xsec;
 
}

double hplus::get_xsec_raw_exact(double mass, double tanb)
{
 	
 if(m_debug) cout<<"get_xsec_raw_exact()"<<endl;

 double xsec=-1.0;
 int found=0;
 for(int m=0;m<xsec_mH.size();m++)
 {
 	if((int)xsec_mH[m]==(int)mass && (int)(xsec_tanb[m]*10.0)==(int)(tanb*10.0))
  {
 	 xsec=xsec_val[m];
 	 found=1;
 	}
 }
 
 if(!found)
 {
 	//cout<<"### Problem: raw xsec cannot be found for mass="<<mass<<" and tanb="<<tanb<<" ###"<<endl;
 	xsec=-1.0;
 }
 
 return xsec;
 
}


double hplus::linear(double x1,double x2,double y1,double y2,double x)
{
 double y=-1;
 
 double eps=0.0000000001;
 if(fabs(x2-x1)<eps)
 {
 	y=y1;
 }
 else
 {
  double m=(y2-y1)/(x2-x1);
  double n=y1-m*x1;
  y=m*x+n;
 }
 
 
 return y;
}


double hplus::interpolate_quad(double x,double y,double x1,double y1,double x2,double y2,
	                             double z11,double z12,double z21,double z22,double y3,double z13,double z23)
{
	
 if(m_debug) cout<<"interpolate_quad"<<endl;
  
 //do linear in direction of mass (x)
 double f_R1=(x2-x)/(x2-x1)*z11+(x-x1)/(x2-x1)*z21;
 double f_R2=(x2-x)/(x2-x1)*z12+(x-x1)/(x2-x1)*z22;
 double f_R3=(x2-x)/(x2-x1)*z13+(x-x1)/(x2-x1)*z23;
 
 //exact quadratic interpolation in direction of tanb:
 //http://isezen.com/2012/01/15/quadratic-interpolation-three-point/
 double A=f_R1/((y1-y2)*(y1-y3))+f_R2/((y2-y1)*(y2-y3))+f_R3/((y3-y1)*(y3-y2));
 double B=-(f_R1*(y2+y3)/((y1-y2)*(y1-y3))+f_R2*(y1+y3)/((y2-y1)*(y2-y3))+f_R3*(y1+y2)/((y3-y1)*(y3-y2)));
 double C=f_R1*y2*y3/((y1-y2)*(y1-y3))+f_R2*y1*y3/((y2-y1)*(y2-y3))+f_R3*y1*y2/((y3-y1)*(y3-y2));
 double z=A*y*y+B*y+C;
 return z;
 
}

void hplus::a_vs_b()
{
 
 cout<<endl;
 cout<<"(a1) BR(H+->..) vs mass"<<endl;
 cout<<"(a2) BR(H+->..) vs tanb"<<endl;
 cout<<"(b1) xsec vs mass"<<endl;
 cout<<"(b2) xsec vs tanb"<<endl;
 cout<<"(c1) BR(t->H+b) vs mass"<<endl;
 cout<<"(c2) BR(t->H+b) vs tanb"<<endl;
 cout<<"(d1) xsec*BR vs mass"<<endl;
 cout<<"(d2) xsec*BR vs tanb"<<endl;
 cout<<"(0) Go back to main menu"<<endl;
 
 //cout<<"(d1) Total width vs mass"<<endl;
 //cout<<"(d2) Total width vs tanb"<<endl;
 string choice="dummy";
 cout<<endl;
 while(choice!="a1" && choice!="a2" && choice!="b1" && choice!="b2" && choice!="c1" && choice!="c2" && choice!="d1" && choice!="d2" && choice!="0")
 {
  cout<<"Enter your choice: "; cin>>choice;
 }
 
 if(choice=="a1") vs_mass("BR");
 if(choice=="a2") vs_tanb("BR");
 
 if(choice=="b1") vs_mass("xsec");
 if(choice=="b2") vs_tanb("xsec");
 
 if(choice=="c1") vs_mass("topBR");
 if(choice=="c2") vs_tanb("topBR");
 
 if(choice=="d1") vs_mass("xsecBR");
 if(choice=="d2") vs_tanb("xsecBR");
 
}

void hplus::vs_tanb(string firstvari)
{
 
 double tanb_start,tanb_end,tanb_step,mass;

 cout<<endl;
 cout<<"-------- "<<firstvari<<" vs tanb ------"<<endl;
 cout<<endl;
 cout<<"Enter: "<<endl;
 mass=-1;
 cout<<" mass [GeV] \t"; cin>>mass;
 cout<<" tanb start (0.5-60) \t"; cin>>tanb_start; 
 cout<<" tanb end (0.5-60) \t"; cin>>tanb_end;
 cout<<" tanb step     \t"; cin>>tanb_step;
 cout<<endl;
 
 int verbose=1;
  
 //mass to string:
 ostringstream strs;
 strs << mass;
 string mass_string = strs.str();
 
 if(firstvari=="BR")
 {
 	for(int m=0;m<model_list.size();m++)
 	{
 	 TCanvas* can=new TCanvas("can","can",0,0,600,600);
 	 can->SetFillColor(0);
 	 TMultiGraph* mg=new TMultiGraph(); mg->SetName("mg");
 	 TLegend* leg=new TLegend(0.6,0.6,0.9,0.9);
 	 leg->SetFillStyle(0);
   vector<TGraph*> vector_graphs;
 	 for(int b=0;b<BR_list.size();b++)
 	 {
 	 	string filename="../output/"+firstvari+"_"+model_list[m]+"_"+BR_list[b]+"_mass"+mass_string+".txt";
 	 	fstream out; out.open(filename.c_str(),ios::out);
 	 	TGraph* graph=new TGraph();
 	 	string graphname="BR_"+BR_list[b];
 	 	graph->SetName(graphname.c_str());
 	  double tanb=tanb_start;
 	  int p=0;
 	  while(tanb<=tanb_end)
 	  {
 	   double BRval=mymodel[m].get_BR(BR_list[b],mass,tanb);
 	   if(verbose) cout<<model_list[m]<<" tanb = "<<tanb<<" -> BR(H+->"<<BR_list[b]<<") = "<<BRval<<endl;
 	   out<<tanb<<" "<<BRval<<endl;
 	   graph->SetPoint(p,tanb,BRval);
 	   p++;
 	   tanb+=tanb_step;
 	  }
 	  out.close();
 	 	cout<<" -> Output safed in "<<filename<<endl;
 	 	graph->SetMarkerStyle(8);
 	 	graph->SetLineWidth(2);
 	 	if((1+b)!=10) {graph->SetLineColor(1+b);graph->SetMarkerColor(1+b);}
 	 	else          {graph->SetLineColor(49);graph->SetMarkerColor(49);}
 	 	leg->AddEntry(graph,BR_list[b].c_str(),"l");
 	 	mg->Add(graph);
 	 	graph->GetXaxis()->SetTitle("tan#beta");
 	 	string ytitle="BR(H^{+}#rightarrow"+BR_list[b]+")";	graph->GetYaxis()->SetTitle(ytitle.c_str());
 	 	vector_graphs.push_back(graph);
 	 }
 	 cout<<endl;
 	 mg->Draw("al");
 	 mg->GetXaxis()->SetTitleOffset(1.4);
 	 mg->GetYaxis()->SetTitleOffset(1.2);
 	 mg->GetXaxis()->SetTitle(((string)(model_list[m]+" m_{H^{+}}="+mass_string+" GeV   tan#beta")).c_str());
 	 mg->GetYaxis()->SetTitle(firstvari.c_str());
 	 leg->Draw();
 	 string safename="../plots/"+firstvari+"_"+model_list[m]+"_mass"+mass_string+".pdf";
 	 can->Print(safename.c_str());
 	 //save the canvas into a root file:
 	 safename="../plots/"+firstvari+"_"+model_list[m]+"_mass"+mass_string+".root";
 	 TFile* file=new TFile(safename.c_str(),"RECREATE");
 	 file->Add(can);
 	 for(int g=0;g<vector_graphs.size();g++)
 	  file->Add(vector_graphs[g]);
 	 file->Write();
 	}
 } // firstvari=="BR"
 
 
 if(firstvari=="xsec")
 {
 	for(int m=0;m<model_list.size();m++)
 	{
 	 TCanvas* can=new TCanvas("can","can",0,0,600,600);
 	 can->SetFillColor(0);
 	 
 	 string filename="../output/"+firstvari+"_"+model_list[m]+"_mass"+mass_string+".txt";
 	 fstream out; out.open(filename.c_str(),ios::out);
 	 TGraph* graph=new TGraph(); graph->SetName("graph");
 	 double tanb=tanb_start;
 	 int p=0;
 	 while(tanb<=tanb_end)
 	 {
    double xsec_mssm=get_xsec_mssm(mass,tanb,m);
 	  if(verbose) cout<<model_list[m]<<" tanb = "<<tanb<<" -> xsec = "<<xsec_mssm<<" pb"<<endl;
 	  out<<tanb<<" "<<xsec_mssm<<endl;
 	  graph->SetPoint(p,tanb,xsec_mssm);
 	  p++;
 	  tanb+=tanb_step;
 	 }
 	 out.close();
   
 	 cout<<" -> Output safed in "<<filename<<endl;
 	 graph->SetMarkerStyle(8);
 	 graph->SetLineWidth(2);
 	 
 	 graph->Draw("apl");
   graph->GetXaxis()->SetTitleOffset(1.4);
 	 graph->GetYaxis()->SetTitleOffset(1.2);
 	 graph->GetXaxis()->SetTitle(((string)(model_list[m]+" m_{H^{+}}="+mass_string+" GeV   tan#beta")).c_str());
 	 graph->GetYaxis()->SetTitle(firstvari.c_str());
 	 
 	 string safename="../plots/"+firstvari+"_"+model_list[m]+"_mass"+mass_string+".pdf";
 	 can->Print(safename.c_str());
 	 
 	 //save the canvas into a root file:
 	 safename="../plots/"+firstvari+"_"+model_list[m]+"_mass"+mass_string+".root";
 	 TFile* file=new TFile(safename.c_str(),"RECREATE");
 	 file->Add(can);
 	 file->Add(graph);
 	 file->Write();
 	} //for model
 } //firstvari xsec
 
 if(firstvari=="topBR")
 {
 	for(int m=0;m<model_list.size();m++)
 	{
 	 TCanvas* can=new TCanvas("can","can",0,0,600,600);
 	 can->SetFillColor(0);
 	 
 	 string filename="../output/"+firstvari+"_"+model_list[m]+"_mass"+mass_string+".txt";
 	 fstream out; out.open(filename.c_str(),ios::out);
 	 TGraph* graph=new TGraph(); graph->SetName("graph");
 	 double tanb=tanb_start;
 	 int p=0;
 	 while(tanb<=tanb_end)
 	 {
 	  double topBRval=mymodel[m].get_topBR(mass,tanb);
 	  if(verbose) cout<<model_list[m]<<" tanb = "<<tanb<<" -> BR(t+->H+b) = "<<topBRval<<endl;
 	  out<<tanb<<" "<<topBRval<<endl;
 	  graph->SetPoint(p,tanb,topBRval);
 	  p++;
 	  tanb+=tanb_step;
 	 }
 	 out.close();
   
 	 cout<<" -> Output safed in "<<filename<<endl;
 	 graph->SetMarkerStyle(8);
 	 graph->SetLineWidth(2);
 	 
 	 graph->Draw("apl");
   graph->GetXaxis()->SetTitleOffset(1.4);
 	 graph->GetYaxis()->SetTitleOffset(1.2);
 	 graph->GetXaxis()->SetTitle(((string)(model_list[m]+" m_{H^{+}}="+mass_string+" GeV   tan#beta")).c_str());
 	 graph->GetYaxis()->SetTitle(firstvari.c_str());
 	 
 	 string safename="../plots/"+firstvari+"_"+model_list[m]+"_mass"+mass_string+".pdf";
 	 can->Print(safename.c_str());
 	 
 	 //save the canvas into a root file:
 	 safename="../plots/"+firstvari+"_"+model_list[m]+"_mass"+mass_string+".root";
 	 TFile* file=new TFile(safename.c_str(),"RECREATE");
 	 file->Add(can);
 	 file->Add(graph);
 	 file->Write();
 	} //for model
 } //firstvari topBR
 
 
 if(firstvari=="xsecBR")
 {
 	for(int m=0;m<model_list.size();m++)
 	{
 	 TCanvas* can=new TCanvas("can","can",0,0,600,600);
 	 can->SetFillColor(0);
 	 TMultiGraph* mg=new TMultiGraph(); mg->SetName("mg");
 	 TLegend* leg=new TLegend(0.6,0.6,0.9,0.9);
 	 leg->SetFillStyle(0);
   vector<TGraph*> vector_graphs;
   
 	 for(int b=0;b<BR_list.size();b++)
 	 {
 	 	string filename="../output/"+firstvari+"_"+model_list[m]+"_"+BR_list[b]+"_mass"+mass_string+".txt";
 	 	fstream out; out.open(filename.c_str(),ios::out);
 	 	TGraph* graph=new TGraph();
 	 	string graphname="xsecBR_"+BR_list[b];
 	 	graph->SetName(graphname.c_str());
 	  double tanb=tanb_start;
 	  int p=0;
 	  while(tanb<=tanb_end)
 	  {
     //get the xsec:
     double xsec_mssm=get_xsec_mssm(mass,tanb,m);
 	   double BRval=mymodel[m].get_BR(BR_list[b],mass,tanb);
 	   double value=xsec_mssm*BRval;
 	   if(xsec_mssm<0 || BRval<0) value=-1.0;
 	   if(verbose) cout<<model_list[m]<<" tanb = "<<tanb<<" -> xsec_mssm*BR(H+->"<<BR_list[b]<<") = "<<value<<endl;
 	   out<<tanb<<" "<<value<<endl;
 	   graph->SetPoint(p,tanb,value);
 	   p++;
 	   tanb+=tanb_step;
 	  }
 	  out.close();
 	 	cout<<" -> Output safed in "<<filename<<endl;
 	 	graph->SetMarkerStyle(8);
 	 	graph->SetLineWidth(2);
 	 	if((1+b)!=10) {graph->SetLineColor(1+b);graph->SetMarkerColor(1+b);}
 	 	else          {graph->SetLineColor(49);graph->SetMarkerColor(49);}
 	 	leg->AddEntry(graph,BR_list[b].c_str(),"lp");
 	 	mg->Add(graph);
 	 	graph->GetXaxis()->SetTitle("tan#beta");
 	 	string ytitle="xsec*BR(H^{+}#rightarrow"+BR_list[b]+")";	graph->GetYaxis()->SetTitle(ytitle.c_str());
 	 	vector_graphs.push_back(graph);
 	 }
 	 cout<<endl;
 	 mg->Draw("alp");
 	 mg->GetXaxis()->SetTitleOffset(1.4);
 	 mg->GetYaxis()->SetTitleOffset(1.2);
 	 mg->GetXaxis()->SetTitle(((string)(model_list[m]+" m_{H^{+}}="+mass_string+" GeV   tan#beta")).c_str());
 	 mg->GetYaxis()->SetTitle(firstvari.c_str());
 	 leg->Draw();
 	 string safename="../plots/"+firstvari+"_"+model_list[m]+"_mass"+mass_string+".pdf";
 	 can->Print(safename.c_str());
 	 //save the canvas into a root file:
 	 safename="../plots/"+firstvari+"_"+model_list[m]+"_mass"+mass_string+".root";
 	 TFile* file=new TFile(safename.c_str(),"RECREATE");
 	 file->Add(can);
 	 for(int g=0;g<vector_graphs.size();g++)
 	  file->Add(vector_graphs[g]);
 	 file->Write();
 	}
 } // firstvari=="xsecBR"
 
 
}


void hplus::vs_mass(string firstvari)
{
 
 double mass_start,mass_end,mass_step,tanb;

 cout<<endl;
 cout<<"-------- "<<firstvari<<" vs mass ------"<<endl;
 cout<<endl;
 cout<<"Enter: "<<endl;
 tanb=-1;
 while(tanb<0.5 || tanb>60)
 {
 	cout<<" tanbeta (0.5-60) \t"; cin>>tanb;
 }
 cout<<" mass start [GeV] \t"; cin>>mass_start; 
 cout<<" mass end [GeV] \t"; cin>>mass_end;
 cout<<" mass step [GeV] \t"; cin>>mass_step;
 cout<<endl;
 
 int verbose=1;
 
 //tanb to string:
 ostringstream strs;
 strs << tanb;
 string tanb_string = strs.str();
 
 if(firstvari=="BR")
 {
 	for(int m=0;m<model_list.size();m++)
 	{
 	 TCanvas* can=new TCanvas("can","can",0,0,600,600);
 	 can->SetFillColor(0);
 	 TMultiGraph* mg=new TMultiGraph(); mg->SetName("mg");
 	 TLegend* leg=new TLegend(0.6,0.6,0.9,0.9);
 	 leg->SetFillStyle(0);
   vector<TGraph*> vector_graphs;
 	 for(int b=0;b<BR_list.size();b++)
 	 {
 	 	string filename="../output/"+firstvari+"_"+model_list[m]+"_"+BR_list[b]+"_tanb"+tanb_string+".txt";
 	 	fstream out; out.open(filename.c_str(),ios::out);
 	 	TGraph* graph=new TGraph();
 	 	string graphname="BR_"+BR_list[b];
 	 	graph->SetName(graphname.c_str());
 	  double mass=mass_start;
 	  int p=0;
 	  while(mass<=mass_end)
 	  {
 	   double BRval=mymodel[m].get_BR(BR_list[b],mass,tanb);
 	   if(verbose) cout<<model_list[m]<<" mass = "<<mass<<" GeV -> BR(H+->"<<BR_list[b]<<") = "<<BRval<<endl;
 	   out<<mass<<" "<<BRval<<endl;
 	   graph->SetPoint(p,mass,BRval);
 	   p++;
 	   mass+=mass_step;
 	  }
 	  out.close();
 	 	cout<<" -> Output safed in "<<filename<<endl;
 	 	graph->SetMarkerStyle(8);
 	 	graph->SetLineWidth(2);
 	 	if((1+b)!=10) {graph->SetLineColor(1+b);graph->SetMarkerColor(1+b);}
 	 	else          {graph->SetLineColor(49);graph->SetMarkerColor(49);}
 	 	leg->AddEntry(graph,BR_list[b].c_str(),"lp");
 	 	mg->Add(graph);
 	 	graph->GetXaxis()->SetTitle("m_H^{+} [GeV]");
 	 	string ytitle="BR(H^{+}#rightarrow"+BR_list[b]+")";	graph->GetYaxis()->SetTitle(ytitle.c_str());
 	 	vector_graphs.push_back(graph);
 	 }
 	 cout<<endl;
 	 mg->Draw("alp");
 	 mg->GetXaxis()->SetTitleOffset(1.4);
 	 mg->GetYaxis()->SetTitleOffset(1.2);
 	 mg->GetXaxis()->SetTitle(((string)(model_list[m]+" tanb="+tanb_string+"   m_{H^{+}} [GeV]")).c_str());
 	 mg->GetYaxis()->SetTitle(firstvari.c_str());
 	 leg->Draw();
 	 string safename="../plots/"+firstvari+"_"+model_list[m]+"_tanb"+tanb_string+".pdf";
 	 can->Print(safename.c_str());
 	 //save the canvas and the tgraphs into a root file:
 	 safename="../plots/"+firstvari+"_"+model_list[m]+"_tanb"+tanb_string+".root";
 	 TFile* file=new TFile(safename.c_str(),"RECREATE");
 	 file->Add(can);
 	 for(int g=0;g<vector_graphs.size();g++)
 	  file->Add(vector_graphs[g]);
 	 file->Write();
 	}
 } // firstvari=="BR"
 
 
 if(firstvari=="xsec")
 {
 	for(int m=0;m<model_list.size();m++)
 	{
 	 TCanvas* can=new TCanvas("can","can",0,0,600,600);
 	 can->SetFillColor(0);
 	 
 	 string filename="../output/"+firstvari+"_"+model_list[m]+"_tanb"+tanb_string+".txt";
 	 fstream out; out.open(filename.c_str(),ios::out);
 	 TGraph* graph=new TGraph(); graph->SetName("graph");
 	 double mass=mass_start;
 	 int p=0;
 	 while(mass<=mass_end)
 	 {
    double xsec_mssm=get_xsec_mssm(mass,tanb,m);
 	  if(verbose) cout<<model_list[m]<<" mass = "<<mass<<" GeV -> xsec = "<<xsec_mssm<<" pb"<<endl;
 	  out<<mass<<" "<<xsec_mssm<<endl;
 	  graph->SetPoint(p,mass,xsec_mssm);
 	  p++;
 	  mass+=mass_step;
 	 }
 	 out.close();
   
 	 cout<<" -> Output safed in "<<filename<<endl;
 	 graph->SetMarkerStyle(8);
 	 graph->SetLineWidth(2);
 	 
 	 graph->Draw("apl");
   graph->GetXaxis()->SetTitleOffset(1.4);
 	 graph->GetYaxis()->SetTitleOffset(1.2);
 	 graph->GetXaxis()->SetTitle(((string)(model_list[m]+" tanb="+tanb_string+"   m_{H^{+}} [GeV]")).c_str());
 	 graph->GetYaxis()->SetTitle(firstvari.c_str());
 	 
 	 string safename="../plots/"+firstvari+"_"+model_list[m]+"_tanb"+tanb_string+".pdf";
 	 can->Print(safename.c_str());
 	 
 	 //save the canvas into a root file:
 	 safename="../plots/"+firstvari+"_"+model_list[m]+"_tanb"+tanb_string+".root";
 	 TFile* file=new TFile(safename.c_str(),"RECREATE");
 	 file->Add(can);
 	 file->Add(graph);
 	 file->Write();
 	} //for model
 } //firstvari xsec
 
 if(firstvari=="topBR")
 {
 	for(int m=0;m<model_list.size();m++)
 	{
 	 TCanvas* can=new TCanvas("can","can",0,0,600,600);
 	 can->SetFillColor(0);
 	 
 	 string filename="../output/"+firstvari+"_"+model_list[m]+"_tanb"+tanb_string+".txt";
 	 fstream out; out.open(filename.c_str(),ios::out);
 	 TGraph* graph=new TGraph(); graph->SetName("graph");
 	 double mass=mass_start;
 	 int p=0;
 	 while(mass<=mass_end)
 	 {
 	  double topBRval=mymodel[m].get_topBR(mass,tanb);
 	  if(verbose) cout<<model_list[m]<<" mass = "<<mass<<" GeV -> BR(t+->H+b) "<<topBRval<<endl;
 	  out<<mass<<" "<<topBRval<<endl;
 	  graph->SetPoint(p,mass,topBRval);
 	  p++;
 	  mass+=mass_step;
 	 }
 	 out.close();
   
 	 cout<<" -> Output safed in "<<filename<<endl;
 	 graph->SetMarkerStyle(8);
 	 graph->SetLineWidth(2);
 	 
 	 graph->Draw("apl");
   graph->GetXaxis()->SetTitleOffset(1.4);
 	 graph->GetYaxis()->SetTitleOffset(1.2);
 	 graph->GetXaxis()->SetTitle(((string)(model_list[m]+" tanb="+tanb_string+"   m_{H^{+}} [GeV]")).c_str());
 	 graph->GetYaxis()->SetTitle(firstvari.c_str());
 	 
 	 string safename="../plots/"+firstvari+"_"+model_list[m]+"_tanb"+tanb_string+".pdf";
 	 can->Print(safename.c_str());
 	 
 	 //save the canvas into a root file:
 	 safename="../plots/"+firstvari+"_"+model_list[m]+"_tanb"+tanb_string+".root";
 	 TFile* file=new TFile(safename.c_str(),"RECREATE");
 	 file->Add(can);
 	 file->Add(graph);
 	 file->Write();
 	} //for model
 } //firstvari topBR
 
 
 if(firstvari=="xsecBR")
 {
 	for(int m=0;m<model_list.size();m++)
 	{
 	 TCanvas* can=new TCanvas("can","can",0,0,600,600);
 	 can->SetFillColor(0);
 	 TMultiGraph* mg=new TMultiGraph(); mg->SetName("mg");
 	 TLegend* leg=new TLegend(0.6,0.6,0.9,0.9);
 	 leg->SetFillStyle(0);
   vector<TGraph*> vector_graphs;
   
 	 for(int b=0;b<BR_list.size();b++)
 	 {
 	 	string filename="../output/"+firstvari+"_"+model_list[m]+"_"+BR_list[b]+"_tanb"+tanb_string+".txt";
 	 	fstream out; out.open(filename.c_str(),ios::out);
 	 	TGraph* graph=new TGraph();
 	 	string graphname="xsecBR_"+BR_list[b];
 	 	graph->SetName(graphname.c_str());
 	  double mass=mass_start;
 	  int p=0;
 	  while(mass<=mass_end)
 	  {
     double xsec_mssm=get_xsec_mssm(mass,tanb,m);
 	   double BRval=mymodel[m].get_BR(BR_list[b],mass,tanb);
 	   double value=xsec_mssm*BRval;
 	   if(xsec_mssm<0 || BRval<0) value=-1.0;
 	   if(verbose) cout<<model_list[m]<<" mass = "<<mass<<" GeV -> xsec*BR(H+->"<<BR_list[b]<<") = "<<value<<endl;
 	   out<<mass<<" "<<value<<endl;
 	   graph->SetPoint(p,mass,value);
 	   p++;
 	   mass+=mass_step;
 	  }
 	  out.close();
 	 	cout<<" -> Output safed in "<<filename<<endl;
 	 	graph->SetMarkerStyle(8);
 	 	graph->SetLineWidth(2);
 	 	if((1+b)!=10) {graph->SetLineColor(1+b);graph->SetMarkerColor(1+b);}
 	 	else          {graph->SetLineColor(49);graph->SetMarkerColor(49);}
 	 	leg->AddEntry(graph,BR_list[b].c_str(),"lp");
 	 	mg->Add(graph);
 	 	graph->GetXaxis()->SetTitle("m_H^{+} [GeV]");
 	 	string ytitle="xsec*BR(H^{+}#rightarrow"+BR_list[b]+")";	graph->GetYaxis()->SetTitle(ytitle.c_str());
 	 	vector_graphs.push_back(graph);
 	 }
 	 cout<<endl;
 	 mg->Draw("alp");
 	 mg->GetXaxis()->SetTitleOffset(1.4);
 	 mg->GetYaxis()->SetTitleOffset(1.2);
 	 mg->GetXaxis()->SetTitle(((string)(model_list[m]+" tanb="+tanb_string+"   m_{H^{+}} [GeV]")).c_str());
 	 mg->GetYaxis()->SetTitle(firstvari.c_str());
 	 leg->Draw();
 	 string safename="../plots/"+firstvari+"_"+model_list[m]+"_tanb"+tanb_string+".pdf";
 	 can->Print(safename.c_str());
 	 //save the canvas and the tgraphs into a root file:
 	 safename="../plots/"+firstvari+"_"+model_list[m]+"_tanb"+tanb_string+".root";
 	 TFile* file=new TFile(safename.c_str(),"RECREATE");
 	 file->Add(can);
 	 for(int g=0;g<vector_graphs.size();g++)
 	  file->Add(vector_graphs[g]);
 	 file->Write();
 	}
 } // firstvari=="BR"
 
 
 
}


bool hplus::fexists(string filename)
{
  ifstream ifile(filename.c_str());
  return ifile;
}


int hplus::exclusions()
{
 
 int do_obs=0;
 
 string limitname;
 cout<<"---- Exclusions ----"<<endl;
 cout<<endl;
 cout<<"Your decay mode:"<<endl;
 for(int b=0;b<BR_list.size();b++)
 {
  cout<<" ("<<b+1<<") "<<BR_list[b]<<endl;
 }
 int decay=0;
 while(decay<1 || decay > BR_list.size())
 {
  cout<<" Enter number: "; cin>>decay;
 }
 decay-=1;
 cout<<endl;
 cout<<"Provide the NAME of your input file (limit_NAME.txt): ";
 cin>>limitname;
 cout<<endl;
 
 string filename="../limit_input/limit_"+limitname+".txt";
 
 ifstream file(filename.c_str());
 
 if(fexists(filename))
  cout<<"loading "<<filename<<" successful"<<endl;
 if(!fexists(filename))
 {
  cout<<"ERROR: Loading "<<filename<<" failed!"<<endl;
  return 0;
 }
 
 cout<<endl;
 
 vector<double> mass;
 vector<double> obs;
 vector<double> exp;
 vector<double> sigma2up;
 vector<double> sigma1up;
 vector<double> sigma1dn;
 vector<double> sigma2dn;
 
 string line;
 while(getline(file, line))
 {
  
  string comment="#";
  size_t found_nocomment = line.find(comment);
  
  string obs_string="doobserved";
  size_t found_obs = line.find(obs_string);
    
  if(found_nocomment)
  {
   
   stringstream linestream(line);

   if(!found_obs)
   {
   	string word1,word2;
    linestream >> word1 >> word2;
    if(word1=="doobserved" && word2=="1")
     do_obs=1;
   }
   
   if(found_obs)
   {
    double this_mass, this_obs, this_exp, this_sigma2up, this_sigma1up, this_sigma1dn, this_sigma2dn;
    this_mass=this_obs=this_exp=this_sigma2up=this_sigma1up=this_sigma1dn=this_sigma2dn=-1;
    
    if(do_obs)  linestream >> this_mass >> this_exp >> this_sigma2up >> this_sigma1up >> this_sigma1dn >> this_sigma2dn  >> this_obs;
    if(!do_obs) linestream >> this_mass >> this_exp >> this_sigma2up >> this_sigma1up >> this_sigma1dn >> this_sigma2dn;
    
    mass.push_back(this_mass);
    obs.push_back(this_obs);
    exp.push_back(this_exp);
    sigma1up.push_back(this_sigma1up);
    sigma2up.push_back(this_sigma2up);
    sigma1dn.push_back(this_sigma1dn);
    sigma2dn.push_back(this_sigma2dn);
   }
   
  }
  
 } // while file
 
 //construct the limit objects:
 limit thislimit_obs;
 thislimit_obs.set_name("observed");
 thislimit_obs.init_names();
 for(int i=0;i<obs.size();i++)
 {
  thislimit_obs.add_val(obs[i]);
 }
 if(do_obs) mylimit.push_back(thislimit_obs);
 
 limit thislimit_exp;
 thislimit_exp.set_name("expected");
 thislimit_exp.init_names();
 for(int i=0;i<exp.size();i++) thislimit_exp.add_val(exp[i]);
 mylimit.push_back(thislimit_exp);
 
 limit thislimit_2up;
 thislimit_2up.set_name("sigma2up");
 thislimit_2up.init_names();
 for(int i=0;i<sigma2up.size();i++) thislimit_2up.add_val(sigma2up[i]);
 mylimit.push_back(thislimit_2up);
 
 limit thislimit_1up;
 thislimit_1up.set_name("sigma1up");
 thislimit_1up.init_names();
 for(int i=0;i<sigma1up.size();i++) thislimit_1up.add_val(sigma1up[i]);
 mylimit.push_back(thislimit_1up);
 
 limit thislimit_1dn;
 thislimit_1dn.set_name("sigma1dn");
 thislimit_1dn.init_names();
 for(int i=0;i<sigma1dn.size();i++) thislimit_1dn.add_val(sigma1dn[i]);
 mylimit.push_back(thislimit_1dn);
 
 limit thislimit_2dn;
 thislimit_2dn.set_name("sigma2dn");
 thislimit_2dn.init_names();
 for(int i=0;i<sigma2dn.size();i++) thislimit_2dn.add_val(sigma2dn[i]);
 mylimit.push_back(thislimit_2dn);
 
 for(int m=0;m<model_list.size();m++)
 {
  
  for(int l=0;l<mylimit.size();l++)
  {
   
   for(int i=0;i<mass.size();i++)
   {
 	  
 	  cout<<endl;
    cout<<"***** Now calculating exclusions of "<<mylimit[l].get_name()<<" for "<<model_list[m]<<" at mass "<<mass[i]<<" GeV *****"<<endl;
    cout<<endl;
 	  double crossing_up,crossing_dn;
 	  
 	  get_crossings(m,decay,mass[i],mylimit[l].get_val(i),crossing_up,crossing_dn, mylimit[l].get_name());
 	  
 	  mylimit[l].add_mass(mass[i]);
 	  mylimit[l].add_up(crossing_up);
 	  mylimit[l].add_dn(crossing_dn);
 	  
 	  //save the results into the limit object, directly as TGraphs
 	  mylimit[l].graph_up->SetPoint(i,mass[i],crossing_up);
 	  mylimit[l].graph_dn->SetPoint(i,mass[i],crossing_dn);
 	  
   } //Mass
   
  } //limit
  
  //save the 10 or 12 TGraphs into one file
  string filename="../output/exclusion_"+limitname+"_"+model_list[m]+"_"+BR_list[decay]+".root";
  cout<<"saving graphs into "<<filename<<endl;
  TFile* file=new TFile(filename.c_str(),"RECREATE");
  for(int l=0;l<mylimit.size();l++)
  {
   file->Add(mylimit[l].graph_up);
   file->Add(mylimit[l].graph_dn);
  }
  file->Write();
  
 } //Model
 
 return 1;
 
}


void hplus::get_crossings(int model_index, int decay_index, double mass, double limit, double &crossing_up, double &crossing_dn, string limitname)
{
 
 //Def up crossing: Ratio changes from <1 tp >1 at this tanb
 //Def dn crossing: Ratio changes from >1 tp -1 at this tanb
 
 //Low mass does not work yet!!!
 int LM=0;
 if(mass<200)
 {
  LM=1;
  cout<<"Low mass exclusions do not work yet. Exiting Program."<<endl;
  exit(1);
 }
 
 //scan fine for tanb<3, coarser in tanb>3
 vector<double> tanb_grid;
 double val=0.5;
 while(val<=1.0)
 {
 	if(!LM) tanb_grid.push_back(val);
  val+=0.01;	
 }
 //LM starts only at tanb=1
 while(val<=3.0)
 {
 	tanb_grid.push_back(val);
  val+=0.1;	
 }
 
 int maxtanb=mymodel[model_index].get_maxtanb();
 
 for(int t=3;t<=maxtanb;t++)
 {
 	tanb_grid.push_back((double)t);
 }
 
 int debug=0;
 
 double prev_ratio=999;
 int found_up=0;
 int found_dn=0;
 int write_rfile=0;
 TGraph* g_ratio=new TGraph(); g_ratio->SetName("g_ratio");
 for(int t=0;t<tanb_grid.size();t++)
 {
 	
 	double BR=mymodel[model_index].get_BR(BR_list[decay_index],mass,tanb_grid[t]);
  
  double prod;
  if(mass>199)
  {
   prod=BR*get_xsec_mssm(mass,tanb_grid[t],model_index);
  }
  else
  {
   prod=BR*mymodel[model_index].get_topBR(mass,tanb_grid[t]);
   //cout<<"tanb="<<tanb_grid[t]<<" BR "<<BR<<" topBR "<<mymodel[model_index].get_topBR(mass,tanb_grid[t])<<" prod "<<prod<<endl;
  }
  
  double ratio=9999999;
  if(t>0) ratio=prev_ratio;
  if(prod>0) ratio=limit/prod;
    
  //cout<<"t "<<t<<" tanb "<<tanb_grid[t]<<" prod "<<prod<<" ratio "<<ratio<<endl;
  
  g_ratio->SetPoint(t,tanb_grid[t],ratio);
  if(t>0)
  {
   if(ratio<1 && prev_ratio>1)
   {
   	//get the crossing from linear extrapolation:
   	double crossing=linear(prev_ratio,ratio,tanb_grid[t-1],tanb_grid[t],1.0);
   	crossing_dn=crossing;
   	cout<<" Down crossing found at tanb="<<crossing<<endl;
   	found_dn=1;
   }
   if(ratio>1 && prev_ratio<1)
   {
   	//get the crossing from linear extrapolation:
   	double crossing=linear(prev_ratio,ratio,tanb_grid[t-1],tanb_grid[t],1.0);
   	crossing_up=crossing;
   	cout<<" Up crossing found at tanb="<<crossing<<endl;
   	found_up=1;
   }
  }
  
  prev_ratio=ratio;
  if(!LM && !write_rfile && found_up && found_dn)
   t=tanb_grid.size();
  
 } //tanb scan
 
 if(write_rfile)
 {
 	string rfilename="rfile_"+limitname+"_"+".root";
  TFile* rfile=new TFile(rfilename.c_str(),"RECREATE");
  rfile->Add(g_ratio);
  rfile->Write();
 }
 
 if(!found_up)
 {
 	//the script requires that for every mas spoint we find an up and a down crossing
 	//the treatment now will dependon the decay and the mass
 	
 	if(BR_list[decay_index]=="taunu")
 	{
 	 if(!LM)
    crossing_up=-1;
 	 if(LM)
    crossing_up=8.0;
 	 cout<<"setting up crosssing to "<<crossing_up<<endl;
 	}
 	
 	if(BR_list[decay_index]=="tb")
 	{
 	 //attempting extrapolation to lower tanb
 	 
 	 cout<<"Attention: No up crossing for mass = "<<mass<<" GeV found. Attempting extrapolation to low tanb."<<endl;
 	 
 	 TGraph *graph=new TGraph(); graph->SetName("graph");
 	 int p=0;
 	 for(int t=0;t<tanb_grid.size();t++)
 	 {
 	  if(tanb_grid[t]<1.0)
 	  {
 	 	 
 	 	 double BR=mymodel[model_index].get_BR(BR_list[decay_index],mass,tanb_grid[t]);
     double prod;
     prod=BR*get_xsec_mssm(mass,tanb_grid[t],model_index);
 	 	 
 	 	 if(prod>0)
 	 	 {
 	 	  double ratio=limit/prod;
 	 	  graph->SetPoint(p,tanb_grid[t],ratio);
 	 	  p++;
 	 	 }
 	  }
 	 }
 	 
 	 graph->SetMarkerStyle(8);
   
 	 TCanvas* can=new TCanvas("can","can",0,0,600,600);
 	 graph->Draw("ap");
 	 TF1* fit=new TF1("fit","pol3",0.5,1);
   graph->Fit("fit","R");
   if(debug)
   {
    can->Print("can.pdf");
    TFile* file=new TFile("file.root","RECREATE");
    file->Add(can);
    file->Write();
   }
   
   double a,b,c,d,e;
   a=fit->GetParameter(0);
   b=fit->GetParameter(1);
   c=fit->GetParameter(2);
   d=fit->GetParameter(3);
   
   double result=-1;
   double prev_val;
   int found_extrap_up=0;
   int count=0;
   int i=0;
   double tb=0.6;
   while(tb>=0)
   {
    double val=a+b*tb+c*pow(tb,2)+d*pow(tb,3);
    i++;
    if(count>0)
    {
   	 if(val<1 && prev_val>=1) //UP crossing, but going backwards from larger to smaller tanb
     {
      cout<<"prev_val "<<prev_val<<" val "<<val<<" tb "<<tb<<endl;
      result=linear(prev_val,val,tb+0.01,tb,1.0);
      found_extrap_up=1;
      tb=-1;
     }
    }
    tb-=0.01;
    count++;
    prev_val=val;
   }
   if(result>0.5 || result<0) result=-1;
 	 cout<<endl;
 	 cout<<"Found up crossing at "<<result<< " (extrapolated)"<<endl;
 	 crossing_up=result;
  } //if tb
  
 	 
 } // if!found_up
 	
 if(!found_dn)
 {
 	
 	if(LM)
 	{
 	 crossing_dn=8;
 	 cout<<"setting dn crosssing to "<<crossing_dn<<endl;
  }
 	
 	if(!LM)
 	{
   
   cout<<"Attention: No dn crossing for mass = "<<mass<<" GeV found. Attempting extrapolation to high tanb."<<endl;
 	 
 	 TGraph *graph2=new TGraph(); graph2->SetName("graph2");
 	 int p=0;
 	 //find the point when prod becomes 0 (thats when the model is no longer well defined)
 	 int tb_negprod=tanb_grid.size()-1;
 	 for(int t=0;t<tanb_grid.size();t++)
 	 {
    double prod;
    double BR=mymodel[model_index].get_BR(BR_list[decay_index],mass,tanb_grid[t]);
    prod=BR*get_xsec_mssm(mass,tanb_grid[t],model_index);
 	 	if(prod<0)
 	 	{
 	   tb_negprod=t;
 	   t=tanb_grid.size();
 	  }
 	 }
 	 
 	 //do a simple linear extrapolation using the last 2 points!
 	 
 	 double tanb1=tanb_grid[tb_negprod-2];
 	 double tanb2=tanb_grid[tb_negprod-1];
 	 double prod1=mymodel[model_index].get_BR(BR_list[decay_index],mass,tanb1)*get_xsec_mssm(mass,tanb1,model_index);
 	 double prod2=mymodel[model_index].get_BR(BR_list[decay_index],mass,tanb2)*get_xsec_mssm(mass,tanb2,model_index);
   double extrap_tanb=linear(prod1/limit,prod2/limit,tanb1,tanb2,1.0);
   
   if(debug) cout<<"-----------> tanb1 "<<tanb1<<" tanb2 "<<tanb2<<" prod1 "<<prod1<<" prod2 "<<prod2<<" extrap_tanb "<<extrap_tanb<<endl;
   
   int found_extrap_dn=0;
   if(extrap_tanb<1000) found_extrap_dn=1;
   double result=extrap_tanb;
   if(!found_extrap_dn) result=1000.0;
   cout<<endl;
 	 cout<<"Found dn crossing at "<<result<< " (extrapolated)"<<endl;
 	 
 	 crossing_dn=result;
 	 
 	 
  } //if !LM
  
 } //if !found_dn
 
}


