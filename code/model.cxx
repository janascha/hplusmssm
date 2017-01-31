#include <iterator>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <math.h>

#include "model.h"

string model::get_name()
{
 
 return m_name;	
 
}

void model::init_topBR()
{
 
 string filename="../WG3_inputs/"+m_name+"-BRt.out";
 
 ifstream file(filename.c_str());
 string line;
 
 while(getline(file, line))
 {
  stringstream linestream(line);
  double mHp,tanb,gamma,mA,topBR;
  linestream >> mHp >> mA >> tanb >> gamma >> topBR;
  
  topBR_val.push_back(topBR);
  topBR_mHp.push_back(mHp);
  topBR_tanb.push_back(tanb);
 
 } // while file
 
 
}


void model::init_deltab()
{
 
 string filename="../WG3_inputs/"+m_name+".dat2";
 
 ifstream file(filename.c_str());
 
 m_dodeltab=0;
 
 if(file)
 {
  m_dodeltab=1;
  
  cout<<" "<<filename<<" read"<<endl;
  string line; 
  
  while(getline(file, line))
  {
   stringstream linestream(line);
   double tanb;
   double deltab;
   linestream >> tanb >> deltab;
   
   deltab_val.push_back(deltab);
   deltab_tanb.push_back(tanb);
  
  } // while file
  
 }
 
 if(!m_dodeltab) cout<<"### INFO: No delta-b corrections will be applied for "<<m_name<<endl;
 
 
}


void model::init_BR()
{
 
 vector<string> filename;
 filename.push_back("../WG3_inputs/BR."+m_name+".LM.Hp.output");
 filename.push_back("../WG3_inputs/BR."+m_name+".Hp.output");
 filename.push_back("../WG3_inputs/BR."+m_name+"-MA2.Hp.output");
 
 BR_mHp.set_name("mHp");
 BR_tanb.set_name("tanb");
 for(int b=0;b<m_BR_list.size();b++)
 {
  BR thisBR;
  thisBR.set_name(m_BR_list[b]);
  allBR.push_back(thisBR);
 }
 BR BR_width;
 BR_width.set_name("width");
 allBR.push_back(BR_width);
 
 BR BR_mA;
 BR_mA.set_name("mA");
 allBR.push_back(BR_mA);

 for(int f=0;f<filename.size();f++)
 {  
  
  ifstream file(filename[f].c_str());
  
  if(file) cout<<" "<<filename[f]<<" read"<<endl;
  
  string line;
  
  while(getline(file, line))
  {
   stringstream   linestream(line);
   double M_A,mue,M_p,totWidth,BRenue,BRmunumu,BRtaunutau,BRtbb,BRtsb,BRtdb,BRcbb,BRcsb,BRcdb,BRubb,BRusb,BRudb,BRH0W,BRHHW,BRA0W,BRSUSY,zero;
   totWidth=BRenue=BRmunumu=BRtaunutau=BRtbb=BRtsb=BRtdb=BRcbb=BRcsb=BRcdb=BRubb=BRusb=BRudb=BRH0W=BRHHW=BRA0W=BRSUSY=zero=-1;
   string         name;
   double tanbeta;
   
   if(m_name!="hmssm")
   {
    linestream >> M_A >> tanbeta >> mue >> M_p >> totWidth >> BRenue >> BRmunumu >> BRtaunutau >>BRtbb >> BRtsb >> BRtdb >> BRcbb >> BRcsb
               >> BRcdb >> BRubb >> BRusb >> BRudb >> BRH0W >> BRHHW >> BRA0W >> BRSUSY;
   }
   
   if(m_name=="hmssm")
   {
    linestream >> tanbeta >> M_A >> M_p >> totWidth >> BRcbb >> BRtaunutau >> BRmunumu >> BRusb >> BRcsb >> BRtbb >> zero >> 
                  BRcdb >> BRubb >> BRtsb >> BRtdb >> BRH0W >> BRA0W;
    
   }
   
   BR_mHp.add_val(M_p);
   BR_tanb.add_val(tanbeta);
   
   for(int b=0;b<allBR.size();b++)
   {
    if(allBR[b].get_name()=="mA")     allBR[b].add_val(M_A);
    if(allBR[b].get_name()=="width")  allBR[b].add_val(totWidth);
    if(allBR[b].get_name()=="tb")  	  allBR[b].add_val(BRtbb);
    if(allBR[b].get_name()=="taunu")	allBR[b].add_val(BRtaunutau);
    if(allBR[b].get_name()=="cs")  	  allBR[b].add_val(BRcsb);
    if(allBR[b].get_name()=="susy") 	allBR[b].add_val(BRSUSY);
    if(allBR[b].get_name()=="cb")   	allBR[b].add_val(BRcbb);
    if(allBR[b].get_name()=="ts")   	allBR[b].add_val(BRtsb);
    if(allBR[b].get_name()=="munu") 	allBR[b].add_val(BRmunumu);
    if(allBR[b].get_name()=="enu")  	allBR[b].add_val(BRenue);
    if(allBR[b].get_name()=="cd")   	allBR[b].add_val(BRcdb);
    if(allBR[b].get_name()=="ub")   	allBR[b].add_val(BRubb);
    if(allBR[b].get_name()=="us")   	allBR[b].add_val(BRusb);
    if(allBR[b].get_name()=="ud")   	allBR[b].add_val(BRudb);
    if(allBR[b].get_name()=="hW")   	allBR[b].add_val(BRH0W);
    if(allBR[b].get_name()=="HW")   	allBR[b].add_val(BRHHW);
    if(allBR[b].get_name()=="AW")   	allBR[b].add_val(BRA0W);
   }
  
  } // while file
 
 } //for filename
 


}


model::model(string name, vector<string> BR_list)
{
 
 m_name=name;
 
 for(int b=0;b<BR_list.size();b++)
  m_BR_list.push_back(BR_list[b]);
 
}

double model::get_topBR(double mass, double tanb)
{
 
 double tanb_dn,tanb_up;
 
 if(tanb<1.0)
 {
  tanb_dn=myround_dn(tanb,1);
  tanb_up=myround_up(tanb,1);
 }
 else
 {
  tanb_dn=myround_dn(tanb,0);
  tanb_up=myround_up(tanb,0);
 }
 
 double val_dn=get_topBRval(mass, tanb_dn);
 double val_up=get_topBRval(mass, tanb_up);
 
 //cout<<"val_dn "<<val_dn<<" val_up "<<val_up<<endl;
 
 double val=-1;
 if(val_dn>0 && val_up>0) val=linear(tanb_dn,tanb_up,val_dn,val_up,tanb);
 
 return val;
 
}

double model::get_topBRval(double mass, double tanb)
{
 
 vector<double> safe_br;
 vector<double> safe_mass;
 
 for(int i=0;i<topBR_val.size();i++)
 {
  int rounded_BR_tanb=(int)(myround(topBR_tanb[i],1)*10.0);
  int rounded_input=(int)(tanb*10.0);
  
 	if(rounded_BR_tanb==rounded_input)
 	{
   safe_br.push_back(topBR_val[i]);
 	 safe_mass.push_back(topBR_mHp[i]);
 	}
 }
 
 double br1,br2,m1,m2,value;
 value=-1;
 
 if(safe_br.size()>0)
 	{
 
 //order the "safe" vectors by mass:
 vector<int> sorted_indices;
 vector<double> sorted_mass;
 vector<double> sorted_br;
 int min=get_min(safe_mass);
 sorted_indices.push_back(min);
 sorted_mass.push_back(safe_mass[min]);
 sorted_br.push_back(safe_br[min]);
 vector<int> temp_indices;
 vector<double> temp_mass;
 vector<double> temp_br;
 for(int a=0;a<safe_mass.size();a++)
 {
 	if(a!=min)
 	{
 	 temp_mass.push_back(safe_mass[a]);
 	 temp_br.push_back(safe_br[a]);
 	 temp_indices.push_back(a);
 	}
 }
 int n=temp_indices.size();
 for(int i=0;i<n;i++)
 {
  int min=get_min(temp_mass); 
  sorted_indices.push_back(temp_indices[min]);
  sorted_mass.push_back(temp_mass[min]);
  sorted_br.push_back(temp_br[min]);
 	temp_indices.erase(temp_indices.begin()+min);
 	temp_mass.erase(temp_mass.begin()+min);
 	temp_br.erase(temp_br.begin()+min);
 }
  
 int found=0;
 int seems_ok=1;
 if(sorted_mass.size()>0)
 {
  if(sorted_mass[0]>mass)
   seems_ok=0;
 }
 
 if(seems_ok)
 {
  for(int i=0;i<sorted_mass.size();i++)
  {
   if(sorted_mass[i]>=mass && i>0)
   {
    m1=sorted_mass[i-1];
    m2=sorted_mass[i];
    br1=sorted_br[i-1];
    br2=sorted_br[i];
    found=1;
   } 
   if(found) i=sorted_mass.size();
  }
  if(found) value=linear(m1,m2,br1,br2,mass);
 }
}

 return value;
 
}

double model::get_BR(string name, double mass, double tanb)
{

 int index=-1;
 for(int b=0;b<allBR.size();b++)
 {
 	if(allBR[b].get_name()==name)
 	 index=b;
 }
 
 if(index<0) cout<<"  Error: BR index for "<<name<<" not found"<<endl;
 
 double tanb_dn,tanb_up;
 
 if(tanb<1.0)
 {
  tanb_dn=myround_dn(tanb,1);
  tanb_up=myround_up(tanb,1);
 }
 else
 {
  tanb_dn=myround_dn(tanb,0);
  tanb_up=myround_up(tanb,0);
 }
 
 //cout<<"tanb "<<tanb<<" up "<<tanb_up<<" dn "<<tanb_dn<<endl;
 
 double val=-1;
 if((m_name!="tauphobic" && !(tanb_up>60 || tanb_dn>60)) || (m_name=="tauphobic" && !(tanb_up>50 || tanb_dn>50)) )
 {
  double val_dn=get_BRval(index, mass, tanb_dn);
  double val_up=get_BRval(index, mass, tanb_up);
  //cout<<"val_up "<<val_up<<" val_dn "<<val_dn<<endl;
  val=linear(tanb_dn,tanb_up,val_dn,val_up,tanb);
 }
 
 return val;
}

double model::get_BRval(int index, double mass, double tanb)
{
	
 vector<double> safe_br;
 vector<double> safe_mass;
 
 for(int i=0;i<allBR[index].get_size();i++)
 {
  int rounded_BR_tanb=0;
  if(BR_tanb.get_val(i)<1.0) rounded_BR_tanb=(int)(myround(BR_tanb.get_val(i),1)*10.0);
  else                       rounded_BR_tanb=(int)(BR_tanb.get_val(i)*10.0);
  
  int rounded_input=(int)(tanb*10.0);
  
 	if(rounded_BR_tanb==rounded_input)
 	{
   //cout<<"tanb "<<tanb<<" rounded_input "<<rounded_input<<" val "<<BR_tanb.get_val(i)<<" rounded "<<rounded_BR_tanb<<endl;
 	 safe_br.push_back(allBR[index].get_val(i));
 	 safe_mass.push_back(BR_mHp.get_val(i));
 	}
 }
 
 //if(safe_br.size()==0) cout<<"   Error: safe BR vector for index "<<index<<" has size 0. mass "<<mass<<" tanb "<<tanb<<endl;
 
 double br1,br2,m1,m2,value;
 value=-1;
  
 //order the "safe" vectors by mass:
 if(safe_br.size()>0)
 {
 vector<int> sorted_indices;
 vector<double> sorted_mass;
 vector<double> sorted_br;
 int min=get_min(safe_mass);
 sorted_indices.push_back(min);
 sorted_mass.push_back(safe_mass[min]);
 sorted_br.push_back(safe_br[min]);
 vector<int> temp_indices;
 vector<double> temp_mass;
 vector<double> temp_br;
 for(int a=0;a<safe_mass.size();a++)
 {
 	if(a!=min)
 	{
 	 temp_mass.push_back(safe_mass[a]);
 	 temp_br.push_back(safe_br[a]);
 	 temp_indices.push_back(a);
 	}
 }
 int n=temp_indices.size();
 for(int i=0;i<n;i++)
 {
  int min=get_min(temp_mass); 
  sorted_indices.push_back(temp_indices[min]);
  sorted_mass.push_back(temp_mass[min]);
  sorted_br.push_back(temp_br[min]);
 	temp_indices.erase(temp_indices.begin()+min);
 	temp_mass.erase(temp_mass.begin()+min);
 	temp_br.erase(temp_br.begin()+min);
 }
  
 int seems_ok=1;
 int found=0;
 if(sorted_mass.size()>0)
 {
 	if(sorted_mass[0]>mass)
 	 seems_ok=0;
 }
 if(seems_ok)
 {
  for(int i=0;i<sorted_mass.size();i++)
  {
   if(sorted_mass[i]>=mass && i>0)
   {
    m1=sorted_mass[i-1];
    m2=sorted_mass[i];
    br1=sorted_br[i-1];
    br2=sorted_br[i];
    found=1;
   } 
   if(found) i=sorted_mass.size();
  }
  if(found) value=linear(m1,m2,br1,br2,mass);
 }
 } //safe size >0
 
 return value;
 
}

double model::linear(double x1,double x2,double y1,double y2,double x)
{
 double y=-1;
 
 double eps=0.0000000001;
 if((x2-x1)<eps) y=y1;
 else
 {
  double m=(y2-y1)/(x2-x1);
  double n=y1-m*x1;
  y=m*x+n;
 }
 
 return y;
}


double model::myround(double x,int n)
{
 
 double wert;
 x*=(double)(pow(10.,(double)n));
 if((x-0.5)<floor(x))
  wert=floor(x)/(double)(pow(10.0,(double)n));
 else
  wert=ceil(x)/(double)(pow(10.0,(double)n));
 
 return wert;
 
}


double model::myround_dn(double x,int n)
{
 
 double wert;
 x*=(double)(pow(10.,(double)n));
 wert=floor(x)/(double)(pow(10.0,(double)n));
 
 return wert;
 
}

double model::myround_up(double x,int n)
{
 
 double wert;
 x*=(double)(pow(10.,(double)n));
 wert=ceil(x)/(double)(pow(10.0,(double)n));
 
// cout<<"x "<<x<<" n "<<n<<" (double)(pow(10.,(double)n)) "<<(double)(pow(10.,(double)n))<<" wert "<<wert<<endl;
 
 return wert;
 
}



int model::get_min(vector<double> vector)
{
 double min=100000.0;
 int index;
 for(int i=0;i<vector.size();i++)
 {
 	if(vector[i]<min)
 	{
 	 min=vector[i];
 	 index=i;
 	}
 }
 return index;
}


double model::get_deltab(double tanb)
{
 
 double deltab=-9999.0;
 
 double t1,t2,d1,d2;
 
 int found=0;
 for(int i=0;i<deltab_tanb.size();i++)
 {
 	if(deltab_tanb[i]>=tanb && i>0)
 	{
 	 t1=deltab_tanb[i-1];
   t2=deltab_tanb[i];
   d1=deltab_val[i-1];
   d2=deltab_val[i];
   found=1;
 	}
 	if(found)
 	 i=deltab_tanb.size();
 }
 
 if(found) deltab=linear(t1,t2,d1,d2,tanb);
 
 //cout<<"t1 "<<t1<<" t2 "<<t2<<" d1 "<<d1<<" d2 "<<d2<<" deltab "<<deltab<<endl;
 
 if(deltab<-9000) cout<<"### Error: Delta-b is not correct ###"<<endl;
 
 return deltab;
 
}

int model::get_maxtanb()
{
 
 int maxtanb=60;
 
 double mass=400;
 
 for(int t=40;t<=60;t++)
 {
 	if(get_BR("width", mass, (double)t)<0)
 	{
 	 maxtanb=t-1;
 	 t=61;
 	}
 }
 
 return maxtanb;
 
}


