#include <iterator>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <math.h>

#include "limit.h"

void limit::set_name(string name)
{
 m_name=name; 
}

string limit::get_name()
{
 return m_name;
}

void limit::add_val(double value)
{
 val.push_back(value);
}


double limit::get_val(int index)
{
 double thisval=-1.0;
 if(index<(int)val.size())
  thisval=val[index];	
 else
  cout<<"### PROBLEM: Limit val array for "<<m_name<<" out of range. ###"<<endl;
 return thisval;
}

void limit::add_mass(double value)
{
 mass.push_back(value);
}

void limit::add_up(double value)
{
 crossings_up.push_back(value);
}

void limit::add_dn(double value)
{
 crossings_dn.push_back(value);
}

double limit::get_mass(int index)
{
 double thisval=-1.0;
 if(index<(int)mass.size())
  thisval=mass[index];	
 else
  cout<<"### PROBLEM: Limit mass array for "<<m_name<<" out of range. ###"<<endl;
 return thisval;
}

double limit::get_up(int index)
{
 double thisval=-1.0;
 if(index<(int)crossings_up.size())
  thisval=crossings_up[index];	
 else
  cout<<"### PROBLEM: Limit crossings_up array for "<<m_name<<" out of range. ###"<<endl;
 return thisval;
}

double limit::get_dn(int index)
{
 double thisval=-1.0;
 if(index<(int)crossings_dn.size())
  thisval=crossings_dn[index];	
 else
  cout<<"### PROBLEM: Limit crossings_dn array for "<<m_name<<" out of range. ###"<<endl;
 return thisval;
}

void limit::init_names()
{
 
 string gname=m_name+"_low";
 graph_up=new TGraph(); graph_up->SetName(gname.c_str());
 gname=m_name+"_high";
 graph_dn=new TGraph(); graph_dn->SetName(gname.c_str());
 
}

