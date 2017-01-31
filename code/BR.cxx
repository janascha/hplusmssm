#include <iterator>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <math.h>

#include "BR.h"

void BR::set_name(string name)
{
 m_name=name; 
}

string BR::get_name()
{
 return m_name;
}

void BR::add_val(double value)
{
 
 val.push_back(value);
 
}

int BR::get_size()
{
 return val.size();		
}

double BR::get_val(int index)
{
 double thisval=-1.0;
 if(index<val.size())
  thisval=val[index];	
 else
  cout<<"   Problem: BR array for "<<m_name<<" out of range."<<endl;
 return thisval;
}


