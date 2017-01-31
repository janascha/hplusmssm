#include <iostream>
#include <iomanip>
#include <fstream>

#include "hplus.h"

using namespace std;

/*
---------------------------------
Hplus MSSM Interpretation Package
---------------------------------
This code serves as an interface to the MSSM H+ inputs from the LHCHXSWG.
Twiki:  https://twiki.cern.ch/twiki/bin/view/AtlasProtected/HplusMSSMInterpretation
Author: jana.schaarschmidt@cern.ch
*/

int main()
{
	
 cout<<endl;
 cout<<"-------------------------"<<endl;
 cout<<"Hplus MSSM interface v0.5"<<endl;
 cout<<"-------------------------"<<endl;
 cout<<endl;
 
 hplus myhplus;
 cout<<"Initialize..."<<endl;
 myhplus.init();
 cout<<"Initialize raw xsec..."<<endl;
 myhplus.init_xsec();
 cout<<"Initialize models..."<<endl;
 myhplus.init_model();
 cout<<endl;
 myhplus.run();
 
 cout<<endl;
 cout<<"Good bye. Live long and prosper!"<<endl;
 cout<<endl;
 
}

