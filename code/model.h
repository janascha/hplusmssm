#include <string>
#include <vector>

#include "BR.h"

using namespace std;

class model
{
	
 private:
 string m_name;
 
 public:
 	
 model(string,vector<string>);
 ~model() {};
 
 vector<string> m_BR_list;
 
 void init_BR();
 void init_deltab();
 void init_topBR();
 
 string get_name();
 void set_name(string);
 int get_min(vector<double> vector);
 
 BR BR_mHp;
 BR BR_tanb;
 BR BR_width;
 vector<BR> allBR;
 
 int m_dodeltab;
 vector<double> deltab_val;
 vector<double> deltab_tanb;
 
 vector<double> topBR_val;
 vector<double> topBR_mHp;
 vector<double> topBR_tanb;
 
 double get_BR(string,double,double);
 double get_BRval(int,double,double);
 double get_deltab(double tanb);
 
 double get_topBR(double mass, double tanb);
 double get_topBRval(double mass, double tanb);
 
 double myround(double,int);
 double myround_up(double,int);
 double myround_dn(double,int);
 double linear(double,double,double,double,double);
 
 int    get_maxtanb();
 
};

