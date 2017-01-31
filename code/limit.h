#include <string>
#include <vector>
#include "TGraph.h"

using namespace std;

class limit
{

 private:
 
 string m_name;

 public:
 	
 limit() { };
 ~limit() {};
 
 void set_name(string);
 string get_name();
 
 void add_val(double);
 double get_val(int);
 vector<double> val;

 vector<double> mass;
 void add_mass(double);
 double get_mass(int);

 void add_up(double);
 vector<double> crossings_up;
 double get_up(int);
 
 void add_dn(double);
 vector<double> crossings_dn;
 double get_dn(int);
 
 TGraph* graph_up;
 TGraph* graph_dn;
 
 void init_names(); 
 
};

