#include <vector>

#include "model.h"
#include "limit.h"

using namespace std;

class hplus
{
 
 public:
 
 hplus() {};
 ~hplus() {};

 vector<model> mymodel;
 
 vector<string> model_list;
 vector<string> BR_list;
 string energy;
 
 void run();
 void init();
 void init_xsec();
 void init_model();
 
 double get_xsec(double,double);
 double get_width(double,double);
 void simple_output();
 //double round(double,int);
 double interpolate_quad(double,double,double,double,double,double,double,double,double,double,double,double,double);
 double get_xsec_raw_exact(double,double);
 
 double get_xsec_mssm(double,double,int);
 
 vector<double> xsec_input_mass;
 vector<double> xsec_input_tanb;
 
 vector<double> xsec_val;
 vector<double> xsec_up;
 vector<double> xsec_dn;
 vector<double> xsec_mH;
 vector<double> xsec_tanb;
 
 vector<double> width;
 vector<double> BR_taunu;
 vector<double> BR_tb;
 
 vector<double> BR_topHplus;
 
 void a_vs_b();
 void vs_mass(string); 
 void vs_tanb(string);
 
 int exclusions();
 bool fexists(string);
 
 vector<limit> mylimit;
 
 void get_crossings(int , int , double , double , double&, double&, string);
 double linear(double,double,double,double,double);
 void draw_exclusions(string,int,int,double,double,double,int);
 
 int m_debug;
 
};

