#include <string>
#include <vector>

using namespace std;

class BR
{

 private:
 
 string m_name;

 public:
 	
 BR() {};
 ~BR() {};
    
 vector<double> val;

 void set_name(string);
 string get_name();
 void add_val(double);
 int get_size();
 double get_val(int);
 
 
};

