#include "mdp_all.h"
#include "string"
#include "map"
#include "vector"
using namespace std;

class Step {
public:
  string algorithm;
  map<string,float> parameters;
};


void parse(int argc, char** argv) {
  vector<Step> steps;
  Step step;
  string s;
  int i,j,k;

  for(i=1; i<argc; i++) {
    step
    s=string(argv[i])
      if(argv[i][0]=='-') {
	j=s.find(":");
	if(j<0) 
      }
      
    cout << argv[i] << endl;
}

int main(int argc, char** argv) {
  parse(argc,argv);
  return 0;
}
