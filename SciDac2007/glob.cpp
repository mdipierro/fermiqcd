#include "glob.h"
#include "string"
#include "iostream"
#include "vector"
using namespace std;

vector<string> glob(string pattern) {
  vector<string> v;
  glob_t pglob;
  pglob.gl_offs=2;
  if(glob(pattern.c_str(),0,0,&pglob)!=0) v.push_back("?");
  else
    for(int i=0; i<pglob.gl_pathc; i++)
      v.push_back(string(pglob.gl_pathv[i]));
  globfree(&pglob);
  return v;
}

int main() { 
  vector<string> v=glob("*.vtk");
  for(int i=0; i<v.size(); i++) cout << v[i] << endl;
  return 0; 
}

