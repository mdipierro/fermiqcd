/////////////////////////////////////////////////////////////////
/// @file mdp_utils.h
/// @version 2009-12-21
/// @author Massimo Di Pierro <mdipierro@cs.depaul.edu>
///
/// Other junk that did not fit anywhere else
///
/// Licensed under GPL2 license
/// Read attached license in file mdp_license.pdf
/// This file cannot be distributed without file mdp_license.pdf
//////////////////////////////////////////////////////////////////

string tostring(int k, int length=5) {
  char buf[128];
  if(length>=5)
    sprintf(buf,"%.5i",k);
  else if(length==4)
    sprintf(buf,"%.4i",k);
  else if(length==3)
    sprintf(buf,"%.3i",k);
  else if(length==2)
    sprintf(buf,"%.2i",k);
  else if(length==1)
    sprintf(buf,"%.1i",k);
  return string(buf);
}

bool file_exists(string filename) {
  if (FILE * file = fopen(filename.c_str(), "rb")) {
    fclose(file);
    return true;
  }
  return false;
}

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

string latest_file(string pattern) {
  vector<string> v=glob(pattern);
  return v[v.size()-1];
}

string next_to_latest_file(string pattern) {
  int i=pattern.find("*");
  if(i<0) return pattern;
  vector<string> v=glob(pattern);
  string latest=v[v.size()-1];
  if(latest=="?") return pattern.replace(i,1,tostring(0));
  int k;
  latest=latest.substr(i,5);
  k=atoi(latest.c_str());
  k++;
  return pattern.replace(i,1,tostring(k).c_str());
}


string tostring(float k) {
  char buf[128];
  sprintf(buf,"%f",k);
  return string(buf);
}

int is_file(string filename, char permission[]="r") {
  FILE *fp=fopen(filename.c_str(), permission);
  if(fp>0) {
    fclose(fp);
    return true;
  }
  return false;
}

mdp_field_file_header get_info(string filename, int proc=0) {
  mdp_field_file_header myheader;
  if(ME==proc) {
    FILE *fp=fopen(filename.c_str(), "r");
    if(fp==0) error("Unable to open file");
    fread(&myheader, sizeof(char),
	  sizeof(mdp_field_file_header)/sizeof(char), fp);
    switch_header_endianess(myheader);
    fclose(fp); // fixed by Lucky [lucky@sfu.ca]
  }
  mpi.broadcast(myheader,proc);
  return myheader;
}

int mail(string email, string message) {
  string s;
  static int ret;
  if(ME==0) {
    s = "echo '"+message+"' | mail -s 'MDP MESSAGE' " +email;
    ret=system(s.c_str());
  }
  mpi.broadcast(ret,0);
  return ret;
}

int mail_file(string email, string filename) {
  string s;
  static int ret;
  if(ME==0) {
    s = "more "+filename+" | mail -s 'MDP MESSAGE' "+email;
    ret=system(s.c_str());
  }
  mpi.broadcast(ret,0);
  return ret;
}

bool startswith(string a, string b) {
  return (a.substr(0,b.length())==b);
}

bool endswith(string a, string b) {
  int i = a.length();
  int j = b.length();
  if(i<j) return false;
  return (a.substr(i-j,j)==b);
}

int parse_int(string a, string b, int value=0) {
  int i=a.find(string(":")+b),j=0;
  if(i<0) return value;
  else {
    i += b.length()+2;
    j = a.find(":",i);
    if(j<0) j=a.length();
    sscanf(a.substr(i,j-i).c_str(),"%i",&value);
    return value;
  }
}

float parse_float(string a, string b, float value=0.0) {
  int i=a.find(string(":")+b),j=0;
  if(i<0) return value;
  else {
    i += b.length()+2;
    j = a.find(":",i);
    if(j<0) j=a.length();
    sscanf(a.substr(i,j-i).c_str(),"%f",&value);
    return value;
  }
}

string parse_string(string a, string b, string value="") {
  int i=a.find(string(":")+b),j=0;
  char cvalue[512];
  if(i<0) return value;
  else {
    i += b.length()+2;
    j = a.find(":",i);
    if(j<0) j=a.length();
    return a.substr(i,j-i);
  }
}

class mdp_args {
 public:
  vector<string> args;
  mdp_args(int argc, char** argv) {
    for(int i=1; i<argc; i++)
      this->args.push_back(argv[i]);
  }
  int length() {
    return args.size();
  }
  bool have(string name) {
    for(int i=0; i<this->args.size(); i++)
      if(this->args[i]==name || startswith(this->args[i],name+":"))
	return true;
    return false;
  }
  float get(string name, string key, float value=0.0) {
    for(int i=0; i<this->args.size(); i++)
      if(startswith(this->args[i],name+":"))
	return parse_float(this->args[i],key,value);
    return value;
  }
  float get(string name, string key, double value=0.0) {
    for(int i=0; i<this->args.size(); i++)
      if(startswith(this->args[i],name+":")) {
	value = parse_float(this->args[i],key,value);
	break;
      }
    mdp << "INPUT " << name << ":" << key << "=" << value << endl;
    return value;
  }
  int get(string name, string key, int value=0) {
    for(int i=0; i<this->args.size(); i++)
      if(startswith(this->args[i],name+":")) {
	value = parse_int(this->args[i],key,value);
	break;
      }
    mdp << "INPUT " << name << ":" << key << "=" << value << endl;
    return value;
  }
  string get(string name, string key, string value="") {
    int i = value.find('|') ;
    if (i>=0)
      value = value.substr(0,i);
    for(int i=0; i<this->args.size(); i++)
      if(startswith(this->args[i],name+":")) {
	value = parse_string(this->args[i],key,value);
	break;
      }
    mdp << "INPUT " << name << ":" << key << "=" << value << endl;
    return value;
  }
};  

