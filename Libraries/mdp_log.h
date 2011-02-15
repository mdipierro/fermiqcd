/////////////////////////////////////////////////////////////////
/// @file mdp_log.h
/// @version 2009-12-21
/// @author Massimo Di Pierro <mdipierro@cs.depaul.edu>
///
/// Contains class mdp_log
///
/// Licensed under GPL2 license
/// Read attached license in file mdp_license.pdf
/// This file cannot be distributed without file mdp_license.pdf
//////////////////////////////////////////////////////////////////

/// @brief base class of class mdp_communicator (DO NOT INSTANTIATE)
/// @see class mdp_communicator
class mdp_log {
private:
  int level;
  int max_level;
  vector<string> level_tag;
  ostream* os;
public:
  bool print;
  void abort() {
    exit(-1);
  }
  void set_level(int i) {
    max_level=i;
  }
  mdp_log() {   
    level=0;
    max_level=100000;
    print=true;
    connect(cout);
  }
  void connect(ostream &os1) {
    os=&os1;
  }
  void connect(ofstream &os2) {
    os=&os2; // is this correct? I think so!
  }
  void error_message(string s, string file, int line) {
    if(print) {
      begin_function("error");
      *os << "In file \"" << file;
      *os << "\", before line " << line;
      *os << ", this error occurred: " << s << '\n';
      for(;level; level--) 
	if(level<max_level) 
	  *os << "</" << level_tag[level-1] << ">" << '\n';
    }
    throw s;
  }
  void begin_function(string s) {
    level_tag.resize(++level);
    level_tag[level-1]=s;  
    if(print && level<max_level) 
      *os << "<" << s << ">" << '\n';
  }
  void end_function(string s) {
    if(level_tag[level-1]==s) {
      if(print && level<max_level) 
	*os << "</" << level_tag[level-1] << ">" << '\n';
      level--;
    } else error_message("missing end_function()", "unkown", 0);
  }
  template<class T>
    mdp_log &operator<< (const T x) {
    if (print && level<max_level) *os << x;
    return (*this);
  }
};

