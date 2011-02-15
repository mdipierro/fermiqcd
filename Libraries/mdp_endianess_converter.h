/////////////////////////////////////////////////////////////////
/// @file mdp_endianess_converter.h
/// @version 2009-12-21
/// @author Massimo Di Pierro <mdipierro@cs.depaul.edu>
///
/// Contains declaration of function swicth_endianess_byte4()
///
/// Licensed under GPL2 license
/// Read attached license in file mdp_license.pdf
/// This file cannot be distributed without file mdp_license.pdf
//////////////////////////////////////////////////////////////////

/// Converts endianess of object passed by reference
template<class T>
void switch_endianess_byte4(T &a) {
  char *p=(char*) &a;
  static char q[4];
  if(sizeof(T)==4) {
    q[0]=p[0];
    q[1]=p[1];
    q[2]=p[2];
    q[3]=p[3];
    p[0]=q[3];
    p[1]=q[2];
    p[2]=q[1];
    p[3]=q[0];  
  } else error("switch_endianess_byte4: sizeof(T)!=4");
}

template<class T>
void switch_endianess_byte8(T &a) {
  char *p=(char*) &a;
  static char q[8];
  if(sizeof(T)==8) {
    q[0]=p[0];
    q[1]=p[1];
    q[2]=p[2];
    q[3]=p[3];
    q[4]=p[4];
    q[5]=p[5];
    q[6]=p[6];
    q[7]=p[7];
    p[0]=q[7];
    p[1]=q[6];
    p[2]=q[5];
    p[3]=q[4];
    p[4]=q[3];
    p[5]=q[2];
    p[6]=q[1];
    p[7]=q[0];

  } else error("switch_endianess_byte8 sizeof(T)!=8");
}

