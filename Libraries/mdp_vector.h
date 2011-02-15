/////////////////////////////////////////////////////////////////
/// @file mdp_vector.h
/// @version 2009-12-21
/// @author Massimo Di Pierro <mdipierro@cs.depaul.edu>
///
/// Contains class mdp_vector
///
/// Licensed under GPL2 license
/// Read attached license in file mdp_license.pdf
/// This file cannot be distributed without file mdp_license.pdf
//////////////////////////////////////////////////////////////////

/// @brief discerete vectors to navigate on a lattice
/// 
class mdp_vector {
public:
  int x[10];
  mdp_vector() {
    x[0]=x[1]=x[2]=x[3]=x[4]=0;
    x[5]=x[6]=x[7]=x[8]=x[9]=0;
  }
  mdp_vector(int x0,   int x1=0, int x2=0, int x3=0, int x4=0,
	 int x5=0, int x6=0, int x7=0, int x8=0, int x9=0) {
    x[0]=x0; x[1]=x1; x[2]=x2; x[3]=x3; x[4]=x4;
    x[5]=x5; x[6]=x6; x[7]=x7; x[8]=x8; x[9]=x9;
  }
};

inline mdp_vector binary2versor(mdp_int a) {
  mdp_vector v((a)      & 0x1, 
	   (a >> 1) & 0x1, 
	   (a >> 2) & 0x1, 
	   (a >> 3) & 0x1, 
	   (a >> 4) & 0x1, 
	   (a >> 5) & 0x1, 
	   (a >> 6) & 0x1, 
	   (a >> 7) & 0x1, 
	   (a >> 8) & 0x1, 
	   (a >> 9) & 0x1);
  return v;
}

inline int versor2binary(int x0,   int x1=0, int x2=0, int x3=0, int x4=0,
			 int x5=0, int x6=0, int x7=0, int x8=0, int x9=0) {
#ifdef CHECK_ALL
  if((fabs(0.5-x0)>1) || (fabs(0.5-x1)>1) || 
     (fabs(0.5-x2)>1) || (fabs(0.5-x3)>1) ||
     (fabs(0.5-x4)>1) || (fabs(0.5-x5)>1) || 
     (fabs(0.5-x6)>1) || (fabs(0.5-x7)>1) ||
     (fabs(0.5-x8)>1) || (fabs(0.5-x9)>1)) error("versor2binary"); 				   
#endif
  return x0+2*x1+4*x2+8*x3+16*x4+32*x5+64*x6+128*x7+256*x8+512*x9;
}

inline mdp_int vector2binary(mdp_vector v) {
#ifdef CHECK_ALL
  if((fabs(0.5-v.x[0])>1) || (fabs(0.5-v.x[1])>1) || 
     (fabs(0.5-v.x[2])>1) || (fabs(0.5-v.x[3])>1) ||
     (fabs(0.5-v.x[4])>1) || (fabs(0.5-v.x[5])>1) || 
     (fabs(0.5-v.x[6])>1) || (fabs(0.5-v.x[7])>1) ||
     (fabs(0.5-v.x[8])>1) || (fabs(0.5-v.x[9])>1) || 
     (fabs(0.5-v.x[2])>1)) error("vector2binary");
#endif
  return v.x[0]+2*v.x[1]+4*v.x[2]+8*v.x[3]+16*v.x[4]+
    32*v.x[5]+64*v.x[6]+128*v.x[7]+256*v.x[8]+512*v.x[9];
}
