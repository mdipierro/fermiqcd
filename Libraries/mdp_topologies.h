/////////////////////////////////////////////////////////////////
/// @file mdp_topologies.h
/// @version 2009-12-21
/// @author Massimo Di Pierro <mdipierro@cs.depaul.edu>
///
/// Examples of lattice topologies
///
/// Licensed under GPL2 license
/// Read attached license in file mdp_license.pdf
/// This file cannot be distributed without file mdp_license.pdf
//////////////////////////////////////////////////////////////////

// //////////////////////////////////////////////////////
// Basic topologies:
// //////////////////////////////////////////////////////

void torus_topology(int mu, 
		    int *x_dw, 
		    int *x, 
		    int *x_up, 
		    int ndim, 
		    int *nx) {
  for(int nu=0; nu<ndim; nu++) if(nu==mu) {
    x_dw[mu]=(x[mu]-1+nx[mu]) % nx[mu];
    x_up[mu]=(x[mu]+1) % nx[mu];
  } else x_up[nu]=x_dw[nu]=x[nu];
}

void box_topology(int mu, 
		  int *x_dw, 
		  int *x, 
		  int *x_up, 
		  int ndim, 
		  int *nx) { 
  torus_topology(mu,x_dw,x,x_up,ndim,nx);
  if(x[mu]==0)        x_dw[mu]=x[mu];
  if(x[mu]==nx[mu]-1) x_up[mu]=x[mu];
} 

void moebious_topolgy(int mu, 
		      int *x_dw, 
		      int *x, 
		      int *x_up, 
		      int ndim, 
		      int *nx) {
  torus_topology(mu,x_dw,x,x_up,ndim,nx);
  if(mu==0) {
    if(x[0]==0)       x_dw[1]=nx[1]-x[1]-1;
    if(x[0]==nx[0]-1) x_up[1]=nx[1]-x[1]-1;
  }
}
