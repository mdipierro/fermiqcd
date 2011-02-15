/////////////////////////////////////////////////////////////////
/// @file fermiqcd_staggered_mesons.h
/// @version 2009-12-21
/// @author Massimo Di Pierro <mdipierro@cs.depaul.edu>
///
/// Convenience functions to make staggered mesons
///
/// Distributed under GPL2 License
/// 
/// Created with support from the US Department of Energy
//////////////////////////////////////////////////////////////////

static const int KS_NDIM=4;

#ifndef TWISTED_BOUNDARY

class phase_field : public mdp_field<int> {
public:
  phase_field(mdp_lattice &a) {
    if(a.ndim!=4) error("fermiqcd_staggered_mesons/phase_field: ndim!=4");
    allocate_field(a,16);
  };
  int component(site x, site y) { 
    int i=0, mu;
    for(mu=0; mu<KS_NDIM; mu++) {
      if(x(mu)/2 != y(mu)/2) 
	return 0;
      i = (i << 1) + (x(mu) % 2);
    }
    return *address(x,i);
  }
  void compute(mdp_matrix GAMMA, mdp_matrix ZETA) {
    int A[KS_NDIM], B[KS_NDIM];
    mdp_matrix G1, G2, ZETA_H;
    site x(lattice());
    int i,mu,nu;
    ZETA_H=hermitian(ZETA);
    forallsites(x) {
      G1=mdp_identity(KS_NDIM);
      for(mu=1; mu<=KS_NDIM; mu++) {
	nu=mu % KS_NDIM;
	A[nu]=x(nu) % 2;
	if(A[nu]!=0) G1=-1.0*Gamma[nu]*G1;
      }
      for(i=0; i<16; i++) {
	B[0]=(i >> 3) & 0x1;
	B[1]=(i >> 2) & 0x1;
	B[2]=(i >> 1) & 0x1;
	B[3]=(i >> 0) & 0x1;
	G2=mdp_identity(KS_NDIM);
	for(mu=1; mu<=KS_NDIM; mu++) {
	  nu=mu % KS_NDIM;
	  if(B[nu]!=0) 
	    G2=G2*Gamma[nu];
	}
	*address(x,i)=(int) (0.25*real(trace(G1*GAMMA*G2*ZETA_H)));
      }
    }
  }
};


void operator_staggered_meson(staggered_field &out,
		    staggered_field &in,
		    phase_field &phases,
		    gauge_field &U) {

  site x(U.lattice());
  site x1(U.lattice());
  site x2(U.lattice());
  site y(U.lattice());
  int i, c, mu, A[KS_NDIM], B[KS_NDIM], P[4][2];
  int d, i0, i1, i2, i3;
  mdp_matrix paths, path;
  path=mdp_identity(U.nc);
  paths.dimension(U.nc,U.nc);
  paths=0;

  forallsites(x) {
    for(c=0; c<U.nc; c++) out(x,c)=0;
    for(mu=0; mu<KS_NDIM; mu++) 
      A[mu]=x(mu) % 2;
    for(i=0; i<16; i++) if(phases(x,i)!=0) {
      B[0]=(i >> 3) & 0x1;
      B[1]=(i >> 2) & 0x1;
      B[2]=(i >> 1) & 0x1;
      B[3]=(i >> 0) & 0x1;
      y.set(x(0)-A[0]+B[0],x(1)-A[1]+B[1], x(2)-A[2]+B[2],x(3)-A[3]+B[3]);
      // only for goldstone pion if(x!=y || phases(x,i)!=1) error("what!");
      if(phases(x,i)!=0) {
	for(d=0, mu=0; mu<KS_NDIM; mu++) 
	  if(B[mu]-A[mu]!=0) {
	    P[d][0]=B[mu]-A[mu];
	    P[d][1]=mu;
	    d++;
	  }
	if(d==0) paths=mdp_identity(U.nc);
	if(d==1) paths=U(x,P[0][0],P[0][1]);
	if(d==2) {
	  paths=
	    U(x,P[0][0],P[0][1])*U(x.hop(P[0][0],P[0][1]),P[1][0],P[1][1])+
	    U(x,P[1][0],P[1][1])*U(x.hop(P[1][0],P[1][1]),P[0][0],P[0][1]);
	  paths/=2;
	}
	if(d==3) {
	  x1=x.hop(P[0][0],P[0][1]);
	  paths=U(x,P[0][0],P[0][1])*U(x1,P[1][0],P[1][1])*U(x1.hop(P[1][0],P[1][1]),P[2][0],P[2][1]);
	  x1=x.hop(P[0][0],P[0][1]);
	  paths+=U(x,P[0][0],P[0][1])*U(x1,P[2][0],P[2][1])*U(x1.hop(P[2][0],P[2][1]),P[1][0],P[1][1]);
	  x1=x.hop(P[1][0],P[1][1]);
	  paths+=U(x,P[1][0],P[1][1])*U(x1,P[0][0],P[0][1])*U(x1.hop(P[0][0],P[0][1]),P[2][0],P[2][1]);
	  x1=x.hop(P[1][0],P[1][1]);
	  paths+=U(x,P[1][0],P[1][1])*U(x1,P[2][0],P[2][1])*U(x1.hop(P[2][0],P[2][1]),P[0][0],P[0][1]);
	  x1=x.hop(P[2][0],P[2][1]);
	  paths+=U(x,P[2][0],P[2][1])*U(x1,P[0][0],P[0][1])*U(x1.hop(P[0][0],P[0][1]),P[1][0],P[1][1]);
	  x1==x.hop(P[2][0],P[2][1]);
	  paths+=U(x,P[2][0],P[2][1])*U(x1,P[1][0],P[1][1])*U(x1.hop(P[1][0],P[1][1]),P[0][0],P[0][1]);
	  paths/=(3*2);
	}
	if(d==4) {
	  paths=0;
	  for(i0=0; i0<d; i0++)
	    for(i1=0; i1<d; i1++) if(i1!=i0)
	      for(i2=0; i2<d; i2++) if(i2!=i0 && i2!=i1)
		for(i3=0; i3<d; i3++) if(i3!=i0 && i3!=i1 && i3!=i2) {
		  x1=x.hop(P[i0][0],P[i0][1]);
		  x2=x1.hop(P[i1][0],P[i1][1]);
		  path=
		    U(x,P[i0][0],P[i0][1])*
		    U(x1,P[i1][0],P[i1][1])*
		    U(x2,P[i2][0],P[i2][1])*
		    U(x2.hop(P[i2][0],P[i2][1]),P[i3][0],P[i3][1]);
		  paths+=path;
		}
	  paths/=(4*3*2);
	}
	out(x)+=phases(x,i)*paths*in(y);
      }
    }
  }
}

const int local_source=1;  // 2^0
const int wall_source =2;  // 2^1

mdp_matrix make_meson(gauge_field &U, gauge_field &V,
		  mdp_matrix GAMMA,
		  mdp_matrix ZETA, 
		  coefficients &coeff1, 
		  coefficients &coeff2,
		  int source1_type=wall_source,
		  int source2_type=wall_source & local_source,
		  mdp_real precision =1e-7) {

  if(source1_type!=wall_source && source1_type!=local_source)
    error("fermiqcd_staggered_mesons/make_meson: source option not implemented");
  if(source2_type!=local_source)
    error("fermiqcd_staggered_mesons/make_meson: source option not implemented");
  
  int t, t_source, nsources=1, sourcestep=1;
  int i, j, nt=U.lattice().size(0);
  int nc=U.nc;
  mdp_complex c;
  site x(U.lattice());
  staggered_field tmp(U.lattice(),nc);
  staggered_field quark_source(U.lattice(),nc);
  staggered_field quark_prop(U.lattice(),nc);
  staggered_field anti_prop(U.lattice(),nc);

  phase_field phases(U.lattice());
  phases.compute(GAMMA, ZETA);

  mdp_matrix prop(1,nt);
  prop=0;

  U.update();
  V.update();

  for(t_source=0; t_source<nsources*sourcestep; t_source+=sourcestep) {
    for(i=0; i<U.nc; i++) {
      forallsites(x)
	for(j=0; j<U.nc; j++)
	  if(i==j && 
	     (x(0)==t_source || x(0)==(t_source+1)) && 
	     source1_type==wall_source) 
	    quark_source(x,j)=mdp_complex(1,0);
	  else if(i==j && 
		  (x(0)==t_source || x(0)==(t_source+1)) && 
		  (x(1)/2==0) && (x(2)/2==0) && (x(3)/2==0) &&
		  source1_type==local_source) 
	    quark_source(x,j)=mdp_complex(1,0);
	  else     
	    quark_source(x,j)=mdp_complex(0,0);

      quark_source.update();
      mul_invQ(quark_prop,quark_source, V,coeff1,precision);
      quark_prop.update();
      operator_staggered_meson(tmp, quark_source, phases, U);
      tmp.update();
      mul_invQ(anti_prop, tmp, V,coeff2,precision);
      quark_prop.update();      
      operator_staggered_meson(tmp, quark_prop, phases, U);
      forallsites(x) { 
	for(j=0, c=mdp_complex(0,0); j<U.nc; j++)
	  c+=conj(anti_prop(x,j))*tmp(x,j);
	t=(x(0)-t_source+nt) % nt; if(t%2==1) t--;
	prop(0,t)+=c;
      }
    }
  }
  return prop;
}


mdp_matrix GoldstonBoson_5x5(gauge_field &U, // input gauge field
			 gauge_field &V, // input fat and link links
			 coefficients &coeff,    // input quark mass
			 float precision=1e-6)        // color source index
{
  // // Local variable /////////////  
  int i,j;                        //
  unsigned int t_source=0;        // location of the source (timeslice)
  unsigned int t, nt              // number of timeslices
    =U.lattice().size(0);         //
  mdp_matrix tmp(nt,U.nc);            // auxiliary var
  mdp_matrix prop(2, nt);             // output vector
  site x(U.lattice());            // 
  // ///////////////////////////////

  // // Local fields ////////////////////////////////
  staggered_field quark_source(U.lattice(), U.nc); //
  staggered_field quark_prop(U.lattice(), U.nc);   //
  // ////////////////////////////////////////////////

  // // Initialization of parameters ////////////////
  prop=0;                                          //
  // ////////////////////////////////////////////////

  printf("Starting to make 5x5 pion...\n");
  for(i=0; i<U.nc; i++) {
    // // Make wall source /////////////////////////////
    forallsites(x)                                    //
      for(j=0; j<U.nc; j++)                           //
	if(i==j && x(0)-t_source==0)                  //
	  quark_source(x,j)=mdp_complex(1,0);             //
	else                                          //
	  quark_source(x,j)=mdp_complex(0,0);             //
    quark_source.update();                            //
    ////////////////////////////////////////////////////
    
    // // Make propagator - antipropagator - sink //////
    mul_invQ(quark_prop,quark_source,V,coeff,1e-5);   //
    // /////////////////////////////////////////////////
    
    // // trace with local sink ////////////////////////
    tmp=0;                                            //
    forallsites(x) {                                  //
      t=(x(0)-t_source+nt) % nt;                      //
      for(j=0; j<U.nc; j++)                           //
	tmp(t,j)+=conj(quark_prop(x,j))*quark_prop(x,j);
    }                                                 //
    mpi.add(tmp);                                     //
    for(t=0; t<nt; t++)                               //
      for(j=0; j<U.nc; j++)                           //
	prop(0,t)+=tmp(t,j);                          //
    // /////////////////////////////////////////////////
    
    // // trace with wall sink /////////////////////////
    tmp=0;                                            //
    forallsites(x) {                                  //
      t=(x(0)-t_source+nt) % nt;                      //
      for(j=0; j<U.nc; j++)                           //
	tmp(t,j)+=quark_prop(x,j);                    //
    }                                                 //
    mpi.add(tmp);                                     //
    for(t=0; t<nt; t++)                               //
      for(j=0; j<U.nc; j++)                           //
	prop(1,t)+=conj(tmp(t,j))*tmp(t,j);           //
    // /////////////////////////////////////////////////
  }
  // // Normalize output and retrun //////////////////
  int volume=U.lattice().size(1)*
    U.lattice().size(2)*
    U.lattice().size(3);
  for(t=0; t<nt; t++) {
    prop(0,t)/=U.nc*volume;
    prop(1,t)/=U.nc*volume*volume;
  }

  return prop;
}

#endif

/* EXAMPLE:

// #define BLOCKSITE 4
// #define TWISTED_BOUNDARY

#include "/home/mdp/mdpqcd2/MDP_Lib2.h"
#include "/home/mdp/mdpqcd2/MDP_MPI.h"
#include "/home/mdp/mdpqcd2/MDP_Prompt.h"
#include "/home/mdp/mdpqcd2/MDP_Gauge.h"
#include "/home/mdp/mdpqcd2/MDP_iGauge.h"
#include "/home/mdp/mdpqcd2/MDP_GaugeFix.h"
#include "/home/mdp/mdpqcd2/MDP_Fermi.h"
#include "/home/mdp/mdpqcd2/MDP_Staggered.h"
#include "/home/mdp/mdpqcd2/MDP_StaggeredMesons.h"

int main(int argc, char **argv) {
  mpi.open_wormholes(argc, argv);
  define_base_matrices("FERMILAB");

 
  //  define_twist_matrices();

  printf("COMPUTING C2(t) FOR STAGGERED MESON: %s\n", argv[1]);
  
  float seed=0;
  int t, Nt=24, Ns=8;
  int Nc=3;
  int box[]={Nt, Ns, Ns, Ns};

  float pU, u0=0.8629;
  float *c;
  char s[256];

  mdp_lattice space_time(4,box,
			     default_partitioning<0>,
			     torus_topology,
			     seed,3);

  int conf, nconf=(int) val(argv[2]);
  
  gauge_field     U(space_time,Nc);
  gauge_field     V(space_time,Nc);
  eta_field       eta(space_time);
  site x(space_time);
  mdp_matrix prop;
  float mass_a, mass_b;

  mass_a=0.03;
  mass_b=0.03;

  default_staggered_action=staggered_mul_Q_TWO;
  default_staggered_inversion_method=staggered_BiCG_UML;

  u0=0.8629;

  JackBoot c2(nconf, Nt);

  for(conf=0; conf<nconf; conf++) {
    c2.conf=conf;
    sprintf(s, "/home/mdp/data/gauge_improved_b7.4_u0.8629/gauge32x08_b7.4_u0.8629_n%.6i", conf+1);
 
    U.load(s);
    sprintf(s, "/home/mdp/data/gauge_improved_b7.4_u0.8629/gauge24x08_b7.4_u0.8629_n%.6i", conf+1);
    U.save(s);

    pU=pow(u0,4);
    c=lepage_coeffiecients(pU, "Full");
    set_antiperiodic_phases(U,0,true);
    lepage_improved_links(V,U,c);

    if(strcmp(argv[1],"5x5")==0) 
      prop= make_meson(U,V, Gamma5, Gamma5,mass_a, mass_b, 
		       wall_source, wall_source, 1e-5);
    if(strcmp(argv[1],"05x05")==0) 
      prop= make_meson(U,V, Gamma[0]*Gamma5, Gamma[0]*Gamma5,mass_a, mass_b, 
		       wall_source, wall_source, 1e-5);
    if(strcmp(argv[1],"5x35")==0) 
      prop= make_meson(U,V, Gamma5, Gamma[3]*Gamma5,mass_a, mass_b, 
		       wall_source, wall_source, 1e-5);
    if(strcmp(argv[1],"05x12")==0) 
      prop= make_meson(U,V, Gamma[0]*Gamma5, Gamma[1]*Gamma[2],mass_a, mass_b, 
		       wall_source, wall_source, 1e-5);
    if(strcmp(argv[1],"5x03")==0) 
      prop= make_meson(U,V, Gamma5, Gamma[0]*Gamma[3],mass_a, mass_b, wall_source, 
		       wall_source, 1e-5);
    if(strcmp(argv[1],"05x3")==0) 
      prop= make_meson(U,V, Gamma[0]*Gamma5, Gamma[3],mass_a, mass_b, 
		       wall_source, wall_source, 1e-5);
    if(strcmp(argv[1],"5x0")==0) 
      prop= make_meson(U,V, Gamma5, Gamma[0],mass_a, mass_b, 
		       wall_source, wall_source, 1e-5);
    if(strcmp(argv[1],"05xI")==0) 
      prop= make_meson(U,V, Gamma[0]*Gamma5, Gamma1,mass_a, mass_b, 
		       wall_source, wall_source, 1e-5);

    for(t=0; t<Nt; t++) c2(t)=real(prop(0,t));

    float x, dx;

    if(ME==0) {
      for(t=0; t<Nt; t=t+2) {
	c2.plain(t);
	x=c2.mean();
	dx=c2.j_err();
	printf("%i,\t%e,\t%e\n", t, log(x), dx/x);
      }
      printf("\n");
    }
  }
  mpi.close_wormholes();
  return 0;
}




NOTES:
make meson dow not work is you enable twist. The rest does work!


*/


