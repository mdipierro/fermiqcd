/////////////////////////////////////////////////////////////////
/// @file fermiqcd_fermilab_action.h
/// @version 2009-12-21
/// @author Massimo Di Pierro <mdipierro@cs.depaul.edu>
///
/// Fermilab action with all dimension 5 and 6 terms (requires SSE)
///
/// Distributed under GPL2 License
/// 
/// Created with support from the US Department of Energy
//////////////////////////////////////////////////////////////////

#if defined(SSE2)

class FermiFermilabActionNew {
 public:

  static void mul_Q(fermi_field &psi_out,
		    fermi_field &psi_in,
		    gauge_field &U,
		    coefficients &coeff, int parity=EVENODD) {
    
    mdp_lattice& lattice=psi_in.lattice();

    if(parity!=EVENODD)
      error("FermiFermilabAction::mul_Q\nparity must be EVENODD");

    if(psi_in.nspin!=4)
      error("FermiFermilabAction::mul_Q\ndoes not work for nspin!=4");
    if(psi_in.nc!=U.nc)
      error("FermiFermilabAction::mul_Q\nincompatible number of colors");

      mdp_real sign, kappa, kappat,kappas,zeta,rs,rt,cSW,cE,cB;
      mdp_real c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,mac1,mac2,mac3,mac4,mac5;

    if(coeff.has_key("sign")) sign=coeff["sign"];
    else sign=1;
    if(sign!=+1)
      error("FermiFermilabAction::mul_Q\nsign must be +1.0");

    if(coeff.has_key("kappa")) kappat=kappas=coeff["kappa"];
    if(coeff.has_key("zeta"))  kappat=kappas/coeff["zeta"];

    if(coeff.has_key("kappa_s")) kappas=coeff["kappa_s"];
    if(coeff.has_key("kappa_t")) kappat=coeff["kappa_t"];

    if(kappat==0 || kappas==0) 
      error("FermiFermilabAction::mul_Q\nparameter kappa not assigned");

    zeta=kappas/kappat;
   
    if(coeff.has_key("r_t")) rt=coeff["r_t"];
    else rt=1;

    if(coeff.has_key("r_s")) rs=coeff["r_s"];
    else rs=1;

    if(coeff.has_key("c_{sw}")) cSW=coeff["c_{sw}"];
    else cSW=0;

    if(coeff.has_key("c_E")) cE=coeff["c_E"];
    else cE=1; // this definition different from Clover

    if(coeff.has_key("c_B")) cB=coeff["c_B"];
    else cB=1; // this definition different from Clover

    if(coeff.has_key("alpha_1")) c1=coeff["alpha_1"];
    else c1=0;
    if(coeff.has_key("alpha_2")) c2=coeff["alpha_2"];
    else c2=0;
    if(coeff.has_key("alpha_3")) c3=coeff["alpha_3"];
    else c3=0;
    if(coeff.has_key("alpha_4")) c4=coeff["alpha_4"];
    else c4=0;
    if(coeff.has_key("alpha_5")) c5=coeff["alpha_5"];
    else c5=0;
    if(coeff.has_key("alpha_6")) c6=coeff["alpha_6"];
    else c6=0;
    if(coeff.has_key("alpha_7")) c7=coeff["alpha_7"];
    else c7=0;
    if(coeff.has_key("alpha_8")) c8=coeff["alpha_8"];
    else c8=0;
    if(coeff.has_key("alpha_9")) c9=coeff["alpha_9"];
    else c9=0;
    if(coeff.has_key("alpha_10")) c10=coeff["alpha_10"];
    else c10=0;
    if(coeff.has_key("mac_1")) mac1=coeff["mac_1"];
    else mac1=0;
    if(coeff.has_key("mac_2")) mac2=coeff["mac_2"];
    else mac2=0;
    if(coeff.has_key("mac_3")) mac3=coeff["mac_3"];
    else mac3=0;
    if(coeff.has_key("mac_4")) mac4=coeff["mac_4"];
    else mac4=0;
    if(coeff.has_key("mac_5")) mac5=coeff["mac_5"];
    else mac5=0;
    
    site x(lattice);

    int i,j,a,b,k;

    mdp_matrix_field Up(lattice,4,3);
    mdp_matrix_field Dw(lattice,4,3);
    mdp_matrix_field Ei(lattice,3,3);
    mdp_matrix_field Bi(lattice,3,3);

    mdp_matrix out(3,1);


static mdp_matrix u0u0(3,1);
static mdp_matrix d0d0(3,1);
static mdp_matrix u1u1(3,1);
static mdp_matrix d1d1(3,1);
static mdp_matrix e1u1(3,1);
static mdp_matrix e1d1(3,1);
static mdp_matrix u1e1(3,1);
static mdp_matrix d1e1(3,1);
static mdp_matrix b1u1(3,1);
static mdp_matrix b1d1(3,1);
static mdp_matrix u1b1(3,1);
static mdp_matrix d1b1(3,1);
static mdp_matrix u1u2(3,1);
static mdp_matrix u1d2(3,1);
static mdp_matrix d1u2(3,1);
static mdp_matrix d1d2(3,1);
static mdp_matrix e1u2(3,1);
static mdp_matrix e1d2(3,1);
static mdp_matrix u1e2(3,1);
static mdp_matrix d1e2(3,1);
static mdp_matrix b1u2(3,1);
static mdp_matrix b1d2(3,1);
static mdp_matrix u1b2(3,1);
static mdp_matrix d1b2(3,1);
static mdp_matrix u1u3(3,1);
static mdp_matrix u1d3(3,1);
static mdp_matrix d1u3(3,1);
static mdp_matrix d1d3(3,1);
static mdp_matrix e1u3(3,1);
static mdp_matrix e1d3(3,1);
static mdp_matrix u1e3(3,1);
static mdp_matrix d1e3(3,1);
static mdp_matrix b1u3(3,1);
static mdp_matrix b1d3(3,1);
static mdp_matrix u1b3(3,1);
static mdp_matrix d1b3(3,1);
static mdp_matrix u2u1(3,1);
static mdp_matrix u2d1(3,1);
static mdp_matrix d2u1(3,1);
static mdp_matrix d2d1(3,1);
static mdp_matrix e2u1(3,1);
static mdp_matrix e2d1(3,1);
static mdp_matrix u2e1(3,1);
static mdp_matrix d2e1(3,1);
static mdp_matrix b2u1(3,1);
static mdp_matrix b2d1(3,1);
static mdp_matrix u2b1(3,1);
static mdp_matrix d2b1(3,1);
static mdp_matrix u2u2(3,1);
static mdp_matrix d2d2(3,1);
static mdp_matrix e2u2(3,1);
static mdp_matrix e2d2(3,1);
static mdp_matrix u2e2(3,1);
static mdp_matrix d2e2(3,1);
static mdp_matrix b2u2(3,1);
static mdp_matrix b2d2(3,1);
static mdp_matrix u2b2(3,1);
static mdp_matrix d2b2(3,1);
static mdp_matrix u2u3(3,1);
static mdp_matrix u2d3(3,1);
static mdp_matrix d2u3(3,1);
static mdp_matrix d2d3(3,1);
static mdp_matrix e2u3(3,1);
static mdp_matrix e2d3(3,1);
static mdp_matrix u2e3(3,1);
static mdp_matrix d2e3(3,1);
static mdp_matrix b2u3(3,1);
static mdp_matrix b2d3(3,1);
static mdp_matrix u2b3(3,1);
static mdp_matrix d2b3(3,1);
static mdp_matrix u3u1(3,1);
static mdp_matrix u3d1(3,1);
static mdp_matrix d3u1(3,1);
static mdp_matrix d3d1(3,1);
static mdp_matrix e3u1(3,1);
static mdp_matrix e3d1(3,1);
static mdp_matrix u3e1(3,1);
static mdp_matrix d3e1(3,1);
static mdp_matrix b3u1(3,1);
static mdp_matrix b3d1(3,1);
static mdp_matrix u3b1(3,1);
static mdp_matrix d3b1(3,1);
static mdp_matrix u3u2(3,1);
static mdp_matrix u3d2(3,1);
static mdp_matrix d3u2(3,1);
static mdp_matrix d3d2(3,1);
static mdp_matrix e3u2(3,1);
static mdp_matrix e3d2(3,1);
static mdp_matrix u3e2(3,1);
static mdp_matrix d3e2(3,1);
static mdp_matrix b3u2(3,1);
static mdp_matrix b3d2(3,1);
static mdp_matrix u3b2(3,1);
static mdp_matrix d3b2(3,1);
static mdp_matrix u3u3(3,1);
static mdp_matrix d3d3(3,1);
static mdp_matrix e3u3(3,1);
static mdp_matrix e3d3(3,1);
static mdp_matrix u3e3(3,1);
static mdp_matrix d3e3(3,1);
static mdp_matrix b3u3(3,1);
static mdp_matrix b3d3(3,1);
static mdp_matrix u3b3(3,1);
static mdp_matrix d3b3(3,1);
mdp_complex c_id[16];     //fermilabread_id(c_id);
mdp_complex c_u[16*4];    //fermilabread_u(c_u);
mdp_complex c_d[16*4];    //fermilabread_d(c_d);
mdp_complex c_e[16*4];    //fermilabread_e(c_e);
mdp_complex c_b[16*4];    //fermilabread_b(c_b);
mdp_complex c_uu[16*4*4]; //fermilabread_uu(c_uu);
mdp_complex c_ud[16*4*4]; //fermilabread_ud(c_ud);
mdp_complex c_du[16*4*4]; //fermilabread_du(c_du);
mdp_complex c_dd[16*4*4]; //fermilabread_dd(c_dd);
mdp_complex c_eu[16*4*4]; //fermilabread_eu(c_eu);
mdp_complex c_ed[16*4*4]; //fermilabread_ed(c_ed);
mdp_complex c_ue[16*4*4]; //fermilabread_ue(c_ue);
mdp_complex c_de[16*4*4]; //fermilabread_de(c_de);
mdp_complex c_bu[16*4*4]; //fermilabread_bu(c_bu);
mdp_complex c_bd[16*4*4]; //fermilabread_bd(c_bd);
mdp_complex c_ub[16*4*4]; //fermilabread_ub(c_ub);
mdp_complex c_db[16*4*4]; //fermilabread_db(c_db);


    #include "fermiqcd_fermilab_coefficients.h"

    psi_out=0;

    
for(a=0; a<4; a++) {

forallsites(x) {

_sse_mulABC_set_331(&U(x,0,0,0),&psi_in(x+0,a,0),&Up(x,0,0));
_sse_mulABC_set_331(&U(x,1,0,0),&psi_in(x+1,a,0),&Up(x,1,0));
_sse_mulABC_set_331(&U(x,2,0,0),&psi_in(x+2,a,0),&Up(x,2,0));
_sse_mulABC_set_331(&U(x,3,0,0),&psi_in(x+3,a,0),&Up(x,3,0));
_sse_mulAHBC_set_331(&U(x-0,0,0,0),&psi_in(x-0,a,0),&Dw(x,0,0));
_sse_mulAHBC_set_331(&U(x-1,1,0,0),&psi_in(x-1,a,0),&Dw(x,1,0));
_sse_mulAHBC_set_331(&U(x-2,2,0,0),&psi_in(x-2,a,0),&Dw(x,2,0));
_sse_mulAHBC_set_331(&U(x-3,3,0,0),&psi_in(x-3,a,0),&Dw(x,3,0));
_sse_mulABC_set_331(&U.em(x,0,1,0,0),&psi_in(x,a,0),&Ei(x,0,0));
_sse_mulABC_set_331(&U.em(x,0,2,0,0),&psi_in(x,a,0),&Ei(x,1,0));
_sse_mulABC_set_331(&U.em(x,0,3,0,0),&psi_in(x,a,0),&Ei(x,2,0));
_sse_mulABC_set_331(&U.em(x,2,3,0,0),&psi_in(x,a,0),&Bi(x,0,0));
_sse_mulABC_set_331(&U.em(x,1,3,0,0),&psi_in(x,a,0),&Bi(x,1,0));
_sse_mulABC_set_331(&U.em(x,1,2,0,0),&psi_in(x,a,0),&Bi(x,2,0));

}

Up.update();
Dw.update();
Ei.update();
Bi.update();

forallsites(x) {

_sse_mulABC_set_331(&U(x,0,0,0),&Up(x+0,0,0),&u0u0(0,0));
_sse_mulAHBC_set_331(&U(x-0,0,0,0),&Dw(x-0,0,0),&d0d0(0,0));
_sse_mulABC_set_331(&U(x,1,0,0),&Up(x+1,1,0),&u1u1(0,0));
_sse_mulAHBC_set_331(&U(x-1,1,0,0),&Dw(x-1,1,0),&d1d1(0,0));
_sse_mulABC_set_331(&U(x,1,0,0),&Bi(x+1,0,0),&u1b1(0,0));
_sse_mulAHBC_set_331(&U(x-1,1,0,0),&Bi(x-1,0,0),&d1b1(0,0));
_sse_mulABC_set_331(&U.em(x,2,3,0,0),&Up(x,1,0),&b1u1(0,0));
_sse_mulABC_set_331(&U.em(x,2,3,0,0),&Dw(x,1,0),&b1d1(0,0));
_sse_mulABC_set_331(&U(x,1,0,0),&Up(x+1,2,0),&u1u2(0,0));
_sse_mulAHBC_set_331(&U(x-1,1,0,0),&Dw(x-1,2,0),&d1d2(0,0));
_sse_mulABC_set_331(&U(x,1,0,0),&Dw(x+1,2,0),&u1d2(0,0));
_sse_mulAHBC_set_331(&U(x-1,1,0,0),&Up(x-1,2,0),&d1d2(0,0));
_sse_mulABC_set_331(&U(x,1,0,0),&Bi(x+1,1,0),&u1b2(0,0));
_sse_mulAHBC_set_331(&U(x-1,1,0,0),&Bi(x-1,1,0),&d1b2(0,0));
_sse_mulABC_set_331(&U.em(x,2,3,0,0),&Up(x,2,0),&b1u2(0,0));
_sse_mulABC_set_331(&U.em(x,2,3,0,0),&Dw(x,2,0),&b1d2(0,0));
_sse_mulABC_set_331(&U(x,1,0,0),&Up(x+1,3,0),&u1u3(0,0));
_sse_mulAHBC_set_331(&U(x-1,1,0,0),&Dw(x-1,3,0),&d1d3(0,0));
_sse_mulABC_set_331(&U(x,1,0,0),&Dw(x+1,3,0),&u1d3(0,0));
_sse_mulAHBC_set_331(&U(x-1,1,0,0),&Up(x-1,3,0),&d1d3(0,0));
_sse_mulABC_set_331(&U(x,1,0,0),&Bi(x+1,2,0),&u1b3(0,0));
_sse_mulAHBC_set_331(&U(x-1,1,0,0),&Bi(x-1,2,0),&d1b3(0,0));
_sse_mulABC_set_331(&U.em(x,2,3,0,0),&Up(x,3,0),&b1u3(0,0));
_sse_mulABC_set_331(&U.em(x,2,3,0,0),&Dw(x,3,0),&b1d3(0,0));
_sse_mulABC_set_331(&U(x,2,0,0),&Up(x+2,1,0),&u2u1(0,0));
_sse_mulAHBC_set_331(&U(x-2,2,0,0),&Dw(x-2,1,0),&d2d1(0,0));
_sse_mulABC_set_331(&U(x,2,0,0),&Dw(x+2,1,0),&u2d1(0,0));
_sse_mulAHBC_set_331(&U(x-2,2,0,0),&Up(x-2,1,0),&d2d1(0,0));
_sse_mulABC_set_331(&U(x,2,0,0),&Bi(x+2,0,0),&u2b1(0,0));
_sse_mulAHBC_set_331(&U(x-2,2,0,0),&Bi(x-2,0,0),&d2b1(0,0));
_sse_mulABC_set_331(&U.em(x,1,3,0,0),&Up(x,1,0),&b2u1(0,0));
_sse_mulABC_set_331(&U.em(x,1,3,0,0),&Dw(x,1,0),&b2d1(0,0));
_sse_mulABC_set_331(&U(x,2,0,0),&Up(x+2,2,0),&u2u2(0,0));
_sse_mulAHBC_set_331(&U(x-2,2,0,0),&Dw(x-2,2,0),&d2d2(0,0));
_sse_mulABC_set_331(&U(x,2,0,0),&Bi(x+2,1,0),&u2b2(0,0));
_sse_mulAHBC_set_331(&U(x-2,2,0,0),&Bi(x-2,1,0),&d2b2(0,0));
_sse_mulABC_set_331(&U.em(x,1,3,0,0),&Up(x,2,0),&b2u2(0,0));
_sse_mulABC_set_331(&U.em(x,1,3,0,0),&Dw(x,2,0),&b2d2(0,0));
_sse_mulABC_set_331(&U(x,2,0,0),&Up(x+2,3,0),&u2u3(0,0));
_sse_mulAHBC_set_331(&U(x-2,2,0,0),&Dw(x-2,3,0),&d2d3(0,0));
_sse_mulABC_set_331(&U(x,2,0,0),&Dw(x+2,3,0),&u2d3(0,0));
_sse_mulAHBC_set_331(&U(x-2,2,0,0),&Up(x-2,3,0),&d2d3(0,0));
_sse_mulABC_set_331(&U(x,2,0,0),&Bi(x+2,2,0),&u2b3(0,0));
_sse_mulAHBC_set_331(&U(x-2,2,0,0),&Bi(x-2,2,0),&d2b3(0,0));
_sse_mulABC_set_331(&U.em(x,1,3,0,0),&Up(x,3,0),&b2u3(0,0));
_sse_mulABC_set_331(&U.em(x,1,3,0,0),&Dw(x,3,0),&b2d3(0,0));
_sse_mulABC_set_331(&U(x,3,0,0),&Up(x+3,1,0),&u3u1(0,0));
_sse_mulAHBC_set_331(&U(x-3,3,0,0),&Dw(x-3,1,0),&d3d1(0,0));
_sse_mulABC_set_331(&U(x,3,0,0),&Dw(x+3,1,0),&u3d1(0,0));
_sse_mulAHBC_set_331(&U(x-3,3,0,0),&Up(x-3,1,0),&d3d1(0,0));
_sse_mulABC_set_331(&U(x,3,0,0),&Bi(x+3,0,0),&u3b1(0,0));
_sse_mulAHBC_set_331(&U(x-3,3,0,0),&Bi(x-3,0,0),&d3b1(0,0));
_sse_mulABC_set_331(&U.em(x,1,2,0,0),&Up(x,1,0),&b3u1(0,0));
_sse_mulABC_set_331(&U.em(x,1,2,0,0),&Dw(x,1,0),&b3d1(0,0));
_sse_mulABC_set_331(&U(x,3,0,0),&Up(x+3,2,0),&u3u2(0,0));
_sse_mulAHBC_set_331(&U(x-3,3,0,0),&Dw(x-3,2,0),&d3d2(0,0));
_sse_mulABC_set_331(&U(x,3,0,0),&Dw(x+3,2,0),&u3d2(0,0));
_sse_mulAHBC_set_331(&U(x-3,3,0,0),&Up(x-3,2,0),&d3d2(0,0));
_sse_mulABC_set_331(&U(x,3,0,0),&Bi(x+3,1,0),&u3b2(0,0));
_sse_mulAHBC_set_331(&U(x-3,3,0,0),&Bi(x-3,1,0),&d3b2(0,0));
_sse_mulABC_set_331(&U.em(x,1,2,0,0),&Up(x,2,0),&b3u2(0,0));
_sse_mulABC_set_331(&U.em(x,1,2,0,0),&Dw(x,2,0),&b3d2(0,0));
_sse_mulABC_set_331(&U(x,3,0,0),&Up(x+3,3,0),&u3u3(0,0));
_sse_mulAHBC_set_331(&U(x-3,3,0,0),&Dw(x-3,3,0),&d3d3(0,0));
_sse_mulABC_set_331(&U(x,3,0,0),&Bi(x+3,2,0),&u3b3(0,0));
_sse_mulAHBC_set_331(&U(x-3,3,0,0),&Bi(x-3,2,0),&d3b3(0,0));
_sse_mulABC_set_331(&U.em(x,1,2,0,0),&Up(x,3,0),&b3u3(0,0));
_sse_mulABC_set_331(&U.em(x,1,2,0,0),&Dw(x,3,0),&b3d3(0,0));
for(b=0; b<4; b++) {
k=4*a+b;
_sse_mulAbC_set_31(&psi_in(x,a,0),c_id[k],&out(0,0));
k=16*a+4*b+0;
if(c_u[k]!=0) _sse_mulAbC_add_31(&Up(x,0,0),c_u[k],&out(0,0));
if(c_d[k]!=0) _sse_mulAbC_add_31(&Dw(x,0,0),c_d[k],&out(0,0));
k=16*a+4*b+1;
if(c_u[k]!=0) _sse_mulAbC_add_31(&Up(x,1,0),c_u[k],&out(0,0));
if(c_d[k]!=0) _sse_mulAbC_add_31(&Dw(x,1,0),c_d[k],&out(0,0));
if(c_e[k]!=0) _sse_mulAbC_add_31(&Ei(x,0,0),c_e[k],&out(0,0));
if(c_b[k]!=0) _sse_mulAbC_add_31(&Bi(x,0,0),c_b[k],&out(0,0));
k=16*a+4*b+2;
if(c_u[k]!=0) _sse_mulAbC_add_31(&Up(x,2,0),c_u[k],&out(0,0));
if(c_d[k]!=0) _sse_mulAbC_add_31(&Dw(x,2,0),c_d[k],&out(0,0));
if(c_e[k]!=0) _sse_mulAbC_add_31(&Ei(x,1,0),c_e[k],&out(0,0));
if(c_b[k]!=0) _sse_mulAbC_add_31(&Bi(x,1,0),c_b[k],&out(0,0));
k=16*a+4*b+3;
if(c_u[k]!=0) _sse_mulAbC_add_31(&Up(x,3,0),c_u[k],&out(0,0));
if(c_d[k]!=0) _sse_mulAbC_add_31(&Dw(x,3,0),c_d[k],&out(0,0));
if(c_e[k]!=0) _sse_mulAbC_add_31(&Ei(x,2,0),c_e[k],&out(0,0));
if(c_b[k]!=0) _sse_mulAbC_add_31(&Bi(x,2,0),c_b[k],&out(0,0));
k=64*a+16*b+4*0+0;
if(c_uu[k]!=0) _sse_mulAbC_add_31(&u0u0(0,0),c_uu[k],&out(0,0));
if(c_dd[k]!=0) _sse_mulAbC_add_31(&d0d0(0,0),c_dd[k],&out(0,0));
k=64*a+16*b+4*1+1;
if(c_uu[k]!=0) _sse_mulAbC_add_31(&u1u1(0,0),c_uu[k],&out(0,0));
if(c_dd[k]!=0) _sse_mulAbC_add_31(&d1d1(0,0),c_dd[k],&out(0,0));
if(c_ub[k]!=0) _sse_mulAbC_add_31(&u1b1(0,0),c_ub[k],&out(0,0));
if(c_db[k]!=0) _sse_mulAbC_add_31(&d1b1(0,0),c_db[k],&out(0,0));
if(c_bu[k]!=0) _sse_mulAbC_add_31(&b1u1(0,0),c_bu[k],&out(0,0));
if(c_bd[k]!=0) _sse_mulAbC_add_31(&b1d1(0,0),c_bd[k],&out(0,0));
k=64*a+16*b+4*1+2;
if(c_uu[k]!=0) _sse_mulAbC_add_31(&u1u2(0,0),c_uu[k],&out(0,0));
if(c_dd[k]!=0) _sse_mulAbC_add_31(&d1d2(0,0),c_dd[k],&out(0,0));
if(c_ud[k]!=0) _sse_mulAbC_add_31(&u1d2(0,0),c_ud[k],&out(0,0));
if(c_du[k]!=0) _sse_mulAbC_add_31(&d1u2(0,0),c_du[k],&out(0,0));
if(c_ub[k]!=0) _sse_mulAbC_add_31(&u1b2(0,0),c_ub[k],&out(0,0));
if(c_db[k]!=0) _sse_mulAbC_add_31(&d1b2(0,0),c_db[k],&out(0,0));
if(c_bu[k]!=0) _sse_mulAbC_add_31(&b1u2(0,0),c_bu[k],&out(0,0));
if(c_bd[k]!=0) _sse_mulAbC_add_31(&b1d2(0,0),c_bd[k],&out(0,0));
k=64*a+16*b+4*1+3;
if(c_uu[k]!=0) _sse_mulAbC_add_31(&u1u3(0,0),c_uu[k],&out(0,0));
if(c_dd[k]!=0) _sse_mulAbC_add_31(&d1d3(0,0),c_dd[k],&out(0,0));
if(c_ud[k]!=0) _sse_mulAbC_add_31(&u1d3(0,0),c_ud[k],&out(0,0));
if(c_du[k]!=0) _sse_mulAbC_add_31(&d1u3(0,0),c_du[k],&out(0,0));
if(c_ub[k]!=0) _sse_mulAbC_add_31(&u1b3(0,0),c_ub[k],&out(0,0));
if(c_db[k]!=0) _sse_mulAbC_add_31(&d1b3(0,0),c_db[k],&out(0,0));
if(c_bu[k]!=0) _sse_mulAbC_add_31(&b1u3(0,0),c_bu[k],&out(0,0));
if(c_bd[k]!=0) _sse_mulAbC_add_31(&b1d3(0,0),c_bd[k],&out(0,0));
k=64*a+16*b+4*2+1;
if(c_uu[k]!=0) _sse_mulAbC_add_31(&u2u1(0,0),c_uu[k],&out(0,0));
if(c_dd[k]!=0) _sse_mulAbC_add_31(&d2d1(0,0),c_dd[k],&out(0,0));
if(c_ud[k]!=0) _sse_mulAbC_add_31(&u2d1(0,0),c_ud[k],&out(0,0));
if(c_du[k]!=0) _sse_mulAbC_add_31(&d2u1(0,0),c_du[k],&out(0,0));
if(c_ub[k]!=0) _sse_mulAbC_add_31(&u2b1(0,0),c_ub[k],&out(0,0));
if(c_db[k]!=0) _sse_mulAbC_add_31(&d2b1(0,0),c_db[k],&out(0,0));
if(c_bu[k]!=0) _sse_mulAbC_add_31(&b2u1(0,0),c_bu[k],&out(0,0));
if(c_bd[k]!=0) _sse_mulAbC_add_31(&b2d1(0,0),c_bd[k],&out(0,0));
k=64*a+16*b+4*2+2;
if(c_uu[k]!=0) _sse_mulAbC_add_31(&u2u2(0,0),c_uu[k],&out(0,0));
if(c_dd[k]!=0) _sse_mulAbC_add_31(&d2d2(0,0),c_dd[k],&out(0,0));
if(c_ub[k]!=0) _sse_mulAbC_add_31(&u2b2(0,0),c_ub[k],&out(0,0));
if(c_db[k]!=0) _sse_mulAbC_add_31(&d2b2(0,0),c_db[k],&out(0,0));
if(c_bu[k]!=0) _sse_mulAbC_add_31(&b2u2(0,0),c_bu[k],&out(0,0));
if(c_bd[k]!=0) _sse_mulAbC_add_31(&b2d2(0,0),c_bd[k],&out(0,0));
k=64*a+16*b+4*2+3;
if(c_uu[k]!=0) _sse_mulAbC_add_31(&u2u3(0,0),c_uu[k],&out(0,0));
if(c_dd[k]!=0) _sse_mulAbC_add_31(&d2d3(0,0),c_dd[k],&out(0,0));
if(c_ud[k]!=0) _sse_mulAbC_add_31(&u2d3(0,0),c_ud[k],&out(0,0));
if(c_du[k]!=0) _sse_mulAbC_add_31(&d2u3(0,0),c_du[k],&out(0,0));
if(c_ub[k]!=0) _sse_mulAbC_add_31(&u2b3(0,0),c_ub[k],&out(0,0));
if(c_db[k]!=0) _sse_mulAbC_add_31(&d2b3(0,0),c_db[k],&out(0,0));
if(c_bu[k]!=0) _sse_mulAbC_add_31(&b2u3(0,0),c_bu[k],&out(0,0));
if(c_bd[k]!=0) _sse_mulAbC_add_31(&b2d3(0,0),c_bd[k],&out(0,0));
k=64*a+16*b+4*3+1;
if(c_uu[k]!=0) _sse_mulAbC_add_31(&u3u1(0,0),c_uu[k],&out(0,0));
if(c_dd[k]!=0) _sse_mulAbC_add_31(&d3d1(0,0),c_dd[k],&out(0,0));
if(c_ud[k]!=0) _sse_mulAbC_add_31(&u3d1(0,0),c_ud[k],&out(0,0));
if(c_du[k]!=0) _sse_mulAbC_add_31(&d3u1(0,0),c_du[k],&out(0,0));
if(c_ub[k]!=0) _sse_mulAbC_add_31(&u3b1(0,0),c_ub[k],&out(0,0));
if(c_db[k]!=0) _sse_mulAbC_add_31(&d3b1(0,0),c_db[k],&out(0,0));
if(c_bu[k]!=0) _sse_mulAbC_add_31(&b3u1(0,0),c_bu[k],&out(0,0));
if(c_bd[k]!=0) _sse_mulAbC_add_31(&b3d1(0,0),c_bd[k],&out(0,0));
k=64*a+16*b+4*3+2;
if(c_uu[k]!=0) _sse_mulAbC_add_31(&u3u2(0,0),c_uu[k],&out(0,0));
if(c_dd[k]!=0) _sse_mulAbC_add_31(&d3d2(0,0),c_dd[k],&out(0,0));
if(c_ud[k]!=0) _sse_mulAbC_add_31(&u3d2(0,0),c_ud[k],&out(0,0));
if(c_du[k]!=0) _sse_mulAbC_add_31(&d3u2(0,0),c_du[k],&out(0,0));
if(c_ub[k]!=0) _sse_mulAbC_add_31(&u3b2(0,0),c_ub[k],&out(0,0));
if(c_db[k]!=0) _sse_mulAbC_add_31(&d3b2(0,0),c_db[k],&out(0,0));
if(c_bu[k]!=0) _sse_mulAbC_add_31(&b3u2(0,0),c_bu[k],&out(0,0));
if(c_bd[k]!=0) _sse_mulAbC_add_31(&b3d2(0,0),c_bd[k],&out(0,0));
k=64*a+16*b+4*3+3;
if(c_uu[k]!=0) _sse_mulAbC_add_31(&u3u3(0,0),c_uu[k],&out(0,0));
if(c_dd[k]!=0) _sse_mulAbC_add_31(&d3d3(0,0),c_dd[k],&out(0,0));
if(c_ub[k]!=0) _sse_mulAbC_add_31(&u3b3(0,0),c_ub[k],&out(0,0));
if(c_db[k]!=0) _sse_mulAbC_add_31(&d3b3(0,0),c_db[k],&out(0,0));
if(c_bu[k]!=0) _sse_mulAbC_add_31(&b3u3(0,0),c_bu[k],&out(0,0));
if(c_bd[k]!=0) _sse_mulAbC_add_31(&b3d3(0,0),c_bd[k],&out(0,0));
_sse_sumAC_add_31(&out(0,0),&psi_out(x,b,0));

}

}

}
}
};
#endif
