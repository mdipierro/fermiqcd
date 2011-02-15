/////////////////////////////////////////////////////////////////
/// @file fermiqcd_fermilab_coefficients.h
/// @version 2009-12-21
/// @author Massimo Di Pierro <mdipierro@cs.depaul.edu>
///
/// Included by fermiqcd_fermilab_action.h
///
/// Distributed under GPL2 License
/// 
/// Created with support from the US Department of Energy
//////////////////////////////////////////////////////////////////

// include clover + c1,c2,c3,c4,c5,c10

 mdp_matrix M(4,4), SE[4], SB[4];


 M=1.0+(2.0*42*c3+2.0*18*c4)*kappat+Gamma[0]*2.0*kappat*(mac1-2.0*mac2+6.0*mac4)+Gamma1*2.0*kappat*mac5;
 for(a=0;a<4; a++)
   for(b=0; b<4; b++) 
     c_id[4*a+b]=M(b,a);
 
 M=kappat*Gamma[0]-rt*kappat+Gamma[0]*2.0*kappat*(mac2-mac3-4.0*mac4);
 for(a=0;a<4; a++)
   for(b=0; b<4; b++) 
     c_u[16*a+4*b]=M(b,a);
 
 M=-kappat*Gamma[0]-rt*kappat+Gamma[0]*2.0*kappat*(mac2+mac3-4.0*mac4);
 for(a=0;a<4; a++)
   for(b=0; b<4; b++) 
     c_d[16*a+4*b]=M(b,a);

 M=Gamma[0]*2.0*kappat*(0.5*mac3+mac4)+Gamma1*2.0*kappat*mac5 *(-0.5);
 for(a=0;a<4; a++)
   for(b=0; b<4; b++) 
     c_uu[64*a+16*b+4*0+0]=M(b,a);

M=Gamma[0]*2.0*kappat*(-0.5*mac3+mac4)+Gamma1*2.0*kappat*mac5 *(-0.5);
 for(a=0;a<4; a++)
   for(b=0; b<4; b++) 
     c_dd[64*a+16*b+4*0+0]=M(b,a);

 for(i=1; i<4; i++) {
   switch(i) {
   case 1: j=2; k=3; break;
   case 2: j=1; k=3; break;
   case 3: j=1; k=2; break;
   }   
   SE[i]=(-I)*Sigma[0][i];
   SB[i]=(-I)*Sigma[j][k];

   M=(zeta-12.0*c1-2.0*c2)*kappat*Gamma[i]+(-rs*zeta+(-24.0)*c3+(-8.0)*c4)*kappat;
   for(a=0;a<4; a++)
     for(b=0; b<4; b++)
       c_u[16*a+4*b+i]=M(b,a);
   
   M=(-zeta+12.0*c1+2.0*c2)*kappat*Gamma[i]+(-rs*zeta+(-24.0)*c3+(-8.0)*c4)*kappat;
   for(a=0;a<4; a++)
     for(b=0; b<4; b++) 
       c_d[16*a+4*b+i]=M(b,a);
       
   M=-cE*cSW*zeta*kappat*SE[i];
     for(a=0;a<4; a++)
       for(b=0; b<4; b++) 
	 c_e[16*a+4*b+i]=M(b,a);
   
   M=-cB*cSW*zeta*kappat*SB[i]-16*c10*kappat*SB[i];
   for(a=0;a<4; a++)
     for(b=0; b<4; b++) 
       c_b[16*a+4*b+i]=M(b,a);
 } 

 for(i=1; i<4; i++)
   for(j=1; j<4; j++) {

     M=c1*kappat*(Gamma[i]+Gamma[j])
       +c2*kappat*Gamma[i]*delta(i,j)
       +2.0*c3*kappat
       +2.0*c4*kappat*delta(i,j);
     for(a=0;a<4; a++)
       for(b=0; b<4; b++) 
	 c_uu[64*a+16*b+4*i+j]=M(b,a);

     M=c1*kappat*(Gamma[i]-Gamma[j])
       +c2*kappat*Gamma[i]*delta(i,j)
       +2.0*c3*kappat
       +2.0*c4*kappat*delta(i,j);
     for(a=0;a<4; a++)
       for(b=0; b<4; b++) 
	 c_ud[64*a+16*b+4*i+j]=M(b,a);

     M=-c1*kappat*(Gamma[i]-Gamma[j])
       -c2*kappat*Gamma[i]*delta(i,j)
       +2.0*c3*kappat
       +2.0*c4*kappat*delta(i,j);
     for(a=0;a<4; a++)
       for(b=0; b<4; b++) 
	 c_du[64*a+16*b+4*i+j]=M(b,a);

     M=-c1*kappat*(Gamma[i]+Gamma[j])
       -c2*kappat*Gamma[i]*delta(i,j)
       +2.0*c3*kappat
       +2.0*c4*kappat*delta(i,j);
     for(a=0;a<4; a++)
       for(b=0; b<4; b++) 
	 c_dd[64*a+16*b+4*i+j]=M(b,a);

     M=c5*kappat*SB[i]*Gamma[j]+(1.0-delta(i,j))*2.0*c10*SB[i];
     for(a=0;a<4; a++)
       for(b=0; b<4; b++) 
	 c_bu[64*a+16*b+4*i+j]=M(b,a);

     M=-c5*kappat*SB[i]*Gamma[j]+(1.0-delta(i,j))*2.0*c10*SB[i];
     for(a=0;a<4; a++)
       for(b=0; b<4; b++) 
	 c_bd[64*a+16*b+4*i+j]=M(b,a);

     M=c5*kappat*Gamma[i]*SB[j]+(1.0-delta(i,j))*2.0*c10*SB[j];
     for(a=0;a<4; a++)
       for(b=0; b<4; b++) 
	 c_ub[64*a+16*b+4*i+j]=M(b,a);

     M=-c5*kappat*Gamma[i]*SB[j]+(1.0-delta(i,j))*2.0*c10*SB[j];
     for(a=0;a<4; a++)
       for(b=0; b<4; b++) 
	 c_db[64*a+16*b+4*i+j]=M(b,a);
   }
 
 
