/////////////////////////////////////////////////////////////////
/// @file mdp_fitting_functions.h
/// @version 2009-12-21
/// @author Massimo Di Pierro <mdipierro@cs.depaul.edu>
///
/// Contains mdp_field::update()
///
/// Licensed under GPL2 license
/// Read attached license in file mdp_license.pdf
/// This file cannot be distributed without file mdp_license.pdf
//////////////////////////////////////////////////////////////////


/// Fits y[i], x[i] for i0<=i<in with y=a[0]*x+a[1]
void linear_fit(float *x, Measure *y, mdp_int i0, mdp_int in, Measure *a) {
  mdp_int i;
  double S=0,Sx=0,Sy=0,Sxx=0,Sxy=0,det;
  double dy2;
  for(i=i0; i<in; i++) {
    dy2=pow(y[i].error,2.0);
    if(dy2<=0) dy2=pow(y[i].mean*1e-5,2.0);
    S=S+1.0/dy2;
    Sx=Sx+1.0/dy2*x[i];
    Sy=Sy+y[i].mean/dy2;
    Sxx=Sxx+x[i]*x[i]/dy2;
    Sxy=Sxy+y[i].mean/dy2*x[i];
  }
  det=(S*Sxx-Sx*Sx);
  // this is m
  a[0].mean=(S*Sxy-Sx*Sy)/det;
  a[0].error=sqrt(S/det);
  // this is q
  a[1].mean=(Sxx*Sy-Sx*Sxy)/det;
  a[1].error=sqrt(Sxx/det);
}


///  finds x=xmin that minimizes
///  (*fp)(&x,1,dummy)
///   must be:
///  (*fp)(&ax) > (*fp)(&bx) && (*fp)(&cx) > (*fp)(&bx) 
float golden_rule(float (*fp)(float*, mdp_int, void*), float &xmin, 
		  float ax, float bx, float cx,
		  float tol=0.001, mdp_int niter=100, void *dummy=0) {
  static const float CGOLD=0.3819660;
  mdp_int iter;
  float a,b,d=0,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm;
  float e=0.0;
  a=(ax<cx ? ax:cx);  
  b=(ax>cx ? ax:cx);
  x=w=v=bx;
  fw=fv=fx=(*fp)(&x,1,dummy);
  for(iter=0; iter<niter; iter++) {
    xm=0.5*(a+b);
    tol2=2.0*(tol1=tol*fabs(x)+PRECISION);
    if(fabs(x-xm)<=(tol2-0.5*(b-a))) {
      xmin=x;
      return fx;
    }
    if(fabs(e)>tol1) {
      r=(x-w)*(fx-fv);
      q=(x-v)*(fx-fw);
      p=(x-v)*q-(x-w)*r;
      q=2.0*(q-r);
      if(q>0.0) p=-p;
      q=fabs(q);
      etemp=e;
      e=d;
      if(fabs(p)>=fabs(0.5*q*etemp) || p<=q*(a-x) || p>=q*(b-x))
	d=CGOLD*(e=((x>=xm) ? (a-x) : (b-x)));
      else {
	d=p/q;
	u=x+d;
	if(u-a < tol2 || (b-u) < tol2) {
	  d=fabs(tol1);
	  if(xm - x < 0.0) d=-d;
	}
      }
    } else {
      d=CGOLD*(e=(x>=xm ? (a-x) : (b-x)));
    }
    if(fabs(d)>=tol1) u=x+d; else u=x+(d>0 ? fabs(tol1): -fabs(tol1));
    fu=(*fp)(&u,1,dummy);
    if(fu <= fx) {
      if(u>x) a=x; else b=x;
      v=w; w=x; x=u;
      fv=fw; fw=fx; fx=fu;
    } else {
      if(u<x) a=u; else b=u;
      if(fu<fw || w==x) {
	v=w; w=u;
	fv=fw; fw=fu;
      } else if(fu <= fv || v==x || v==w) {
	v=u;
	fv=fu;
      }
    }
  }
  printf("Too many iterations in golden_rune\n");
  xmin=x;
  return fx;
}

typedef float (*BLM_function)(float, float*, mdp_int, void*);

/////////////////////////////////////////////////////////////////////
/// This function is used by the BayesianLevenbergMarquardt
/// It computes the chi_square (including the Baesyan term)
/// and fills alpha and beta
///
/// alpha(j,k)=\sum_i (Dy(x[i],a)/Da[j])*(Dy(x[i],a)/Da[k])/dy[i]^2
/// beta(j)=sum_i     (y[i]-y(x[i],a))*(dy(x[i],a)/da[j])/dy[i]^2
///
/// chi_square=\sum_i (y[i]-y(x[i],a))*(y[i]-y(x[i],a))/dy[i]^2
///           +\sum_{j,k} (a[j]-a0[j])*(a[k]-a0[k])*sigma(j,k)
///
/// This function take into account multipliticty factors
/// y[i].num, i.e. the numbers of measures used to determine y[i].mean
/// This is used as a weight factor!
/////////////////////////////////////////////////////////////////////

float BLMaux(float *x, Measure *y, 
mdp_int i_min, mdp_int i_max, 
float *a, float *a0, mdp_matrix &sigma, int ma,
	     mdp_matrix &alpha,
	     mdp_matrix &beta,
	     BLM_function func,
	     float h,
	     void *junk) {
  mdp_int i,j,k;
  double sig2i, chi_square, dy, ymod, wt;
  double *dyda=new double[ma];

  for(i=0; i<ma; i++) {
    for(j=0; j<ma; j++) 
      alpha(i,j)=0.0;
    beta(i,0)=0.0;
  }
  chi_square=0.0;
  for(i=i_min; i<i_max; i++) {
    ymod=(*func)(x[i],a,ma,junk);
    for(j=0; j<ma; j++) { 
      a[j]+=h;
      dyda[j]=(*func)(x[i],a,ma,junk);
      a[j]-=2*h;
      dyda[j]-=(*func)(x[i],a,ma,junk);
      a[j]+=h;
      dyda[j]/=(2.0*h);
    }
    sig2i=1.0/pow(y[i].error,2.0);
    dy=y[i].mean-ymod;
    for(j=0; j<ma; j++) {
      wt=dyda[j]*sig2i;
      for(k=0; k<ma; k++)
	alpha(j,k)+=dyda[k]*wt*y[i].num;
      beta(j,0)+=dy*wt*y[i].num;
    }
    chi_square+=dy*dy*sig2i;
  }
  
  for(i=0; i<ma; i++) 
    for(j=0; j<ma; j++)
      chi_square+=(a0[i]-a[i])*sigma(i,j).re*(a0[j]-a[j]);
  
  delete[] dyda;
  return chi_square;
}

/////////////////////////////////////////////////////////////////////
/// This implements the BaesyanLevenbergMarquardt
/// It uses mdp_matrix.
/// Arguments are:
///
/// x[i]          : an array of float
/// y[i]          : an array of Measures
/// i_min, i_max  : range to be used in the fit
///                 points within the range that have y[i].num=0
///                 are ignored
/// a[i], ma      : vector of paramters for the fit and number of parameters
///                 they are all used in the fit
///                 the initial values are used as preons
/// covar(i,j)    : covariance matrix for the preons
/// func(x,a,ma,junk) : the function to be used in the fit
/// h             : a float used to evaluate derivatives 
/// nmax          : max number of iterations
/// junk          : junk to be passed to func
///
/// Return the Baesyan ChiSquare. To obtain the correct chi_square
///  rerun it with same ftting values and nmax=1;
/////////////////////////////////////////////////////////////////////

float BaesyanLevenbergMarquardt(float *x, Measure *y, 
				mdp_int i_min, mdp_int i_max, 
				float *a, int ma,
				mdp_matrix &covar, 
				BLM_function func,
				float h=0.001, 
				mdp_int nmax=1000,
				void *junk=0) {
  
  double lambda=0.1;
  double scale=2;
  double chi_square=0;
  double old_chi_square=0;
  mdp_matrix alpha(ma,ma);
  float *atry=new float[ma];
  float *a0=new float[ma];
  mdp_matrix sigma(ma,ma);
  mdp_matrix beta(ma,1);
  mdp_matrix oneda(ma,1);
  mdp_int i,j,k,n;

  if(max(covar)<=1e-6) 
    for(i=0; i<ma; i++)
      sigma(i,i)=0;
  else 
    sigma=inv(covar);

  for(i=0; i<ma; i++) {
    atry[i]=a[i];
    a0[i]=a[i];
  }
  chi_square=BLMaux(x,y,i_min,i_max,a,a0,sigma,
		    ma,alpha,beta,func,h,junk);
  old_chi_square=chi_square;
  
  for(n=0; n<nmax; n++) {
    for(i=0; i<ma; i++) {
      for(j=0; j<ma; j++) covar(i,j)=alpha(i,j);
      covar(i,i)*=(1.0+lambda);
      oneda(i,0)=beta(i,0);
    }
    
    // this is gaussj 
    covar=inv(covar);
    oneda=covar*oneda;

    for(i=0; i<ma; i++) atry[i]=a[i]+real(oneda(i,0));
    chi_square=BLMaux(x,y,i_min,i_max,atry,a0,sigma,
		      ma,covar,oneda,func,h,junk);


    if(fabs(chi_square - old_chi_square)<1e-4 || lambda==0) {
      covar=inv(covar);
      
      delete[] atry;
      delete[] a0;
      return chi_square;
    } else if(chi_square <= old_chi_square) {
      lambda/=scale;
      old_chi_square=chi_square;
      for(i=0; i<ma; i++) {
	for(j=0; j<ma; j++)
	  alpha(i,j)=covar(i,j);
	beta(i,0)=oneda(i,0);
      }
      for(j=0; j<ma; j++) a[j]=atry[j];
    } else {
      lambda*=scale;
      chi_square=old_chi_square;
    } 
    //    printf("%f %f %f\n", a[0], a[1], lambda);
  }
}

/**

// Example of a program that uses the Leverberger-Marquardt

// fitting function  
float f(float x, float *a, mdp_int ma, void *junk) {
  return a[0]*(exp(-a[1]*x)+exp(-a[1]*(10.0-x)))*sin(a[2]*x);
}

int main() {
  float   x[10]={0,1,2,3,4,5,6,7,8,9};
  Measure  y[10];
  int i,ma=3;
  float a[3]={1000.0, 0.1, 0.1};

  mdp_matrix covar(ma,ma);
  for(i=0; i<10; i++) {
    x[i]=i;
    y[i].set(f(x[i],a,ma,0),0.1);
  }

  a[0]=1; a[1]=1.0, a[2]=0.2;

  BaesyanLevenbergMarquardt(x,y, 0,10, a, ma, covar, f);

  for(i=0; i<10; i++) 
    printf("%f %f (%f) %f\n", 
	   x[i], y[i].mean, y[i].error, f(x[i], a, ma,0));


  
  printf("%f %f %f\n", a[0], a[1], a[2]);
  print(covar);

}

*/
