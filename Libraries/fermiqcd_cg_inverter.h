/////////////////////////////////////////////////////////////////
/// @file fermiqcd_cg_inverter.h
/// @version 5-8-2007
/// @author Joseph Schneible <>
///
//////////////////////////////////////////////////////////////////

/// @brief the conjugate gradient inverter
///
/// It inverts mul_Q(psi_out,psi_in,U,coeff) ///not really
/// @param psi_out the output field passed by reference
/// @param psi_in the input field passed by reference
/// @param U the gauge field to be passed to mul_Q
/// @param coeff the gauge parameters to be passed to mul_Q
/// @param absolute_precision the target absolute precision
/// @param relative_precision the target relative precision
/// @param max_steps the maximum number of steps
///
/// Example:
/// @verbatim
/// gauge_field U(lattice,nc);
/// fermi_field psi(lattice,nc);
/// fermi_field chi(lattice,nc);
/// coefficinets coeff;
/// coeff["kappa"]=1.12;
/// U.load("myfield");
/// psi.load("myfield_psi");
/// CG2::inverter(chi,psi,U,coeff);
/// chi.save("myfield_chi");
/// @endverbatim
/// Note that mul_invQ(chi,psi,U,coeff) reads 
/// \f$ \chi=(/\!\!\!D[U]+m)^{-1}\psi \f$ 

class CG2 {
 public:  
  template <class fieldT, class fieldG>
  static inversion_stats inverter(fieldT &psi_out, 
				  fieldT &psi_in, 
				  fieldG &U, 
				  coefficients &coeff, 
				  mdp_real absolute_precision=mdp_precision,
				  mdp_real relative_precision=0,
				  int max_steps=2000,
				  bool qdaggerq=false) {
    mpi.begin_function("ConjugateGradientInverter");
    int             step=0;
    double          residue, rresidue=-1, old_rresidue;
    double          time=mpi.time();
    inversion_stats stats;

    fieldT r(psi_in);
    fieldT rnew(psi_in);
    fieldT p(psi_in);
    fieldT pnew(psi_in);
    fieldT t(psi_in);
    fieldT s(psi_in);
    fieldT solnew(psi_in);
    fieldT b(psi_in);
	
    double alpha, beta, rrtmp, rrtmp2, psdot;

    //mpi << "\tstep\tresidue\t\ttime (sec)\n";


    int mulQ=1, guesszero=1,debug=0;

    psi_out=0; 
    r=psi_in;
    p=r;    

    //mpi << "\t<target>\n" 
	//    << "\t\t<max_steps>" << max_steps << "</max_steps>\n"
	//    << "\t\t<absolute_precision>" << absolute_precision << "</absolute_precision>\n"
	//    << "\t\t<relative_precision>" << relative_precision << "</relative_precision>\n"
	//    << "\t</target>\n";
	 
    do {
      p.update();
      mul_Q(t,p,U,coeff); t/=2.0*coeff["kappa"]; t.update();
      coeff["sign"]=-1; 
      mul_Q(s,t,U,coeff);
      coeff["sign"]=1;
      s/=2.0*coeff["kappa"];
      rrtmp=real_scalar_product(r,r);		
      psdot=real_scalar_product(s,p);	       
      alpha=rrtmp/psdot;
      rnew=r;
      mdp_add_scaled_field(rnew, -alpha, s);	
      solnew=psi_out;
      mdp_add_scaled_field(solnew, alpha, p);
      rrtmp2=real_scalar_product(rnew,rnew);	  
      beta=rrtmp2/rrtmp;	  
      pnew=rnew;
      mdp_add_scaled_field(pnew, beta, p);
      old_rresidue=rresidue;
      residue=sqrt(rrtmp2/(psi_in.global_size()));	  
      rresidue=relative_residue(rnew,psi_out);         
      p=pnew;
      r=rnew;
      psi_out=solnew;	  
      step++;
      
      //mpi << "\t\t<step>" << step << "</step>\n"
      //    << "\t\t<residue>" << residue << "</residue>\n"
      //    << "\t\t<relative_residue>" << rresidue << "</relative_residue>\n"
      //    << "\t\t<time>" << mpi.time()-time << "</time>\n\n";

      // debug: cout << residue << endl;

      //if((step>10) && (rresidue==old_rresidue))
      //    error("not converging"); 
      step++;
      
    } while (residue>absolute_precision && 
	     rresidue>relative_precision && 
	     step<max_steps);
    
    if(!qdaggerq) {
      solnew.update();
      mul_Q(psi_out,solnew,U,coeff);
    }

    psi_out.update();
    
    stats.target_absolute_precision=absolute_precision;
    stats.target_relative_precision=relative_precision;
    stats.max_steps=max_steps;
    stats.absolute_precision=residue;
    stats.relative_precision=rresidue;
    stats.residue=residue;
    stats.steps=step;
    stats.mul_Q_steps=2*step+1;
    stats.time=mpi.time()-time;

    mpi << "\t<stats>\n" 
	<< "\t\t<max_steps>" << step << "</max_steps>\n"
	<< "\t\t<absolute_precision>" << residue << "</absolute_precision>\n"
	<< "\t\t<relative_precision>" << rresidue << "</relative_precision>\n"
	<< "\t\t<time>" << stats.time << "</time>\n"
	<< "\t</stats>\n";
    
    mpi.end_function("ConjugateGradientInverter");
    return stats;
  }
};
