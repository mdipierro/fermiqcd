/////////////////////////////////////////////////////////////////
/// @file fermiqcd_coefficients.h
/// @version 2009-12-21
/// @author Massimo Di Pierro <mdipierro@cs.depaul.edu>
///
/// Container for action parameters
///
/// Distributed under GPL2 License
/// 
/// Created with support from the US Department of Energy
//////////////////////////////////////////////////////////////////

/// @brief container for action parameters
///
/// All FermiQCD actions are classe and share the same prototype.
/// Parameters are passed to the action via coefficients objects which
/// are nothing more than hash tables.
/// 
/// Example:
/// @verbatim
///    gauge_field U(lattice,nc);
///    coefficients gauge;
///    gauge["beta"]=6.0;
///    WilsonGaugeAction::heatbath(U,gauge);
/// @endverbatim
/// Please check the spalling of the variables you store into the 
/// coefficients object (each action has its own coefficients).
///
/// Why? This allows the creating of new actions while reusing inverters
/// and simplify passing parameters to the action.
class coefficients : public map<string,mdp_real> {
 public:
  bool has_key(const string s) const {
    if(this->find(s)==this->end()) return false;
    return true;
  }
  /* for some reason the const declaration does not do anything. removing
  const mdp_real &operator[] (const string s) const {
    cout << "HERE\n";
    if(!has_key(s)) {
      mdp << "coefficient " << s << " is undefined but required" << endl;
      exit(1);
    }
    return static_cast<map<string,mdp_real> >(*this)[s];
  }
  mdp_real &operator[] (const string s) {
    cout << "THERE\n";
    return static_cast<map<string,mdp_real> * >(this)->operator[](s);
  }
  */
};

void dagger(coefficients &coeff) {
  if(!coeff.has_key("sign")) coeff["sign"]=-1;
  else coeff["sign"]=-coeff["sign"];
}

ostream& operator<< (ostream& os, const coefficients& coeff ) {
  begin_function("print_coefficients");
  coefficients::const_iterator iter;
  for (iter=coeff.begin(); iter != coeff.end(); iter++) {
    cout << "<coefficient name=\"" << iter->first << "\">" << iter->second 
	 << "</coefficient>" << '\n';
  }
  end_function("print_coefficients");
  return os;
}
