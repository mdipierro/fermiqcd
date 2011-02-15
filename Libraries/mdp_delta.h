/////////////////////////////////////////////////////////////////
/// @file mdp_delta.h
/// @version 2009-12-21
/// @author Massimo Di Pierro <mdipierro@cs.depaul.edu>
///
/// Contains declaration delta function
///
/// Licensed under GPL2 license
/// Read attached license in file mdp_license.pdf
/// This file cannot be distributed without file mdp_license.pdf
//////////////////////////////////////////////////////////////////

/// True if i==j, false otherwise
template<class T>
const bool delta(const T& i, const T& j) {
  return (i==j)?true:false;
}
   
