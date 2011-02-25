/////////////////////////////////////////////////////////////////
/// @file fermiqcd.h
/// @version 2009-12-21
/// @author Massimo Di Pierro <mdipierro@cs.depaul.edu>
///
/// Main header file for FermiQCD libraries
///
/// Distributed under GPL2 License
/// 
/// Created with support from the US Department of Energy
//////////////////////////////////////////////////////////////////

#define FERMIQCD 1
#include "mdp.h"

#include "fermiqcd_su_generators.h"
#include "fermiqcd_so_generators.h"
#include "fermiqcd_global_vars.h"
#include "fermiqcd_gamma_matrices.h"
#include "fermiqcd_default_parameters.h"
#include "fermiqcd_check_differences.h"
#include "fermiqcd_set_random.h"
#include "fermiqcd_coefficients.h"

#include "fermiqcd_sse_su3.h"
#include "fermiqcd_gauge_field.h"
#include "fermiqcd_gauge_routines.h"
#include "fermiqcd_gauge_actions.h"
#include "fermiqcd_gauge_algorithms.h"
#include "fermiqcd_gauge_fixing.h" 
#include "fermiqcd_topological_charge.h"

#include "fermiqcd_minres_inverter.h"
#include "fermiqcd_cg_inverter.h"
#include "fermiqcd_bicgstab_inverter.h"

#include "fermiqcd_fermi_field.h"
#include "fermiqcd_fermi_actions.h"
#include "fermiqcd_fermi_actions_sse2.h"
#include "fermiqcd_fermi_algorithms.h"
#include "fermiqcd_fermi_propagator.h"
#include "fermiqcd_fermi_rotation.h"
#include "fermiqcd_fermi_smearing.h"
#include "fermiqcd_lanczos.h"

#include "fermiqcd_staggered_field.h"
#include "fermiqcd_staggered_actions.h"
#include "fermiqcd_staggered_actions_sse2.h"
#include "fermiqcd_staggered_algorithms.h"
#include "fermiqcd_staggered_uml_inverter.h"
#include "fermiqcd_staggered_propagator.h"
#include "fermiqcd_staggered_mesons.h"

#include "fermiqcd_dwfermi_field.h"
#include "fermiqcd_dwfermi_actions.h"
#include "fermiqcd_dwfermi_algorithms.h"
#include "fermiqcd_hmc.h"

// #include "fermiqcd_instanton4d.h"
// #include "fermiqcd_ffts.h"



