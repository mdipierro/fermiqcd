void select_action_and_inverter(string quark_action, string inverter) {
  if (quark_action == "clover_fast")
    default_fermi_action=FermiCloverActionFast::mul_Q;
  else if (quark_action == "clover_slow")
    default_fermi_action=FermiCloverActionSlow::mul_Q;
#if defined(SSE2)
  else if (quark_action == "clover_sse2")
    default_fermi_action=FermiCloverActionSSE2::mul_Q;
#endif
  else if (quark_action == "staggered_fast")
    default_staggered_action=StaggeredAsqtadActionFast::mul_Q;
  else if (quark_action == "staggered_slow")
    default_staggered_action=StaggeredAsqtadActionSlow::mul_Q;
#if defined(SSE2)
  else if (quark_action == "staggered_sse2")
    default_staggered_action=StaggeredAsqtadActionSSE2::mul_Q;
#endif
  else
    mdp.error_message("quark action not supported");
  if (inverter == "minres")
    default_fermi_inverter=MinRes::inverter<fermi_field,gauge_field>;
  else if (inverter == "bicgstab")
    default_fermi_inverter=BiCGStab::inverter<fermi_field,gauge_field>;
  else if (inverter == "minres-vtk")
    default_fermi_inverter=MinResVtk::inverter<fermi_field,gauge_field>;
  else if (inverter == "bicgstab-vtk")
    default_fermi_inverter=BiCGStabVtk::inverter<fermi_field,gauge_field>;
  else
    mdp.error_message("quark inverter not supported");
}
