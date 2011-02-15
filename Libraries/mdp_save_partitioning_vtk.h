void save_partitioning_vtk(mdp_lattice &lattice, string filename) {
  mdp_field<int> where(lattice);
  mdp_site x(lattice);
  forallsites(x) {
    where(x)=ME;
  }
  where.save_vtk(filename);
}
