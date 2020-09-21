#ifndef VORVQ_VOROBJ_HPP
#include "VorVQ_types.hpp"
#endif


RCPP_MODULE(vor_module){
  using namespace Rcpp; // Added (if not done globally)

  class_<VOR>("VOR")

    .constructor()

    // Initializer
    .field_readonly("W", &VOR::W, "Prototype matrix")
    .field_readonly("nW", &VOR::nW, "Number of prototypes in the tessellation")
    .field_readonly("d", &VOR::d, "Dimenson of prototypes")
    .field_readonly("lb", &VOR::lb, "Lower bound of data dimension")
    .field_readonly("ub", &VOR::ub, "Upper bound of data dimension")
    .field_readonly("CADJ", &VOR::CADJ, "CADJ adjacency matrix")
    .field_readonly("GADJ", &VOR::GADJ, "Gabriel graph adjacency matrix")

    //.method("set_W", &VOR::set_W, "Set the prototypes generating the Voronoi tessellation")
    .method("set_bounds", &VOR::set_bounds, "Set the dimension-wise bounds of the Voronoi tessellation")
    //.method("set_CADJ", &VOR::set_CADJ, "Set the CADJ adjacency matrix")
    .method("calc_GADJ", &VOR::calc_GADJ, "Calculate the Gabriel Graph Adjacency")
    .method("initialize_VOR", &VOR::initialize_VOR, "Initialize the Voronoi object")

    
    // Control 
    .method("set_params", &VOR::set_params, "Set a VOR object parameter")
    .method("get_params", &VOR::get_params, "Get a VOR object parameters")

    // Active indices
    .field_readonly("vor1_active", &VOR::vor1_active, "List of active Vor1 cells")
    .field_readonly("vor2_active", &VOR::vor2_active, "List of active Vor2 cells")
    .method("set_vor1_active", &VOR::set_vor1_active, "Set the list of active Vor1 cells")
    .method("set_vor2_active", &VOR::set_vor2_active, "Set the list of active Vor2 cells")

    // Polydef 
    .method("get_vor1_polytope", &VOR::get_vor1_polytope, "Get a vor1 polytope definition")
    .method("get_vor2_polytope", &VOR::get_vor2_polytope, "Get a vor2 polytope definition")

    // DADJ
    .field_readonly("DADJ", &VOR::DADJ, "Delaunay graph of Voronoi Tessellation")
    .method("calc_DADJ", &VOR::calc_DADJ, "Calculate the Delaunay graph")
    .method("set_DADJ", &VOR::set_DADJ, "Set the Delaunay graph")


    // Interior points
    .field_readonly("vor1_centers", &VOR::vor1_centers, "Chebyshev centers of vor1 cells")
    .field_readonly("vor2_centers", &VOR::vor2_centers, "Chebyshev centers of vor2 cells")

    .method("calc_vor1_centers", &VOR::calc_vor1_centers, "Calc Chebyshev centers of vor1 cells")
    .method("calc_vor2_centers", &VOR::calc_vor2_centers, "Calc Chebyshev centers of vor2 cells")
    .method("set_vor1_centers", &VOR::set_vor1_centers, "Set the centers of vor1 cells")
    .method("set_vor2_centers", &VOR::set_vor2_centers, "Set the centers of vor2 cells")
    .method("clear_vor1_centers", &VOR::clear_vor1_centers, "Clear the centers of vor1 cells")
    .method("clear_vor2_centers", &VOR::clear_vor2_centers, "Clear the centers of vor2 cells")
    .method("get_vor1_centers", &VOR::get_vor1_centers, "Get the centers of vor1 cells")
    .method("get_vor2_centers", &VOR::get_vor2_centers, "Get the centers of vor2 cells")



    
    // MVIE
    .field_readonly("vor1_MVIE_c", &VOR::vor1_MVIE_c, "Centers of active vor1 MVIEs")
    .field_readonly("vor1_MVIE_E", &VOR::vor1_MVIE_E, "Rotations of active vor1 MVIEs")
    .field_readonly("vor1_MVIE_logdet", &VOR::vor1_MVIE_logdet, "Log determinant of rotation matrix defining active vor1 MVIEs")
    .field_readonly("vor1_MVIE_volratio", &VOR::vor1_MVIE_volratio, "Ratio of volumes of rotation matrices defining active vor1 MVIEs")
    .field_readonly("vor1_MVIE_status", &VOR::vor1_MVIE_status, "Return flags from computing active vor1 MVIEs")
    .field_readonly("vor2_MVIE_c", &VOR::vor2_MVIE_c, "Centers of active vor2 MVIEs")
    .field_readonly("vor2_MVIE_E", &VOR::vor2_MVIE_E, "Rotations of active vor2 MVIEs")
    .field_readonly("vor2_MVIE_logdet", &VOR::vor2_MVIE_logdet, "Log determinant of rotation matrix defining active vor2 MVIEs")
    .field_readonly("vor2_MVIE_volratio", &VOR::vor2_MVIE_volratio, "Ratio of volumes of rotation matrices defining active vor2 MVIEs")
    .field_readonly("vor2_MVIE_status", &VOR::vor2_MVIE_status, "Return flags from computing active vor2 MVIEs")

    .method("calc_vor1_MVIE", &VOR::calc_vor1_MVIE, "Calc the active vor1 MVIEs")
    .method("calc_vor2_MVIE", &VOR::calc_vor2_MVIE, "Calc the active vor2 MVIEs")


    // Dikin
    .field_readonly("vor1_Dikin_E", &VOR::vor1_Dikin_E, "The Dikin ellipsoid rotation matrices in vor1 cells")
    .field_readonly("vor2_Dikin_E", &VOR::vor2_Dikin_E, "The Dikin ellipsoid rotation matrices in vor2 cells")

    .method("calc_vor1_Dikin", &VOR::calc_vor1_Dikin, "Calc the active vor1 Dikin ellipsoids")
    .method("calc_vor2_Dikin", &VOR::calc_vor2_Dikin, "Calc the active vor2 Dikin ellipsoids")

    // Calc all 
    .method("calc_all", &VOR::calc_all, "Calc all available Voronoi quantities")

    // Vis 
    .field_readonly("vis_par", &VOR::vis_par, "Plotting par parameters")
    .method("set_vis_par", &VOR::set_vis_par, "Set the plotting par parameters")
    
    // Save & Load 
    .method("save", &VOR::save, "Save a VOR object to disk")
    .method("load", &VOR::load, "Load a VOR object from disk")

    ;
}




