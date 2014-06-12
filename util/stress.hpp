#ifndef STRESS_HPP
#define STRESS_HPP

#include "FSI.hpp"

#include <libmesh/libmesh.h>
#include <libmesh/equation_systems.h>


// monolithic block for displacement
void compute_stress(libMesh::EquationSystems& es, std::string sys_name);

// different systems for each disp component
void compute_stress(libMesh::EquationSystems& es);

#endif // STRESS_HPP
