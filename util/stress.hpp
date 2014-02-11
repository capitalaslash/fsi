#ifndef STRESS_HPP
#define STRESS_HPP

#include "FSIconfig.h"

#include <libmesh/libmesh.h>
#include <libmesh/equation_systems.h>
#include <libmesh/linear_implicit_system.h>
#include <libmesh/fe_type.h>
#include <libmesh/dof_map.h>
#include <libmesh/quadrature_gauss.h>
#include <libmesh/numeric_vector.h>
#include <libmesh/utility.h>


// monolithic block for displacement
void compute_stress(EquationSystems& es, std::string sys_name);

// different systems for each disp component
void compute_stress(libMesh::EquationSystems& es);

#endif // STRESS_HPP
