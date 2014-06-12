#ifndef MAT_HPP
#define MAT_HPP

#include "FSI.hpp"

#include <libmesh/equation_systems.h>

void assemble_mat( libMesh::EquationSystems & es, std::string const & name );

#endif // MAT_HPP
