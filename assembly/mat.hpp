#ifndef MAT_HPP
#define MAT_HPP

#include "FSIconfig.h"

#include <libmesh/libmesh.h>
#include <libmesh/equation_systems.h>

void assemble_mat( libMesh::EquationSystems & es, std::string const & name );

#endif // MAT_HPP
