#ifndef INTERFACE_HPP
#define INTERFACE_HPP

#include "FSIconfig.h"

#include <libmesh/libmesh.h>
#include <libmesh/equation_systems.h>

void assemble_interface( libMesh::EquationSystems & es, std::string const & name );

#endif // INTERFACE_HPP
