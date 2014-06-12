#ifndef ASSEMBLE_DISP_HPP
#define ASSEMBLE_DISP_HPP

#include "FSI.hpp"

#include <libmesh/libmesh.h>
#include <libmesh/equation_systems.h>

void assemble_disp( libMesh::EquationSystems & es, std::string const & name );

#endif // ASSEMBLE_DISP_HPP
