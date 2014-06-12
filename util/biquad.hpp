#ifndef BIQUAD_HPP
#define BIQUAD_HPP

#include "FSI.hpp"

#include <libmesh/mesh.h>
#include <libmesh/mesh_generation.h>
#include <libmesh/getpot.h>
#include <libmesh/elem.h>

void generate_biquad( libMesh::Mesh & mesh, GetPot & param );

void generate_barrel( libMesh::Mesh & mesh, GetPot & param );

#endif // BIQUAD_HPP
