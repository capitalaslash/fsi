#ifndef BIQUAD_HPP
#define BIQUAD_HPP

#include "FSIconfig.h"

#include <libmesh/libmesh.h>
#include <libmesh/mesh.h>
#include <libmesh/mesh_generation.h>
#include <libmesh/getpot.h>
#include <libmesh/elem.h>

using namespace libMesh;

void generate_biquad( Mesh & mesh, GetPot & param );

#endif // BIQUAD_HPP
