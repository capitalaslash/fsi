#ifndef INIT_HPP
#define INIT_HPP

#include "FSIconfig.h"

#include <libmesh/libmesh.h>
#include <libmesh/point.h>
#include <libmesh/parameters.h>
#include <libmesh/equation_systems.h>
#include <libmesh/explicit_system.h>

using namespace libMesh;

inline Number f_zero (const Point& /*p*/,
                      const Parameters& /*parameters*/,
                      const std::string&,
                      const std::string&)
{
    return 0.;
}

void init_zero (EquationSystems& es,
                const std::string& system_name);

#endif // INIT_HPP
