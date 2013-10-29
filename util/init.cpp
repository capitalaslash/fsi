#include "init.hpp"

void init_zero (EquationSystems& es,
               const std::string& system_name)
{
    // Get a reference to the Convection-Diffusion system object.
    ExplicitSystem & system =
            es.get_system<ExplicitSystem>(system_name);

    // Project initial conditions at time t_in
    es.parameters.set<Real> ("time") = system.time = es.parameters.get<Real>("t_in");

    system.project_solution(f_zero, NULL, es.parameters);
}

