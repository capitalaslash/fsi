#ifndef STRESS_HPP
#define STRESS_HPP

#include "FSI.hpp"

#include <libmesh/libmesh.h>
#include <libmesh/equation_systems.h>
#include <libmesh/explicit_system.h>
#include <libmesh/linear_implicit_system.h>

// monolithic block for displacement
void compute_stress(libMesh::EquationSystems& es, std::string sys_name);

// different systems for each disp component
void compute_stress(libMesh::EquationSystems& es);

enum CoordType
{
    TWOD,
    TWOD_AXI,
    THREED
};

class StressSystem: public libMesh::ExplicitSystem
{
public:
    StressSystem (libMesh::EquationSystems& es,
                  const std::string& name,
                  const uint number):
        libMesh::ExplicitSystem (es, name, number),
        M_CoordType(TWOD_AXI),
        M_systemDx(es.get_system<libMesh::LinearImplicitSystem>("dx")),
        M_systemDy(es.get_system<libMesh::LinearImplicitSystem>("dy")),
        M_systemDz(es.get_system<libMesh::LinearImplicitSystem>("dx")),
        M_sVar( 3, std::vector<uint>(3) )
    {
        for (uint i=0; i<3; i++)
            for (uint j=0; j<3; j++)
            {
                std::string varname = "s_" + std::to_string(i) + std::to_string(j);
                M_sVar[i][j] = this->add_variable(varname.c_str(), libMesh::CONSTANT, libMesh::MONOMIAL);
            }
        M_vonMisesVar = this->add_variable("vonMises", libMesh::CONSTANT, libMesh::MONOMIAL);

    }

    // build interface position in solution vector
    void assemble();

    template <CoordType>
    void localStress (libMesh::DenseMatrix<Real> & s);

protected:
    CoordType M_CoordType;
    libMesh::LinearImplicitSystem const& M_systemDx;
    libMesh::LinearImplicitSystem const& M_systemDy;
    libMesh::LinearImplicitSystem const& M_systemDz;

    std::vector<std::vector<uint> > M_sVar;
    uint M_vonMisesVar;

};

#endif // STRESS_HPP
