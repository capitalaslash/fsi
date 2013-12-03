#ifndef LAPLACIAN_HPP
#define LAPLACIAN_HPP

#include <libmesh/libmesh.h>
#include <libmesh/equation_systems.h>

using namespace libMesh;

class Laplacian : public System::Assembly
{
public:
    Laplacian(EquationSystems &es):_es(es) {}

    void assemble();

    Real rhs( Point const & p ) const
    {
        return p(0) * (1.-p(0));
    }

protected:
    EquationSystems& _es;
};

#endif // LAPLACIAN_HPP
