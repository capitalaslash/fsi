#ifndef EXTVTKIO_HPP
#define EXTVTKIO_HPP

#include "FSIconfig.h"

#include <libmesh/vtk_io.h>
#include <libmesh/parameters.h>

#ifdef FSI_HAS_LIBXML2
#include <libxml/xmlwriter.h>
#endif

namespace libMesh
{

class ExtVTKIO : public VTKIO
{
public:
    explicit ExtVTKIO (MeshBase const & mesh, Parameters const & par);

    void write_solution(const EquationSystems & es,
                        const std::set<std::string> *system_names = NULL);
protected:

#ifdef FSI_HAS_LIBXML2
    void write_pvd(const Real t_in, const Real t_out, const Real dt);
#endif

    std::string _dirname;
    std::string _basename;
};

}

#endif // EXTVTKIO_HPP
