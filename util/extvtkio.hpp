#ifndef EXTVTKIO_HPP
#define EXTVTKIO_HPP

#include "FSI.hpp"

#include <libmesh/vtk_io.h>
#include <libmesh/parameters.h>

#ifdef FSI_HAS_LIBXML2
#include <libxml/xmlwriter.h>
#endif

class ExtVTKIO : public libMesh::VTKIO
{
public:
    explicit ExtVTKIO (libMesh::MeshBase const & mesh, libMesh::Parameters const & par);

    void write_solution(const libMesh::EquationSystems & es,
                        const std::set<std::string> *system_names = NULL);
protected:

#ifdef FSI_HAS_LIBXML2
    void write_pvd(const Real t_in, const Real t_out, const Real dt);
#endif

    std::string _dirname;
    std::string const _basename;
    uint const _print_step;
};

#endif // EXTVTKIO_HPP
