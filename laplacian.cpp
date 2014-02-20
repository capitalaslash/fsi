#include <libmesh/libmesh.h>
#include <libmesh/linear_implicit_system.h>
#include <libmesh/dof_map.h>
#include <libmesh/quadrature_gauss.h>
#include <libmesh/sparse_matrix.h>
#include <libmesh/numeric_vector.h>

#include "laplacian.hpp"

void Laplacian::assemble()
{
    PerfLog perf_log("Laplacian::assemble");

    const MeshBase& mesh = _es.get_mesh();

    const uint dim = mesh.mesh_dimension();

    LinearImplicitSystem& system = _es.get_system<LinearImplicitSystem>("lap");

    DofMap const& dof_map = system.get_dof_map();

    FEType fe_type = dof_map.variable_type(0);

    AutoPtr<FEBase> fe(FEBase::build(dim, fe_type));

    QGauss qrule(dim, FIFTH);

    fe->attach_quadrature_rule(&qrule);

    const std::vector<Real>& JxW = fe->get_JxW();
    const std::vector<Point>& q_point = fe->get_xyz();
    const std::vector<std::vector<Real> >& v = fe->get_phi();
    const std::vector<std::vector<RealGradient> >& grad_v = fe->get_dphi();

    DenseMatrix<Number> Ke;
    DenseVector<Number> Fe;

    std::vector<dof_id_type> dof_indices;

    MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
    const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();

    for ( ; el != end_el; ++el)
    {
        perf_log.push("elem init");

        const Elem* elem = *el;

        dof_map.dof_indices (elem, dof_indices);

        const uint n_dofs   = dof_indices.size();

        fe->reinit (elem);

        Ke.resize (n_dofs, n_dofs);
        Fe.resize (n_dofs);

        perf_log.pop("elem init");

        perf_log.push("Ke");
        // Now we will build the element matrix.
        for (uint qp=0; qp<qrule.n_points(); qp++)
        {
            for (uint i=0; i<n_dofs; i++)
            {
                for (uint j=0; j<n_dofs; j++)
                {
                    Ke(i,j) += JxW[qp]*(grad_v[i][qp]*grad_v[j][qp]);
                }
            }
        } // end of the quadrature point qp-loop
        perf_log.pop("Ke");

        perf_log.push("Fe");
        // Now we will build the element matrix.
        for (uint qp=0; qp<qrule.n_points(); qp++)
        {
            for (uint i=0; i<n_dofs; i++)
            {
                Fe(i) += JxW[qp]*this->rhs(q_point[qp])*v[i][qp];
            }
        } // end of the quadrature point qp-loop
        perf_log.pop("Fe");

        // If this assembly program were to be used on an adaptive mesh,
        // we would have to apply any hanging node constraint equations.
        dof_map.constrain_element_matrix_and_vector (Ke, Fe, dof_indices);

        perf_log.push("insertion");
        // The element matrix and right-hand-side are now built
        // for this element.  Add them to the global matrix and
        // right-hand-side vector.  The \p NumericMatrix::add_matrix()
        // and \p NumericVector::add_vector() members do this for us.
        system.matrix->add_matrix (Ke, dof_indices);
        system.rhs->add_vector    (Fe, dof_indices);
        perf_log.pop("insertion");
    } // end of element loop
}
