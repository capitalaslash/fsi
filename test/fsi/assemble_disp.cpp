
#include <libmesh/libmesh.h>
#include <libmesh/equation_systems.h>
#include <libmesh/transient_system.h>
#include <libmesh/linear_implicit_system.h>
#include <libmesh/quadrature_gauss.h>
#include <libmesh/sparse_matrix.h>

#include "assembly/interface.hpp"

using namespace libMesh;

void assemble_disp (EquationSystems& es,
                    const std::string& system_name)
{
    // It is a good idea to make sure we are assembling
    // the proper system.
    if (! ((system_name == "dx") || (system_name == "dy")))
    {
        libmesh_error();
    }

    // Get a constant reference to the mesh object.
    const MeshBase& mesh = es.get_mesh();

    // The dimension that we are running
    const uint dim = mesh.mesh_dimension();

    // Get a reference to the disp system object.
    TransientLinearImplicitSystem & system =
            es.get_system<TransientLinearImplicitSystem> (system_name);
    TransientLinearImplicitSystem & system_vel =
            es.get_system<TransientLinearImplicitSystem> ("fsi");
    InterfaceSystem & system_i = es.get_system<InterfaceSystem>("int");

    // Numeric ids corresponding to each variable in the system
    const uint d_var = system.variable_number (system_name);

    std::string vel_name;
    if (system_name == "dx")
        vel_name = "ux";
    else if (system_name == "dy")
        vel_name = "uy";

    uint const u_var = system_vel.variable_number(vel_name);

    // Get the Finite Element type for "u".  Note this will be
    // the same as the type for "v".
    FEType fe_d_type = system.variable_type(d_var);
    FEType fe_u_type = system_vel.variable_type(u_var);

    // Build a Finite Element object of the specified type for
    // the velocity variables.
    AutoPtr<FEBase> fe_d (FEBase::build(dim, fe_d_type));
    AutoPtr<FEBase> fe_u (FEBase::build(dim, fe_u_type));

    // A Gauss quadrature rule for numerical integration.
    // Let the \p FEType object decide what order rule is appropriate.
    QGauss qrule (dim, fe_d_type.default_quadrature_order());

    // Tell the finite element objects to use our quadrature rule.
    fe_d->attach_quadrature_rule (&qrule);
    fe_u->attach_quadrature_rule (&qrule);

    // Here we define some references to cell-specific data that
    // will be used to assemble the linear system.
    //
    // The element Jacobian * quadrature weight at each integration point.
    const std::vector<Real>& JxW = fe_d->get_JxW();

    // The element shape function gradients for the velocity
    // variables evaluated at the quadrature points.
    const std::vector<std::vector<Real> >& b = fe_d->get_phi();
    const std::vector<std::vector<Gradient> >& grad_b = fe_d->get_dphi();
    const std::vector<std::vector<Real> >& v = fe_u->get_phi();

    // A reference to the \p DofMap object for this system.  The \p DofMap
    // object handles the index translation from node and element numbers
    // to degree of freedom numbers.  We will talk more about the \p DofMap
    // in future examples.
    const DofMap & dof_map = system.get_dof_map();
    const DofMap & dof_map_vel = system_vel.get_dof_map();

    // Define data structures to contain the element matrix
    // and right-hand-side vector contribution.  Following
    // basic finite element terminology we will denote these
    // "Ke" and "Fe".
    DenseMatrix<Number> Ke;
    DenseVector<Number> Fe;

    // This vector will hold the degree of freedom indices for
    // the element.  These define where in the global system
    // the element degrees of freedom get mapped.
    std::vector<dof_id_type> dof_indices;
    std::vector<dof_id_type> dof_indices_d;
    std::vector<dof_id_type> dof_indices_u;
    std::vector<bool> on_interface;

    // Now we will loop over all the elements in the mesh that
    // live on the local processor. We will compute the element
    // matrix and right-hand-side contribution.  In case users later
    // modify this program to include refinement, we will be safe and
    // will only consider the active elements; hence we use a variant of
    // the \p active_elem_iterator.

    Real dt = es.parameters.get<Real>("dt");

    subdomain_id_type const flag_s = es.parameters.get<subdomain_id_type>("flag_s");
    subdomain_id_type const flag_f = es.parameters.get<subdomain_id_type>("flag_f");

    MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
    const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();

    for ( ; el != end_el; ++el)
    {
        // Store a pointer to the element we are currently
        // working on.  This allows for nicer syntax later.
        const Elem* elem = *el;

        // Get the degree of freedom indices for the
        // current element.  These define where in the global
        // matrix and right-hand-side this element will
        // contribute to.
        dof_map.dof_indices (elem, dof_indices);
        dof_map.dof_indices (elem, dof_indices_d, d_var);
        dof_map_vel.dof_indices (elem, dof_indices_u, u_var);
        system_i.on_interface(elem, on_interface);

        const uint n_dofs   = dof_indices.size();
        const uint n_d_dofs = dof_indices_d.size();
        const uint n_u_dofs = dof_indices_u.size();

        // Compute the element-specific data for the current
        // element.  This involves computing the location of the
        // quadrature points (q_point) and the shape functions
        // (phi, dphi) for the current element.
        fe_d->reinit (elem);
        fe_u->reinit (elem);

        // Zero the element matrix and right-hand side before
        // summing them.  We use the resize member here because
        // the number of degrees of freedom might have changed from
        // the last element.  Note that this will be the case if the
        // element type is different (i.e. the last element was a
        // triangle, now we are on a quadrilateral).
        Ke.resize (n_dofs, n_dofs);
        Fe.resize (n_dofs);

        // solid domain
        if (elem->subdomain_id() == flag_s)
        {

            // Now we will build the element matrix.
            for (uint qp=0; qp<qrule.n_points(); qp++)
            {
                // compute old solution
                Number d_old = 0.;
                for (uint l=0; l<n_d_dofs; l++)
                {
                    d_old += b[l][qp]*system.old_solution (dof_indices_d[l]);
                }
                Number u_old = 0.;
                for (uint l=0; l<n_u_dofs; l++)
                {
                    u_old += v[l][qp]*system_vel.current_solution (dof_indices_u[l]);
                }

                for (uint i=0; i<n_d_dofs; i++)
                {
                    // Fe(i) = system.old_solution(dof_indices_d[i])
                    //         + dt*system_vel.current_solution(dof_indices_u[i]);
                    // Ke(i,i) = 1.;
                    Fe(i) += JxW[qp]*(d_old + dt*u_old)*b[i][qp];
                    for (uint j=0; j<n_d_dofs; j++)
                    {
                        Ke(i,j) += JxW[qp]*b[i][qp]*b[j][qp];
                    }
                }
            } // end of the quadrature point qp-loop
        }

        // fluid domain
        else if (elem->subdomain_id() == flag_f)
        {
            // Now we will build the element matrix.
            for (uint qp=0; qp<qrule.n_points(); qp++)
            {
                for (uint i=0; i<n_d_dofs; i++)
                {
                    if(!on_interface[i])
                    {
                        for (uint j=0; j<n_d_dofs; j++)
                        {
                            Ke(i,j) += JxW[qp]*dt*grad_b[i][qp]*grad_b[j][qp];
                        }
                    }
                }
            } // end of the quadrature point qp-loop
        }

        // If this assembly program were to be used on an adaptive mesh,
        // we would have to apply any hanging node constraint equations.
        dof_map.constrain_element_matrix_and_vector (Ke, Fe, dof_indices);

        // The element matrix and right-hand-side are now built
        // for this element.  Add them to the global matrix and
        // right-hand-side vector.  The \p NumericMatrix::add_matrix()
        // and \p NumericVector::add_vector() members do this for us.
        system.matrix->add_matrix (Ke, dof_indices);
        system.rhs->add_vector    (Fe, dof_indices);
    } // end of element loop

    // std::stringstream ss;
    // ss << "mat_" << system_name << ".m";
    // system.matrix->close();
    // system.matrix->print_matlab(ss.str().c_str());

    // ss.str("");
    // ss << "rhs_" << system_name << ".m";
    // system.rhs->close();
    // system.rhs->print_matlab(ss.str().c_str());

    // That's it.
    return;
}
