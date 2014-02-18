/* The libMesh Finite Element Library. */
/* Copyright (C) 2003  Benjamin S. Kirk */

/* This library is free software; you can redistribute it and/or */
/* modify it under the terms of the GNU Lesser General Public */
/* License as published by the Free Software Foundation; either */
/* version 2.1 of the License, or (at your option) any later version. */

/* This library is distributed in the hope that it will be useful, */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU */
/* Lesser General Public License for more details. */

/* You should have received a copy of the GNU Lesser General Public */
/* License along with this library; if not, write to the Free Software */
/* Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA */

// <h1>Systems Example 1 - Stokes Equations</h1>
//
// This example shows how a simple, linear system of equations
// can be solved in parallel.  The system of equations are the familiar
// Stokes equations for low-speed incompressible fluid flow.

// C++ include files that we need
#include <iostream>
#include <algorithm>
#include <cmath>

// Basic include file needed for the mesh functionality.
#include "libmesh/libmesh.h"
#include "libmesh/getpot.h"
#include "libmesh/mesh.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/vtk_io.h"
#include "libmesh/xdr_io.h"
#include "libmesh/equation_systems.h"
#include "libmesh/fe.h"
#include "libmesh/quadrature_gauss.h"
#include "libmesh/dof_map.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/dense_matrix.h"
#include "libmesh/dense_vector.h"
#include "libmesh/transient_system.h"
#include "libmesh/linear_implicit_system.h"
#include "libmesh/zero_function.h"
#include "libmesh/dirichlet_boundaries.h"
#include <libmesh/analytic_function.h>

// For systems of equations the \p DenseSubMatrix
// and \p DenseSubVector provide convenient ways for
// assembling the element matrix and vector on a
// component-by-component basis.
#include "libmesh/dense_submatrix.h"
#include "libmesh/dense_subvector.h"

// The definition of a geometric element
#include "libmesh/elem.h"

#include <petsc.h>

#include "bc/pressureramp.hpp"
#include "util/init.hpp"

// Bring in everything from the libMesh namespace
using namespace libMesh;

void assemble_vel (EquationSystems& es,
                   const std::string& system_name);

Real external_pressure( Point const & /*point*/, Parameters const & /*param*/ )
{
    return 1.;
}

void inlet_vel_x (DenseVector<Number>& output, const Point& p, const Real /*t*/)
{
    //output(0) = p(0)*p(1)*p(1);
    output(0) = p(0);
}

void inlet_vel_y (DenseVector<Number>& output, const Point& p, const Real /*t*/)
{
    //output(1) = - p(0)*p(0)*p(1);
    output(1) = - 2*p(1);
}

// The main program.
int main (int argc, char** argv)
{
    // Initialize libMesh.
    LibMeshInit init (argc, argv);

    GetPot param_file("param.dat");

    Mesh mesh(init.comm(), 2);

    std::string mesh_file = param_file("mesh_file", "structured");

    if( mesh_file == "structured" )
    {
        uint const nx = param_file("nx", 10);
        uint const ny = param_file("ny", 20);
        Real const ox = param_file("ox", 0.0);
        Real const oy = param_file("oy", 0.0);
        Real const lx = param_file("lx", 1.0);
        Real const ly = param_file("ly", 2.0);

        MeshTools::Generation::build_square (mesh,
                                             nx, ny,
                                             ox, lx,
                                             oy, ly,
                                             QUAD9);
    }
    else
    {
        std::cout << "reading mesh from file..." << std::endl;
        mesh.read(mesh_file);
    }

    // Print information about the mesh to the screen.
    mesh.print_info();

    // Create an equation systems object.
    EquationSystems equation_systems (mesh);

    // Declare the system and its variables.
    TransientLinearImplicitSystem & system_vel =
            equation_systems.add_system<TransientLinearImplicitSystem> ("vel");

    const uint u_var = system_vel.add_variable ("u", SECOND);
    const uint v_var = system_vel.add_variable ("v", SECOND);
    const uint p_var = system_vel.add_variable ("p", FIRST);

    const std::set<boundary_id_type> bc_sym = {3}; // left side
    //const std::set<boundary_id_type> bc_wall = {0}; // right side
    const std::set<boundary_id_type> bc_in = {0,1,2}; // bottom
    //const std::set<boundary_id_type> bc_out = {1}; // top

    //system_vel.get_dof_map().add_dirichlet_boundary( DirichletBoundary( bc_wall, {u_var, v_var}, ZeroFunction<Real>() ) );
    system_vel.get_dof_map().add_dirichlet_boundary( DirichletBoundary( bc_sym, {u_var}, ZeroFunction<Real>() ) );
//    system_vel.get_dof_map().add_dirichlet_boundary( DirichletBoundary( bc_sym, {p_var}, ConstFunction<Real>(1.0) ) );
//    system_vel.get_dof_map().add_dirichlet_boundary( DirichletBoundary( bc_sym, {p_var}, PressureRamp(0.1,0.9,2.0) ) );
    system_vel.get_dof_map().add_dirichlet_boundary( DirichletBoundary( bc_in, {u_var}, AnalyticFunction<>(&inlet_vel_x) ) );
    system_vel.get_dof_map().add_dirichlet_boundary( DirichletBoundary( bc_in, {v_var}, AnalyticFunction<>(&inlet_vel_y) ) );
//    system_vel.get_dof_map().add_dirichlet_boundary( DirichletBoundary( bc_in, {p_var}, ConstFunction<>(1.0) ) );
  //  system_vel.get_dof_map().add_dirichlet_boundary( DirichletBoundary( bc_out, {u_var},        ZeroFunction<Real>() ) );

    // Give the system a pointer to the matrix assembly
    // function.
    system_vel.attach_assemble_function (assemble_vel);
    system_vel.attach_init_function (init_zero);

    std::string const output_file = param_file("output_file", "axisym_test.e");

    equation_systems.parameters.set<Real>("t_in") = param_file("t_in", 0.);
    equation_systems.parameters.set<Real>("t_out") = param_file("t_out", 5.0);
    equation_systems.parameters.set<Real>("dt") = param_file("dt", 0.1);

    equation_systems.parameters.set<bool>("axisym") = param_file("axisym", false);

    equation_systems.parameters.set<Real>("f_ux") = param_file("f_u", 0.);
    equation_systems.parameters.set<Real>("f_uy") = param_file("f_v", 0.);

    equation_systems.parameters.set<Real>("mu") = param_file("mu_f", 1e-1);

    equation_systems.parameters.set<uint>("linear solver maximum iterations") = 250;
    equation_systems.parameters.set<Real>("linear solver tolerance") = TOLERANCE;

//    PetscOptionsSetValue("-ksp_monitor_true_residual",PETSC_NULL);

    // Initialize the data structures for the equation system.
    equation_systems.init ();

    // Prints information about the system to the screen.
    equation_systems.print_info();

    system_vel.time = equation_systems.parameters.get<Real>("t_in");
    uint timestep = 0;

#ifdef LIBMESH_HAVE_EXODUS_API
    ExodusII_IO (mesh).write_equation_systems (output_file, equation_systems);
#endif

    Real dt = equation_systems.parameters.get<Real>("dt");
    Real t_out = equation_systems.parameters.get<Real>("t_out");
    while (system_vel.time + dt < t_out + 1e-12)
    {
        timestep++;
        // Incremenet the time counter, set the time and the
        // time step size as parameters in the EquationSystem.
        system_vel.time += dt;

        equation_systems.parameters.set<Real> ("time") = system_vel.time;
        equation_systems.parameters.set<Real> ("dt")   = dt;

        // Recreate any hanging node constraints
        system_vel.get_dof_map().create_dof_constraints(mesh, system_vel.time);

        // Apply any user-defined constraints
        system_vel.user_constrain();

        // Expand any recursive constraints
        system_vel.get_dof_map().process_constraints(mesh);

        // A pretty update message
        std::cout << " Solving time step ";
        {
            std::ostringstream out;

            out << std::setw(2)
                << std::right
                << timestep
                << ", time="
                << std::fixed
                << std::setw(6)
                << std::setprecision(6)
                << std::setfill('0')
                << std::left
                << system_vel.time
                <<  "...";

            std::cout << out.str() << std::endl;
        }

        // At this point we need to update the old
        // solution vector.  The old solution vector
        // will be the current solution vector from the
        // previous time step.  We will do this by extracting the
        // system from the \p EquationSystems object and using
        // vector assignment.  Since only \p TransientSystems
        // (and systems derived from them) contain old solutions
        // we need to specify the system type when we ask for it.
        *system_vel.old_local_solution = *system_vel.current_local_solution;

        // Assemble & solve the linear system
        system_vel.solve();

        // Output evey 1 timesteps to file.
        if ((timestep)%1 == 0)
        {
#ifdef LIBMESH_HAVE_EXODUS_API
            ExodusII_IO io (mesh);
            io.append(true);
            io.write_timestep (output_file, equation_systems, timestep, system_vel.time);
#endif
        }
    }

    // All done.
    return 0;
}

void assemble_vel (EquationSystems& es,
                   const std::string& system_name)
{
    // It is a good idea to make sure we are assembling
    // the proper system.
    libmesh_assert_equal_to (system_name, "vel");

    // Get a constant reference to the mesh object.
    const MeshBase& mesh = es.get_mesh();

    // The dimension that we are running
    const uint dim = mesh.mesh_dimension();

    // Get a reference to the Convection-Diffusion system object.
    TransientLinearImplicitSystem & system =
            es.get_system<TransientLinearImplicitSystem> ("vel");

    // Numeric ids corresponding to each variable in the system
    const uint u_var = system.variable_number ("u");
    const uint v_var = system.variable_number ("v");
    const uint p_var = system.variable_number ("p");

    // Get the Finite Element type for "u".  Note this will be
    // the same as the type for "v".
    FEType fe_u_type = system.variable_type(u_var);
    FEType fe_p_type = system.variable_type(p_var);

    // Build a Finite Element object of the specified type for
    // the velocity variables.
    AutoPtr<FEBase> fe_u (FEBase::build(dim, fe_u_type));
    AutoPtr<FEBase> fe_p (FEBase::build(dim, fe_p_type));
    AutoPtr<FEBase> fe_face (FEBase::build(dim, fe_u_type));

    // A Gauss quadrature rule for numerical integration.
    // Let the \p FEType object decide what order rule is appropriate.
    QGauss qrule (dim, fe_u_type.default_quadrature_order());
    QGauss qface(dim-1, fe_u_type.default_quadrature_order());

    // Tell the finite element objects to use our quadrature rule.
    fe_u->attach_quadrature_rule (&qrule);
    fe_p->attach_quadrature_rule (&qrule);    
    fe_face->attach_quadrature_rule (&qface);

    // Here we define some references to cell-specific data that
    // will be used to assemble the linear system.
    //
    // The element Jacobian * quadrature weight at each integration point.
    const std::vector<Real>& JxW = fe_u->get_JxW();

    // The element shape function gradients for the velocity
    // variables evaluated at the quadrature points.
    const std::vector<std::vector<Real> >& phi = fe_u->get_phi();
    const std::vector<std::vector<RealGradient> >& grad_phi = fe_u->get_dphi();
    const std::vector<std::vector<Real> >& psi = fe_p->get_phi();
    const std::vector<Point >& q_point = fe_u->get_xyz();

    // A reference to the \p DofMap object for this system.  The \p DofMap
    // object handles the index translation from node and element numbers
    // to degree of freedom numbers.  We will talk more about the \p DofMap
    // in future examples.
    const DofMap & dof_map = system.get_dof_map();

    // Define data structures to contain the element matrix
    // and right-hand-side vector contribution.  Following
    // basic finite element terminology we will denote these
    // "Ke" and "Fe".
    DenseMatrix<Number> Ke;
    DenseVector<Number> Fe;

    DenseSubMatrix<Number>
            Kuu(Ke), Kuv(Ke), Kup(Ke),
            Kvu(Ke), Kvv(Ke), Kvp(Ke),
            Kpu(Ke), Kpv(Ke), Kpp(Ke);

    DenseSubVector<Number>
            Fu(Fe),
            Fv(Fe),
            Fp(Fe);

    // This vector will hold the degree of freedom indices for
    // the element.  These define where in the global system
    // the element degrees of freedom get mapped.
    std::vector<dof_id_type> dof_indices;
    std::vector<dof_id_type> dof_indices_u;
    std::vector<dof_id_type> dof_indices_v;
    std::vector<dof_id_type> dof_indices_p;

    // Now we will loop over all the elements in the mesh that
    // live on the local processor. We will compute the element
    // matrix and right-hand-side contribution.  In case users later
    // modify this program to include refinement, we will be safe and
    // will only consider the active elements; hence we use a variant of
    // the \p active_elem_iterator.

    Real dt = es.parameters.get<Real>("dt");
    bool axisym = es.parameters.get<bool>("axisym");
    Real f_u = es.parameters.get<Real>("f_ux");
    Real f_v = es.parameters.get<Real>("f_uy");
    Real mu = es.parameters.get<Real>("mu");

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
        dof_map.dof_indices (elem, dof_indices_u, u_var);
        dof_map.dof_indices (elem, dof_indices_v, v_var);
        dof_map.dof_indices (elem, dof_indices_p, p_var);

        const uint n_dofs   = dof_indices.size();
        const uint n_u_dofs = dof_indices_u.size();
        const uint n_p_dofs = dof_indices_p.size();

        // Compute the element-specific data for the current
        // element.  This involves computing the location of the
        // quadrature points (q_point) and the shape functions
        // (phi, dphi) for the current element.
        fe_u->reinit (elem);
        fe_p->reinit (elem);

        // Zero the element matrix and right-hand side before
        // summing them.  We use the resize member here because
        // the number of degrees of freedom might have changed from
        // the last element.  Note that this will be the case if the
        // element type is different (i.e. the last element was a
        // triangle, now we are on a quadrilateral).
        Ke.resize (n_dofs, n_dofs);
        Fe.resize (n_dofs);

        Kuu.reposition (u_var*n_u_dofs, u_var*n_u_dofs, n_u_dofs, n_u_dofs);
        Kuv.reposition (u_var*n_u_dofs, v_var*n_u_dofs, n_u_dofs, n_u_dofs);
        Kup.reposition (u_var*n_u_dofs, p_var*n_u_dofs, n_u_dofs, n_p_dofs);

        Kvu.reposition (v_var*n_u_dofs, u_var*n_u_dofs, n_u_dofs, n_u_dofs);
        Kvv.reposition (v_var*n_u_dofs, v_var*n_u_dofs, n_u_dofs, n_u_dofs);
        Kvp.reposition (v_var*n_u_dofs, p_var*n_u_dofs, n_u_dofs, n_p_dofs);

        Kpu.reposition (p_var*n_u_dofs, u_var*n_u_dofs, n_p_dofs, n_u_dofs);
        Kpv.reposition (p_var*n_u_dofs, v_var*n_u_dofs, n_p_dofs, n_u_dofs);
        Kpp.reposition (p_var*n_u_dofs, p_var*n_u_dofs, n_p_dofs, n_p_dofs);

        Fu.reposition (u_var*n_u_dofs, n_u_dofs);
        Fv.reposition (v_var*n_u_dofs, n_u_dofs);
        Fp.reposition (p_var*n_u_dofs, n_p_dofs);

        // Now we will build the element matrix.
        for (uint qp=0; qp<qrule.n_points(); qp++)
        {
            Real const invR = 1./q_point[qp](0);
            Real const invR2 = invR * invR;

            Real const JxWqp = (axisym)? JxW[qp]*q_point[qp](0):JxW[qp];

            // compute old solution
            Number u_old = 0.;
            Number v_old = 0.;
            for (uint l=0; l<n_u_dofs; l++)
            {
                u_old += phi[l][qp]*system.old_solution (dof_indices_u[l]);
                v_old += phi[l][qp]*system.old_solution (dof_indices_v[l]);
            }

            for (uint i=0; i<n_u_dofs; i++)
            {
                Fu(i) += JxWqp*(u_old+dt*f_u)*phi[i][qp];
                for (uint j=0; j<n_u_dofs; j++)
                {
                    Kuu(i,j) += JxWqp*( phi[i][qp]*phi[j][qp] +
                                          dt*mu*(
                                              grad_phi[i][qp]*grad_phi[j][qp]        // grad(dt*u) : grad(v)
                                              + grad_phi[i][qp](0)*grad_phi[j][qp](0) // grad(dt*u)^T : grad(v)
                                              + axisym*2.*invR2* phi[i][qp]*phi[j][qp]
                                              )
                                          );
                    Kuv(i,j) += JxWqp*dt*mu*grad_phi[i][qp](1)*grad_phi[j][qp](0); // stress1 \eps(d):\eps(v)
                }
                for (uint j=0; j<n_p_dofs; j++ )
                {
                    Kup(i,j) += -JxWqp*dt*psi[j][qp]*(
                                grad_phi[i][qp](0)
                                + axisym*invR * phi[i][qp]
                                );
                }

                Fv(i) += JxWqp*(v_old+dt*f_v)*phi[i][qp];
                for (uint j=0; j<n_u_dofs; j++)
                {
                    Kvu(i,j) += JxWqp*dt*mu*grad_phi[i][qp](0)*grad_phi[j][qp](1); // stress1 \eps(d):\eps(v)
                    Kvv(i,j) += JxWqp*( phi[i][qp]*phi[j][qp] +
                                       dt*mu*(
                                              grad_phi[i][qp]*grad_phi[j][qp]  // stress1 \eps(d):\eps(v)
                                              + grad_phi[i][qp](1)*grad_phi[j][qp](1)
                                       )
                                  );
                }
                for (uint j=0; j<n_p_dofs; j++)
                {
                    Kvp(i,j) += -JxWqp*dt*psi[j][qp]*grad_phi[i][qp](1);
                }
            }
            for (uint i=0; i<n_p_dofs; i++)
            {
                for (uint j=0; j<n_u_dofs; j++)
                {
                    Kpu(i,j) += -JxWqp*dt*psi[i][qp]*grad_phi[j][qp](0);
                    Kpv(i,j) += -JxWqp*dt*psi[i][qp]*grad_phi[j][qp](1);
                }
            }

        } // end of the quadrature point qp-loop

        // loop on element sides
        for (uint s=0; s<elem->n_sides(); s++)
        {
            // check if side in on boundary
            if (elem->neighbor(s) == NULL)
            {
                if (mesh.boundary_info->has_boundary_id(elem, s, 0)) // inlet
                {
                    // AutoPtr<Elem> side (elem->build_side(s));

                    const std::vector<std::vector<Real> >&  phi_face = fe_face->get_phi();
                    const std::vector<Real>& JxWface = fe_face->get_JxW();
                    const std::vector<Point >& qface_point = fe_face->get_xyz();
                    const std::vector<Point>& normals = fe_face->get_normals();

                    fe_face->reinit(elem, s);

                    for (uint qp=0; qp<qface.n_points(); qp++)
                    {
                        Real const JxWfaceqp = (axisym)? JxWface[qp]*q_point[qp](0):JxWface[qp];

                        // const Real xf = qface_point[qp](0);
                        // const Real yf = qface_point[qp](1);

                        const Real value = external_pressure( qface_point[qp], es.parameters );

                        for (uint i=0; i<phi_face.size(); i++)
                        {
                            Fu(i) -= JxWfaceqp*dt*value*normals[qp](0)*phi_face[i][qp];
                            Fv(i) -= JxWfaceqp*dt*value*normals[qp](1)*phi_face[i][qp];
                        }
                    }
                }
            }
        }

        // If this assembly program were to be used on an adaptive mesh,
        // we would have to apply any hanging node constraint equations.
        dof_map.heterogenously_constrain_element_matrix_and_vector (Ke, Fe, dof_indices);

        // The element matrix and right-hand-side are now built
        // for this element.  Add them to the global matrix and
        // right-hand-side vector.  The \p NumericMatrix::add_matrix()
        // and \p NumericVector::add_vector() members do this for us.

        system.matrix->add_matrix (Ke, dof_indices);
        system.rhs->add_vector    (Fe, dof_indices);
    } // end of element loop

//    system.matrix->close();
//    system.matrix->print_matlab("mat.m");

//    system.rhs->close();
//    system.rhs->print_matlab("rhs.m");

    // That's it.
    return;
}
