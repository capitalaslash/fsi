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
#include <math.h>

// Basic include file needed for the mesh functionality.
#include "libmesh/libmesh.h"
#include "libmesh/mesh.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/xdr_io.h"
#include "libmesh/gmsh_io.h"
#include "libmesh/equation_systems.h"
#include "libmesh/fe.h"
#include "libmesh/quadrature_gauss.h"
#include "libmesh/dof_map.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/dense_matrix.h"
#include "libmesh/dense_vector.h"
#include "libmesh/linear_implicit_system.h"
#include "libmesh/transient_system.h"
#include "libmesh/zero_function.h"
#include "libmesh/dirichlet_boundaries.h"

// For systems of equations the \p DenseSubMatrix
// and \p DenseSubVector provide convenient ways for
// assembling the element matrix and vector on a
// component-by-component basis.
#include "libmesh/dense_submatrix.h"
#include "libmesh/dense_subvector.h"

// The definition of a geometric element
#include "libmesh/elem.h"

// Bring in everything from the libMesh namespace
using namespace libMesh;

// Function prototype.  This function will assemble the system
// matrix and right-hand-side.
void assemble_fsi (EquationSystems& es,
                   const std::string& system_name);

void init_fsi (EquationSystems& es,
                   const std::string& system_name);

Number init_zero (const Point& /*p*/,
                   const Parameters& /*parameters*/,
                   const std::string&,
                   const std::string&)
{
    return 0.;
}

struct F
{
    virtual Real operator() (Point const & p) = 0;
};

struct F_f: public F
{
    Real operator() (Point const & p)
    {
        return std::exp(p(0));
    }
};

struct F_s: public F
{
    Real operator() (Point const & p)
    {
        return std::exp(1.);
    }
};

// The main program.
int main (int argc, char** argv)
{
    LibMeshInit init (argc, argv);

    Mesh mesh(init.comm(), 2);
    mesh.read("quad-2.msh");

    mesh.print_info();

    EquationSystems equation_systems (mesh);

    TransientLinearImplicitSystem & system =
            equation_systems.add_system<TransientLinearImplicitSystem> ("fsi");

    const uint dx_var = system.add_variable ("dx", FIRST);
    const uint dy_var = system.add_variable ("dy", FIRST);
    const uint ux_var = system.add_variable ("ux", FIRST);
    const uint uy_var = system.add_variable ("uy", FIRST);
    system.add_variable ("p", FIRST);

    std::set<boundary_id_type> dirichlet_bc;
    dirichlet_bc.insert(2); // right side

    std::vector<uint> vars (4);
    vars[0] = dx_var;
    vars[1] = dy_var;
    vars[2] = ux_var;
    vars[3] = uy_var;

    ZeroFunction<Real> zero;

    system.get_dof_map().add_dirichlet_boundary( libMesh::DirichletBoundary( dirichlet_bc, vars, &zero ) );

    system.attach_assemble_function (assemble_fsi);
    system.attach_init_function (init_fsi);

    std::string const output_file = "fsi_test.e";

    equation_systems.parameters.set<uint>("linear solver maximum iterations") = 250;
    equation_systems.parameters.set<Real>("linear solver tolerance") = TOLERANCE;

    equation_systems.parameters.set<uint>("flag_s") = 12;
    equation_systems.parameters.set<uint>("flag_f") = 11;

    equation_systems.parameters.set<Real>("t_in") = 0.;
    equation_systems.parameters.set<Real>("t_out") = 0.2;
    equation_systems.parameters.set<Real>("dt") = 1.e-3;

    equation_systems.parameters.set<Real>("f_ux") = 0.;
    equation_systems.parameters.set<Real>("f_uy") = 1.;

    Real const E = 1e8;
    Real const ni = 0.3;

    equation_systems.parameters.set<Real>("mu_s") = E / ( 2.*(1.+ni) );
    equation_systems.parameters.set<Real>("lambda") = E*ni / ( (1.+ni)*(1.-2.*ni) );

    equation_systems.parameters.set<uint>("linear solver maximum iterations") = 250;
    equation_systems.parameters.set<Real>("linear solver tolerance") = TOLERANCE;

//    PetscOptionsSetValue("-ksp_monitor_true_residual",PETSC_NULL);

    // Initialize the data structures for the equation system.
    equation_systems.init ();

    // Prints information about the system to the screen.
    equation_systems.print_info();

    system.time   = equation_systems.parameters.get<Real>("t_in");
    uint timestep = 0;

#ifdef LIBMESH_HAVE_EXODUS_API
    ExodusII_IO io (mesh);
    io.write_equation_systems (output_file, equation_systems);
    //io.write_timestep (output_file, equation_systems, 0, system.time);
    io.append(true);
#endif

    Real dt = equation_systems.parameters.get<Real>("dt");
    Real t_out = equation_systems.parameters.get<Real>("t_out");
    while (system.time - dt < t_out + 1e-12)
    {
        timestep++;
        // Incremenet the time counter, set the time and the
        // time step size as parameters in the EquationSystem.
        system.time += dt;

        equation_systems.parameters.set<Real> ("time") = system.time;
        equation_systems.parameters.set<Real> ("dt")   = dt;

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
                << system.time
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
        *system.old_local_solution = *system.current_local_solution;

        // Assemble & solve the linear system
        equation_systems.get_system("fsi").solve();

        // Output evey 1 timesteps to file.
        if ((timestep)%1 == 0)
        {
#ifdef LIBMESH_HAVE_EXODUS_API
            io.write_timestep (output_file, equation_systems, timestep, system.time);
#endif
        }
    }

    // All done.
    return 0;
}

void init_fsi (EquationSystems& es,
               const std::string& system_name)
{
    // It is a good idea to make sure we are initializing
    // the proper system.
    libmesh_assert_equal_to (system_name, "fsi");

    // Get a reference to the Convection-Diffusion system object.
    TransientLinearImplicitSystem & system =
            es.get_system<TransientLinearImplicitSystem>("fsi");

    // Project initial conditions at time t_in
    es.parameters.set<Real> ("time") = system.time = es.parameters.get<Real>("t_in");

    system.project_solution(init_zero, NULL, es.parameters);
}

void assemble_fsi (EquationSystems& es,
                   const std::string& system_name)
{
    // It is a good idea to make sure we are assembling
    // the proper system.
    libmesh_assert_equal_to (system_name, "fsi");

    // Get a constant reference to the mesh object.
    const MeshBase& mesh = es.get_mesh();

    // The dimension that we are running
    const uint dim = mesh.mesh_dimension();

    // Get a reference to the Convection-Diffusion system object.
    TransientLinearImplicitSystem & system =
            es.get_system<TransientLinearImplicitSystem> ("fsi");

    // Numeric ids corresponding to each variable in the system
    const uint dx_var = system.variable_number ("dx");
    const uint dy_var = system.variable_number ("dy");
    const uint ux_var = system.variable_number ("ux");
    const uint uy_var = system.variable_number ("uy");
    const uint p_var = system.variable_number ("p");

    // Get the Finite Element type for "u".  Note this will be
    // the same as the type for "v".
    FEType fe_d_type = system.variable_type(dx_var);
    FEType fe_u_type = system.variable_type(ux_var);
    FEType fe_p_type = system.variable_type(p_var);

    // Build a Finite Element object of the specified type for
    // the velocity variables.
    AutoPtr<FEBase> fe_d (FEBase::build(dim, fe_d_type));
    AutoPtr<FEBase> fe_u (FEBase::build(dim, fe_u_type));
    AutoPtr<FEBase> fe_p (FEBase::build(dim, fe_p_type));

    // A Gauss quadrature rule for numerical integration.
    // Let the \p FEType object decide what order rule is appropriate.
    QGauss qrule (dim, fe_d_type.default_quadrature_order());

    // Tell the finite element objects to use our quadrature rule.
    fe_d->attach_quadrature_rule (&qrule);
    fe_u->attach_quadrature_rule (&qrule);
    fe_p->attach_quadrature_rule (&qrule);

    // Here we define some references to cell-specific data that
    // will be used to assemble the linear system.
    //
    // The element Jacobian * quadrature weight at each integration point.
    const std::vector<Real>& JxW = fe_d->get_JxW();

    // The element shape function gradients for the velocity
    // variables evaluated at the quadrature points.
    const std::vector<std::vector<Real> >& b = fe_d->get_phi();
    const std::vector<std::vector<RealGradient> >& grad_b = fe_d->get_dphi();
    const std::vector<std::vector<Real> >& v = fe_u->get_phi();
    const std::vector<std::vector<RealGradient> >& grad_v = fe_u->get_dphi();
    const std::vector<std::vector<Real> >& q = fe_p->get_phi();

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
            Kdxdx(Ke), Kdxdy(Ke), Kdxux(Ke), Kdxuy(Ke), Kdxp(Ke),
            Kdydx(Ke), Kdydy(Ke), Kdyux(Ke), Kdyuy(Ke), Kdyp(Ke),
            Kuxdx(Ke), Kuxdy(Ke), Kuxux(Ke), Kuxuy(Ke), Kuxp(Ke),
            Kuydx(Ke), Kuydy(Ke), Kuyux(Ke), Kuyuy(Ke), Kuyp(Ke),
            Kpdx(Ke),  Kpdy(Ke),  Kpux(Ke),  Kpuy(Ke),  Kpp(Ke);

    DenseSubVector<Number>
            Fdx(Fe),
            Fdy(Fe),
            Fux(Fe),
            Fuy(Fe),
            Fp(Fe);

    // This vector will hold the degree of freedom indices for
    // the element.  These define where in the global system
    // the element degrees of freedom get mapped.
    std::vector<dof_id_type> dof_indices;
    std::vector<dof_id_type> dof_indices_dx;
    std::vector<dof_id_type> dof_indices_dy;
    std::vector<dof_id_type> dof_indices_ux;
    std::vector<dof_id_type> dof_indices_uy;
    std::vector<dof_id_type> dof_indices_p;

    // Now we will loop over all the elements in the mesh that
    // live on the local processor. We will compute the element
    // matrix and right-hand-side contribution.  In case users later
    // modify this program to include refinement, we will be safe and
    // will only consider the active elements; hence we use a variant of
    // the \p active_elem_iterator.

    Real dt = es.parameters.get<Real>("dt");
    Real f_ux = es.parameters.get<Real>("f_ux");
    Real f_uy = es.parameters.get<Real>("f_uy");
    Real mu_s = es.parameters.get<Real>("mu_s");
    Real lambda = es.parameters.get<Real>("lambda");

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
        dof_map.dof_indices (elem, dof_indices_dx, dx_var);
        dof_map.dof_indices (elem, dof_indices_dy, dy_var);
        dof_map.dof_indices (elem, dof_indices_ux, ux_var);
        dof_map.dof_indices (elem, dof_indices_uy, uy_var);
        dof_map.dof_indices (elem, dof_indices_p, p_var);

        const uint n_dofs   = dof_indices.size();
        const uint n_d_dofs = dof_indices_dx.size();
        const uint n_u_dofs = dof_indices_ux.size();
        const uint n_p_dofs = dof_indices_p.size();

        // Compute the element-specific data for the current
        // element.  This involves computing the location of the
        // quadrature points (q_point) and the shape functions
        // (phi, dphi) for the current element.
        fe_d->reinit (elem);
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

        // Reposition the submatrices...  The idea is this:
        //
        //         -           -          -  -
        //        | Kuu Kuv Kup |        | Fu |
        //   Ke = | Kvu Kvv Kvp |;  Fe = | Fv |
        //        | Kpu Kpv Kpp |        | Fp |
        //         -           -          -  -
        //
        // The \p DenseSubMatrix.repostition () member takes the
        // (row_offset, column_offset, row_size, column_size).
        //
        // Similarly, the \p DenseSubVector.reposition () member
        // takes the (row_offset, row_size)
        Kdxdx.reposition (dx_var*n_d_dofs, dx_var*n_d_dofs, n_d_dofs, n_d_dofs);
        Kdxdy.reposition (dx_var*n_d_dofs, dy_var*n_d_dofs, n_d_dofs, n_d_dofs);
        Kdxux.reposition (dx_var*n_d_dofs, ux_var*n_d_dofs, n_d_dofs, n_u_dofs);
        Kdxuy.reposition (dx_var*n_d_dofs, uy_var*n_d_dofs, n_d_dofs, n_u_dofs);
        Kdxp .reposition (dx_var*n_d_dofs,  p_var*n_d_dofs, n_d_dofs, n_p_dofs);

        Kdydx.reposition (dy_var*n_d_dofs, dx_var*n_d_dofs, n_d_dofs, n_d_dofs);
        Kdydy.reposition (dy_var*n_d_dofs, dy_var*n_d_dofs, n_d_dofs, n_d_dofs);
        Kdyux.reposition (dy_var*n_d_dofs, ux_var*n_d_dofs, n_d_dofs, n_u_dofs);
        Kdyuy.reposition (dy_var*n_d_dofs, uy_var*n_d_dofs, n_d_dofs, n_u_dofs);
        Kdyp .reposition (dy_var*n_d_dofs,  p_var*n_d_dofs, n_d_dofs, n_p_dofs);

        Kuxdx.reposition (ux_var*n_u_dofs, dx_var*n_u_dofs, n_d_dofs, n_d_dofs);
        Kuxdy.reposition (ux_var*n_u_dofs, dy_var*n_u_dofs, n_d_dofs, n_d_dofs);
        Kuxux.reposition (ux_var*n_u_dofs, ux_var*n_u_dofs, n_d_dofs, n_u_dofs);
        Kuxuy.reposition (ux_var*n_u_dofs, uy_var*n_u_dofs, n_d_dofs, n_u_dofs);
        Kuxp .reposition (ux_var*n_u_dofs,  p_var*n_u_dofs, n_d_dofs, n_p_dofs);

        Kuydx.reposition (uy_var*n_u_dofs, dx_var*n_u_dofs, n_d_dofs, n_d_dofs);
        Kuydy.reposition (uy_var*n_u_dofs, dy_var*n_u_dofs, n_d_dofs, n_d_dofs);
        Kuyux.reposition (uy_var*n_u_dofs, ux_var*n_u_dofs, n_d_dofs, n_u_dofs);
        Kuyuy.reposition (uy_var*n_u_dofs, uy_var*n_u_dofs, n_d_dofs, n_u_dofs);
        Kuyp .reposition (uy_var*n_u_dofs,  p_var*n_u_dofs, n_d_dofs, n_p_dofs);

        Kpdx .reposition ( p_var*n_d_dofs, dx_var*n_d_dofs, n_p_dofs, n_d_dofs);
        Kpdy .reposition ( p_var*n_d_dofs, dy_var*n_d_dofs, n_p_dofs, n_d_dofs);
        Kpux .reposition ( p_var*n_d_dofs, ux_var*n_d_dofs, n_p_dofs, n_u_dofs);
        Kpuy .reposition ( p_var*n_d_dofs, uy_var*n_d_dofs, n_p_dofs, n_u_dofs);
        Kpp  .reposition ( p_var*n_d_dofs,  p_var*n_d_dofs, n_p_dofs, n_p_dofs);

        Fdx.reposition (dx_var*n_d_dofs, n_d_dofs);
        Fdy.reposition (dy_var*n_d_dofs, n_d_dofs);
        Fux.reposition (ux_var*n_d_dofs, n_u_dofs);
        Fuy.reposition (uy_var*n_d_dofs, n_u_dofs);
        Fp .reposition ( p_var*n_d_dofs, n_p_dofs);

        // Now we will build the element matrix.
        for (uint qp=0; qp<qrule.n_points(); qp++)
        {
            // compute old solution
            Number dx_old = 0.;
            Number dy_old = 0.;
            Gradient grad_dx_old;
            Gradient grad_dy_old;
            for (uint l=0; l<n_d_dofs; l++)
            {
                dx_old += b[l][qp]*system.old_solution (dof_indices_dx[l]);
                dy_old += b[l][qp]*system.old_solution (dof_indices_dy[l]);
                grad_dx_old.add_scaled (grad_b[l][qp],system.old_solution (dof_indices_dx[l]));
                grad_dy_old.add_scaled (grad_b[l][qp],system.old_solution (dof_indices_dy[l]));
            }
            Number ux_old = 0.;
            Number uy_old = 0.;
            Gradient grad_ux_old;
            Gradient grad_uy_old;
            for (uint l=0; l<n_u_dofs; l++)
            {
                ux_old += v[l][qp]*system.old_solution (dof_indices_ux[l]);
                uy_old += v[l][qp]*system.old_solution (dof_indices_uy[l]);
                grad_ux_old.add_scaled (grad_v[l][qp],system.old_solution (dof_indices_ux[l]));
                grad_uy_old.add_scaled (grad_v[l][qp],system.old_solution (dof_indices_uy[l]));
            }

            if (elem->subdomain_id() == es.parameters.get<uint>("flag_s"))
            {

                for (uint i=0; i<n_d_dofs; i++)
                {
                    Fdx(i) += JxW[qp]*dx_old*b[i][qp];
                    for (uint j=0; j<n_d_dofs; j++)
                    {
                        Kdxdx(i,j) += JxW[qp]*b[i][qp]*b[j][qp];
                    }
                    for (uint j=0; j<n_u_dofs; j++)
                    {
                        Kdxux(i,j) += -JxW[qp]*dt*b[i][qp]*v[j][qp];
                    }

                    Fdy(i) += JxW[qp]*dy_old*b[i][qp];
                    for (uint j=0; j<n_d_dofs; j++)
                    {
                        Kdydy(i,j) += JxW[qp]*b[i][qp]*b[j][qp];
                    }
                    for (uint j=0; j<n_u_dofs; j++)
                    {
                        Kdyuy(i,j) += -JxW[qp]*dt*b[i][qp]*v[j][qp];
                    }

                    Fux(i) += JxW[qp]*( (ux_old+dt*f_ux)*v[i][qp]
                                        - dt*(mu_s*(2.*grad_v[i][qp](0)*grad_dx_old(0)
                                                  +grad_v[i][qp](1)*grad_dx_old(1)
                                                  +grad_v[i][qp](1)*grad_dy_old(0))
                                                  + lambda*(grad_v[i][qp](0)*grad_dx_old(0)+grad_v[i][qp](0)*grad_dy_old(1)))
                                                            );
                                                  for (uint j=0; j<n_d_dofs; j++)
                    {
                                            Kuxdx(i,j) += JxW[qp]*dt*(
                                            mu_s*(2.*grad_v[i][qp](0)*grad_b[j][qp](0)+grad_v[i][qp](1)*grad_b[j][qp](1)) // stress1 \eps(d):\eps(v)
                                            + lambda*grad_v[i][qp](0)*grad_b[j][qp](0) // stress2 tr(\eps(d))*tr(\eps(v))
                                            );
                                            Kuxdy(i,j) += JxW[qp]*dt*(
                                            mu_s*grad_v[i][qp](1)*grad_b[j][qp](0) // stress1 \eps(d):\eps(v)
                                            + lambda*grad_v[i][qp](0)*grad_b[j][qp](1) // stress2 tr(\eps(d))*tr(\eps(v))
                                            );
                                        }
                                        for (uint j=0; j<n_u_dofs; j++)
                    {
                        Kuxux(i,j) += JxW[qp]*(
                                    v[i][qp]*v[j][qp] // mass matrix u*v/dt
                                    );
                    }
                    //                for (uint j=0; j<n_p_dofs; j++ )
                    //                {
                    //                    Kuxp(i,j) += -JxW[qp]*dt*q[j][qp]*grad_v[i][qp](1);
                    //                }

                    Fuy(i) += JxW[qp]*( (uy_old+dt*f_uy)*v[i][qp]
                                        - dt*(mu_s*(grad_v[i][qp](0)*grad_dx_old(1)
                                                  +grad_v[i][qp](0)*grad_dy_old(0)
                                                  +2.*grad_v[i][qp](1)*grad_dy_old(1))
                                                  + lambda*(grad_v[i][qp](1)*grad_dx_old(0)+grad_v[i][qp](1)*grad_dy_old(1)))
                                                            );
                                                  for (uint j=0; j<n_d_dofs; j++)
                    {
                                            Kuydx(i,j) += JxW[qp]*dt*(
                                            mu_s*grad_v[i][qp](0)*grad_b[j][qp](1) // stress1 \eps(d):\eps(v)
                                            + lambda*grad_v[i][qp](1)*grad_b[j][qp](0) // stress2 tr(\eps(d))*tr(\eps(v))
                                            );
                                            Kuydy(i,j) += JxW[qp]*dt*(
                                            mu_s*(grad_v[i][qp](0)*grad_b[j][qp](0) + 2.*grad_v[i][qp](1)*grad_b[j][qp](1)) // stress1 \eps(d):\eps(v)
                                            + lambda*grad_v[i][qp](1)*grad_b[j][qp](1) // stress2 tr(\eps(d))*tr(\eps(v))
                                            );
                                        }
                                        for (uint j=0; j<n_u_dofs; j++)
                    {
                        Kuyuy(i,j) += JxW[qp]*v[i][qp]*v[j][qp];
                    }
                    //                for (uint j=0; j<n_p_dofs; j++)
                    //                {
                    //                    Kuyp(i,j) += -JxW[qp]*dt*q[j][qp]*grad_v[i][qp](1);
                    //                }
                }
                for (uint i=0; i<n_p_dofs; i++)
                {
                    //                Kpp(i,i) = 1.;
                    //                for (uint j=0; j<n_u_dofs; j++)
                    //                {
                    //                    Kpux(i,j) += -JxW[qp]*dt*q[i][qp]*grad_v[j][qp](0);
                    //                    Kpuy(i,j) += -JxW[qp]*dt*q[i][qp]*grad_v[j][qp](1);
                    //                }
                    for(uint j=0; j<n_p_dofs; j++)
                    {
                        Kpp(i,j) += JxW[qp]*dt*q[i][qp]*q[j][qp];
                    }
                }
            }
            else if (elem->subdomain_id() == es.parameters.get<uint>("flag_f"))
            {
                for (uint i=0; i<n_d_dofs; i++)
                {
                    Fdx(i) = 0.;
                    Kdxdx(i,i) = 1.;
                    Fdy(i) = 0.;
                    Kdydy(i,i) = 1.;
                    Fdx(i) = 0.;
                    Kuxux(i,i) = 1.;
                    Fdx(i) = 0.;
                    Kuyuy(i,i) = 1.;
                }
                for (uint i=0; i<n_p_dofs; i++)
                {
                    Fp(i) = 0.;
                    Kpp(i,i) = 1.;
                }
            }
            else
            {
                abort();
            }

        } // end of the quadrature point qp-loop

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

//    system.matrix->close();
//    system.matrix->print_matlab("mat.m");

    //    system.rhs->close();
    //    system.rhs->print();

    // That's it.
    return;
}
