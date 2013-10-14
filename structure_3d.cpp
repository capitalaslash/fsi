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

// For systems of equations the \p DenseSubMatrix
// and \p DenseSubVector provide convenient ways for
// assembling the element matrix and vector on a
// component-by-component basis.
#include "libmesh/dense_submatrix.h"
#include "libmesh/dense_subvector.h"

// The definition of a geometric element
#include "libmesh/elem.h"

#include <petsc.h>

// Bring in everything from the libMesh namespace
using namespace libMesh;

Number init_zero (const Point& /*p*/,
                   const Parameters& /*parameters*/,
                   const std::string&,
                   const std::string&)
{
    return 0.;
}

void init_structure (EquationSystems& es,
                     const std::string& system_name);

// Function prototype.  This function will assemble the system
// matrix and right-hand-side.
void assemble_structure (EquationSystems& es,
                         const std::string& system_name);

void compute_stresses(EquationSystems& es);

// The main program.
int main (int argc, char** argv)
{
    // Initialize libMesh.
    LibMeshInit init (argc, argv);

    Mesh mesh(init.comm());

    MeshTools::Generation::build_cube (mesh,
                                       60, 4, 4,
                                       0., 3.0,
                                       0., 0.2,
                                       0., 0.2,
                                       HEX20);

    //    XdrIO mesh_io(mesh);
    //    mesh_io.read("one_tri.xda");

    // Print information about the mesh to the screen.
    mesh.print_info();

    // Create an equation systems object.
    EquationSystems equation_systems (mesh);

    // Declare the system and its variables.
    // Create a transient system named "Convection-Diffusion"
    TransientLinearImplicitSystem & system =
            equation_systems.add_system<TransientLinearImplicitSystem> ("Structure");

    // Add the variables "u" & "v" to "Stokes".  They
    // will be approximated using second-order approximation.
    const uint dx_var = system.add_variable ("dx", SECOND);
    const uint dy_var = system.add_variable ("dy", SECOND);
    const uint dz_var = system.add_variable ("dz", SECOND);
    const uint ux_var = system.add_variable ("ux", SECOND);
    const uint uy_var = system.add_variable ("uy", SECOND);
    const uint uz_var = system.add_variable ("uz", SECOND);
    system.add_variable ("p", FIRST);

    std::set<boundary_id_type> zero_bc;
    zero_bc.insert(4); // left side

//    std::vector<uint> vars (4);
//    vars[0] = dx_var;
//    vars[1] = dy_var;
//    vars[2] = ux_var;
//    vars[3] = uy_var;

    std::vector<uint> vars (6);
    vars[0] = dx_var;
    vars[1] = dy_var;
    vars[2] = dz_var;
    vars[3] = ux_var;
    vars[4] = uy_var;
    vars[5] = uz_var;

    ZeroFunction<Real> zero;

    system.get_dof_map().add_dirichlet_boundary( libMesh::DirichletBoundary( zero_bc, vars, &zero ) );

    // Give the system a pointer to the matrix assembly
    // function.
    system.attach_assemble_function (assemble_structure);
    system.attach_init_function (init_structure);

    ExplicitSystem & stress =
            equation_systems.add_system<ExplicitSystem> ("Stress");
    stress.add_variable ("s_00", CONSTANT, MONOMIAL);
    stress.add_variable ("s_01", CONSTANT, MONOMIAL);
    stress.add_variable ("s_02", CONSTANT, MONOMIAL);
    stress.add_variable ("s_10", CONSTANT, MONOMIAL);
    stress.add_variable ("s_11", CONSTANT, MONOMIAL);
    stress.add_variable ("s_12", CONSTANT, MONOMIAL);
    stress.add_variable ("s_20", CONSTANT, MONOMIAL);
    stress.add_variable ("s_21", CONSTANT, MONOMIAL);
    stress.add_variable ("s_22", CONSTANT, MONOMIAL);
    stress.add_variable ("vonMises", CONSTANT, MONOMIAL);

    std::string const output_file = "structure3d.e";

    equation_systems.parameters.set<uint>("linear solver maximum iterations") = 250;
    equation_systems.parameters.set<Real>("linear solver tolerance") = TOLERANCE;

    equation_systems.parameters.set<Real>("t_in") = 0.;
    equation_systems.parameters.set<Real>("t_out") = 0.2;
    equation_systems.parameters.set<Real>("dt") = 1.e-3;

    equation_systems.parameters.set<Real>("f_ux") = 0.;
    equation_systems.parameters.set<Real>("f_uy") = 0.;
    equation_systems.parameters.set<Real>("f_uz") = 1.;

    Real const E = 1e8;
    Real const ni = 0.3;

    equation_systems.parameters.set<Real>("mu") = E / ( 2.*(1.+ni) );
    equation_systems.parameters.set<Real>("lambda") = E*ni / ( (1.+ni)*(1.-2.*ni) );

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
    io.append(true);
    io.write_timestep (output_file, equation_systems, 0, system.time);
#endif

    Real dt = equation_systems.parameters.get<Real>("dt");
    Real t_out = equation_systems.parameters.get<Real>("t_out");
    while (system.time + dt < t_out + 1e-12)
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
        equation_systems.get_system("Structure").solve();

        // Post-process the solution to compute the stresses
        compute_stresses(equation_systems);

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

void init_structure (EquationSystems& es,
                     const std::string& system_name)
{
    // It is a good idea to make sure we are initializing
    // the proper system.
    libmesh_assert_equal_to (system_name, "Structure");

    // Get a reference to the Convection-Diffusion system object.
    TransientLinearImplicitSystem & system =
            es.get_system<TransientLinearImplicitSystem>("Structure");

    // Project initial conditions at time t_in
    es.parameters.set<Real> ("time") = system.time = es.parameters.get<Real>("t_in");

    system.project_solution(init_zero, NULL, es.parameters);
}

void assemble_structure (EquationSystems& es,
                         const std::string& system_name)
{
    // It is a good idea to make sure we are assembling
    // the proper system.
    libmesh_assert_equal_to (system_name, "Structure");

    // Get a constant reference to the mesh object.
    const MeshBase& mesh = es.get_mesh();

    // The dimension that we are running
    const uint dim = mesh.mesh_dimension();

    // Get a reference to the Convection-Diffusion system object.
    TransientLinearImplicitSystem & system =
            es.get_system<TransientLinearImplicitSystem> ("Structure");

    // Numeric ids corresponding to each variable in the system
    const uint dx_var = system.variable_number ("dx");
    const uint dy_var = system.variable_number ("dy");
    const uint dz_var = system.variable_number ("dz");
    const uint ux_var = system.variable_number ("ux");
    const uint uy_var = system.variable_number ("uy");
    const uint uz_var = system.variable_number ("uz");
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
            Kdxdx(Ke), Kdxdy(Ke), Kdxdz(Ke), Kdxux(Ke), Kdxuy(Ke), Kdxuz(Ke), Kdxp(Ke),
            Kdydx(Ke), Kdydy(Ke), Kdydz(Ke), Kdyux(Ke), Kdyuy(Ke), Kdyuz(Ke), Kdyp(Ke),
            Kdzdx(Ke), Kdzdy(Ke), Kdzdz(Ke), Kdzux(Ke), Kdzuy(Ke), Kdzuz(Ke), Kdzp(Ke),
            Kuxdx(Ke), Kuxdy(Ke), Kuxdz(Ke), Kuxux(Ke), Kuxuy(Ke), Kuxuz(Ke), Kuxp(Ke),
            Kuydx(Ke), Kuydy(Ke), Kuydz(Ke), Kuyux(Ke), Kuyuy(Ke), Kuyuz(Ke), Kuyp(Ke),
            Kuzdx(Ke), Kuzdy(Ke), Kuzdz(Ke), Kuzux(Ke), Kuzuy(Ke), Kuzuz(Ke), Kuzp(Ke),
            Kpdx(Ke),  Kpdy(Ke),  Kpdz(Ke),  Kpux(Ke),  Kpuy(Ke),  Kpuz(Ke),  Kpp(Ke);

    DenseSubVector<Number>
            Fdx(Fe),
            Fdy(Fe),
            Fdz(Fe),
            Fux(Fe),
            Fuy(Fe),
            Fuz(Fe),
            Fp(Fe);

    // This vector will hold the degree of freedom indices for
    // the element.  These define where in the global system
    // the element degrees of freedom get mapped.
    std::vector<dof_id_type> dof_indices;
    std::vector<dof_id_type> dof_indices_dx;
    std::vector<dof_id_type> dof_indices_dy;
    std::vector<dof_id_type> dof_indices_dz;
    std::vector<dof_id_type> dof_indices_ux;
    std::vector<dof_id_type> dof_indices_uy;
    std::vector<dof_id_type> dof_indices_uz;
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
    Real f_uz = es.parameters.get<Real>("f_uz");
    Real mu = es.parameters.get<Real>("mu");
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
        dof_map.dof_indices (elem, dof_indices_dz, dz_var);
        dof_map.dof_indices (elem, dof_indices_ux, ux_var);
        dof_map.dof_indices (elem, dof_indices_uy, uy_var);
        dof_map.dof_indices (elem, dof_indices_uz, uz_var);
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
        Kdxdz.reposition (dx_var*n_d_dofs, dz_var*n_d_dofs, n_d_dofs, n_d_dofs);
        Kdxux.reposition (dx_var*n_d_dofs, ux_var*n_d_dofs, n_d_dofs, n_u_dofs);
        Kdxuy.reposition (dx_var*n_d_dofs, uy_var*n_d_dofs, n_d_dofs, n_u_dofs);
        Kdxuz.reposition (dx_var*n_d_dofs, uz_var*n_d_dofs, n_d_dofs, n_u_dofs);
        Kdxp .reposition (dx_var*n_d_dofs,  p_var*n_d_dofs, n_d_dofs, n_p_dofs);

        Kdydx.reposition (dy_var*n_d_dofs, dx_var*n_d_dofs, n_d_dofs, n_d_dofs);
        Kdydy.reposition (dy_var*n_d_dofs, dy_var*n_d_dofs, n_d_dofs, n_d_dofs);
        Kdydz.reposition (dy_var*n_d_dofs, dz_var*n_d_dofs, n_d_dofs, n_d_dofs);
        Kdyux.reposition (dy_var*n_d_dofs, ux_var*n_d_dofs, n_d_dofs, n_u_dofs);
        Kdyuy.reposition (dy_var*n_d_dofs, uy_var*n_d_dofs, n_d_dofs, n_u_dofs);
        Kdyuz.reposition (dy_var*n_d_dofs, uz_var*n_d_dofs, n_d_dofs, n_u_dofs);
        Kdyp .reposition (dy_var*n_d_dofs,  p_var*n_d_dofs, n_d_dofs, n_p_dofs);

        Kdzdx.reposition (dz_var*n_d_dofs, dx_var*n_d_dofs, n_d_dofs, n_d_dofs);
        Kdzdy.reposition (dz_var*n_d_dofs, dy_var*n_d_dofs, n_d_dofs, n_d_dofs);
        Kdzdz.reposition (dz_var*n_d_dofs, dz_var*n_d_dofs, n_d_dofs, n_d_dofs);
        Kdzux.reposition (dz_var*n_d_dofs, ux_var*n_d_dofs, n_d_dofs, n_u_dofs);
        Kdzuy.reposition (dz_var*n_d_dofs, uy_var*n_d_dofs, n_d_dofs, n_u_dofs);
        Kdzuz.reposition (dz_var*n_d_dofs, uz_var*n_d_dofs, n_d_dofs, n_u_dofs);
        Kdzp .reposition (dz_var*n_d_dofs,  p_var*n_d_dofs, n_d_dofs, n_p_dofs);

        Kuxdx.reposition (ux_var*n_u_dofs, dx_var*n_u_dofs, n_u_dofs, n_d_dofs);
        Kuxdy.reposition (ux_var*n_u_dofs, dy_var*n_u_dofs, n_u_dofs, n_d_dofs);
        Kuxdz.reposition (ux_var*n_u_dofs, dz_var*n_u_dofs, n_u_dofs, n_d_dofs);
        Kuxux.reposition (ux_var*n_u_dofs, ux_var*n_u_dofs, n_u_dofs, n_u_dofs);
        Kuxuy.reposition (ux_var*n_u_dofs, uy_var*n_u_dofs, n_u_dofs, n_u_dofs);
        Kuxuz.reposition (ux_var*n_u_dofs, uz_var*n_u_dofs, n_u_dofs, n_u_dofs);
        Kuxp .reposition (ux_var*n_u_dofs,  p_var*n_u_dofs, n_u_dofs, n_p_dofs);

        Kuydx.reposition (uy_var*n_u_dofs, dx_var*n_u_dofs, n_u_dofs, n_d_dofs);
        Kuydy.reposition (uy_var*n_u_dofs, dy_var*n_u_dofs, n_u_dofs, n_d_dofs);
        Kuydz.reposition (uy_var*n_u_dofs, dz_var*n_u_dofs, n_u_dofs, n_d_dofs);
        Kuyux.reposition (uy_var*n_u_dofs, ux_var*n_u_dofs, n_u_dofs, n_u_dofs);
        Kuyuy.reposition (uy_var*n_u_dofs, uy_var*n_u_dofs, n_u_dofs, n_u_dofs);
        Kuyuz.reposition (uy_var*n_u_dofs, uz_var*n_u_dofs, n_u_dofs, n_u_dofs);
        Kuyp .reposition (uy_var*n_u_dofs,  p_var*n_u_dofs, n_u_dofs, n_p_dofs);

        Kuzdx.reposition (uz_var*n_u_dofs, dx_var*n_u_dofs, n_u_dofs, n_d_dofs);
        Kuzdy.reposition (uz_var*n_u_dofs, dy_var*n_u_dofs, n_u_dofs, n_d_dofs);
        Kuzdz.reposition (uz_var*n_u_dofs, dz_var*n_u_dofs, n_u_dofs, n_d_dofs);
        Kuzux.reposition (uz_var*n_u_dofs, ux_var*n_u_dofs, n_u_dofs, n_u_dofs);
        Kuzuy.reposition (uz_var*n_u_dofs, uy_var*n_u_dofs, n_u_dofs, n_u_dofs);
        Kuzuz.reposition (uz_var*n_u_dofs, uz_var*n_u_dofs, n_u_dofs, n_u_dofs);
        Kuzp .reposition (uz_var*n_u_dofs,  p_var*n_u_dofs, n_u_dofs, n_p_dofs);

        Kpdx .reposition ( p_var*n_d_dofs, dx_var*n_d_dofs, n_p_dofs, n_d_dofs);
        Kpdy .reposition ( p_var*n_d_dofs, dy_var*n_d_dofs, n_p_dofs, n_d_dofs);
        Kpdz .reposition ( p_var*n_d_dofs, dz_var*n_d_dofs, n_p_dofs, n_d_dofs);
        Kpux .reposition ( p_var*n_d_dofs, ux_var*n_d_dofs, n_p_dofs, n_u_dofs);
        Kpuy .reposition ( p_var*n_d_dofs, uy_var*n_d_dofs, n_p_dofs, n_u_dofs);
        Kpuz .reposition ( p_var*n_d_dofs, uz_var*n_d_dofs, n_p_dofs, n_u_dofs);
        Kpp  .reposition ( p_var*n_d_dofs,  p_var*n_d_dofs, n_p_dofs, n_p_dofs);

        Fdx.reposition (dx_var*n_d_dofs, n_d_dofs);
        Fdy.reposition (dy_var*n_d_dofs, n_d_dofs);
        Fdy.reposition (dz_var*n_d_dofs, n_d_dofs);
        Fux.reposition (ux_var*n_d_dofs, n_u_dofs);
        Fuy.reposition (uy_var*n_d_dofs, n_u_dofs);
        Fuz.reposition (uz_var*n_d_dofs, n_u_dofs);
        Fp .reposition ( p_var*n_d_dofs, n_p_dofs);

        // Now we will build the element matrix.
        for (uint qp=0; qp<qrule.n_points(); qp++)
        {
            // compute old solution
            Number dx_old = 0.;
            Number dy_old = 0.;
            Number dz_old = 0.;
            Gradient grad_dx_old;
            Gradient grad_dy_old;
            Gradient grad_dz_old;
            for (uint l=0; l<n_d_dofs; l++)
            {
                dx_old += b[l][qp]*system.old_solution (dof_indices_dx[l]);
                dy_old += b[l][qp]*system.old_solution (dof_indices_dy[l]);
                dz_old += b[l][qp]*system.old_solution (dof_indices_dz[l]);
                grad_dx_old.add_scaled (grad_b[l][qp],system.old_solution (dof_indices_dx[l]));
                grad_dy_old.add_scaled (grad_b[l][qp],system.old_solution (dof_indices_dy[l]));
                grad_dz_old.add_scaled (grad_b[l][qp],system.old_solution (dof_indices_dz[l]));
            }
            Number ux_old = 0.;
            Number uy_old = 0.;
            Number uz_old = 0.;
//            Gradient grad_ux_old;
//            Gradient grad_uy_old;
//            Gradient grad_uz_old;
            for (uint l=0; l<n_u_dofs; l++)
            {
                ux_old += v[l][qp]*system.old_solution (dof_indices_ux[l]);
                uy_old += v[l][qp]*system.old_solution (dof_indices_uy[l]);
                uz_old += v[l][qp]*system.old_solution (dof_indices_uz[l]);
//                grad_ux_old.add_scaled (grad_v[l][qp],system.old_solution (dof_indices_ux[l]));
//                grad_uy_old.add_scaled (grad_v[l][qp],system.old_solution (dof_indices_uy[l]));
//                grad_uz_old.add_scaled (grad_v[l][qp],system.old_solution (dof_indices_uz[l]));
            }

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

                Fdz(i) += JxW[qp]*dz_old*b[i][qp];
                for (uint j=0; j<n_d_dofs; j++)
                {
                    Kdzdz(i,j) += JxW[qp]*b[i][qp]*b[j][qp];
                }
                for (uint j=0; j<n_u_dofs; j++)
                {
                    Kdzuz(i,j) += -JxW[qp]*dt*b[i][qp]*v[j][qp];
                }

                Fux(i) += JxW[qp]*
                        (
                            ( ux_old + dt*f_ux )*v[i][qp]
                            - dt*
                            (
                                mu*
                                (
                                    grad_v[i][qp]*grad_dx_old          // grad(d_old) : grad(v)
                                    + grad_v[i][qp](0)*grad_dx_old(0)  // grad(d_old)^T : grad(v)
                                    + grad_v[i][qp](1)*grad_dy_old(0)  // |
                                    + grad_v[i][qp](2)*grad_dz_old(0)  // |
                                 )
                                 + lambda*
                                 (
                                    grad_v[i][qp](0)*grad_dx_old(0)    // tr(eps(d_old)) * grad(v)_x
                                    + grad_v[i][qp](0)*grad_dy_old(1)  // |
                                    + grad_v[i][qp](0)*grad_dz_old(2)  // |
                                 )
                             )
                        );
                for (uint j=0; j<n_d_dofs; j++)
                {
                    Kuxdx(i,j) += JxW[qp]*dt*(
                                        mu*( 2.*grad_v[i][qp](0)*grad_b[j][qp](0)
                                        + grad_v[i][qp](1)*grad_b[j][qp](1)
                                        + grad_v[i][qp](2)*grad_b[j][qp](2)
                                        ) // stress1 \eps(d):\eps(v)
                                        + lambda*grad_v[i][qp](0)*grad_b[j][qp](0) // stress2 tr(\eps(d))*tr(\eps(v))
                    );
                    Kuxdy(i,j) += JxW[qp]*dt*(
                                mu*grad_v[i][qp](1)*grad_b[j][qp](0) // stress1 \eps(d):\eps(v)
                                + lambda*grad_v[i][qp](0)*grad_b[j][qp](1) // stress2 tr(\eps(d))*tr(\eps(v))
                                );
                    Kuxdz(i,j) += JxW[qp]*dt*(
                                mu*grad_v[i][qp](2)*grad_b[j][qp](0) // stress1 \eps(d):\eps(v)
                                + lambda*grad_v[i][qp](0)*grad_b[j][qp](2) // stress2 tr(\eps(d))*tr(\eps(v))
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

                Fuy(i) += JxW[qp]*( ( uy_old + dt*f_uy )*v[i][qp]
                                    - dt*( mu*( grad_v[i][qp]*grad_dy_old
                                               + grad_v[i][qp](0)*grad_dx_old(1)
                                               + grad_v[i][qp](1)*grad_dy_old(1)
                                               + grad_v[i][qp](2)*grad_dz_old(1))
                                              + lambda*(
                                                    grad_v[i][qp](1)*grad_dx_old(0)
                                                    + grad_v[i][qp](1)*grad_dy_old(1)
                                                    + grad_v[i][qp](1)*grad_dz_old(2)))
                                                       );
                for (uint j=0; j<n_d_dofs; j++)
                {
                    Kuydx(i,j) += JxW[qp]*dt*(
                                        mu*grad_v[i][qp](0)*grad_b[j][qp](1) // stress1 \eps(d):\eps(v)
                                        + lambda*grad_v[i][qp](1)*grad_b[j][qp](0) // stress2 tr(\eps(d))*tr(\eps(v))
                                        );
                    Kuydy(i,j) += JxW[qp]*dt*(
                                mu*(
                                    grad_v[i][qp](0)*grad_b[j][qp](0)
                                    + 2.*grad_v[i][qp](1)*grad_b[j][qp](1)
                                    + grad_v[i][qp](2)*grad_b[j][qp](2) ) // stress1 \eps(d):\eps(v)
                                    + lambda*grad_v[i][qp](1)*grad_b[j][qp](1) // stress2 tr(\eps(d))*tr(\eps(v))
                                    );
                    Kuydz(i,j) += JxW[qp]*dt*(
                                        mu*grad_v[i][qp](2)*grad_b[j][qp](1) // stress1 \eps(d):\eps(v)
                                        + lambda*grad_v[i][qp](1)*grad_b[j][qp](2) // stress2 tr(\eps(d))*tr(\eps(v))
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

                Fuz(i) += JxW[qp]*( ( uz_old + dt*f_uz )*v[i][qp]
                                    - dt*( mu*( grad_v[i][qp]*grad_dz_old
                                               + grad_v[i][qp](0)*grad_dx_old(2)
                                               + grad_v[i][qp](1)*grad_dy_old(2)
                                               + grad_v[i][qp](2)*grad_dz_old(2))
                                              + lambda*(
                                                    grad_v[i][qp](2)*grad_dx_old(0)
                                                    + grad_v[i][qp](2)*grad_dy_old(1)
                                                    + grad_v[i][qp](2)*grad_dz_old(2)))
                                                       );
                for (uint j=0; j<n_d_dofs; j++)
                {
                    Kuzdx(i,j) += JxW[qp]*dt*(
                                        mu*grad_v[i][qp](0)*grad_b[j][qp](2) // stress1 \eps(d):\eps(v)
                                        + lambda*grad_v[i][qp](2)*grad_b[j][qp](0) // stress2 tr(\eps(d))*tr(\eps(v))
                                        );
                    Kuzdy(i,j) += JxW[qp]*dt*(
                                        mu*grad_v[i][qp](1)*grad_b[j][qp](2) // stress1 \eps(d):\eps(v)
                                        + lambda*grad_v[i][qp](2)*grad_b[j][qp](1) // stress2 tr(\eps(d))*tr(\eps(v))
                                        );
                    Kuzdz(i,j) += JxW[qp]*dt*(
                                mu*(
                                    grad_v[i][qp](0)*grad_b[j][qp](0)
                                    + grad_v[i][qp](1)*grad_b[j][qp](1)
                                    + 2.*grad_v[i][qp](2)*grad_b[j][qp](2) ) // stress1 \eps(d):\eps(v)
                                    + lambda*grad_v[i][qp](2)*grad_b[j][qp](2) // stress2 tr(\eps(d))*tr(\eps(v))
                                    );
                }
                for (uint j=0; j<n_u_dofs; j++)
                {
                    Kuzuz(i,j) += JxW[qp]*v[i][qp]*v[j][qp];
                }
//                for (uint j=0; j<n_p_dofs; j++)
//                {
//                    Kuzp(i,j) += -JxW[qp]*dt*q[j][qp]*grad_v[i][qp](2);
//                }
            }
            for (uint i=0; i<n_p_dofs; i++)
            {
//                Kpp(i,i) = 1.;
//                for (uint j=0; j<n_u_dofs; j++)
//                {
//                    Kpux(i,j) += -JxW[qp]*dt*q[i][qp]*grad_v[j][qp](0);
//                    Kpuy(i,j) += -JxW[qp]*dt*q[i][qp]*grad_v[j][qp](1);
//                    Kpuy(i,j) += -JxW[qp]*dt*q[i][qp]*grad_v[j][qp](2);
//                }
                for(uint j=0; j<n_p_dofs; j++)
                {
                    Kpp(i,j) += JxW[qp]*dt*q[i][qp]*q[j][qp];
                }
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

uint kron_d (uint const i, uint const j)
{
    return i == j;
}

void compute_stresses(EquationSystems& es)
{
  const MeshBase& mesh = es.get_mesh();

  const uint dim = mesh.mesh_dimension();

  LinearImplicitSystem& system = es.get_system<LinearImplicitSystem>("Structure");

  std::vector<uint> d_vars(dim);
  d_vars[0] = system.variable_number ("dx");
  d_vars[1] = system.variable_number ("dy");
  if (dim > 2) d_vars[2] = system.variable_number ("dz");

  const uint dx_var = system.variable_number ("dx");

  const DofMap& dof_map = system.get_dof_map();
  FEType fe_type = dof_map.variable_type(dx_var);
  AutoPtr<FEBase> fe (FEBase::build(dim, fe_type));
  QGauss qrule (dim, fe_type.default_quadrature_order());
  fe->attach_quadrature_rule (&qrule);

  const std::vector<Real>& JxW = fe->get_JxW();
  const std::vector<std::vector<RealGradient> >& dphi = fe->get_dphi();

  // Also, get a reference to the ExplicitSystem
  ExplicitSystem& stress = es.get_system<ExplicitSystem>("Stress");
  const DofMap& stress_dof_map = stress.get_dof_map();
  std::vector<std::vector<uint> > sigma_vars( 3, std::vector<uint>(3) );
  sigma_vars[0][0] = stress.variable_number ("s_00");
  sigma_vars[0][1] = stress.variable_number ("s_01");
  sigma_vars[1][0] = stress.variable_number ("s_10");
  sigma_vars[1][1] = stress.variable_number ("s_11");
  if (dim>2)
  {
      sigma_vars[0][2] = stress.variable_number ("s_02");
      sigma_vars[1][2] = stress.variable_number ("s_12");
      sigma_vars[2][0] = stress.variable_number ("s_20");
      sigma_vars[2][1] = stress.variable_number ("s_21");
      sigma_vars[2][2] = stress.variable_number ("s_22");
  }
  uint vonMises_var = stress.variable_number ("vonMises");

  // Storage for the stress dof indices on each element
  std::vector< std::vector<dof_id_type> > dof_indices_var(system.n_vars());
  std::vector<dof_id_type> stress_dof_indices_var;

  // To store the stress tensor on each element
  DenseMatrix<Number> elem_sigma;

  Real mu = es.parameters.get<Real>("mu");
  Real lambda = es.parameters.get<Real>("lambda");

  MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();

  for ( ; el != end_el; ++el)
  {
    const Elem* elem = *el;

    for(uint var=0; var<dim; var++)
    {
      dof_map.dof_indices (elem, dof_indices_var[var], d_vars[var]);
    }

    fe->reinit (elem);

    elem_sigma.resize(3,3);

    for (uint qp=0; qp<qrule.n_points(); qp++)
    {
        for (uint k=0; k<dim; k++)
        {
            const uint n_var_dofs = dof_indices_var[k].size();

            // Get the gradient at this quadrature point
            Gradient grad_dk;
            for(uint l=0; l<n_var_dofs; l++)
            {
                grad_dk.add_scaled(dphi[l][qp], system.current_solution(dof_indices_var[k][l]));
            }

            for (uint j=0; j<dim; j++)
            {
                elem_sigma(k,j) += JxW[qp]*mu*grad_dk(j);
                elem_sigma(j,k) += JxW[qp]*mu*grad_dk(j);
            }
            elem_sigma(k,k) += JxW[qp]*lambda*grad_dk(k);
        }
    }

    // Get the average stresses by dividing by the element volume
    elem_sigma.scale(1./elem->volume());

    // load elem_sigma data into stress_system
    for(uint i=0; i<dim; i++)
      for(uint j=0; j<dim; j++)
      {
        stress_dof_map.dof_indices (elem, stress_dof_indices_var, sigma_vars[i][j]);

        // We are using CONSTANT MONOMIAL basis functions, hence we only need to get
        // one dof index per variable
        dof_id_type dof_index = stress_dof_indices_var[0];

        if( (stress.solution->first_local_index() <= dof_index) &&
            (dof_index < stress.solution->last_local_index()) )
        {
          stress.solution->set(dof_index, elem_sigma(i,j));
        }

      }

    // Also, the von Mises stress
    Number vonMises_value = std::sqrt( 0.5*( pow(elem_sigma(0,0) - elem_sigma(1,1),2.) +
                                             pow(elem_sigma(1,1) - elem_sigma(2,2),2.) +
                                             pow(elem_sigma(2,2) - elem_sigma(0,0),2.) +
                                             6.*(pow(elem_sigma(0,1),2.) + pow(elem_sigma(1,2),2.) + pow(elem_sigma(2,0),2.))
                                           ) );
    stress_dof_map.dof_indices (elem, stress_dof_indices_var, vonMises_var);
    dof_id_type dof_index = stress_dof_indices_var[0];
    if( (stress.solution->first_local_index() <= dof_index) &&
        (dof_index < stress.solution->last_local_index()) )
    {
      stress.solution->set(dof_index, vonMises_value);
    }

  }

  // Should call close and update when we set vector entries directly
  stress.solution->close();
  stress.update();
}
