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
#include "libmesh/equation_systems.h"
#include "libmesh/fe.h"
#include "libmesh/quadrature_gauss.h"
#include "libmesh/dof_map.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/dense_matrix.h"
#include "libmesh/dense_vector.h"
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

// Bring in everything from the libMesh namespace
using namespace libMesh;

Real const f_x = 0.;
Real const f_y = 0.;
Real const f_z = 1.;

Real const E = 1e8;
Real const ni = 0.3;

Real const mu = E / ( 2.*(1.+ni) );
Real const lmbda = E*ni / ( (1.+ni)*(1.-2.*ni) );

// Function prototype.  This function will assemble the system
// matrix and right-hand-side.
void assemble_structure (EquationSystems& es,
                         const std::string& system_name);

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
    LinearImplicitSystem & system =
            equation_systems.add_system<LinearImplicitSystem> ("Structure");

    // Add the variables "u" & "v" to "Stokes".  They
    // will be approximated using second-order approximation.
    const uint dx_var = system.add_variable ("dx", SECOND);
    const uint dy_var = system.add_variable ("dy", SECOND);
    const uint dz_var = system.add_variable ("dz", SECOND);

    std::set<boundary_id_type> dirichlet_bc;
    dirichlet_bc.insert(4); // left side

    std::vector<uint> vars (3);
    vars[0] = dx_var;
    vars[1] = dy_var;
    vars[2] = dz_var;

    ZeroFunction<Real> zero;

    system.get_dof_map().add_dirichlet_boundary( libMesh::DirichletBoundary( dirichlet_bc, vars, &zero ) );

    // Give the system a pointer to the matrix assembly
    // function.
    system.attach_assemble_function (assemble_structure);

    // Initialize the data structures for the equation system.
    equation_systems.init ();

    equation_systems.parameters.set<uint>("linear solver maximum iterations") = 250;
    equation_systems.parameters.set<Real>        ("linear solver tolerance") = TOLERANCE;

    // Prints information about the system to the screen.
    equation_systems.print_info();

    equation_systems.parameters.print();

    // Assemble & solve the linear system,
    // then write the solution.
    equation_systems.get_system("Structure").solve();

#ifdef LIBMESH_HAVE_EXODUS_API
    ExodusII_IO(mesh).write_equation_systems ("structure_static_3d.e",
                                              equation_systems);
#endif // #ifdef LIBMESH_HAVE_EXODUS_API

    // All done.
    return 0;
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
    LinearImplicitSystem & system =
            es.get_system<LinearImplicitSystem> ("Structure");

    // Numeric ids corresponding to each variable in the system
    const uint dx_var = system.variable_number ("dx");
    const uint dy_var = system.variable_number ("dy");
    const uint dz_var = system.variable_number ("dz");

    // Get the Finite Element type for "u".  Note this will be
    // the same as the type for "v".
    FEType fe_d_type = system.variable_type(dx_var);

    // Build a Finite Element object of the specified type for
    // the velocity variables.
    AutoPtr<FEBase> fe_d (FEBase::build(dim, fe_d_type));

    // A Gauss quadrature rule for numerical integration.
    // Let the \p FEType object decide what order rule is appropriate.
    QGauss qrule (dim, fe_d_type.default_quadrature_order());

    // Tell the finite element objects to use our quadrature rule.
    fe_d->attach_quadrature_rule (&qrule);

    // Here we define some references to cell-specific data that
    // will be used to assemble the linear system.
    //
    // The element Jacobian * quadrature weight at each integration point.
    const std::vector<Real>& JxW = fe_d->get_JxW();

    // The element shape function gradients for the velocity
    // variables evaluated at the quadrature points.
    const std::vector<std::vector<Real> >& phi = fe_d->get_phi();
    const std::vector<std::vector<RealGradient> >& dphi = fe_d->get_dphi();

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
            Kdxdx(Ke), Kdxdy(Ke), Kdxdz(Ke),
            Kdydx(Ke), Kdydy(Ke), Kdydz(Ke),
            Kdzdx(Ke), Kdzdy(Ke), Kdzdz(Ke);

    DenseSubVector<Number>
            Fdx(Fe),
            Fdy(Fe),
            Fdz(Fe);

    // This vector will hold the degree of freedom indices for
    // the element.  These define where in the global system
    // the element degrees of freedom get mapped.
    std::vector<dof_id_type> dof_indices;
    std::vector<dof_id_type> dof_indices_dx;
    std::vector<dof_id_type> dof_indices_dy;
    std::vector<dof_id_type> dof_indices_dz;

    // Now we will loop over all the elements in the mesh that
    // live on the local processor. We will compute the element
    // matrix and right-hand-side contribution.  In case users later
    // modify this program to include refinement, we will be safe and
    // will only consider the active elements; hence we use a variant of
    // the \p active_elem_iterator.

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

        const uint n_dofs   = dof_indices.size();
        const uint n_dx_dofs = dof_indices_dx.size();
        const uint n_dy_dofs = dof_indices_dy.size();
        const uint n_dz_dofs = dof_indices_dz.size();

        // Compute the element-specific data for the current
        // element.  This involves computing the location of the
        // quadrature points (q_point) and the shape functions
        // (phi, dphi) for the current element.
        fe_d->reinit  (elem);

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
        Kdxdx.reposition (dx_var*n_dx_dofs, dx_var*n_dx_dofs, n_dx_dofs, n_dx_dofs);
        Kdxdy.reposition (dx_var*n_dx_dofs, dy_var*n_dx_dofs, n_dx_dofs, n_dy_dofs);
        Kdxdz.reposition (dx_var*n_dx_dofs, dz_var*n_dx_dofs, n_dx_dofs, n_dz_dofs);

        Kdydx.reposition (dy_var*n_dy_dofs, dx_var*n_dy_dofs, n_dy_dofs, n_dx_dofs);
        Kdydy.reposition (dy_var*n_dy_dofs, dy_var*n_dy_dofs, n_dy_dofs, n_dy_dofs);
        Kdydz.reposition (dy_var*n_dy_dofs, dz_var*n_dy_dofs, n_dy_dofs, n_dz_dofs);

        Kdzdx.reposition (dz_var*n_dy_dofs, dx_var*n_dz_dofs, n_dz_dofs, n_dx_dofs);
        Kdzdy.reposition (dz_var*n_dy_dofs, dy_var*n_dz_dofs, n_dz_dofs, n_dy_dofs);
        Kdzdz.reposition (dz_var*n_dy_dofs, dz_var*n_dz_dofs, n_dz_dofs, n_dz_dofs);

        Fdx.reposition (dx_var*n_dx_dofs, n_dx_dofs);
        Fdy.reposition (dy_var*n_dx_dofs, n_dy_dofs);
        Fdz.reposition (dz_var*n_dx_dofs, n_dz_dofs);

        // Now we will build the element matrix.
        for (uint qp=0; qp<qrule.n_points(); qp++)
        {
            for (uint i=0; i<n_dx_dofs; i++)
            {
                Fdx(i) += JxW[qp]*f_x*phi[i][qp];
                for (uint j=0; j<n_dx_dofs; j++)
                {
                    Kdxdx(i,j) += JxW[qp]*(
                                mu*( 2.*dphi[i][qp](0)*dphi[j][qp](0)
                                     + dphi[i][qp](1)*dphi[j][qp](1)
                                     + dphi[i][qp](2)*dphi[j][qp](2) )
                                     + lmbda*dphi[i][qp](0)*dphi[j][qp](0)
                                     );
                }
                for (uint j=0; j<n_dy_dofs; j++)
                {
                    Kdxdy(i,j) += JxW[qp]*(
                                mu*dphi[i][qp](1)*dphi[j][qp](0)
                                + lmbda*dphi[i][qp](0)*dphi[j][qp](1)
                                );
                }
                for (uint j=0; j<n_dz_dofs; j++)
                {
                    Kdxdz(i,j) += JxW[qp]*(
                                mu*dphi[i][qp](2)*dphi[j][qp](0)
                                + lmbda*dphi[i][qp](0)*dphi[j][qp](2)
                                );
                }
            }

            for (uint i=0; i<n_dy_dofs; i++)
            {
                Fdy(i) += JxW[qp]*f_y*phi[i][qp];
                for (uint j=0; j<n_dx_dofs; j++)
                {
                    Kdydx(i,j) += JxW[qp]*(
                                mu*dphi[i][qp](0)*dphi[j][qp](1)
                                + lmbda*dphi[i][qp](1)*dphi[j][qp](0)
                                );
                }
                for (uint j=0; j<n_dy_dofs; j++)
                {
                    Kdydy(i,j) += JxW[qp]*(
                                mu*( dphi[i][qp](0)*dphi[j][qp](0)
                                     + 2.*dphi[i][qp](1)*dphi[j][qp](1)
                                     + dphi[i][qp](2)*dphi[j][qp](2) )
                                     + lmbda*dphi[i][qp](1)*dphi[j][qp](1)
                                     );
                }
                for (uint j=0; j<n_dz_dofs; j++)
                {
                    Kdydz(i,j) += JxW[qp]*(
                                mu*dphi[i][qp](2)*dphi[j][qp](1)
                                + lmbda*dphi[i][qp](1)*dphi[j][qp](2)
                                );
                }
            }

            for (uint i=0; i<n_dz_dofs; i++)
            {
                Fdz(i) += JxW[qp]*f_z*phi[i][qp];
                for (uint j=0; j<n_dx_dofs; j++)
                {
                    Kdzdx(i,j) += JxW[qp]*(
                                mu*dphi[i][qp](0)*dphi[j][qp](2)
                                + lmbda*dphi[i][qp](2)*dphi[j][qp](0)
                                );
                }
                for (uint j=0; j<n_dy_dofs; j++)
                {
                    Kdzdy(i,j) += JxW[qp]*(
                                mu*dphi[i][qp](1)*dphi[j][qp](2)
                                + lmbda*dphi[i][qp](2)*dphi[j][qp](1)
                                );
                }
                for (uint j=0; j<n_dz_dofs; j++)
                {
                    Kdzdz(i,j) += JxW[qp]*(
                                mu*( dphi[i][qp](0)*dphi[j][qp](0)
                                     + dphi[i][qp](1)*dphi[j][qp](1)
                                     + 2.*dphi[i][qp](2)*dphi[j][qp](2) )
                                     + lmbda*dphi[i][qp](2)*dphi[j][qp](2)
                                     );
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
//    system.matrix->print();

//    system.rhs->close();
//    system.rhs->print();

    // That's it.
    return;
}
