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

// Function prototype.  This function will assemble the system
// matrix and right-hand-side.
void assemble_system (EquationSystems& es,
                      const std::string& system_name);

Real f( Point const & p )
{
    return std::exp(p(0));
}

// The main program.
int main (int argc, char** argv)
{
    LibMeshInit init (argc, argv);

    Mesh mesh(init.comm());

    MeshTools::Generation::build_square (mesh,
                                         4, 4,
                                         0.0, 1.0,
                                         0.0, 1.0,
                                         QUAD4);

//    XdrIO mesh_io(mesh);
//    mesh_io.read("one_tri.xda");

    mesh.print_info();

    EquationSystems equation_systems (mesh);

    LinearImplicitSystem & system =
            equation_systems.add_system<LinearImplicitSystem> ("system");

    /*const uint u_var =*/ system.add_variable ("u", CONSTANT, MONOMIAL);

//    std::set<boundary_id_type> dirichlet_bc;
//    dirichlet_bc.insert(3); // left side

//    std::vector<uint> vars;
//    vars.push_back (u_var);

//    ZeroFunction<Real> zero;

//    system.get_dof_map().add_dirichlet_boundary( libMesh::DirichletBoundary( dirichlet_bc, vars, &zero ) );

    system.attach_assemble_function (assemble_system);

    // Initialize the data structures for the equation system.
    equation_systems.init ();

    equation_systems.parameters.set<uint>("linear solver maximum iterations") = 250;
    equation_systems.parameters.set<Real>("linear solver tolerance") = TOLERANCE;

    // Prints information about the system to the screen.
    equation_systems.print_info();

    // Assemble & solve the linear system,
    // then write the solution.
    equation_systems.get_system("system").solve();

#ifdef LIBMESH_HAVE_EXODUS_API
    ExodusII_IO(mesh).write_equation_systems ("monomial2d.e",
                                              equation_systems);
#endif // #ifdef LIBMESH_HAVE_EXODUS_API

    // All done.
    return 0;
}

void assemble_system (EquationSystems& es,
                         const std::string& system_name)
{
    libmesh_assert_equal_to (system_name, "system");

    const MeshBase& mesh = es.get_mesh();

    const uint dim = mesh.mesh_dimension();

    LinearImplicitSystem & system =
            es.get_system<LinearImplicitSystem> ("system");

    const uint u_var = system.variable_number ("u");

    FEType fe_type = system.variable_type(u_var);

    AutoPtr<FEBase> fe (FEBase::build(dim, fe_type));

    QGauss qrule (dim, fe_type.default_quadrature_order());

    // Tell the finite element objects to use our quadrature rule.
    fe->attach_quadrature_rule (&qrule);

    const std::vector<Real>& JxW = fe->get_JxW();

    //const std::vector<std::vector<RealGradient> >& grad_v = fe->get_dphi();
    const std::vector<std::vector<Real> >& v = fe->get_phi();

    const DofMap & dof_map = system.get_dof_map();

    DenseMatrix<Number> Ke;
    DenseVector<Number> Fe;

    std::vector<dof_id_type> dof_indices;

    MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
    const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();

    for ( ; el != end_el; ++el)
    {
        const Elem* elem = *el;

        dof_map.dof_indices (elem, dof_indices);

        const uint n_dofs   = dof_indices.size();

        fe->reinit (elem);

        Ke.resize (n_dofs, n_dofs);
        Fe.resize (n_dofs);

        for (uint qp=0; qp<qrule.n_points(); qp++)
        {
            for (uint i=0; i<n_dofs; i++)
            {
                Real const fp = f(fe->get_xyz()[qp]);
                Fe(i) += JxW[qp]*fp*v[i][qp];
                for (uint j=0; j<n_dofs; j++)
                {
                    Ke(i,j) += JxW[qp]*v[i][qp]*v[j][qp];
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
