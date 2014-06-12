#include "FSI.hpp"

#include <iostream>
#include <algorithm>
#include <cmath>

#include <libmesh/getpot.h>
#include <libmesh/mesh.h>
#include <libmesh/mesh_generation.h>
#include <libmesh/exodusII_io.h>
#include <libmesh/gmv_io.h>
#include <libmesh/equation_systems.h>
#include <libmesh/quadrature_gauss.h>
#include <libmesh/dof_map.h>
#include <libmesh/sparse_matrix.h>
#include <libmesh/numeric_vector.h>
#include <libmesh/system.h>
#include <libmesh/linear_implicit_system.h>
#include <libmesh/transient_system.h>
#include <libmesh/zero_function.h>
#include <libmesh/dirichlet_boundaries.h>
#include <libmesh/dense_submatrix.h>
#include <libmesh/dense_subvector.h>

// extended VTK IO
#include "util/extvtkio.hpp"

#include "util/init.hpp"

// Bring in everything from the libMesh namespace
using namespace libMesh;

// Function prototype.  This function will assemble the system
// matrix and right-hand-side.
void assemble_disp (EquationSystems& es,
                    const std::string& system_name);

void assemble_vel (EquationSystems& es,
                   const std::string& system_name);

// The main program.
int main (int argc, char** argv)
{
    // Initialize libMesh.
    LibMeshInit init (argc, argv);

    Mesh mesh(init.comm());

    GetPot param_file("../param.dat");

    std::string mesh_file = param_file ("mesh_file", "structured");

    if (mesh_file == "structured")
    {
        uint const nx = param_file("nx", 15);
        uint const ny = param_file("ny", 1);
        uint const nz = param_file("nz", 1);
        Real const ox = param_file("ox", 0.0);
        Real const oy = param_file("oy", 0.0);
        Real const oz = param_file("oz", 0.0);
        Real const lx = param_file("lx", 3.0);
        Real const ly = param_file("ly", 0.2);
        Real const lz = param_file("lz", 0.2);

        MeshTools::Generation::build_cube (mesh,
                                           nx, ny, nz,
                                           ox, lx,
                                           oy, ly,
                                           oz, lz,
                                           HEX20);
    }
    else
    {
        mesh.read (mesh_file);
    }

    //    XdrIO mesh_io(mesh);
    //    mesh_io.read("one_tri.xda");

    // Print information about the mesh to the screen.
    mesh.print_info();

    // Create an equation systems object.
    EquationSystems equation_systems (mesh);

    // Declare the system and its variables.
    TransientLinearImplicitSystem & system_dx =
            equation_systems.add_system<TransientLinearImplicitSystem> ("dx");

    const uint dx_var = system_dx.add_variable ("dx", SECOND);

    TransientLinearImplicitSystem & system_dy =
            equation_systems.add_system<TransientLinearImplicitSystem> ("dy");

    const uint dy_var = system_dy.add_variable ("dy", SECOND);

    TransientLinearImplicitSystem & system_dz =
            equation_systems.add_system<TransientLinearImplicitSystem> ("dz");

    const uint dz_var = system_dz.add_variable ("dz", SECOND);

    TransientLinearImplicitSystem & system_vel =
            equation_systems.add_system<TransientLinearImplicitSystem> ("vel");

    const uint u_var = system_vel.add_variable ("ux", SECOND);
    const uint v_var = system_vel.add_variable ("uy", SECOND);
    const uint w_var = system_vel.add_variable ("uz", SECOND);
    system_vel.add_variable ("p", FIRST);

    std::set<boundary_id_type> zero_bc;
    zero_bc.insert (4); // left side

    std::vector<uint> vars_dx (1, dx_var);
    std::vector<uint> vars_dy (1, dy_var);
    std::vector<uint> vars_dz (1, dz_var);

    std::vector<uint> vars_vel (3);
    vars_vel[0] = u_var;
    vars_vel[1] = v_var;
    vars_vel[2] = w_var;

    ZeroFunction<Real> zero;

    system_dx .get_dof_map().add_dirichlet_boundary( libMesh::DirichletBoundary( zero_bc, vars_dx,  &zero ) );
    system_dy .get_dof_map().add_dirichlet_boundary( libMesh::DirichletBoundary( zero_bc, vars_dy,  &zero ) );
    system_dz .get_dof_map().add_dirichlet_boundary( libMesh::DirichletBoundary( zero_bc, vars_dz,  &zero ) );
    system_vel.get_dof_map().add_dirichlet_boundary( libMesh::DirichletBoundary( zero_bc, vars_vel, &zero ) );

    // Give the system a pointer to the matrix assembly
    // function.
    system_dx.attach_assemble_function (assemble_disp);
    system_dx.attach_init_function (init_zero);

    system_dy.attach_assemble_function (assemble_disp);
    system_dy.attach_init_function (init_zero);

    system_dz.attach_assemble_function (assemble_disp);
    system_dz.attach_init_function (init_zero);

    system_vel.attach_assemble_function (assemble_vel);
    system_vel.attach_init_function (init_zero);

    std::string const output_file = param_file("output_file", "structure_seg3d.e");

    equation_systems.parameters.set<Real>("t_in") = param_file("t_in", 0.);
    equation_systems.parameters.set<Real>("t_out") = param_file("t_out", 2.0);
    equation_systems.parameters.set<Real>("dt") = param_file("dt", 1.e-2);

    equation_systems.parameters.set<Real>("f_ux") = param_file("f_u", 0.);
    equation_systems.parameters.set<Real>("f_uy") = param_file("f_v", 0.);
    equation_systems.parameters.set<Real>("f_uz") = param_file("f_w", 1.);

    Real const E = param_file("E", 1e8);
    Real const ni = param_file("ni", 0.3);

    equation_systems.parameters.set<Real>("mu") = E / ( 2.*(1.+ni) );
    equation_systems.parameters.set<Real>("lambda") = E*ni / ( (1.+ni)*(1.-2.*ni) );

    equation_systems.parameters.set<uint>("linear solver maximum iterations") = 250;
    equation_systems.parameters.set<Real>("linear solver tolerance") = TOLERANCE;

//    PetscOptionsSetValue("-ksp_monitor_true_residual",PETSC_NULL);

    // Initialize the data structures for the equation system.
    equation_systems.init ();

    // Prints information about the system to the screen.
    equation_systems.print_info();

    system_vel.time = equation_systems.parameters.get<Real>("t_in");
    system_dx. time = equation_systems.parameters.get<Real>("t_in");
    system_dy. time = equation_systems.parameters.get<Real>("t_in");
    system_dz. time = equation_systems.parameters.get<Real>("t_in");
    uint timestep = 0;

#ifdef LIBMESH_HAVE_EXODUS_API
    ExodusII_IO io (mesh);
    io.write_equation_systems (output_file, equation_systems);
    io.append(true);
    io.write_timestep (output_file, equation_systems, 0, system_vel.time);
#endif

    Real dt = equation_systems.parameters.get<Real>("dt");
    Real t_out = equation_systems.parameters.get<Real>("t_out");
    while (system_vel.time + dt < t_out + 1e-12)
    {
        timestep++;
        // Incremenet the time counter, set the time and the
        // time step size as parameters in the EquationSystem.
        system_vel .time += dt;
        system_dx.time += dt;
        system_dy.time += dt;
        system_dz.time += dt;

        equation_systems.parameters.set<Real> ("time") = system_vel.time;
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
        *system_dx .old_local_solution = *system_dx .current_local_solution;
        *system_dy .old_local_solution = *system_dy .current_local_solution;
        *system_dz .old_local_solution = *system_dz .current_local_solution;

        // Assemble & solve the linear system
        equation_systems.get_system("vel").solve();
        equation_systems.get_system("dx").solve();
        equation_systems.get_system("dy").solve();
        equation_systems.get_system("dz").solve();

        // Output evey 1 timesteps to file.
        if ((timestep)%1 == 0)
        {
#ifdef LIBMESH_HAVE_EXODUS_API
            io.write_timestep (output_file, equation_systems, timestep, system_vel.time);
#endif
        }
    }

    // All done.
    return 0;
}

void assemble_disp (EquationSystems& es,
                         const std::string& system_name)
{
    // It is a good idea to make sure we are assembling
    // the proper system.
    if (! ((system_name == "dx") || (system_name == "dy") || (system_name == "dz")))
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

    // Get a reference to the disp system object.
    TransientLinearImplicitSystem & system_vel =
            es.get_system<TransientLinearImplicitSystem> ("vel");

    // Numeric ids corresponding to each variable in the system
    const uint d_var = system.variable_number (system_name);

    std::string vel_name;
    if (system_name == "dx")
        vel_name = "ux";
    else if (system_name == "dy")
        vel_name = "uy";
    else if (system_name == "dz")
        vel_name = "uz";

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

    // Now we will loop over all the elements in the mesh that
    // live on the local processor. We will compute the element
    // matrix and right-hand-side contribution.  In case users later
    // modify this program to include refinement, we will be safe and
    // will only consider the active elements; hence we use a variant of
    // the \p active_elem_iterator.

    Real dt = es.parameters.get<Real>("dt");

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
                Fe(i) += JxW[qp]*( d_old*b[i][qp] + dt*u_old*v[i][qp]);
                for (uint j=0; j<n_d_dofs; j++)
                {
                    Ke(i,j) += JxW[qp]*b[i][qp]*b[j][qp];
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

    // That's it.
    return;
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
    TransientLinearImplicitSystem & system_d =
            es.get_system<TransientLinearImplicitSystem> ("dx");
    TransientLinearImplicitSystem & system_e =
            es.get_system<TransientLinearImplicitSystem> ("dy");
    TransientLinearImplicitSystem & system_f =
            es.get_system<TransientLinearImplicitSystem> ("dz");

    // Numeric ids corresponding to each variable in the system
    const uint u_var = system.variable_number ("ux");
    const uint v_var = system.variable_number ("uy");
    const uint w_var = system.variable_number ("uz");
    const uint p_var = system.variable_number ("p");

    const uint d_var = system_d.variable_number ("dx");
    const uint e_var = system_e.variable_number ("dy");
    const uint f_var = system_f.variable_number ("dz");

    // Get the Finite Element type for "u".  Note this will be
    // the same as the type for "v".
    FEType fe_u_type = system.variable_type(u_var);
    FEType fe_p_type = system.variable_type(p_var);
    FEType fe_d_type = system_d.variable_type(d_var);

    // Build a Finite Element object of the specified type for
    // the velocity variables.
    AutoPtr<FEBase> fe_u (FEBase::build(dim, fe_u_type));
    AutoPtr<FEBase> fe_p (FEBase::build(dim, fe_p_type));
    AutoPtr<FEBase> fe_d (FEBase::build(dim, fe_d_type));

    // A Gauss quadrature rule for numerical integration.
    // Let the \p FEType object decide what order rule is appropriate.
    QGauss qrule (dim, fe_u_type.default_quadrature_order());

    // Tell the finite element objects to use our quadrature rule.
    fe_u->attach_quadrature_rule (&qrule);
    fe_p->attach_quadrature_rule (&qrule);
    fe_d->attach_quadrature_rule (&qrule);

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
    const std::vector<std::vector<RealGradient> >& grad_b = fe_d->get_dphi();

    // A reference to the \p DofMap object for this system.  The \p DofMap
    // object handles the index translation from node and element numbers
    // to degree of freedom numbers.  We will talk more about the \p DofMap
    // in future examples.
    const DofMap & dof_map = system.get_dof_map();
    const DofMap & dof_map_d = system_d.get_dof_map();
    const DofMap & dof_map_e = system_e.get_dof_map();
    const DofMap & dof_map_f = system_f.get_dof_map();

    // Define data structures to contain the element matrix
    // and right-hand-side vector contribution.  Following
    // basic finite element terminology we will denote these
    // "Ke" and "Fe".
    DenseMatrix<Number> Ke;
    DenseVector<Number> Fe;

    DenseSubMatrix<Number>
            Kuu(Ke), Kuv(Ke), Kuw(Ke), Kup(Ke),
            Kvu(Ke), Kvv(Ke), Kvw(Ke), Kvp(Ke),
            Kwu(Ke), Kwv(Ke), Kww(Ke), Kwp(Ke),
            Kpu(Ke), Kpv(Ke), Kpw(Ke), Kpp(Ke);

    DenseSubVector<Number>
            Fu(Fe),
            Fv(Fe),
            Fw(Fe),
            Fp(Fe);

    // This vector will hold the degree of freedom indices for
    // the element.  These define where in the global system
    // the element degrees of freedom get mapped.
    std::vector<dof_id_type> dof_indices;
    std::vector<dof_id_type> dof_indices_u;
    std::vector<dof_id_type> dof_indices_v;
    std::vector<dof_id_type> dof_indices_w;
    std::vector<dof_id_type> dof_indices_p;
    std::vector<dof_id_type> dof_indices_d;
    std::vector<dof_id_type> dof_indices_e;
    std::vector<dof_id_type> dof_indices_f;

    // Now we will loop over all the elements in the mesh that
    // live on the local processor. We will compute the element
    // matrix and right-hand-side contribution.  In case users later
    // modify this program to include refinement, we will be safe and
    // will only consider the active elements; hence we use a variant of
    // the \p active_elem_iterator.

    Real dt = es.parameters.get<Real>("dt");
    Real f_u = es.parameters.get<Real>("f_ux");
    Real f_v = es.parameters.get<Real>("f_uy");
    Real f_w = es.parameters.get<Real>("f_uz");
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
        dof_map.dof_indices (elem, dof_indices_u, u_var);
        dof_map.dof_indices (elem, dof_indices_v, v_var);
        dof_map.dof_indices (elem, dof_indices_w, w_var);
        dof_map.dof_indices (elem, dof_indices_p, p_var);
        dof_map_d.dof_indices (elem, dof_indices_d, d_var);
        dof_map_e.dof_indices (elem, dof_indices_e, e_var);
        dof_map_f.dof_indices (elem, dof_indices_f, f_var);

        const uint n_dofs   = dof_indices.size();
        const uint n_u_dofs = dof_indices_u.size();
        const uint n_p_dofs = dof_indices_p.size();
        const uint n_d_dofs = dof_indices_d.size();

        // Compute the element-specific data for the current
        // element.  This involves computing the location of the
        // quadrature points (q_point) and the shape functions
        // (phi, dphi) for the current element.
        fe_u->reinit (elem);
        fe_p->reinit (elem);
        fe_d->reinit (elem);

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

        Kuu.reposition (u_var*n_u_dofs, u_var*n_u_dofs, n_u_dofs, n_u_dofs);
        Kuv.reposition (u_var*n_u_dofs, v_var*n_u_dofs, n_u_dofs, n_u_dofs);
        Kuw.reposition (u_var*n_u_dofs, w_var*n_u_dofs, n_u_dofs, n_u_dofs);
        Kup.reposition (u_var*n_u_dofs, p_var*n_u_dofs, n_u_dofs, n_p_dofs);

        Kvu.reposition (v_var*n_u_dofs, u_var*n_u_dofs, n_u_dofs, n_u_dofs);
        Kvv.reposition (v_var*n_u_dofs, v_var*n_u_dofs, n_u_dofs, n_u_dofs);
        Kvw.reposition (v_var*n_u_dofs, w_var*n_u_dofs, n_u_dofs, n_u_dofs);
        Kvp.reposition (v_var*n_u_dofs, p_var*n_u_dofs, n_u_dofs, n_p_dofs);

        Kwu.reposition (w_var*n_u_dofs, u_var*n_u_dofs, n_u_dofs, n_u_dofs);
        Kwv.reposition (w_var*n_u_dofs, v_var*n_u_dofs, n_u_dofs, n_u_dofs);
        Kww.reposition (w_var*n_u_dofs, w_var*n_u_dofs, n_u_dofs, n_u_dofs);
        Kwp.reposition (w_var*n_u_dofs, p_var*n_u_dofs, n_u_dofs, n_p_dofs);

        Kpu.reposition ( p_var*n_u_dofs, u_var*n_u_dofs, n_p_dofs, n_u_dofs);
        Kpv.reposition ( p_var*n_u_dofs, v_var*n_u_dofs, n_p_dofs, n_u_dofs);
        Kpw.reposition ( p_var*n_u_dofs, w_var*n_u_dofs, n_p_dofs, n_u_dofs);
        Kpp.reposition ( p_var*n_u_dofs, p_var*n_u_dofs, n_p_dofs, n_p_dofs);

        Fu.reposition (u_var*n_u_dofs, n_u_dofs);
        Fv.reposition (v_var*n_u_dofs, n_u_dofs);
        Fw.reposition (w_var*n_u_dofs, n_u_dofs);
        Fp.reposition (p_var*n_u_dofs, n_p_dofs);

        // Now we will build the element matrix.
        for (uint qp=0; qp<qrule.n_points(); qp++)
        {
            // compute old solution
            Gradient grad_d_old;
            Gradient grad_e_old;
            Gradient grad_f_old;
            for (uint l=0; l<n_d_dofs; l++)
            {
                grad_d_old.add_scaled (grad_b[l][qp],system_d.old_solution (dof_indices_d[l]));
                grad_e_old.add_scaled (grad_b[l][qp],system_e.old_solution (dof_indices_e[l]));
                grad_f_old.add_scaled (grad_b[l][qp],system_f.old_solution (dof_indices_f[l]));
            }

            Number u_old = 0.;
            Number v_old = 0.;
            Number w_old = 0.;
            for (uint l=0; l<n_u_dofs; l++)
            {
                u_old += phi[l][qp]*system.old_solution (dof_indices_u[l]);
                v_old += phi[l][qp]*system.old_solution (dof_indices_v[l]);
                w_old += phi[l][qp]*system.old_solution (dof_indices_w[l]);
            }

            for (uint i=0; i<n_u_dofs; i++)
            {
                Fu(i) += JxW[qp]*( (u_old+dt*f_u)*phi[i][qp]
                                   - dt*(
                                       mu*(
                                           grad_phi[i][qp]*grad_d_old        // grad(d_old) : grad(v)
                                           +grad_phi[i][qp](0)*grad_d_old(0) // grad(d_old)^T : grad(v)
                                           +grad_phi[i][qp](1)*grad_e_old(0) // |
                                           +grad_phi[i][qp](2)*grad_f_old(0) // |
                                          )
                                       + lambda*(grad_phi[i][qp](0)*grad_d_old(0)
                                                 +grad_phi[i][qp](0)*grad_e_old(1)
                                                 +grad_phi[i][qp](0)*grad_f_old(2))) // stress2 tr(\eps(d_old))*tr(\eps(v))
                                         );
                for (uint j=0; j<n_u_dofs; j++)
                {
                    Kuu(i,j) += JxW[qp]*( phi[i][qp]*phi[j][qp] +
                                       dt*dt*(
                                       mu*(
                                           grad_phi[i][qp]*grad_phi[j][qp]        // grad(dt*u) : grad(v)
                                           + grad_phi[i][qp](0)*grad_phi[j][qp](0) // grad(dt*u)^T : grad(v)
                                       )
                                       + lambda*grad_phi[i][qp](0)*grad_phi[j][qp](0) // stress2 tr(\eps(dt*u))*tr(\eps(v))
                                       )
                                       );
                    Kuv(i,j) += JxW[qp]*dt*dt*(
                                               mu*grad_phi[i][qp](1)*grad_phi[j][qp](0) // stress1 \eps(d):\eps(v)
                                               + lambda*grad_phi[i][qp](0)*grad_phi[j][qp](1) // stress2 tr(\eps(d))*tr(\eps(v))
                                );
                    Kuw(i,j) += JxW[qp]*dt*dt*(
                                               mu*grad_phi[i][qp](2)*grad_phi[j][qp](0) // stress1 \eps(d):\eps(v)
                                               + lambda*grad_phi[i][qp](0)*grad_phi[j][qp](2) // stress2 tr(\eps(d))*tr(\eps(v))
                                );
                }
//                for (uint j=0; j<n_p_dofs; j++ )
//                {
//                    Kup(i,j) += -JxW[qp]*dt*psi[j][qp]*grad_phi[i][qp](1);
//                }

                Fv(i) += JxW[qp]*( (v_old+dt*f_v)*phi[i][qp]
                                  - dt*(
                                       mu*(
                                           grad_phi[i][qp]*grad_e_old
                                           + grad_phi[i][qp](0)*grad_d_old(1)
                                           + grad_phi[i][qp](1)*grad_e_old(1)
                                           + grad_phi[i][qp](2)*grad_f_old(1)
                                          )
                                       + lambda*(
                                               grad_phi[i][qp](1)*grad_d_old(0)
                                               +grad_phi[i][qp](1)*grad_e_old(1)
                                               +grad_phi[i][qp](1)*grad_f_old(2)
                                               )
                                               )
                                    );
                for (uint j=0; j<n_u_dofs; j++)
                {
                    Kvu(i,j) += JxW[qp]*dt*dt*(
                        mu*grad_phi[i][qp](0)*grad_phi[j][qp](1) // stress1 \eps(d):\eps(v)
                        + lambda*grad_phi[i][qp](1)*grad_phi[j][qp](0) // stress2 tr(\eps(d))*tr(\eps(v))
                                       );
                    Kvv(i,j) += JxW[qp]*( phi[i][qp]*phi[j][qp] +
                                       dt*dt*(
                                              mu*(
                                                  grad_phi[i][qp]*grad_phi[j][qp]  // stress1 \eps(d):\eps(v)
                                                  + grad_phi[i][qp](1)*grad_phi[j][qp](1)
                                       )
                                    + lambda*grad_phi[i][qp](1)*grad_phi[j][qp](1) // stress2 tr(\eps(d))*tr(\eps(v))
                                    )
                                  );
                    Kvw(i,j) += JxW[qp]*dt*dt*(
                                           mu*grad_phi[i][qp](2)*grad_phi[j][qp](1) // stress1 \eps(d):\eps(v)
                                           + lambda*grad_phi[i][qp](1)*grad_phi[j][qp](2) // stress2 tr(\eps(d))*tr(\eps(v))
                                                          );
                }

                Fw(i) += JxW[qp]*( (w_old+dt*f_w)*phi[i][qp]
                                   - dt*(
                                       mu*(
                                           grad_phi[i][qp]*grad_f_old
                                           + grad_phi[i][qp](0)*grad_d_old(2)
                                           + grad_phi[i][qp](1)*grad_e_old(2)
                                           + grad_phi[i][qp](2)*grad_f_old(2)
                                           )
                                        + lambda*(
                                               grad_phi[i][qp](2)*grad_d_old(0)
                                               +grad_phi[i][qp](2)*grad_e_old(1)
                                               +grad_phi[i][qp](2)*grad_f_old(2)
                                               )
                                               )
                                               );
                for (uint j=0; j<n_u_dofs; j++)
                {
                    Kwu(i,j) += JxW[qp]*dt*dt*(
                                       mu*grad_phi[i][qp](0)*grad_phi[j][qp](2) // stress1 \eps(d):\eps(v)
                                       + lambda*grad_phi[i][qp](2)*grad_phi[j][qp](0) // stress2 tr(\eps(d))*tr(\eps(v))
                                       );
                    Kww(i,j) += JxW[qp]*( phi[i][qp]*phi[j][qp] +
                                       dt*dt*(
                                       mu*(
                                       grad_phi[i][qp]*grad_phi[j][qp]  // stress1 \eps(d):\eps(v)
                                       + grad_phi[i][qp](2)*grad_phi[j][qp](2)
                                       )
                                       + lambda*grad_phi[i][qp](2)*grad_phi[j][qp](2) // stress2 tr(\eps(d))*tr(\eps(v))
                                       )
                                       );
                    Kwv(i,j) += JxW[qp]*dt*dt*(
                                       mu*grad_phi[i][qp](1)*grad_phi[j][qp](2) // stress1 \eps(d):\eps(v)
                                       + lambda*grad_phi[i][qp](2)*grad_phi[j][qp](1) // stress2 tr(\eps(d))*tr(\eps(v))
                                       );
                }
//                for (uint j=0; j<n_p_dofs; j++)
//                {
//                    Kvp(i,j) += -JxW[qp]*dt*psi[j][qp]*grad_phi[i][qp](1);
//                }
            }
            for (uint i=0; i<n_p_dofs; i++)
            {
//                for (uint j=0; j<n_u_dofs; j++)
//                {
//                    Kpu(i,j) += -JxW[qp]*dt*psi[i][qp]*grad_phi[j][qp](0);
//                    Kpv(i,j) += -JxW[qp]*dt*psi[i][qp]*grad_phi[j][qp](1);
//                }
                for(uint j=0; j<n_p_dofs; j++)
                {
                    Kpp(i,j) += JxW[qp]*dt*psi[i][qp]*psi[j][qp];
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
