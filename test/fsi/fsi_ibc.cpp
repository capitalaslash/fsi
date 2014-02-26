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
#include "libmesh/xdr_io.h"
#include "libmesh/gmsh_io.h"
#include "libmesh/gmv_io.h"
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
#include <libmesh/boundary_mesh.h>

// For systems of equations the \p DenseSubMatrix
// and \p DenseSubVector provide convenient ways for
// assembling the element matrix and vector on a
// component-by-component basis.
#include "libmesh/dense_submatrix.h"
#include "libmesh/dense_subvector.h"

// The definition of a geometric element
#include "libmesh/elem.h"

// extended VTK IO
#include "util/extvtkio.hpp"

#include "assembly/mat.hpp"
#include "assembly/interface.hpp"
#include "util/biquad.hpp"
#include "util/init.hpp"
#include "util/stress.hpp"
#include "bc/ext_pressure.hpp"
#include "assemble_disp.hpp"

// Bring in everything from the libMesh namespace
using namespace libMesh;

// Function prototype.  This function will assemble the system
// matrix and right-hand-side.
void assemble_fsi (EquationSystems& es,
                   const std::string& system_name);

//Real external_pressure( Point const & /*point*/, Parameters const & /*param*/ )
//{
//    return 1.;
//}

// The main program.
int main (int argc, char** argv)
{
    LibMeshInit init (argc, argv);

    GetPot param_file("param_fsi_ibc.dat");

    Mesh mesh(init.comm(), 2);

    std::string mesh_file = param_file("mesh_file", "structured");

    if( mesh_file == "structured" )
    {
        generate_biquad( mesh, param_file );
    }
    else
    {
        mesh.read(mesh_file);
    }

    mesh.print_info();

//    BoundaryMesh mesh_b(mesh.comm(), 1);
//    mesh.boundary_info->sync(mesh_b);
//    EquationSystems es_b(mesh_b);

//    BInfoSystem & system_b = es_b.add_system<BInfoSystem>("binfo");

    EquationSystems es (mesh);

    // Declare the system and its variables.
    TransientLinearImplicitSystem & system_dx =
            es.add_system<TransientLinearImplicitSystem> ("dx");

    const uint dx_var = system_dx.add_variable ("dx", SECOND);

    TransientLinearImplicitSystem & system_dy =
            es.add_system<TransientLinearImplicitSystem> ("dy");

    const uint dy_var = system_dy.add_variable ("dy", SECOND);

    TransientLinearImplicitSystem & system_vel =
            es.add_system<TransientLinearImplicitSystem> ("fsi");

    const uint u_var = system_vel.add_variable ("ux", SECOND);
    const uint v_var = system_vel.add_variable ("uy", SECOND);
    system_vel.add_variable ("p", FIRST);

    ExplicitSystem & system_mat = es.add_system<ExplicitSystem>("mat");
    system_mat.add_variable("mat", CONSTANT, MONOMIAL);

    InterfaceSystem & interface = es.add_system<InterfaceSystem>("int");

    ExplicitSystem & stress = es.add_system<ExplicitSystem> ("stress");
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

    ExplicitSystem & ext_pressure = es.add_system<ExplicitSystem>("ext_pressure");
    ext_pressure.add_variable("p0", FIRST);

    system_dy .get_dof_map().add_dirichlet_boundary( DirichletBoundary( {1,5}, {dy_var}, ZeroFunction<Real>() ) );
    system_vel.get_dof_map().add_dirichlet_boundary( DirichletBoundary( {1,5}, {v_var},  ZeroFunction<Real>() ) );

    system_dx .get_dof_map().add_dirichlet_boundary( DirichletBoundary( {3,6}, {dx_var}, ZeroFunction<Real>() ) );
    system_vel.get_dof_map().add_dirichlet_boundary( DirichletBoundary( {3}, {u_var},  ZeroFunction<Real>() ) );

//    system_dy .get_dof_map().add_dirichlet_boundary( DirichletBoundary( {3}, {dy_var}, ZeroFunction<Real>() ) );
//    system_vel.get_dof_map().add_dirichlet_boundary( DirichletBoundary( {3}, {v_var},  ZeroFunction<Real>() ) );

    system_dx.attach_assemble_function (assemble_disp);
    system_dx.attach_init_function (init_zero);

    system_dy.attach_assemble_function (assemble_disp);
    system_dy.attach_init_function (init_zero);

    system_vel.attach_assemble_function (assemble_fsi);
    system_vel.attach_init_function (init_zero);

    es.parameters.set<subdomain_id_type>("flag_s") = param_file("flag_s", 12);
    es.parameters.set<subdomain_id_type>("flag_f") = param_file("flag_f", 11);
    es.parameters.set<boundary_id_type>("flag_int") = param_file("flag_int", 10);

    es.parameters.set<Real>("t_in") = param_file("t_in", 0.);
    es.parameters.set<Real>("t_out") = param_file("t_out", 0.2);
    es.parameters.set<Real>("dt") = param_file("dt", 1.e-3);

    es.parameters.set<Real>("f_ux") = param_file("f_u", 0.);
    es.parameters.set<Real>("f_uy") = param_file("f_v", 0.);

    es.parameters.set<bool>("axisym") = param_file("axisym", false);

    Real const E = param_file("E", 1e8);
    Real const ni = param_file("ni", 0.3);

    es.parameters.set<Real>("mu_s") = E / ( 2.*(1.+ni) );
    es.parameters.set<Real>("lambda") = E*ni / ( (1.+ni)*(1.-2.*ni) );
    es.parameters.set<Real>("rho_s") = param_file("rho_s", 1.0);
    es.parameters.set<Real>("damp") = param_file("damp", 1e-1);
    es.parameters.set<Real>("mu_f") = param_file("mu_f", 1e-1);
    es.parameters.set<Real>("rho_f") = param_file("rho_f", 1.0);

    es.parameters.set<Real>("ale_factor") = param_file ("ale_factor", 1e7);

    es.parameters.set<uint>("linear solver maximum iterations") = 250;
    es.parameters.set<Real>("linear solver tolerance") = TOLERANCE;

    es.parameters.set<Real>("pressure_max") = param_file("pressure_max", 1.0);
    es.parameters.set<Real>("pressure_tin") = param_file("pressure_tin", 0.0);
    es.parameters.set<Real>("pressure_tout") = param_file("pressure_tout", 1.0);
    es.parameters.set<Real>("pressure_sigma") = param_file("pressure_sigma", 1.0);

    es.parameters.set<std::string>("output_dir") = param_file("output_dir", "output_ibc/");
    es.parameters.set<std::string>("basename") = param_file("basename", "fsiibc");
    es.parameters.set<uint>("print_step") = param_file("print_step", 1);

    // PetscOptionsSetValue("-ksp_monitor_true_residual",PETSC_NULL);

    // Initialize the data structures for the equation system.
    es.init ();

    assemble_mat( es, "mat" );

    interface.assemble();

    // Prints information about the system to the screen.
    es.print_info();

    system_vel.time = es.parameters.get<Real>("t_in");
    system_dx. time = es.parameters.get<Real>("t_in");
    system_dy. time = es.parameters.get<Real>("t_in");
    uint timestep = 0;

    AutoPtr<ExtVTKIO> io_vtk = AutoPtr<ExtVTKIO>(new ExtVTKIO(mesh,es.parameters));

    AutoPtr<GMVIO> io_gmv = AutoPtr<GMVIO>( new GMVIO(mesh) );

    es.parameters.set<uint>("timestep") = timestep;

    io_vtk->write_solution(es);

    Real const dt = es.parameters.get<Real>("dt");
    Real const t_out = es.parameters.get<Real>("t_out");
    uint const print_step = es.parameters.get<uint>("print_step");
    while (system_vel.time + dt < t_out + 1e-12)
    {
        timestep++;
        es.parameters.set<uint>("timestep") = timestep;
        // Incremenet the time counter, set the time and the
        // time step size as parameters in the EquationSystem.
        system_vel.time += dt;
        system_dx .time += dt;
        system_dy .time += dt;
        ext_pressure.time += dt;

        es.parameters.set<Real> ("time") = system_vel.time;
        es.parameters.set<Real> ("dt")   = dt;

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
        *system_dx .old_local_solution = *system_dx .current_local_solution;
        *system_dy .old_local_solution = *system_dy .current_local_solution;

        // Assemble & solve the linear system
        es.get_system("fsi").solve();
        es.get_system("dx").solve();
        es.get_system("dy").solve();

        // move_mesh( es );

        // Output evey 1 timesteps to file.
        if ((timestep)%print_step == 0)
        {
            // Post-process the solution to compute the stresses
            compute_stress(es);
            compute_ext_pressure(es);

            io_vtk->write_solution(es);
        }
    }

    io_gmv->write_equation_systems("output.gmv", es);

    // All done.
    return 0;
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
    TransientLinearImplicitSystem & system_d =
            es.get_system<TransientLinearImplicitSystem> ("dx");
    TransientLinearImplicitSystem & system_e =
            es.get_system<TransientLinearImplicitSystem> ("dy");
    InterfaceSystem & system_i = es.get_system<InterfaceSystem>("int");

    // Numeric ids corresponding to each variable in the system
    const uint u_var = system.variable_number ("ux");
    const uint v_var = system.variable_number ("uy");
    const uint p_var = system.variable_number ("p");

    // Get the Finite Element type for "u".  Note this will be
    // the same as the type for "v".
    FEType fe_u_type = system.variable_type(u_var);
    FEType fe_p_type = system.variable_type(p_var);
    FEType fe_d_type = system_d.variable_type(0);

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

    AutoPtr<FEBase> fe_face (FEBase::build(dim, fe_u_type));

    QGauss qface(dim-1, fe_u_type.default_quadrature_order());

    fe_face->attach_quadrature_rule (&qface);

    // Here we define some references to cell-specific data that
    // will be used to assemble the linear system.
    //
    // The element Jacobian * quadrature weight at each integration point.
    const std::vector<Real>& JxW = fe_u->get_JxW();

    // The element shape function gradients for the velocity
    // variables evaluated at the quadrature points.
    const std::vector<std::vector<Real> >& phi = fe_u->get_phi();
    const std::vector<std::vector<RealGradient> >& dphi = fe_u->get_dphi();
    const std::vector<std::vector<Real> >& psi = fe_p->get_phi();
    const std::vector<std::vector<Real> >& b = fe_d->get_phi();
    const std::vector<std::vector<RealGradient> >& grad_b = fe_d->get_dphi();
    const std::vector<Point >& q_point = fe_u->get_xyz();

    // A reference to the \p DofMap object for this system.  The \p DofMap
    // object handles the index translation from node and element numbers
    // to degree of freedom numbers.  We will talk more about the \p DofMap
    // in future examples.
    const DofMap & dof_map = system.get_dof_map();
    const DofMap & dof_map_d = system_d.get_dof_map();
    const DofMap & dof_map_e = system_e.get_dof_map();

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
    std::vector<dof_id_type> dof_indices_d;
    std::vector<dof_id_type> dof_indices_e;
    std::vector<bool> on_interface;

    // Now we will loop over all the elements in the mesh that
    // live on the local processor. We will compute the element
    // matrix and right-hand-side contribution.  In case users later
    // modify this program to include refinement, we will be safe and
    // will only consider the active elements; hence we use a variant of
    // the \p active_elem_iterator.

    Real const dt = es.parameters.get<Real>("dt");
    Real const f_u = es.parameters.get<Real>("f_ux");
    Real const f_v = es.parameters.get<Real>("f_uy");
    Real const rho_s = es.parameters.get<Real>("rho_s");
    Real const mu_s = es.parameters.get<Real>("mu_s");
    Real const lambda = es.parameters.get<Real>("lambda");
    // Real const ilambda = 1. / es.parameters.get<Real>("lambda");
    Real const rho_f = es.parameters.get<Real>("rho_f");
    Real const mu_f = es.parameters.get<Real>("mu_f");
    Real const damp = es.parameters.get<Real>("damp");
    bool const axisym = es.parameters.get<bool>("axisym");

    subdomain_id_type const flag_s = es.parameters.get<subdomain_id_type>("flag_s");
    subdomain_id_type const flag_f = es.parameters.get<subdomain_id_type>("flag_f");

    ExtPressure p0(es.parameters);

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
        dof_map_d.dof_indices (elem, dof_indices_d, 0);
        dof_map_e.dof_indices (elem, dof_indices_e, 0);
        system_i.on_interface(elem, on_interface);

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
        //       -           -          -  -
        //      | Kuu Kuv Kup |        | Fu |
        // Ke = | Kvu Kvv Kvp |;  Fe = | Fv |
        //      | Kpu Kpv Kpp |        | Fp |
        //       -           -          -  -
        //
        // The \p DenseSubMatrix.repostition () member takes the
        // (row_offset, column_offset, row_size, column_size).
        //
        // Similarly, the \p DenseSubVector.reposition () member
        // takes the (row_offset, row_size)

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

        // solid domain
        if (elem->subdomain_id() == flag_s)
        {
            for (uint qp=0; qp<qrule.n_points(); qp++)
            {
                Real const invR = 1./q_point[qp](0);
                Real const invR2 = invR * invR;
                Real const JxWqp = (axisym)? JxW[qp]*q_point[qp](0):JxW[qp];

                // compute old solution
                Number d_old = 0.;
                Number e_old = 0.;
                Gradient grad_d_old;
                Gradient grad_e_old;
                for (uint l=0; l<n_d_dofs; l++)
                {
                    d_old += b[l][qp]*system_d.old_solution (dof_indices_d[l]);
                    e_old += b[l][qp]*system_e.old_solution (dof_indices_e[l]);
                    grad_d_old.add_scaled (grad_b[l][qp],system_d.old_solution (dof_indices_d[l]));
                    grad_e_old.add_scaled (grad_b[l][qp],system_e.old_solution (dof_indices_e[l]));
                }

                Number u_old = 0.;
                Number v_old = 0.;
                for (uint l=0; l<n_u_dofs; l++)
                {
                    u_old += phi[l][qp]*system.old_solution (dof_indices_u[l]);
                    v_old += phi[l][qp]*system.old_solution (dof_indices_v[l]);
                }

                for (uint i=0; i<n_u_dofs; i++)
                {
                    Fu(i) += JxWqp*(
                                rho_s*(u_old+f_u*dt)*phi[i][qp]
                                - dt*(
                                    mu_s*(
                                        2.*dphi[i][qp](0)*grad_d_old(0)
                                        + dphi[i][qp](1)*grad_d_old(1)
                                        + dphi[i][qp](1)*grad_e_old(0)
                                        + axisym*2.*invR2 * phi[i][qp]*d_old
                                        )
                                    + lambda*(
                                        dphi[i][qp](0)*grad_d_old(0)
                                        + dphi[i][qp](0)*grad_e_old(1)
                                        +axisym*invR*(
                                            phi[i][qp]*grad_d_old(0)
                                            + dphi[i][qp](0)*d_old
                                            + invR* phi[i][qp]*d_old
                                            + phi[i][qp]*grad_e_old(1)
                                            )
                                        )
                                    )
                                );
                    for (uint j=0; j<n_u_dofs; j++)
                    {
                        Kuu(i,j) += JxWqp*(
                                        rho_s*(1.+damp)*phi[i][qp]*phi[j][qp] // mass matrix u*v/dt
                                        + dt*dt*(
                                            mu_s*(
                                                2.*dphi[i][qp](0)*dphi[j][qp](0)
                                                + dphi[i][qp](1)*dphi[j][qp](1)
                                                + axisym*2.*invR2* phi[i][qp]*phi[j][qp]
                                                )
                                        + lambda*(
                                            dphi[i][qp](0)*dphi[j][qp](0)
                                            + axisym*invR*(
                                                phi[i][qp]*dphi[j][qp](0)
                                                + dphi[i][qp](0)*phi[j][qp]
                                                + invR* phi[i][qp]*phi[j][qp]
                                                )
                                            )
                                        )
                                    );
                        Kuv(i,j) += JxWqp*dt*dt*(
                                        mu_s*(
                                            dphi[i][qp](1)*dphi[j][qp](0)
                                            )
                                        + lambda*(
                                            dphi[i][qp](0)*dphi[j][qp](1)
                                            + axisym*invR* phi[i][qp]*dphi[j][qp](1)
                                            )
                                        );
                    }
                    // for (uint j=0; j<n_p_dofs; j++ )
                    // {
                    //     Kup(i,j) += -JxW[qp]*dt*psi[j][qp]*dphi[i][qp](0);
                    // }

                    Fv(i) += JxWqp*(
                                rho_s*(v_old+f_v*dt)*phi[i][qp]
                                - dt*(
                                    mu_s*(
                                        dphi[i][qp](0)*grad_d_old(1)
                                        + dphi[i][qp](0)*grad_e_old(0)
                                        + 2.*dphi[i][qp](1)*grad_e_old(1)
                                        )
                                    + lambda*(
                                        dphi[i][qp](1)*grad_d_old(0)
                                        + dphi[i][qp](1)*grad_e_old(1)
                                        + axisym*invR* dphi[i][qp](1)*d_old
                                        )
                                    )
                                );
                    for (uint j=0; j<n_u_dofs; j++)
                    {
                        Kvu(i,j) += JxWqp*dt*dt*(
                                    mu_s*(
                                        dphi[i][qp](0)*dphi[j][qp](1)
                                        )
                                    + lambda*(
                                        dphi[i][qp](1)*dphi[j][qp](0)
                                        + axisym*invR* dphi[i][qp](1)*phi[j][qp]
                                        )
                                    );
                        Kvv(i,j) += JxWqp*(
                                        rho_s*(1.+damp)*phi[i][qp]*phi[j][qp]
                                        + dt*dt*(
                                            mu_s*(
                                                dphi[i][qp](0)*dphi[j][qp](0)
                                                + 2.*dphi[i][qp](1)*dphi[j][qp](1)
                                                )
                                            + lambda*(
                                                dphi[i][qp](1)*dphi[j][qp](1)
                                                )
                                            )
                                        );
                    }
                    // for (uint j=0; j<n_p_dofs; j++)
                    // {
                    //     Kvp(i,j) += -JxW[qp]*dt*psi[j][qp]*dphi[i][qp](1);
                    // }
                }

                for (uint i=0; i<n_p_dofs; i++)
                {
                    if(!on_interface[i])
                    {
                        Kpp(i,i) = 1.;
                        //Fp(i) = 1.;
                    }
                    // for (uint j=0; j<n_u_dofs; j++)
                    // {
                    //     Kpu(i,j) += -JxW[qp]*dt*psi[i][qp]*dphi[j][qp](0);
                    //     Kpv(i,j) += -JxW[qp]*dt*psi[i][qp]*dphi[j][qp](1);
                    // }
                    // for(uint j=0; j<n_p_dofs; j++)
                    // {
                    //     Kpp(i,j) += -JxW[qp]*ilambda*psi[i][qp]*psi[j][qp];
                    // }
                }
            }

            // loop on element sides
//            for (uint s=0; s<elem->n_sides(); s++)
//            {
//                AutoPtr<Elem const> side (elem->build_side(s));
//                Point const& centroid = side->centroid();

//                // internal bc
//                if((centroid - point_ibc).size() < 1e-6)
//                {
//                    {
//                        const std::vector<std::vector<Real> >&  phi_face = fe_face->get_phi();
//                        const std::vector<Real>& JxWface = fe_face->get_JxW();
//                        const std::vector<Point >& qface_point = fe_face->get_xyz();
//                        const std::vector<Point>& normals = fe_face->get_normals();

//                        fe_face->reinit(elem, s);

//                        for (uint qp=0; qp<qface.n_points(); qp++)
//                        {
//                            Real const JxWfaceqp = (axisym)? JxWface[qp]*q_point[qp](0):JxWface[qp];

//                            // const Real xf = qface_point[qp](0);
//                            // const Real yf = qface_point[qp](1);

//                            // const Real penalty = 1.e10;

//                            //const Real value = external_pressure( qface_point[qp], es.parameters );
//                            const Real value = pressure_bc(qface_point[qp], system.time);

//                            // for (unsigned int i=0; i<phi_face.size(); i++)
//                            // for (unsigned int j=0; j<phi_face.size(); j++)
//                            // Kuu(i,j) += JxW_face[qp]*penalty*phi_face[i][qp]*phi_face[j][qp];

//                            for (uint i=0; i<phi_face.size(); i++)
//                            {
//                                Fu(i) -= JxWfaceqp*dt*value*normals[qp](0)*phi_face[i][qp];
//                                Fv(i) -= JxWfaceqp*dt*value*normals[qp](1)*phi_face[i][qp];
//                            }
//                        }
//                    }
//                }
//            }
        }
        // fluid domain
        else if (elem->subdomain_id() == flag_f)
        {
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
                    Fu(i) += JxWqp*rho_f*(u_old+(f_u - p0.grad(q_point[qp], system.time))*dt)*phi[i][qp];
                    for (uint j=0; j<n_u_dofs; j++)
                    {
                        Kuu(i,j) += JxWqp*( rho_f*phi[j][qp]*phi[i][qp]
                                              + dt*mu_f*(
                                                  dphi[j][qp]*dphi[i][qp]
                                                  + dphi[j][qp](0)*dphi[i][qp](0)
                                                  + axisym*2.*invR2* phi[i][qp]*phi[j][qp]
                                                  )
                                              );
                        Kuv(i,j) += JxWqp*dt*mu_f*dphi[i][qp](1)*dphi[j][qp](0);
                    }
                    for (uint j=0; j<n_p_dofs; j++ )
                    {
                        Kup(i,j) += -JxWqp*dt*psi[j][qp]*(
                                    dphi[i][qp](0)
                                    + axisym*invR * phi[i][qp]
                                    );
                    }

                    Fv(i) += JxWqp*rho_f*(v_old+f_v*dt)*phi[i][qp];
                    for (uint j=0; j<n_u_dofs; j++)
                    {
                        Kvu(i,j) += JxWqp*dt*mu_f*dphi[i][qp](0)*dphi[j][qp](1);
                        Kvv(i,j) += JxWqp*( rho_f*phi[j][qp]*phi[i][qp]
                                              +dt*mu_f*(
                                                  dphi[j][qp]*dphi[i][qp]
                                                  + dphi[j][qp](1)*dphi[i][qp](1)
                                                  )
                                                  );
                    }
                    for (uint j=0; j<n_p_dofs; j++)
                    {
                        Kvp(i,j) += -JxWqp*dt*psi[j][qp]*dphi[i][qp](1);
                    }
                }
                for (uint i=0; i<n_p_dofs; i++)
                {
                    // Kpp(i,i) = 1.;
                    for (uint j=0; j<n_u_dofs; j++)
                    {
                        Kpu(i,j) += -JxWqp*dt*psi[i][qp]*dphi[j][qp](0);
                        Kpv(i,j) += -JxWqp*dt*psi[i][qp]*dphi[j][qp](1);
                    }
                }
            }

            // loop on element sides
//            for (uint s=0; s<elem->n_sides(); s++)
//            {
//                // AutoPtr<Elem> side (elem->build_side(s));
//                // check if side in on boundary
//                if (elem->neighbor(s) == NULL)
//                {
//                    if (mesh.boundary_info->has_boundary_id(elem, s, 6) )
//                    {
//                        const std::vector<std::vector<Real> >&  phi_face = fe_face->get_phi();
//                        const std::vector<Real>& JxWface = fe_face->get_JxW();
//                        const std::vector<Point >& qface_point = fe_face->get_xyz();
//                        const std::vector<Point>& normals = fe_face->get_normals();

//                        fe_face->reinit(elem, s);

//                        for (uint qp=0; qp<qface.n_points(); qp++)
//                        {
//                            Real const JxWfaceqp = (axisym)? JxWface[qp]*q_point[qp](0):JxWface[qp];
//                            // const Real xf = qface_point[qp](0);
//                            // const Real yf = qface_point[qp](1);

//                            // const Real penalty = 1.e10;

//                            const Real value = 1.;

//                            // for (unsigned int i=0; i<phi_face.size(); i++)
//                            // for (unsigned int j=0; j<phi_face.size(); j++)
//                            // Kuu(i,j) += JxW_face[qp]*penalty*phi_face[i][qp]*phi_face[j][qp];

//                            for (uint i=0; i<phi_face.size(); i++)
//                            {
//                                Fu(i) -= JxWfaceqp*dt*value*normals[qp](0)*phi_face[i][qp];
//                                Fv(i) -= JxWfaceqp*dt*value*normals[qp](1)*phi_face[i][qp];
//                            }
//                        }
//                    }
//                }
//            }
        }
        else
        {
            std::cerr << "elem subdomain not recognized!" << std::endl;
            abort();
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

    // system.matrix->close();
    // system.matrix->print_matlab("mat_fsi.m");

    // system.rhs->close();
    // system.rhs->print_matlab("rhs_fsi.m");

    // That's it.
    return;
}
