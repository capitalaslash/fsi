#include "stress.hpp"

#include <libmesh/linear_implicit_system.h>
#include <libmesh/dof_map.h>
#include <libmesh/quadrature_gauss.h>
#include <libmesh/numeric_vector.h>
#include <libmesh/utility.h>

using namespace libMesh;
using namespace libMesh::Utility;

// monolithic block for displacement
void compute_stress(EquationSystems& /*es*/, std::string /*sys_name*/)
{

}

// different systems for each disp component
void compute_stress(EquationSystems& es)
{
  const MeshBase& mesh = es.get_mesh();

  const uint dim = mesh.mesh_dimension();

  LinearImplicitSystem const& system_d = es.get_system<LinearImplicitSystem>("dx");
  LinearImplicitSystem const& system_e = es.get_system<LinearImplicitSystem>("dy");

  DofMap const& dof_map_d = system_d.get_dof_map();
  DofMap const& dof_map_e = system_e.get_dof_map();
  FEType fe_type = dof_map_d.variable_type(0);
  AutoPtr<FEBase> fe (FEBase::build(dim, fe_type));
  QGauss qrule (dim, fe_type.default_quadrature_order());
  fe->attach_quadrature_rule (&qrule);

  const std::vector<Real>& JxW = fe->get_JxW();
  const std::vector<std::vector<RealGradient> >& dphi = fe->get_dphi();

  // Also, get a reference to the ExplicitSystem
  ExplicitSystem& stress = es.get_system<ExplicitSystem>("stress");
  const DofMap& dof_map_s = stress.get_dof_map();
  std::vector<std::vector<uint> > sigma_vars( dim, std::vector<uint>(dim) );
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
  std::vector< std::vector<dof_id_type> > dof_indices_var(2);
  std::vector<dof_id_type> dof_indices_d;
  std::vector<dof_id_type> dof_indices_e;
  std::vector<dof_id_type> dof_indices_s;

  // To store the stress tensor on each element
  DenseMatrix<Number> elem_sigma;

  Real const mu = es.parameters.get<Real>("mu_s");
  Real const lambda = es.parameters.get<Real>("lambda");
  subdomain_id_type const flag_s = es.parameters.get<subdomain_id_type>("flag_s");

  MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();

  for ( ; el != end_el; ++el)
  {
      const Elem* elem = *el;

      if (elem->subdomain_id() == flag_s)
      {
          dof_map_d.dof_indices (elem, dof_indices_var[0], 0);
          dof_map_e.dof_indices (elem, dof_indices_var[1], 0);
          dof_map_d.dof_indices (elem, dof_indices_d, 0);
          dof_map_e.dof_indices (elem, dof_indices_e, 0);

          const uint n_dofs = dof_indices_d.size();

          fe->reinit (elem);

          elem_sigma.resize(3, 3);

          for (uint qp=0; qp<qrule.n_points(); qp++)
          {
              // Get the gradient at this quadrature point
              Gradient grad_d;
              Gradient grad_e;
              for(uint l=0; l<n_dofs; l++)
              {
                  grad_d.add_scaled(dphi[l][qp], system_d.current_solution(dof_indices_d[l]));
                  grad_e.add_scaled(dphi[l][qp], system_e.current_solution(dof_indices_e[l]));
              }

              for (uint j=0; j<dim; j++)
              {
                  elem_sigma(0,j) += JxW[qp]*mu*grad_d(j);
                  elem_sigma(j,0) += JxW[qp]*mu*grad_d(j);
                  elem_sigma(1,j) += JxW[qp]*mu*grad_e(j);
                  elem_sigma(j,1) += JxW[qp]*mu*grad_e(j);
              }
              const Real div = (grad_d(0) + grad_e(1));
              elem_sigma(0,0) += JxW[qp]*lambda*div;
              elem_sigma(1,1) += JxW[qp]*lambda*div;
          }

          // Get the average stresses by dividing by the element volume
          elem_sigma.scale(1./elem->volume());

          // load elem_sigma data into stress_system
          for(uint i=0; i<dim; i++)
              for(uint j=0; j<dim; j++)
              {
                  dof_map_s.dof_indices (elem, dof_indices_s, sigma_vars[i][j]);

                  // We are using CONSTANT MONOMIAL basis functions, hence we only need to get
                  // one dof index per variable
                  dof_id_type dof_index = dof_indices_s[0];

                  if( (stress.solution->first_local_index() <= dof_index) &&
                          (dof_index < stress.solution->last_local_index()) )
                  {
                      stress.solution->set(dof_index, elem_sigma(i,j));
                  }

              }

          // Also, the von Mises stress
          Number vonMises_value = std::sqrt( 0.5*( pow<2>(elem_sigma(0,0) - elem_sigma(1,1)) +
                                                   pow<2>(elem_sigma(1,1) - elem_sigma(2,2)) +
                                                   pow<2>(elem_sigma(2,2) - elem_sigma(0,0)) +
                                                   6.*(pow<2>(elem_sigma(0,1)) +
                                                       pow<2>(elem_sigma(1,2)) +
                                                       pow<2>(elem_sigma(2,0)))
                                                   ) );
          dof_map_s.dof_indices (elem, dof_indices_s, vonMises_var);
          dof_id_type dof_index = dof_indices_s[0];
          if( (stress.solution->first_local_index() <= dof_index) &&
                  (dof_index < stress.solution->last_local_index()) )
          {
              stress.solution->set(dof_index, vonMises_value);
          }
      }
  }

  // Should call close and update when we set vector entries directly
  stress.solution->close();
  stress.update();
}

void StressSystem::assemble()
{
    EquationSystems const & es = get_equation_systems();

    const MeshBase& mesh = get_mesh();

    const uint dim = mesh.mesh_dimension();

    const DofMap& dof_map = this->get_dof_map();

    DofMap const& dof_map_d = M_systemDx.get_dof_map();
    //DofMap const& dof_map_e = M_systemDy.get_dof_map();
    FEType fe_type_d = dof_map_d.variable_type(0);
    AutoPtr<FEBase> fe_d (FEBase::build(dim, fe_type_d));
    QGauss qrule_d (dim, fe_type_d.default_quadrature_order());
    fe_d->attach_quadrature_rule (&qrule_d);

    const std::vector<Real>& JxW = fe_d->get_JxW();
    const std::vector<std::vector<Real> >& b = fe_d->get_phi();
    const std::vector<std::vector<RealGradient> >& grad_b = fe_d->get_dphi();
    const std::vector<Point >& q_point = fe_d->get_xyz();

    std::vector<dof_id_type> dof_indices;
    std::vector<dof_id_type> dof_indices_d;

    DenseMatrix<Real> elem_sigma;

    Real const mu = es.parameters.get<Real>("mu_s");
    Real const lambda = es.parameters.get<Real>("lambda");
    subdomain_id_type const flag_s = es.parameters.get<subdomain_id_type>("flag_s");

    MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
    const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();

    for ( ; el != end_el; ++el)
    {
        const Elem* elem = *el;

        if (elem->subdomain_id() == flag_s)
        {
            dof_map_d.dof_indices (elem, dof_indices_d, 0);

            const uint n_dofs = dof_indices_d.size();

            fe_d->reinit (elem);

            elem_sigma.resize(3, 3);

            for (uint qp=0; qp<qrule_d.n_points(); qp++)
            {
                Real const invR = 1./q_point[qp](0);

                // Get the gradient at this quadrature point
                Real d_old = 0.;
                Gradient grad_d;
                Gradient grad_e;
                for(uint l=0; l<n_dofs; l++)
                {
                    d_old += b[l][qp] * M_systemDx.current_solution (dof_indices_d[l]);
                    grad_d.add_scaled(grad_b[l][qp], M_systemDx.current_solution(dof_indices_d[l]));
                    grad_e.add_scaled(grad_b[l][qp], M_systemDy.current_solution(dof_indices_d[l]));
                }

                // AXISYM!!!
                Real trace_eps = grad_d(0) + grad_e(1) + d_old*invR;
                elem_sigma(0,0) += JxW[qp] * (2.*mu*grad_d(0)  + lambda*trace_eps);
                elem_sigma(1,1) += JxW[qp] * (2.*mu*grad_e(1)  + lambda*trace_eps);
                elem_sigma(2,2) += JxW[qp] * (2.*mu*d_old*invR + lambda*trace_eps);
                elem_sigma(0,1) += JxW[qp] * (mu*(grad_d(1)+grad_e(0)));
                elem_sigma(1,0) += JxW[qp] * (mu*(grad_d(1)+grad_e(0)));
            }

            // Get the average stresses by dividing by the element volume
            elem_sigma.scale(1./elem->volume());

            // load elem_sigma data into stress_system
            for(uint i=0; i<3; i++)
                for(uint j=0; j<3; j++)
                {
                    dof_map.dof_indices (elem, dof_indices, M_sVar[i][j]);

                    // We are using CONSTANT MONOMIAL basis functions, hence we only need to get
                    // one dof index per variable
                    dof_id_type dof_index = dof_indices[0];

                    if( (this->solution->first_local_index() <= dof_index) &&
                            (dof_index < this->solution->last_local_index()) )
                    {
                        this->solution->set(dof_index, elem_sigma(i,j));
                    }

                }

            // Also, the von Mises stress
            Number vonMises_value = std::sqrt( 0.5*( pow<2>(elem_sigma(0,0) - elem_sigma(1,1)) +
                                                     pow<2>(elem_sigma(1,1) - elem_sigma(2,2)) +
                                                     pow<2>(elem_sigma(2,2) - elem_sigma(0,0)) +
                                                     6.*(pow<2>(elem_sigma(0,1)) +
                                                         pow<2>(elem_sigma(1,2)) +
                                                         pow<2>(elem_sigma(2,0)))
                                                     ) );
            dof_map.dof_indices (elem, dof_indices, M_vonMisesVar);
            dof_id_type dof_index = dof_indices[0];
            if( (this->solution->first_local_index() <= dof_index) &&
                    (dof_index < this->solution->last_local_index()) )
            {
                this->solution->set(dof_index, vonMises_value);
            }
        }
    }

    this->solution->close();
    this->update();
}

template <>
void StressSystem::localStress<TWOD_AXI> (DenseMatrix<Real> & /*s*/)
{

}
