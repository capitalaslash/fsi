#ifndef PRESSURERAMP_HPP
#define PRESSURERAMP_HPP

#include "FSI.hpp"

#include <libmesh/function_base.h>
#include <libmesh/getpot.h>

#include "bc/damp.hpp"

template <typename DampType>
struct PressureRamp: public libMesh::FunctionBase<Real>
{
    PressureRamp(Real const tin, Real const tout, Real const pmax, DampType const & damp):
        M_tin(tin),
        M_tout(tout),
        M_pmax(pmax),
        M_damp(damp)
    {}

    PressureRamp(GetPot const& param):
        M_tin (param("pressure_tin",  0.1)),
        M_tout(param("pressure_tout", 0.9)),
        M_pmax(param("pressure_max",  1.0)),
        M_damp(DampType(param))
    {}

    void operator()(const libMesh::Point& p, const Real t, libMesh::DenseVector<Real>& output)
    {
        if( t < M_tin )       output(2) = 0.0;
        else if( t > M_tout ) output(2) = M_pmax;
        else                  output(2) = M_pmax*(t - M_tin)/(M_tout - M_tin);
        output(2) *= M_damp(p(1));
    }

    Real operator()(const libMesh::Point& p, const Real t)
    {
        if( t < M_tin )       return 0.0;
        else if( t > M_tout ) return M_damp(p(1))*M_pmax;
        else                  return M_damp(p(1))*M_pmax*(t - M_tin)/(M_tout - M_tin);
    }

    libMesh::AutoPtr<libMesh::FunctionBase<Real> > clone () const
    {
      return libMesh::AutoPtr<FunctionBase<Real> >
        ( new PressureRamp(M_tin,M_tout,M_pmax,M_damp) );
    }

    Real const M_tin;
    Real const M_tout;
    Real const M_pmax;
    DampType const M_damp;
};

#endif // PRESSURERAMP_HPP
