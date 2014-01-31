#ifndef PRESSURERAMP_HPP
#define PRESSURERAMP_HPP

#include <libmesh/function_base.h>

using libMesh::Real;

struct PressureRamp: public libMesh::FunctionBase<Real>
{
    PressureRamp( Real const tin, Real const tout, Real const pmax):
        M_tin(tin),
        M_tout(tout),
        M_pmax(pmax)
    {}

    void operator()(const Point& /*p*/, const Real t, libMesh::DenseVector<Real>& output)
    {
        if( t < M_tin )       output(2) = 0.0;
        else if( t > M_tout ) output(2) = M_pmax;
        else                  output(2) = M_pmax*(t - M_tin)/(M_tout - M_tin);
    }

    Real operator()(const libMesh::Point& /*p*/, const Real t)
    {
        if( t < M_tin )       return 0.0;
        else if( t > M_tout ) return M_pmax;
        else                  return M_pmax*(t - M_tin)/(M_tout - M_tin);
    }

    libMesh::AutoPtr<libMesh::FunctionBase<Real> > clone () const
    {
      return libMesh::AutoPtr<FunctionBase<Real> >
        ( new PressureRamp(M_tin,M_tout,M_pmax) );
    }

    Real const M_tin;
    Real const M_tout;
    Real const M_pmax;
};

#endif // PRESSURERAMP_HPP

