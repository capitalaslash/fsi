#ifndef DAMP_HPP
#define DAMP_HPP

#include <libmesh/getpot.h>

using libMesh::Real;

struct Damp
{
    Damp() {}
    explicit Damp(GetPot const & /*param*/) {}

    virtual Real operator()( Real const /*x*/ ) const
    {
        return 1.0;
    }
};

struct DampNull: public Damp
{
    explicit DampNull(GetPot const & param):
        Damp(param)
    {}

    Real operator()( Real const ) const
    {
        return 1.0;
    }
};

struct DampZero: public Damp
{
    explicit DampZero(GetPot const & param):
        Damp(param)
    {}

    Real operator()( Real const ) const
    {
        return 0.0;
    }
};

struct DampLogistic: public Damp
{
    DampLogistic( Real const xmin, Real const xmax, Real const sigma):
        Damp(),
        M_mean(.5*(xmax+xmin)),
        M_delta(xmax-xmin),
        M_sigma(sigma)
    {}

    explicit DampLogistic(GetPot const & param):
        Damp(param),
        M_mean(.5*(param("damp_xmax",1.0)+param("damp_xmin",0.0))),
        M_delta(param("damp_xmax",1.0)-param("damp_xmin",0.0)),
        M_sigma(param("damp_sigma",10))
    {}

    Real operator()( Real const x) const
    {
        return 1.0/(1.0+exp(M_sigma*(x-M_mean)/M_delta));
    }

    Real const M_mean;
    Real const M_delta;
    Real const M_sigma;
};

#endif // DAMP_HPP
