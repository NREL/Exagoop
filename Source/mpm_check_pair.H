#ifndef MPM_CHECK_PAIR_H_
#define MPM_CHECK_PAIR_H_

#include <mpm_specs.H>

struct CheckPair
{
    template <class P>
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    bool operator()(const P& p1, const P& p2) const
    {
        amrex::Real d0 = (p1.pos(0) - p2.pos(0));
        amrex::Real d1 = (p1.pos(1) - p2.pos(1));
        amrex::Real d2 = (p1.pos(2) - p2.pos(2));    
        amrex::Real dsquared = d0*d0 + d1*d1 + d2*d2;

        amrex::Real search_radius=0.5*(p1.rdata(realData::radius)+p2.rdata(realData::radius));
       
        return (dsquared <= 9.0*search_radius*search_radius);
    }
};

#endif
