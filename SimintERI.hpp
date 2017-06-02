#pragma once

#include <pulsar/modulebase/FourCenterIntegral.hpp>
#include <pulsar/system/BasisSet.hpp>
#include "simint/simint.h"

class SimintERI : public pulsar::FourCenterIntegral
{
public:
    using pulsar::FourCenterIntegral::FourCenterIntegral;

    HashType my_hash_(unsigned int deriv,
                      const pulsar::Wavefunction&,
                      const pulsar::BasisSet& bs1,
                      const pulsar::BasisSet& bs2,
                      const pulsar::BasisSet& bs3,
                      const pulsar::BasisSet& bs4)
    {
        return bphash::hash_to_string(bs1.my_hash())+
                bphash::hash_to_string(bs2.my_hash())+
                bphash::hash_to_string(bs3.my_hash())+
                bphash::hash_to_string(bs4.my_hash())+
                std::to_string(deriv);
    }

    virtual void initialize_(unsigned int deriv,
                             const pulsar::Wavefunction & wfn,
                             const pulsar::BasisSet & bs1,
                             const pulsar::BasisSet & bs2,
                             const pulsar::BasisSet & bs3,
                             const pulsar::BasisSet & bs4);

    virtual double* calculate_(size_t shell1, size_t shell2,
                               size_t shell3, size_t shell4);

    typedef std::vector<simint_multi_shellpair> ShellPairVec;
    typedef std::vector<simint_shell> ShellVec;

private:


    std::array<pulsar::BasisSet,4> bs_;
    std::array<std::shared_ptr<const ShellVec>,4> shells_;
    std::array<std::shared_ptr<const ShellPairVec>,2> single_spairs_;
    std::array<std::shared_ptr<const ShellPairVec>,2> multi_spairs_;

    double * buffer_;
    double * sharedwork_;
    double * allwork_;
    double * source_full_;
    double * tformbuf_full_;
};
