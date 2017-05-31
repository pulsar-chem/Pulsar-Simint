#ifndef _GUARD_SIMINTRI_HPP_
#define _GUARD_SIMINTRI_HPP_

#include <pulsar/modulebase/FourCenterIntegral.hpp>
#include <pulsar/system/BasisSet.hpp>

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

private:
    pulsar::BasisSet bs1_, bs2_, bs3_, bs4_;
    std::vector<double> buffer_;
    std::vector<double> work_;
    double * sourcework_;
    double * transformwork_;

};


#endif
