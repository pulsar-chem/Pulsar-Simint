#ifndef _GUARD_SIMINTRI_HPP_
#define _GUARD_SIMINTRI_HPP_

#include <pulsar/modulebase/TwoElectronIntegral.hpp>
#include <pulsar/system/BasisSet.hpp>

class SimintERI : public pulsar::modulebase::TwoElectronIntegral
{
public:
    SimintERI(ID_t id);

    virtual void SetBases_(const pulsar::system::System & sys,
                           const std::string & bs1, const std::string & bs2,
                           const std::string & bs3, const std::string & bs4);


    virtual uint64_t Calculate_(size_t deriv,
                                size_t shell1, size_t shell2,
                                size_t shell3, size_t shell4,
                                double * outbuffer, size_t bufsize);

    virtual ~SimintERI();


private:
    std::shared_ptr<pulsar::system::BasisSet> bs1_, bs2_, bs3_, bs4_;

    std::vector<double> work_;
    double * sourcework_;
    double * transformwork_;

};


#endif
