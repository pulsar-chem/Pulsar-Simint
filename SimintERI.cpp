#include <pulsar/output/OutputStream.hpp>
#include <pulsar/system/AOOrdering.hpp>
#include <pulsar/system/SphericalTransformIntegral.hpp>

#include "SimintERI.hpp"

#include "simint/simint_init.h"
#include "simint/eri/eri.h"


using namespace pulsar::output;
using namespace pulsar::exception;
using namespace pulsar::system;



SimintERI::SimintERI(ID_t id)
    : TwoElectronIntegral(id)
{ }



SimintERI::~SimintERI()
{ }



static
gaussian_shell ToSimintShell(const BasisSetShell & shell, int igen)
{
    size_t nprim = shell.NPrim();

    gaussian_shell gs;
    allocate_gaussian_shell(nprim, &gs);
    gs.nprim = nprim;
    gs.am = shell.GeneralAM(igen);
    gs.x = shell.GetCoord(0);
    gs.y = shell.GetCoord(1);
    gs.z = shell.GetCoord(2);

    const double * alphas = shell.AlphaPtr();
    const double * coefs = shell.CoefPtr(igen);
    std::copy(alphas, alphas + nprim, gs.alpha);
    std::copy(coefs, coefs + nprim, gs.coef);

    normalize_gaussian_shells(1, &gs);

    return gs;
}


uint64_t SimintERI::Calculate_(size_t deriv,
                               size_t shell1, size_t shell2,
                               size_t shell3, size_t shell4,
                               double * outbuffer, size_t bufsize)
{
    if(deriv != 0)
        throw NotYetImplementedException("Not Yet Implemented: Overlap integral with deriv != 0");

    const BasisSetShell & sh1 = bs1_->Shell(shell1);
    const BasisSetShell & sh2 = bs2_->Shell(shell2);
    const BasisSetShell & sh3 = bs3_->Shell(shell3);
    const BasisSetShell & sh4 = bs4_->Shell(shell4);

    size_t nfunc = sh1.NFunctions() * sh2.NFunctions() * sh3.NFunctions() * sh4.NFunctions();

    if(bufsize < nfunc)
        throw GeneralException("Buffer to small for ERI", "bufsize", bufsize, "nfunc", nfunc);

    double * srcptr = sourcework_;

    for(size_t ng1 = 0; ng1 < sh1.NGeneral(); ng1++)
    {
        gaussian_shell gs1 = ToSimintShell(sh1, ng1);
        const size_t ncart1 = NCartesianGaussian(sh1.GeneralAM(ng1));

        for(size_t ng2 = 0; ng2 < sh2.NGeneral(); ng2++)
        {
            gaussian_shell gs2 = ToSimintShell(sh2, ng2);
            const size_t ncart2 = NCartesianGaussian(sh2.GeneralAM(ng2));
            multishell_pair bra_pair = create_multishell_pair(1, &gs1, 1, &gs2); 

            for(size_t ng3 = 0; ng3 < sh3.NGeneral(); ng3++)
            {
                gaussian_shell gs3 = ToSimintShell(sh3, ng3);
                const size_t ncart3 = NCartesianGaussian(sh3.GeneralAM(ng3));

                for(size_t ng4 = 0; ng4 < sh4.NGeneral(); ng4++)
                {
                    // form the shell pair info
                    gaussian_shell gs4 = ToSimintShell(sh4, ng4);
                    multishell_pair ket_pair = create_multishell_pair(1, &gs3, 1, &gs4); 
                    const size_t ncart4 = NCartesianGaussian(sh4.GeneralAM(ng4));
                    const size_t ncart = ncart1 * ncart2 * ncart3 * ncart4;

                    size_t ncalc = simint_compute_eri(bra_pair, ket_pair, srcptr);

                    free_gaussian_shell(gs4);
                    free_multishell_pair(ket_pair);

                    srcptr += ncalc * ncart;
                }

                free_gaussian_shell(gs3);
            }

            free_gaussian_shell(gs2);
            free_multishell_pair(bra_pair);
        }

        free_gaussian_shell(gs1);
    }


    CartesianToSpherical_4Center(sh1, sh2, sh3, sh4, sourcework_, outbuffer, transformwork_);

    return nfunc;
}


void SimintERI::SetBases_(const System & sys,
                          const std::string & bs1, const std::string & bs2,
                          const std::string & bs3, const std::string & bs4)
{
    // initialize the simint library
    simint_init();

    // get the basis sets from the system
    // Note - storing un-normalized
    bs1_ = std::make_shared<BasisSet>(sys.GetBasisSet(bs1));
    bs2_ = std::make_shared<BasisSet>(sys.GetBasisSet(bs2));
    bs3_ = std::make_shared<BasisSet>(sys.GetBasisSet(bs3));
    bs4_ = std::make_shared<BasisSet>(sys.GetBasisSet(bs4));


    // determine work sizes
    size_t maxsize1 = bs1_->MaxProperty(NCartesianGaussianForShellAM);
    size_t maxsize2 = bs2_->MaxProperty(NCartesianGaussianForShellAM);
    size_t maxsize3 = bs3_->MaxProperty(NCartesianGaussianForShellAM);
    size_t maxsize4 = bs4_->MaxProperty(NCartesianGaussianForShellAM);
    size_t transformwork_size = maxsize1*maxsize2*maxsize3*maxsize4;
    
    maxsize1 =  bs1_->MaxProperty(NCartesianGaussianInShell);
    maxsize2 =  bs2_->MaxProperty(NCartesianGaussianInShell);
    maxsize3 =  bs3_->MaxProperty(NCartesianGaussianInShell);
    maxsize4 =  bs4_->MaxProperty(NCartesianGaussianInShell);
    size_t sourcework_size = maxsize1*maxsize2*maxsize3*maxsize4;

    work_.resize(sourcework_size+transformwork_size);
    sourcework_ = work_.data();
    transformwork_ = sourcework_ + sourcework_size;   
}

