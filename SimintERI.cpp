#include <pulsar/output/OutputStream.hpp>
#include <pulsar/system/AOOrdering.hpp>
#include <pulsar/system/SphericalTransformIntegral.hpp>

#include "SimintERI.hpp"

#include "simint/simint_init.h"
#include "simint/eri/eri.h"


using namespace pulsar::output;
using namespace pulsar::exception;
using namespace pulsar::system;
using namespace pulsar::datastore;



static
gaussian_shell ToSimintShell(const BasisSetShell & shell, int igen)
{
    size_t nprim = shell.n_primitives();

    gaussian_shell gs;
    allocate_gaussian_shell(nprim, &gs);
    gs.nprim = nprim;
    gs.am = shell.general_am(igen);
    gs.x = shell.get_coord(0);
    gs.y = shell.get_coord(1);
    gs.z = shell.get_coord(2);

    const double * alphas = shell.alpha_ptr();
    const double * coefs = shell.coef_ptr(igen);
    std::copy(alphas, alphas + nprim, gs.alpha);
    std::copy(coefs, coefs + nprim, gs.coef);

    normalize_gaussian_shells(1, &gs);

    return gs;
}


uint64_t SimintERI::calculate_(size_t shell1, size_t shell2,
                               size_t shell3, size_t shell4,
                               double * outbuffer, size_t bufsize)
{
    const BasisSetShell & sh1 = bs1_.shell(shell1);
    const BasisSetShell & sh2 = bs2_.shell(shell2);
    const BasisSetShell & sh3 = bs3_.shell(shell3);
    const BasisSetShell & sh4 = bs4_.shell(shell4);

    size_t nfunc = sh1.n_functions() * sh2.n_functions() * sh3.n_functions() * sh4.n_functions();

    if(bufsize < nfunc)
        throw GeneralException("Buffer to small for ERI", "bufsize", bufsize, "nfunc", nfunc);

    double * srcptr = sourcework_;

    for(size_t ng1 = 0; ng1 < sh1.n_general_contractions(); ng1++)
    {
        gaussian_shell gs1 = ToSimintShell(sh1, ng1);
        const size_t ncart1 = n_cartesian_gaussian(sh1.general_am(ng1));

        for(size_t ng2 = 0; ng2 < sh2.n_general_contractions(); ng2++)
        {
            gaussian_shell gs2 = ToSimintShell(sh2, ng2);
            const size_t ncart2 = n_cartesian_gaussian(sh2.general_am(ng2));
            multishell_pair bra_pair = create_multishell_pair(1, &gs1, 1, &gs2); 

            for(size_t ng3 = 0; ng3 < sh3.n_general_contractions(); ng3++)
            {
                gaussian_shell gs3 = ToSimintShell(sh3, ng3);
                const size_t ncart3 = n_cartesian_gaussian(sh3.general_am(ng3));

                for(size_t ng4 = 0; ng4 < sh4.n_general_contractions(); ng4++)
                {
                    // form the shell pair info
                    gaussian_shell gs4 = ToSimintShell(sh4, ng4);
                    multishell_pair ket_pair = create_multishell_pair(1, &gs3, 1, &gs4); 
                    const size_t ncart4 = n_cartesian_gaussian(sh4.general_am(ng4));
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


    CartesianToSpherical_4Center(sh1, sh2, sh3, sh4, sourcework_, outbuffer, transformwork_, 1);

    return nfunc;
}


void SimintERI::initialize_(unsigned int deriv,
                            const Wavefunction & wfn,
                            const BasisSet & bs1, const BasisSet & bs2,
                            const BasisSet & bs3, const BasisSet & bs4)
{
    if(deriv != 0)
        throw NotYetImplementedException("Not Yet Implemented: Simint with deriv != 0");

    // initialize the simint library
    simint_init();

    // get the basis sets from the system
    // Note - storing un-normalized
    bs1_ = bs1;
    bs2_ = bs2;
    bs3_ = bs3;
    bs4_ = bs4;


    // determine work sizes
    size_t maxsize1 = bs1_.max_property(n_cartesian_gaussian_for_shell_am);
    size_t maxsize2 = bs2_.max_property(n_cartesian_gaussian_for_shell_am);
    size_t maxsize3 = bs3_.max_property(n_cartesian_gaussian_for_shell_am);
    size_t maxsize4 = bs4_.max_property(n_cartesian_gaussian_for_shell_am);
    size_t transformwork_size = maxsize1*maxsize2*maxsize3*maxsize4;
    
    maxsize1 =  bs1_.max_property(n_cartesian_gaussian_in_shell);
    maxsize2 =  bs2_.max_property(n_cartesian_gaussian_in_shell);
    maxsize3 =  bs3_.max_property(n_cartesian_gaussian_in_shell);
    maxsize4 =  bs4_.max_property(n_cartesian_gaussian_in_shell);
    size_t sourcework_size = maxsize1*maxsize2*maxsize3*maxsize4;

    work_.resize(sourcework_size+transformwork_size);
    sourcework_ = work_.data();
    transformwork_ = sourcework_ + sourcework_size;   
}

