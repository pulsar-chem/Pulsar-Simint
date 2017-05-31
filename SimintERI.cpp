#include <pulsar/output/OutputStream.hpp>
#include <pulsar/system/AOOrdering.hpp>
#include <pulsar/system/SphericalTransformIntegral.hpp>

#include "SimintERI.hpp"

#include "simint/simint.h"

using namespace pulsar;

//This is what gaussian_shell is replaced by
using SimintShell_t=simint_shell;
//This is what multishell_pair is replaced by
using SimintPair=simint_multi_shellpair;

static
SimintShell_t ToSimintShell(const BasisSetShell & shell, int igen)
{
    size_t nprim = shell.n_primitives();

    SimintShell_t gs;
    simint_allocate_shell(nprim, &gs);
    gs.nprim = nprim;
    gs.am = shell.general_am(igen);
    gs.x = shell.get_coord(0);
    gs.y = shell.get_coord(1);
    gs.z = shell.get_coord(2);

    const double * alphas = shell.alpha_ptr();
    const double * coefs = shell.coef_ptr(igen);
    std::copy(alphas, alphas + nprim, gs.alpha);
    std::copy(coefs, coefs + nprim, gs.coef);

    simint_normalize_shells(1, &gs);

    return gs;
}


double *SimintERI::calculate_(size_t shell1, size_t shell2,
                               size_t shell3, size_t shell4)
{
    const BasisSetShell & sh1 = bs1_.shell(shell1);
    const BasisSetShell & sh2 = bs2_.shell(shell2);
    const BasisSetShell & sh3 = bs3_.shell(shell3);
    const BasisSetShell & sh4 = bs4_.shell(shell4);

    double * srcptr = sourcework_;

    for(size_t ng1 = 0; ng1 < sh1.n_general_contractions(); ng1++)
    {
        SimintShell_t gs1 = ToSimintShell(sh1, ng1);
        const size_t ncart1 = n_cartesian_gaussian(sh1.general_am(ng1));

        for(size_t ng2 = 0; ng2 < sh2.n_general_contractions(); ng2++)
        {
            SimintShell_t gs2 = ToSimintShell(sh2, ng2);
            const size_t ncart2 = n_cartesian_gaussian(sh2.general_am(ng2));
            SimintPair bra_pair;
            simint_initialize_multi_shellpair(&bra_pair);
            simint_create_multi_shellpair(1, &gs1, 1, &gs2, &bra_pair, 0);

            for(size_t ng3 = 0; ng3 < sh3.n_general_contractions(); ng3++)
            {
                SimintShell_t gs3 = ToSimintShell(sh3, ng3);
                const size_t ncart3 = n_cartesian_gaussian(sh3.general_am(ng3));

                for(size_t ng4 = 0; ng4 < sh4.n_general_contractions(); ng4++)
                {
                    // form the shell pair info
                    SimintShell_t gs4 = ToSimintShell(sh4, ng4);
                    SimintPair ket_pair;
                    simint_initialize_multi_shellpair(&ket_pair);
                    simint_create_multi_shellpair(1, &gs3, 1, &gs4,&ket_pair,0);
                    const size_t ncart4 = n_cartesian_gaussian(sh4.general_am(ng4));
                    const size_t ncart = ncart1 * ncart2 * ncart3 * ncart4;

                    size_t ncalc = simint_compute_eri(&bra_pair, &ket_pair, 0.0, srcptr);

                    simint_free_shell(&gs4);
                    simint_free_multi_shellpair(&ket_pair);

                    srcptr += ncalc * ncart;
                }

                simint_free_shell(&gs3);
            }

            simint_free_shell(&gs2);
            simint_free_multi_shellpair(&bra_pair);
        }

        simint_free_shell(&gs1);
    }


    CartesianToSpherical_4Center(sh1, sh2, sh3, sh4, sourcework_, buffer_.data(), transformwork_, 1);

    return buffer_.data();
}

//Returns 3 multichoose l, i.e. the number of Cartesian Gaussians for angular
//momentum l
static size_t multichoose(size_t l)
{
    if(l==0)return 1;
    else if(l==1)return 3;
    else if(l==2)return 6;
    else if(l==3)return 10;
    else if(l==4)return 15;
    else if(l==5)return 21;
    else if(l==6)return 28;
    else if(l==7)return 36;
    else
      throw PulsarException("Simint only supports up to angular momentum l==7");
}


void SimintERI::initialize_(unsigned int deriv,
                            const Wavefunction &,
                            const BasisSet & bs1, const BasisSet & bs2,
                            const BasisSet & bs3, const BasisSet & bs4)
{
    if(deriv != 0)
        throw PulsarException("Not Yet Implemented: Simint with deriv != 0");

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
    
    maxsize1 =  multichoose(*bs1_.all_am().rbegin());
    maxsize2 =  multichoose(*bs2_.all_am().rbegin());
    maxsize3 =  multichoose(*bs3_.all_am().rbegin());
    maxsize4 =  multichoose(*bs4_.all_am().rbegin());
    size_t sourcework_size = maxsize1*maxsize2*maxsize3*maxsize4;

    work_.resize(sourcework_size+transformwork_size);
    sourcework_ = work_.data();
    transformwork_ = sourcework_ + sourcework_size;
    buffer_.resize(transformwork_size);
}

