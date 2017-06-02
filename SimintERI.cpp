#include <pulsar/output/OutputStream.hpp>
#include <pulsar/system/AOOrdering.hpp>
#include <pulsar/system/SphericalTransformIntegral.hpp>

#include "SimintERI.hpp"

using namespace pulsar;

using SimintERI::ShellVec;
using SimintERI::ShellPairVec;

//Transforms the igen-th general contraction of shell into a simint shell
static simint_shell psr_to_simint_(const BasisSetShell & shell, size_t igen)
{
    size_t nprim = shell.n_primitives();
    simint_shell gs;
    simint_initialize_shell(&gs);
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

    //simint_normalize_shells(1, &gs);

    return gs;
}

static simint_multi_shellpair make_pair_(const simint_shell& si,
                                         const simint_shell& sj,
                                         double screen_thresh)
{
    simint_multi_shellpair da_pair;
    simint_initialize_multi_shellpair(&da_pair);
    simint_create_multi_shellpair(1, &si, 1, &sj, &da_pair, screen_thresh);
    return da_pair;
}


// Some helpers
static void multishell_vector_deleter_(ShellPairVec * mpv)
{
    for(auto & mp : *mpv)
        simint_free_multi_shellpair(&mp);
    mpv->clear();
}

static std::shared_ptr<ShellVec>
create_shell_vec_(const BasisSet & bs)
{
    auto shells = std::shared_ptr<ShellVec>(new ShellVec, [](ShellVec * mpv){
            for(auto & mp : *mpv)
            simint_free_shell(&mp);
            mpv->clear();
});
    for(auto sh : bs)
        for(size_t ng = 0; ng < sh.n_general_contractions(); ng++)
            shells->push_back(psr_to_simint_(sh, ng));
    return shells;
}

static std::shared_ptr<ShellPairVec>
create_shell_pair_(const ShellVec & bs1, const ShellVec & bs2)
{
    auto spairs = std::shared_ptr<ShellPairVec>(new ShellPairVec, &multishell_vector_deleter_);

    size_t nshell1 = bs1.size();
    size_t nshell2 = bs2.size();
    spairs->resize(nshell1 * nshell2);

    for(size_t i = 0; i < nshell1; i++)
    {
        const simint_shell & sh_i = bs1[i];
        for(size_t j = 0; j < nshell2; j++)
        {
            const simint_shell & sh_j = bs2[j];
            simint_multi_shellpair P;
            simint_initialize_multi_shellpair(&P);
            simint_create_multi_shellpair(1, &sh_i, 1, &sh_j, &P, SIMINT_SCREEN);
            (*spairs)[i*nshell2+j] = P;
        }
    }

    return spairs;
}

static ShellPairVec
create_multi_shellpair_(const std::vector<ShellPairBlock> & vsh,
                        const ShellVec & shell1,
                        const ShellVec & shell2)
{
    ShellPairVec ret;

    // copying simint shells will copy the pointers. That is ok - we won't
    // do anything special to delete them
    for(const auto & spairs : vsh)
    {
        std::vector<simint_shell> simint_shells;

        for(const auto & s : spairs)
        {
            simint_shells.push_back(shell1[s.first]);
            simint_shells.push_back(shell2[s.second]);
        }

        // spairs.size() should be simint_shells.size()/2
        simint_multi_shellpair P;
        simint_initialize_multi_shellpair(&P);
        simint_create_multi_shellpair2(spairs.size(), simint_shells.data(), &P, SIMINT_SCREEN);
        ret.push_back(P);
    }

    return ret;
}

static std::shared_ptr<ShellPairVec>
create_shared_multi_shellpair_(const std::vector<ShellPairBlock> & vsh,
                               const ShellVec & shell1,
                               const ShellVec & shell2)
{
    auto vec = create_multi_shellpair_(vsh, shell1, shell2);
    return std::shared_ptr<ShellPairVec>(new ShellPairVec(std::move(vec)));
}

double *SimintERI::calculate_(size_t shell1, size_t shell2,
                              size_t shell3, size_t shell4)
{
    buffer_ = allwork_;
    double* source =  source_full_;
    double* tformbuf =  tformbuf_full_;

    const size_t nsh2 = shells_[1]->size();
    const size_t nsh4 = shells_[3]->size();

    // get the precomputed simint_multi_shellpair
    const auto * P = &(*single_spairs_[0])[shell1*nsh2 + shell2];
    const auto * Q = &(*single_spairs_[1])[shell3*nsh4 + shell4];

    const bool do_cart = false; //TODO:

    // actually compute
    // if we are doing cartesian, put directly in target. Otherwise, put in source
    // and let pure_transform put it in target
    simint_compute_eri(P, Q, SIMINT_SCREEN_TOL, sharedwork_, do_cart ? buffer_ : source);
    if(do_cart)
        CartesianToSpherical_4Center(shell1, shell2, shell3, shell4, source, buffer_, tformbuf, 1);
    return buffer_;
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

    // get the basis sets from the system
    // Note - storing un-normalized
    bs_[0]=bs1; bs_[1]=bs2; bs_[2]=bs3; bs_[3]=bs4;

    // initialize the simint library
    simint_init();

    batchsize_ = 32;
    std::array<size_t,4> max_ams;
    size_t size=1.0;

    for(size_t i=0; i<4; ++i)
    {
        //all_am returns std::set of AMs sorted by value, rbegin gives us highest
        max_ams[i] = *bs_[i].all_am().rbegin();
        size*=max_ams[i];
    }
    const size_t fullsize = size * batchsize_;
    const size_t allwork_size_ = sizeof(double) * (fullsize * 2 + size);
    const size_t simint_workmem = simint_ostei_workmem(deriv, maxam_);
    sharedwork_ = (double *)SIMINT_ALLOC(simint_workmem);
    allwork_ = (double *)SIMINT_ALLOC(allwork_size_);
    source_full_ = allwork_ + fullsize;
    tformbuf_full_ = allwork_ + 2*fullsize;

    // build plain shells
    for(size_t i=0; i<4; ++i)
    {
        for(size_t j=0;j<i;++j)
        {
            if(bs_[i]==bs_[j])
            {
                shells_[i] = shells_[j];
                break;
            }
        }
        shells_[i] = create_shell_vec_(bs_[i]);
    }

    const bool braket_same = (shells_[0] == shells_[2] &&
                               shells_[1] == shells_[3]);

    single_spairs_[0] = create_shell_pair_(*shells_[0],*shells_[1]);
    single_spairs_[1] = braket_same ? single_spairs_[0] : create_shell_pair_(*shells_[2],*shells_[3]);

//    blocks12_.clear();
//    blocks34_.clear();

//    // sort the basis set AM
//    std::array<std::vector<std::vector<int>>,4> sorted_shells_;
//    for(size_t i=0; i<4; ++i)
//        sorted_shells_[i]=std::vector<std::vector<int>>(max_ams[i]+1);
//    for(size_t i=0;i<4; ++i)
//        for(const auto shelli : *shells_[i])
//            sorted_shells_[i][shelli.am].push_back(j);

//    // form pairs for the bra
//    // these aren't batched

//    for(int iam = 0; iam <= am1; iam++)
//        for(int jam = 0; jam <= am2; jam++)
//        {
//            for(int ishell : sorted_shells1[iam])
//                for(int jshell : sorted_shells2[jam])
//                {
//                    if(!bra_same_ || (bra_same_ && jshell <= ishell) )
//                        blocks12_.push_back( {{ishell, jshell}} );
//                }
//        }

//    // form pairs for the ket
//    for(int iam = 0; iam <= am3; iam++)
//        for(int jam = 0; jam <= am4; jam++)
//        {
//            ShellPairBlock curblock;

//            for(int ishell : sorted_shells3[iam])
//                for(int jshell : sorted_shells4[jam])
//                {
//                    if(!ket_same_ || (ket_same_ && jshell <= ishell) )
//                    {
//                        curblock.push_back({ishell, jshell});
//                        if(curblock.size() == batchsize_)
//                        {
//                            blocks34_.push_back(curblock);
//                            curblock.clear();
//                        }
//                    }
//                }

//            if(curblock.size())
//                blocks34_.push_back(std::move(curblock));
//        }

//    multi_spairs_[0] = create_shared_multi_shellpair_(blocks12_, *shells_[0], *shells_[4]);
//    multi_spairs_[1] = create_shared_multi_shellpair_(blocks34_, *shells_[2], *shells_[3]);
}

