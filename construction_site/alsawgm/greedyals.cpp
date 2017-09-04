#include <iostream>
#include <iomanip>
#include <htucker/htucker.h>
#include <lawa/lawa.h>
#include <vector>
#include <cmath>
#include <cstdlib>

typedef double                                  T;
typedef htucker::HTuckerTree<T>                 HTTree;
typedef htucker::HTuckerTreeNode<T>             HNode;
typedef htucker::GeneralTreeIterator
        <HNode>                                 TreeIt;
typedef htucker::GeneralTreeNode<
        htucker::HTuckerTreeNode<T> >           GNode;
typedef lawa::Basis<T, lawa::Orthogonal,
        lawa::Interval, lawa::Multi>            Basis;
typedef lawa::Index1D                           Index1D;
typedef lawa::IndexSet<Index1D>                 IndexSet;
typedef std::vector<IndexSet>                   IndexSetVec;
typedef lawa::SepCoefficients<
        lawa::Lexicographical, T, Index1D>      SepCoeff;
typedef flens::GeMatrix<flens::FullStorage
        <T, cxxblas::ColMajor> >                GeMat;
typedef flens::GeMatrix<flens::FullStorage
        <int, cxxblas::ColMajor> >              MatInt;
typedef flens::DenseVector<flens::Array<T> >    DenseVector;
typedef lawa::Coefficients<
        lawa::Lexicographical, T,
        lawa::Index1D>                          Coeff1D;
typedef flens::DenseVector<
        flens::Array<unsigned> >                IVector;
typedef lawa::Function<T>                       Function;
typedef lawa::HTCoefficients<T, Basis>          HTCoeffTree;

typedef     lawa::AdaptiveLaplaceOperator1D<T,
            lawa::Orthogonal,
            lawa::Interval,
            lawa::Multi>                        Laplace1D;
typedef     lawa::AdaptiveLaplaceOperator1D<T,
            lawa::Orthogonal,
            lawa::Interval,
            lawa::MultiRefinement>              RefLaplace1D;
typedef     lawa::LocalOperator1D<Basis,
            Basis, RefLaplace1D, Laplace1D>     LOp_Lapl1D;
typedef     lawa::Sepop<LOp_Lapl1D>             Sepop;


double
x2(double x)
{
    return x*x;
}

double
x3(double x)
{
    return x*x*x;
}

double
mycos(double x)
{
    return std::cos(x);
}

    double
mysin(double x)
{
    return std::sin(x);
}

double
myexp(double x)
{
    return std::exp(x);
}


double
zero(double)
{
    return 0.;
}


double
one(double)
{
    return 1.;
}


double
g(double x)
{
    return 9.-std::pow(2.*x-1., 2.);
    //return 1.-std::pow(2.*x-1., 2.);
}


int hashtablelength = 193;
Coeff1D
apply(LOp_Lapl1D& A, const Coeff1D& u, const IndexSet& rows)
{
    Coeff1D ret;
    lawa::TreeCoefficients1D<T> input(hashtablelength,
                                A.getTrialBasis().j0);
    lawa::TreeCoefficients1D<T> output(hashtablelength,
                                A.getTestBasis().j0);
    input = u;
    FillWithZeros(rows, ret);
    output = ret;

    A.eval(input, output, "A");

    fromTreeCoefficientsToCoefficients(output, ret);
    return ret;
}


class
DiffReaction
{
public:
    DiffReaction(const LOp_Lapl1D& _A, const T _alpha, const T _beta):
        A(_A),
        alpha(_alpha),
        beta(_beta)
    {
        assert(alpha>0.);
        assert(beta>0.);
    };

    void
    setRows(const IndexSet& _rows)
    {
        rows = _rows;
    }

    void
    setAlpha(const T _alpha)
    {
        assert(_alpha>0.);
        alpha = _alpha;
    }

    void
    setBeta(const T _beta)
    {
        assert(_beta>0.);
        beta = _beta;
    }

    LOp_Lapl1D
    getDiff()
    {
        return A;
    }

    Coeff1D
    operator()(const Coeff1D& u)
    {
        return alpha*apply(A, u, rows)+beta*u;
    }

private:
    LOp_Lapl1D A;
    IndexSet   rows;
    T          alpha;
    T          beta;
};


double
computeEnergyProduct(const Coeff1D& u, LOp_Lapl1D& A, const Coeff1D& v)
{
    IndexSet rows = supp(u);
    Coeff1D Av = apply(A, v, rows);
    return u*Av;
}


GeMat
computeR1Norms(const SepCoeff& u, LOp_Lapl1D& A, const unsigned j)
{
    assert(u.rank()==1);
    assert(j>=1 && j<=u.dim());

    unsigned d = u.dim();
    GeMat ret(2, d);
    for (unsigned i=1; i<=d; ++i) {
        if (i!=j) {
            ret(1, i) = u(1, i)*u(1, i);
            ret(2, i) = computeEnergyProduct(u(1, i), A, u(1, i));
        }
    }

    return ret;
}


void
updateR1Norms(GeMat& nrms, const Coeff1D& u, LOp_Lapl1D& A, const unsigned j)
{
    assert(j>=1 && j<=(unsigned) nrms.numCols());

    nrms(1, j) = u*u;
    nrms(2, j) = computeEnergyProduct(u, A, u);
}


double
computeAlpha(const GeMat& nrms, const unsigned j)
{
    assert(j>=1 && j<=(unsigned) nrms.numCols());
    assert(nrms.numRows()==2);

    T alpha = 1.;
    for (unsigned i=1; i<=(unsigned) nrms.numCols(); ++i) {
        if (i!=j) {
            alpha *= nrms(1, i);
        }
    }

    return alpha;
}


double
computeBeta(const GeMat& nrms, const unsigned j)
{
    assert(j>=1 && j<=(unsigned) nrms.numCols());
    assert(nrms.numRows()==2);

    T alpha = 1.;
    for (unsigned i=1; i<=(unsigned) nrms.numCols(); ++i) {
        if (i!=j) {
            alpha *= nrms(1, i);
        }
    }

    T sum = 0.;
    for (unsigned k=1; k<=(unsigned) nrms.numCols(); ++k) {
        if (k!=j) {
            sum += nrms(2, k)/nrms(1, k);
        }
    }

    return alpha*sum;
}


GeMat
computeR1RhsIntegrals(const SepCoeff& f, const SepCoeff& u,
                      const unsigned j)
{
    assert(u.rank()==1);
    assert(j>=1 && j<=u.dim());
    assert(f.dim()==u.dim());

    GeMat ret(f.rank(), u.dim());
    for (unsigned i=1; i<=u.dim(); ++i) {
        if (i!=j) {
            for (unsigned k=1; k<=f.rank(); ++k) {
                ret(k, i) = f(k, i)*u(1, i);
            }
        }
    }

    return ret;
}


void
updateR1RhsIntegrals(GeMat& ints, const SepCoeff& f, Coeff1D& u,
                     const unsigned j)
{
    assert(j>=1 && j<=f.dim());

    for (unsigned k=1; k<=f.rank(); ++k) {
        ints(k, j) = f(k, j)*u;
    }
}


Coeff1D
rhsAls(lawa::SeparableRHSD<T, Basis>& f, const IndexSet& Lambda,
       const DenseVector& factors, const unsigned j)
{
    assert(f.rank()==(unsigned) factors.length());
    assert(j>=1 && j<=f.dim());

    Coeff1D ret;
    for (unsigned k=1; k<=f.rank(); ++k) {
        ret += factors(k)*f(k, j, Lambda);
    }

    return ret;
}


class
AlsRhs
{
public:
    AlsRhs(const lawa::SeparableRHSD<T, Basis>& _f):f(_f){};

    void
    setCoeffs(const GeMat& _ints)
    {
        assert(f.rank()==(unsigned) _ints.numRows());
        assert(f.dim() ==(unsigned) _ints.numCols());
        ints = _ints;

    }

    void
    setDim(const unsigned _j)
    {
        assert(_j>=1 && _j<=f.dim());
        j = _j;

    }

    void
    updateFactors()
    {
        if (factors.length()!=ints.numRows()) factors.resize(ints.numRows());
        for (unsigned k=1; k<=f.rank(); ++k) {
            T alpha     = 1.;
            for (unsigned i=1; i<=f.dim(); ++i) {
                if (i!=j) {
                    alpha *= ints(k, i);
                }
            }
            factors(k) = alpha;
        }
        std::cout << "factors=>\n" << factors << std::endl;
    }

    Coeff1D
    operator()(const IndexSet& Lambda)
    {
        return rhsAls(f, Lambda, factors, j);
    }

private:
    lawa::SeparableRHSD<T, Basis>   f;
    GeMat                           ints;
    DenseVector                     factors;
    unsigned                        j;
};


unsigned
solve1D(      DiffReaction&                                A,
              lawa::DiagonalPreconditioner
              <lawa::DiagonalLevelPreconditioner1D<T>, T>& Prec,
              Coeff1D&                                     u,
              AlsRhs&                                      f,
              T&                                           residual,
              IndexSet&                                    Lambda,
        const Basis&                                       basis,
        const lawa::AdaptiveLeafParams&                    params)
{
    assert(Lambda.size()>0);

    // Aux
    IndexSet total = Lambda;
    IndexSet sweep = Lambda;

    // Initial residual
    A.setRows(Lambda);
    Coeff1D b = f(Lambda);
    Coeff1D r = b - A(u);
    Prec(b);
    Prec(r);
    residual  = r.norm(2.);
    T rel     = residual/b.norm(2.);

    if (params.verbose) {
        //std::cout << "solve1D: Iteration 0, r = " << residual << std::endl;
        std::cout << "solve1D: Iteration 0, r = " << rel << std::endl;
    }
    //if (residual<=params.tol) return 0.;
    if (rel<=params.tol) return 0.;

    // AWGM iterations
    T res_cg = 0.;
    for (unsigned k=1; k<=params.maxit; ++k) {
        // Galerkin solve
        Coeff1D b = f(Lambda);
        A.setRows(Lambda);
        unsigned numit = wavPcg(A, Prec, u, b, res_cg,
                                params.gamma*residual,
                                params.cg_maxit,
                                params.cg_verbose);

        if (params.verbose) {
            std::cout << "solve1D: Iteration " << k
                      << ", wavPcg required " << numit
                      << " iterations to reach " << res_cg
                      << std::endl;
        }

        // Extended index set
        extendMultiTree(basis, sweep, total, "standard", false);

        // Evaluate residual
        A.setRows(total);
        b         = f(total);
        r         = b - A(u);
        Prec(r);
        Prec(b);
        residual  = r.norm(2.);
        rel       = residual/b.norm(2.);

        // Check residual
        if (params.verbose) {
            std::cout << "solve1D: Iteration " << k
                      //<< ", r = " << residual << std::endl;
                      << ", r = " << rel << std::endl;
        }
        //if (residual<=params.tol) return k;
        if (rel<=params.tol) return k;
        if (k==params.maxit) break;

        // Bulk chasing
        sweep = bulk(Lambda, r, params.alpha, residual, res_cg, basis,
                     params.bulk_verbose);

        // Additional information
        if (params.verbose) {
            std::cout << "solve1D: Iteration " << k
                      << ", final size of active set = "  << Lambda.size()
                      << std::endl;
            std::cout << "solve1D: Iteration " << k
                      << ", final max active level   = "  << maxlevel(Lambda)
                      << std::endl;
        }
    }

    std::cerr << "solve1D: Reached max iterations " << params.maxit
              << std::endl;

    return params.maxit;
}


double
computeL2NormsProduct(const GeMat& norms, const unsigned j)
{
    assert(j>=1 && j<=(unsigned) norms.numCols());
    assert(norms.numRows()==2);

    T prod = 1.;
    for (unsigned k=1; k<=(unsigned) norms.numCols(); ++k) {
        if (k!=j) {
            prod *= norms(1, k);
        }
    }
    return std::sqrt(prod);
}


double
computeH1NormsProduct(const GeMat& norms, const unsigned j)
{
    assert(j>=1 && j<=(unsigned) norms.numCols());
    assert(norms.numRows()==2);

    T prod = 1.;
    for (unsigned k=1; k<=(unsigned) norms.numCols(); ++k) {
        if (k!=j) {
            T nrm = std::sqrt(norms(1, k) + norms(2, k));
            prod *= nrm;
        }
    }
    return prod;
}


unsigned
rank1Solve(   DiffReaction&                                A,
              lawa::DiagonalPreconditioner
              <lawa::DiagonalLevelPreconditioner1D<T>, T>& Prec,
              SepCoeff&                                    u,
              lawa::SeparableRHSD<T, Basis>                f,
              IndexSetVec&                                 Lambda,
        const Basis&                                       basis,
              lawa::Rank1AdaptiveAlsParams&                params)
{
    assert(f.dim()==u.dim());
    assert(u.rank()==1);
    assert(u.dim()==Lambda.size());

    // Initialize Rhs
    SepCoeff fcp(f.rank(), f.dim());
    for (unsigned k=1; k<=f.rank(); ++k) {
        for (unsigned j=1; j<=f.dim(); ++j) {
            fcp(k, j) = f(k, j, Lambda[j-1]);
        }
    }

    // Initialize constants
    auto  nabla   = A.getDiff();
    GeMat r1Norms = computeR1Norms(u, nabla, 1);
    GeMat rhsInts = computeR1RhsIntegrals(fcp, u, 1);
    T     alpha   = computeAlpha(r1Norms, 1);
    T     beta    = computeBeta(r1Norms, 1);
    std::cout << "eps = " << alpha/beta << std::endl;
    A.setAlpha(alpha);
    A.setBeta(beta);
    AlsRhs fals(f);
    std::cout << "rhsints =>\n" << rhsInts << std::endl;
    fals.setCoeffs(rhsInts);
    fals.setDim(1);
    fals.updateFactors();

    // ALS sweeps
    T tol = params.adaptiveLeaf.tol;
    IndexSetVec start = Lambda;
    SepCoeff    u0    = u;
    for (unsigned sweep=1; sweep<=params.max_sweep; ++sweep) {
        for (unsigned j=1; j<=u.dim(); ++j) {
            // Update constants
            if (sweep!=1 || j!=1) {
                unsigned l = j-1;
                if (j==1) l = u.dim();

                for (unsigned k=1; k<=fcp.rank(); ++k) {
                    fcp(k, l) = f(k, l, Lambda[l-1]);
                }

                updateR1Norms(r1Norms, u(1, l), nabla, l);
                updateR1RhsIntegrals(rhsInts, fcp, u(1, l), l);
                alpha = computeAlpha(r1Norms, j);
                beta  = computeBeta(r1Norms, j);
    std::cout << "eps = " << alpha/beta << std::endl;
                A.setAlpha(alpha);
                A.setBeta(beta);
                fals.setCoeffs(rhsInts);
                fals.setDim(j);
                fals.updateFactors();
            }

            // Approximate factor
            T H1  = computeH1NormsProduct(r1Norms, j);
            T eps = alpha/beta;
            //params.adaptiveLeaf.tol = tol*eps*eps/H1;
            params.adaptiveLeaf.tol = tol;
            std::cout << "Tolerance => " << params.adaptiveLeaf.tol << std::endl;
            IndexSet init = start[j-1];
            Coeff1D  v    = u0(1, j);
            T resleaf;
            unsigned numit = solve1D(A, Prec, v, fals, resleaf, init, basis,
                                     params.adaptiveLeaf);
            if (params.verbose) {
                std::cout << "rank1Solve: Sweep " << sweep
                          << ", leaf " << j
                          << ", solve1D required " << numit
                          << " iterations to reach "
                          << resleaf << std::endl;
            }

            // Update solution and index set
            u(1, j)     = v;
            Lambda[j-1] = init;
        }
    }

    std::cerr << "rank1Solve: Reached max sweeps " << params.max_sweep
              << std::endl;
    return params.max_sweep;
}


int
main(int argc, char* argv[])
{
    if (argc!=3) {
        std::cerr << "Usage: " << argv[0] << " dim level\n";
        return 1;
    }

    int dim  = atoi(argv[1]);
    int lev  = atoi(argv[2]);

    if (dim<=0) {
        std::cerr << "Dimension must be a positive integer\n";
        return 1;
    }
    int rank = 1;

    /* Generate example */
    Basis                       basis(2);
    basis.template              enforceBoundaryCondition<lawa::DirichletBC>();
    lawa::Mapwavind<Index1D>    map(dim);

    IndexSet                    indexset;
    IndexSet                    indexset2;
    IndexSet                    diff;
    getFullIndexSet(basis, indexset,  lev);
    getFullIndexSet(basis, indexset2,  lev-1);
    diff = indexset;
    for (auto& it : indexset2) {
        diff.erase(it);
    }
    IndexSetVec     indexsetvec(dim);
    IndexSetVec     indexsetvec2(dim);
    IndexSetVec     diffvec(dim);
    for (int l=0; (unsigned)l<indexsetvec.size(); ++l) {
        indexsetvec[l] = indexset;
        indexsetvec2[l] = indexset2;
        diffvec[l] = diff;
    }

    double sp = 1.;
    HTCoeffTree                     f(dim, sp, basis, map);
    HTCoeffTree                     u(dim, sp, basis, map);
    SepCoeff                        coeffs(rank, dim);
    SepCoeff                        coeffs2(1, dim);

    DenseVector             sings;
    Function                onef(one, sings);
    Function                gf(g, sings);
    Function                cosf(mycos, sings);
    Function                sinf(mysin, sings);
    Function                expf(myexp, sings);
    std::vector<Function>   fvec;

    for (int i=0; i<dim; ++i) {
        fvec.push_back(onef);
    }

    lawa::SeparableFunctionD<T>     F(fvec, rank, dim);
    GeMat                           deltas(3,2);
    std::vector<GeMat>              _deltas;
    MatInt                          derivs(rank, dim);
    for (int i=1; i<=rank; ++i) {
        for (int j=1; j<=dim; ++j) {
            derivs(i,j) = 0;
            _deltas.push_back(deltas);
        }
    }
    lawa::SeparableRHSD<T, Basis>   Fint(basis, F, _deltas, derivs);
    genCoefficients(coeffs, Fint, indexsetvec);
    set(f, coeffs);

    RefLaplace1D    RefLaplaceBil(basis.refinementbasis);
    Laplace1D       LaplaceBil(basis);
    LOp_Lapl1D      lapl(basis, basis, RefLaplaceBil, LaplaceBil);
    Sepop           A(lapl, dim, dim);
    /* ---------------------------------------------------------- */


    /* Test greedy solver */
    genCoefficientsRnd(coeffs2, indexsetvec, 1., 1);
    set(u, coeffs2);
    lawa::DiagonalLevelPreconditioner1D<T>      p;
    lawa::NoPreconditioner<T, Index1D>          pnull;
    lawa::DiagonalPreconditioner
    <lawa::DiagonalLevelPreconditioner1D<T>, T> dp(p);
    //<lawa::NoPreconditioner<T, Index1D>, T> dp(pnull);

    for (int l=0; (unsigned)l<indexsetvec.size(); ++l) {
        std::cout << "Indexset " << l+1 << " : "
                  << indexsetvec[l].size() << std::endl;
    }

    double delta = 0.5;
    lawa::Sepdiagscal<Basis>    S(u.dim(), u.basis());
    setScaling(S, delta);
    S.set_nu(1e-01);

    lawa::AdaptiveLeafParams par;
    par.tol          = 1e-03;
    par.maxit        = 50;
    par.gamma        = 1e-01;
    par.cg_maxit     = 100;
    par.cg_verbose   = false;
    par.verbose      = true;
    par.alpha        = 0.5;
    par.bulk_verbose = false;
    lawa::Rank1AdaptiveAlsParams pam;
    pam.max_sweep    = 5;
    pam.stag         = 1e-03;
    pam.verbose      = true;
    pam.adaptiveLeaf = par;

    lawa::AdaptiveGreedyParams agparams;
    agparams.r1Als   = pam;
    agparams.verbose = true;
    agparams.maxit   = 2;

    lawa::Rank1UP_Params                    p1;
    lawa::OptTTCoreParams                   p2;
    lawa::GreedyALSParams                   p3;

    /* Start MATLAB session */
//    Engine *ep;
//    if (!(ep = engOpen("matlab -nojvm"))) {
//        std::cerr << "\nCan't start MATLAB engine\n" << std::endl;
//        exit(1);
//    }

    auto start  = std::chrono::system_clock::now();
//    std::cout << "adaptive greedy required "
//              << adaptiveGreedy(ep, A, S, dp, u, indexsetvec, Fint, agparams)
//              << std::endl;

    setCoefficientsJ0(coeffs, u.basis());

    auto rhss = computeR1RhsIntegrals(coeffs, coeffs, 2);
    rhss(1, 2) = 1.;
    DiffReaction nablaR(lapl, 1., 1.);
    nablaR.setRows(indexset);
    Coeff1D w = coeffs(1, 1);
    //AlsRhs ff(Fint);
    //ff.setCoeffs(rhss);
    //ff.setDim(1);
    //(void) solve1D(nablaR, dp, w, ff, indexset, basis, par);
    (void) rank1Solve(nablaR, dp, coeffs, Fint, indexsetvec, basis, pam);
    //w = ff(indexset);

//    lawa::AgALSParams   params;
//    params.maxit              = 15;
//    params.gamma              = 1e-01;
//    params.r1update.update    = false;
//    params.r1update.sw        = true;
//    params.r1update.balance   = 500.;
//    params.r1update.orthog    = true;
//    params.r1update.tol_als   = 5e-02;
//    params.r1update.tol_cg    = 1e-08;
//    params.r1update.check_res = false;
//    params.r1update.max_sweep = 50;
//    params.r1update.verbose   = true;
//    params.r1update.maxit_cg  = 500;
//    params.greedyals.maxit    = 1;
//    params.bulk               = 0.95;
//
//    std::cout << "Solver parameters\n" << params << std::endl;
//    double residual;
//    unsigned it = agals_laplace(ep, A, S, u, Fint, indexsetvec, residual, params);
//    params.r1update.sw        = true;
//    params.r1update.balance   = 500.;
//    params.r1update.orthog    = true;
//    params.r1update.tol_als   = 1e-02;
//    params.r1update.max_sweep = 50;
//    params.r1update.maxit_cg  = 150;
//    params.coreopt.tol        = 1e-08;
//    params.coreopt.stag       = 1e-07;
//    params.greedyals.tol      = 1e-07;
//    params.greedyals.stag     = 1e-07;
//    params.greedyals.maxit    = 30;
//    auto start  = std::chrono::system_clock::now();
//    unsigned it = greedyALS_laplace(ep, A, u, f, indexsetvec, residual,
//                                    params.r1update,
//                                    params.coreopt,
//                                    params.greedyals);
    auto end     = std::chrono::system_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(end-start);
//    std::cout << "AGALS took " << it << " iterations to reach relative residual "
//            << residual << std::endl;
    std::cout << "It took " << elapsed.count() << " secs\n";
    for (unsigned k=1; k<=3; ++k) {
        std::string name = "plot_";
        name            += std::to_string(k);
        name            += ".dat";
        plot(basis, coeffs(1, k), pnull, zero, zero, 0., 1., 1./100.,
             name.c_str());
    }

//    htucker::DimensionIndex idx(1);
//    idx[0] = 1;
//    auto toplot = extract(u, idx);
//    for (unsigned k=1; k<=toplot.rank(); ++k) {
//        std::cout << "Plotting function k=" << k << std::endl;
//        //writeCoefficientsToFile(toplot(k, 1), k, "data_basisd1");
//        std::string name = "basis_functiond1_";
//        name            += std::to_string(k);
//        name            += ".dat";
//        plot(basis, toplot(k, 1), p, zero, zero, 0., 1.1, 1e-03, name.c_str());
//    }
//    std::cout << "Done...\n";

//    engClose(ep);

    return 0;
}
