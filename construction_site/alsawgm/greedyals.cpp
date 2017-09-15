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

typedef flens::SyMatrix<flens::FullStorage
        <T, cxxblas::ColMajor> >                SyMat;
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
typedef     lawa::DiagonalPreconditioner
            <lawa::DiagonalLevelPreconditioner1D<T>, T> DiagonalPrec;


void
saveCoeffVector1D(const Coeff1D &coeff, const Basis &basis, const char* filename)
{
    typedef typename Coeff1D::const_iterator const_coeff_it;

    std::ofstream data(filename);

    if(!data.good()){
        std::cerr << "File "
                  << filename
                  << " could not be opened for writing" << std::endl;
        exit(1);
    }
    data.precision(40);
    data << "# Center_x Value Xtype j k" << std::endl;

    for (const_coeff_it it = coeff.begin(); it != coeff.end(); ++it) {
      FLENS_DEFAULT_INDEXTYPE j=(*it).first.j,
                              k=(*it).first.k;
      lawa::XType type=(*it).first.xtype;

      //center of the support
      double x = 0.5*(basis.generator(type).support(j,k).l2 +
                      basis.generator(type).support(j,k).l1);

      data << x << " " << std::scientific <<  (*it).second << " "
           << type << " " << j << " " << k << std::endl;
    }
    data.close();
}


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

lawa::NoPreconditioner<T, Index1D>          pnull;

double
one(double x)
{
//    if (x<=0.5) return 0.;
    (void) x
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
        assert(beta>=0.);
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
        assert(_beta>=0.);
        beta = _beta;
    }

    LOp_Lapl1D
    getDiff()
    {
        return A;
    }

    T
    getAlpha() const
    {
        return alpha;
    }

    T
    getBeta() const
    {
        return beta;
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


class
PrecDiffReaction
{
public:
    PrecDiffReaction(const T _alpha, const T _beta):
        alpha(_alpha),
        beta(_beta)
    {
        assert(alpha>0.);
        assert(beta>=0.);
    };

    void
    setAlpha(const T _alpha)
    {
        assert(_alpha>0.);
        alpha = _alpha;
    }

    void
    setBeta(const T _beta)
    {
        assert(_beta>=0.);
        beta = _beta;
    }

    T
    operator()(const Index1D& lambda) const
    {
        int j = lambda.j;
        if (lambda.xtype == lawa::XWavelet) ++j;
        T p      = lawa::pow2i<T>(j);
        T peclet = 1./(alpha*p);
        T tau;
        peclet = 0.;
        if (peclet>1.) {
            tau = 1./p;
        } else {
            tau = 1./(alpha*p*p/beta);
        }
        T ret = std::sqrt(1./(1.+alpha*p*p/beta)+tau);
        return ret;
    }

    void
    operator()(Coeff1D& v) const
    {
        for (auto& it : v) {
            it.second *= this->operator()(it.first);
        }
    }

    void
    remove(Coeff1D& v) const
    {
        for (auto& it : v) {
            it.second /= this->operator()(it.first);
        }
    }
private:
    T alpha;
    T beta;
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
updateR1RhsIntegrals_(GeMat& ints, const SepCoeff& f, const Coeff1D& u,
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
        upToDate = false;
    }

    void
    setDim(const unsigned _j)
    {
        assert(_j>=1 && _j<=f.dim());
        j = _j;
        upToDate = false;
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
        upToDate = true;
    }

    void
    updateR1RhsIntegrals(const SepCoeff& fcp, const Coeff1D& v,
                         const unsigned j)
    {
        assert(fcp.rank()==f.rank());
        assert(fcp.dim() ==f.dim());
        assert(j>=1 && j<=f.dim());

        updateR1RhsIntegrals_(ints, fcp, v, j);
        upToDate = false;
    }

    Coeff1D
    operator()(const IndexSet& Lambda, const unsigned t = 0)
    {
        (void) t;
        if (!upToDate) updateFactors();
        return rhsAls(f, Lambda, factors, j);
    }

private:
    lawa::SeparableRHSD<T, Basis>   f;
    GeMat                           ints;
    DenseVector                     factors;
    unsigned                        j;
    bool                            upToDate = false;
};


T
artDiffusion(const T alpha, const T beta,
             const IndexSet& Lambda, const Basis& basis)
{
    int level = maxlevel(Lambda)-1;
    int k     = 1;
    lawa::XType type = lawa::XWavelet;
    T h       = basis.generator(type).support(level, k).length();
    T q       = beta*h/(2.*alpha);
    return q/std::tanh(q);
}


template <typename PrecT, typename Residual, typename MapType>
unsigned
solve1D(      DiffReaction&                 A,
              PrecT&                        Prec,
              Coeff1D&                      u,
              Residual&                     f,
              T&                            rel,
              IndexSet&                     Lambda,
        const Basis&                        basis,
              MapType&                      map,
        const unsigned                      swp,
        const unsigned                      j,
        const unsigned                      up,
        const lawa::AdaptiveLeafParams&     params)
{
    assert(Lambda.size()>0);

    // Temp
//    if (up==4 && swp==1) {
//        auto l = maxlevel(Lambda);
//        IndexSet newset;
//        getFullIndexSet(basis, newset, l+1);
//        Lambda = newset;
//    }

    // Aux
    IndexSet total = Lambda;
    IndexSet sweep = Lambda;

    // Initial residual
    A.setRows(Lambda);
    Coeff1D b = f(Lambda);
    Coeff1D r = b - A(u);
    Prec(b);
    Prec(r);
    T residual = r.norm(2.);
    rel        = residual/b.norm(2.);

    if (params.verbose) {
        std::cout << "solve1D: Iteration 0, r = " << rel << std::endl;
    }
    if (rel<=params.tol) return 0.;

    // AWGM iterations
    T res_cg = 0.;
//    std::string name1 = "suppRes_swp_";
//    name1 += std::to_string(swp);
//    name1 += "_dim_";
//    name1 += std::to_string(j);
//    name1 += "_it_0.dat";
//    saveCoeffVector1D(r, basis, name1.c_str());

//  if (up==4 && swp==1) {
//        auto l = maxlevel(Lambda);
//    for (int i=0; i<=0; ++i) {
//    std::string name2 = "f_up_";
//    name2 += std::to_string(up);
//    name2 += "_swp_";
//    name2 += std::to_string(swp);
//    name2 += "_dim_";
//    name2 += std::to_string(j);
//    name2 += "_it_0.dat";
//    Coeff1D b1 = f(Lambda);
//    plot(basis, b1, pnull, zero, zero, 0., 1.02, 1e-02,
//         name2.c_str());
//
//  std::string name3 = "Au1_up_";
//  name3 += std::to_string(up);
//  name3 += "_swp_";
//  name3 += std::to_string(swp);
//  name3 += "_dim_";
//  name3 += std::to_string(j);
//  name3 += "_it_";
//  name3 += std::to_string(0);
//  name3 += ".dat";
//
//  Coeff1D b2 = f(Lambda, 2);
//  plot(basis, b2, pnull, zero, zero, 0., 1.02, 1e-02,
//       name3.c_str());
//  }
//  }
//    std::string name4 = "Au2_up_";
//    name4 += std::to_string(up);
//    name4 += "_swp_";
//    name4 += std::to_string(swp);
//    name4 += "_dim_";
//    name4 += std::to_string(j);
//    name4 += "_it_0.dat";
//    Coeff1D b3 = f(Lambda, 3);
//    plot(basis, b3, pnull, zero, zero, 0., 1.02, 1e-02,
//         name4.c_str());

    for (unsigned k=1; k<=params.maxit; ++k) {
        // Galerkin solve
        Coeff1D b = f(Lambda);
        A.setRows(Lambda);
        T tol = params.gamma*residual;
        if (up >= 1 && up<=4) tol = b.norm(2.)*1e-08;
        unsigned numit = wavPcg(A, Prec, u, b, res_cg,
                                //params.gamma*residual,
                                tol,
                                params.cg_maxit,
                                params.cg_verbose);
        // Solve exactly once
//        std::cout << "Assembling problem...\n";
//        auto nabla = A.getDiff();
//        SyMat Astiff = assemble_laplace(nabla, map,
//                                        Lambda, j);
//        Coeff1D rhs   = f(Lambda);
//        DenseVector x = convert(rhs, j, map);
//        flens::blas::scal(A.getAlpha(), Astiff);
//        for (int k=1; k<=Astiff.numRows(); ++k) {
//            Astiff(k, k) += A.getBeta();
//        }
//        std::cout << "Assembled!\n";
//        std::cout << "Solving...\n";
//        int info = flens::lapack::posv(Astiff, x);
//        std::cout << "info = " << info << std::endl;
//        std::cout << "Solved!\n";
//        std::cout << "Converting...\n";
//        u = convert(x, Lambda, j, map);
//        std::cout << "Converted!\n";

       // unsigned numit = 1;
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

 //   name1 = "suppRes_swp_";
 //   name1 += std::to_string(swp);
 //   name1 += "_dim_";
 //   name1 += std::to_string(j);
 //   name1 += "_it_";
 //   name1 += std::to_string(k);
 //   name1 += ".dat";
 //   saveCoeffVector1D(r, basis, name1.c_str());

//    name2 = "f_up_";
//    name2 += std::to_string(up);
//    name2 += "_swp_";
//    name2 += std::to_string(swp);
//    name2 += "_dim_";
//    name2 += std::to_string(j);
//    name2 += "_it_";
//    name2 += std::to_string(k);
//    name2 += ".dat";
//    b1 = f(total, 1);
//    plot(basis, b1, pnull, zero, zero, 0., 1.02, 1e-02,
//         name2.c_str());
//
//    name3 = "Au1_up_";
//    name3 += std::to_string(up);
//    name3 += "_swp_";
//    name3 += std::to_string(swp);
//    name3 += "_dim_";
//    name3 += std::to_string(j);
//    name3 += "_it_";
//    name3 += std::to_string(k);
//    name3 += ".dat";
//    b2 = f(total, 2);
//    plot(basis, b2, pnull, zero, zero, 0., 1.02, 1e-02,
//         name3.c_str());
//
//    name4 = "Au2_up_";
//    name4 += std::to_string(up);
//    name4 += "_swp_";
//    name4 += std::to_string(swp);
//    name4 += "_dim_";
//    name4 += std::to_string(j);
//    name4 += "_it_";
//    name4 += std::to_string(k);
//    name4 += ".dat";
//    b3 = f(total, 3);
//    plot(basis, b3, pnull, zero, zero, 0., 1.02, 1e-02,
//         name4.c_str());

        // Check residual
        if (params.verbose) {
            std::cout << "solve1D: Iteration " << k
                      << ", r = " << rel << std::endl;
        }
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


DenseVector
computeR1L2Prods(const SepCoeff& u, const SepCoeff& v)
{
    assert(u.dim()==v.dim());
    assert(u.rank()==1);
    assert(v.rank()==1);

    DenseVector L2Prods(u.dim());
    for (unsigned j=1; j<=u.dim(); ++j) {
        L2Prods(j) = u(1, j)*v(1, j);
    }

    return L2Prods;
}


DenseVector
computeR1H1Prods(const SepCoeff& u, LOp_Lapl1D& A, const SepCoeff& v)
{
    assert(u.dim()==v.dim());
    assert(u.rank()==1);
    assert(v.rank()==1);

    DenseVector H1Prods(u.dim());
    for (unsigned j=1; j<=u.dim(); ++j) {
        H1Prods(j) = computeEnergyProduct(u(1, j), A, v(1, j));
    }

    return H1Prods;
}


double
computeR1H1InnerProduct(const SepCoeff& u, LOp_Lapl1D& A, const SepCoeff& v)
{
    assert(u.dim()==v.dim());
    assert(u.rank()==1);
    assert(v.rank()==1);

    DenseVector L2Prods = computeR1L2Prods(u, v);
    DenseVector H1Prods = computeR1H1Prods(u, A, v);
    T L2Prod = 1.;
    for (unsigned j=1; j<=u.dim(); ++j) {
        L2Prod    *= L2Prods(j);
    }

    T sum = 0.;
    for (unsigned k=1; k<=u.dim(); ++k) {
        T prod = 1.;
        for (unsigned j=1; j<=u.dim(); ++j) {
            if (j!=k) {
                prod *= L2Prods(j);
            }
        }

        sum += H1Prods(k)*prod;
    }

    return L2Prod + sum;
}


double
computeR1H1Difference(const SepCoeff& u, LOp_Lapl1D& A, const SepCoeff& v)
{
    assert(u.dim()==v.dim());
    assert(u.rank()==1);
    assert(v.rank()==1);

    T uu, uv, vv, err2;
    uu   = computeR1H1InnerProduct(u, A, u);
    uv   = computeR1H1InnerProduct(u, A, v);
    vv   = computeR1H1InnerProduct(v, A, v);
    err2 = uu+vv-2.*uv;
    if (err2<0.) {
        std::cerr << "computeR1H1Difference: Warning, squared error = "
                  << err2 << ", setting to 0\n";
        return 0.;
    }

    return std::sqrt(err2);
}


class
AlsOperatorRhs
{
public:
    AlsOperatorRhs(const Sepop&       _A,
                   const HTCoeffTree& _u,
                   const HTCoeffTree& _v):
        A(_A),
        u(_u),
        v(_v)
    {
        assert(A.dim()==(unsigned) u.dim());
        assert(A.dim()==(unsigned) v.dim());
    }

    void
    setRows(const IndexSetVec& _rows)
    {
        assert(_rows.size()==(unsigned) u.dim());

        rows     = _rows;
        upToDate = false;
    }

    void
    setCols(const IndexSetVec& _cols)
    {
        assert(_cols.size()==(unsigned) u.dim());

        cols     = _cols;
        upToDate = false;
    }

    void
    setDim(const unsigned _dim)
    {
        assert(_dim>=1 && _dim<=(unsigned) u.dim());

        dim      = _dim;
        upToDate = false;
    }

    void
    computeProjection()
    {
        assert(rows.size()>0);
        assert(cols.size()>0);
        assert(rows.size()==cols.size());
        assert(rows.size()==(unsigned) u.dim());

        HTCoeffTree Au = evalD_1Dim(A, u, rows, cols, dim);
        Pj             = projection(Au.tree(), v.tree(), dim);
        assert(Pj.numRows()==1);
        upToDate       = true;
    }

    void
    update(const Coeff1D& vj, const unsigned j)
    {
        assert(j>=1 && j<=(unsigned) u.dim());
        assert(cols.size()==(unsigned) u.dim());

        rows[j-1] = supp(vj);
        insert(v, vj, rows[j-1], j);
        upToDate = false;
    }

    Coeff1D
    operator()(const IndexSet& Lambda, const unsigned t = 0)
    {
        using flens::_;
        if (!upToDate) computeProjection();

        GeMat AU = eval1Dim(A, u, Lambda, cols[dim-1], dim);

        assert(AU.numCols()==Pj.numCols());
        assert(Pj.numRows()==1);

            GeMat tmp;
            int r = AU.numCols()/2;
            flens::blas::mm(cxxblas::NoTrans, cxxblas::Trans, 1.,
                            AU(_, _(r+1, 2*r)),
                            Pj(_, _(r+1, 2*r)), 0., tmp);
            Coeff1D ret1 = convert(tmp, u, Lambda, dim);
       //     P(ret1, cols[dim-1]);
        if (t==1) {
            return ret1;
        }

            r = AU.numCols()/2;
            flens::blas::mm(cxxblas::NoTrans, cxxblas::Trans, 1.,
                            AU(_, _(1, r)),
                            Pj(_, _(1, r)), 0., tmp);
            Coeff1D ret2 = convert(tmp, u, Lambda, dim);
        if (t==2) {
            return ret2;
        }

     //   GeMat tmp;
     //   flens::blas::mm(cxxblas::NoTrans, cxxblas::Trans, 1., AU, Pj, 0., tmp);
     //   Coeff1D ret = convert(tmp, u, Lambda, dim);

        return ret1+ret2;
    }

private:
    Sepop       A;
    HTCoeffTree u;
    HTCoeffTree v;
    IndexSetVec rows;
    IndexSetVec cols;
    unsigned    dim;
    GeMat       Pj;
    bool        upToDate = false;
};


class
AlsResidual
{
public:
    AlsResidual(const AlsRhs&         _f,
                const AlsOperatorRhs& _A):
        f(_f),
        A(_A)
    {}

    void
    setCoeffs(const GeMat& ints)
    {
        f.setCoeffs(ints);
    }

    void
    setDim(const unsigned j)
    {
        assert(j>=1);

        f.setDim(j);
        f.updateFactors();

        A.setDim(j);
    }

    void
    updateR1RhsIntegrals(const SepCoeff& fcp, const Coeff1D& v,
                         const unsigned j)
    {
        assert(j>=1 && j<=fcp.dim());

        f.updateR1RhsIntegrals(fcp, v, j);
        A.update(v, j);
    }

    Coeff1D
    operator()(const IndexSet& Lambda, const unsigned t = 0)
    {
        if (t==1) return f(Lambda);
        if (t==2) return A(Lambda, 1);
        if (t==3) return A(Lambda, 2);
        return f(Lambda)-A(Lambda);
    }

private:
    AlsRhs         f;
    AlsOperatorRhs A;
};


template <typename Residual, typename MapType>
unsigned
rank1PoissonUpdate(      DiffReaction&                                A,
                         SepCoeff&                                    u,
                         lawa::SeparableRHSD<T, Basis>&               f,
                         Residual&                                    r,
                         IndexSetVec&                                 Lambda,
                   const Basis&                                       basis,
                         MapType&                                     map,
                   const unsigned                                     up,
                   const lawa::Rank1AdaptiveAlsParams&                params)
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
    A.setAlpha(alpha);
    A.setBeta(beta);
    PrecDiffReaction Prec(alpha, beta);
    std::cout << "epsilon = " << alpha/beta << std::endl;
    r.setCoeffs(rhsInts);
    r.setDim(1);

    // ALS sweeps
    for (unsigned sweep=1; sweep<=params.max_sweep; ++sweep) {
        SepCoeff uold = u;
        for (unsigned j=1; j<=u.dim(); ++j) {
            // Update constants
            if (sweep!=1 || j!=1) {
                unsigned l = j-1;
                if (j==1) l = u.dim();

                for (unsigned k=1; k<=fcp.rank(); ++k) {
                    fcp(k, l) = f(k, l, Lambda[l-1]);
                }

                updateR1Norms(r1Norms, u(1, l), nabla, l);
                alpha = computeAlpha(r1Norms, j);
                beta  = computeBeta(r1Norms, j);
                A.setAlpha(alpha);
                A.setBeta(beta);
                Prec.setAlpha(alpha);
                Prec.setBeta(beta);
                std::cout << "epsilon = " << alpha/beta << std::endl;
                r.updateR1RhsIntegrals(fcp, u(1, l), l);
                r.setDim(j);
            }

            // Add artificial diffusion
            //T sigma = artDiffusion(alpha, beta, Lambda[j-1], basis);
            T sigma = 1.;
            T adiff = alpha*sigma;
            A.setAlpha(adiff);
            Prec.setAlpha(adiff);
            std::cout << "New epsilon => " << adiff/beta << std::endl;

            // Approximate factor
            T resleaf;
            T tol = params.adaptiveLeaf.tol;
            //if (sweep==1) tol = 1e-02;
            auto par = params.adaptiveLeaf;
            par.tol = tol;
            unsigned numit = solve1D(A, Prec, u(1, j), r, resleaf, Lambda[j-1],
                                     basis, map, sweep, j, up,
                                     par);
          //  Prec.remove(u(1, j));
          //  u(1, j) = THRESH(u(1, j), tol*u(1, j).norm(2.), false, false);
          //  Prec(u(1, j));
          //  Lambda[j-1] = supp(u(1, j));
            if (params.verbose) {
                std::cout << "rank1PoissonUpdate: Sweep " << sweep
                          << ", leaf " << j
                          << ", solve1D required " << numit
                          << " iterations to reach "
                          << resleaf << std::endl;
            }

            // Solve exactly once
           // if (sweep == params.max_sweep-1) {
           //     std::cout << "Assembling problem...\n";
           //     SyMat Astiff = assemble_laplace(nabla, map,
           //                                     Lambda[j-1], j);
           //     Coeff1D b     = r(Lambda[j-1]);
           //     DenseVector x = convert(b, j, map);
           //     flens::blas::scal(alpha, Astiff);
           //     for (int k=1; k<=Astiff.numRows(); ++k) {
           //         Astiff(k, k) += beta;
           //     }
           //     std::cout << "Assembled!\n";
           //     std::cout << "Solving...\n";
           //     int info = flens::lapack::posv(Astiff, x);
           //     std::cout << "info = " << info << std::endl;
           //     std::cout << "Solved!\n";
           //     std::cout << "Converting...\n";
           //     u(1, j) = convert(x, Lambda[j-1], j, map);
           //     std::cout << "Converted!\n";
           // }

        }

        // Check for stagnation
        T h1norm = std::sqrt(computeR1H1InnerProduct(u, nabla, u));
        T h1diff = computeR1H1Difference(u, nabla, uold);
        T stag;
        if (h1diff==0.) stag = 0.;
        else            stag = h1diff/h1norm;
        if (params.verbose) {
            std::cout << "rank1PoissonUpdate: Sweep " << sweep
                      << ", stagnation = " << stag << std::endl;
        }
        if (stag<=params.stag) return sweep;
    }

    std::cerr << "rank1PoissonUpdate: Reached max sweeps " << params.max_sweep
              << std::endl;
    return params.max_sweep;
}


void
R1L2Normalize(SepCoeff& u)
{
    assert(u.rank()==1);

    for (unsigned j=1; j<=u.dim(); ++j) {
        T scale  = u(1, j).norm(2.);
        u(1, j) *= 1./scale;
    }
}



template <typename Residual, typename MapType>
unsigned
rank1PoissonUpdate2(      DiffReaction&                                A,
                         SepCoeff&                                    u,
                         lawa::SeparableRHSD<T, Basis>&               f,
                         Residual&                                    r,
                         IndexSetVec&                                 Lambda,
                   const Basis&                                       basis,
                         MapType&                                     map,
                   const unsigned                                     up,
                   const lawa::Rank1AdaptiveAlsParams&                params)
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
    A.setAlpha(alpha);
    A.setBeta(beta);
    PrecDiffReaction Prec(alpha, beta);
    std::cout << "epsilon = " << alpha/beta << std::endl;
    r.setCoeffs(rhsInts);
    r.setDim(1);

    // ALS sweeps
    for (unsigned sweep=1; sweep<=params.max_sweep; ++sweep) {
        SepCoeff uold = u;
        for (unsigned j=1; j<=u.dim(); ++j) {
            // Update constants
            if (sweep!=1 || j!=1) {
                unsigned l = j-1;
                if (j==1) l = u.dim();

                for (unsigned k=1; k<=fcp.rank(); ++k) {
                    fcp(k, l) = f(k, l, Lambda[l-1]);
                }

                updateR1Norms(r1Norms, u(1, l), nabla, l);
                alpha = computeAlpha(r1Norms, j);
                beta  = computeBeta(r1Norms, j);
                A.setAlpha(alpha);
                A.setBeta(beta);
                Prec.setAlpha(alpha);
                Prec.setBeta(beta);
                std::cout << "epsilon = " << alpha/beta << std::endl;
                r.updateR1RhsIntegrals(fcp, u(1, l), l);
                r.setDim(j);
            }

            // Add artificial diffusion
            //T sigma = artDiffusion(alpha, beta, Lambda[j-1], basis);
            T sigma = 1.;
            T adiff = alpha*sigma;
            A.setAlpha(adiff);
            Prec.setAlpha(adiff);
            std::cout << "New epsilon => " << adiff/beta << std::endl;

            // Approximate factor
            T resleaf;
            T tol     = params.adaptiveLeaf.tol;
            auto par  = params.adaptiveLeaf;
            par.tol   = 0.;
            par.maxit = 2.;
            unsigned numit = solve1D(A, Prec, u(1, j), r, resleaf, Lambda[j-1],
                                     basis, map, sweep, j, up,
                                     par);
           if (params.verbose) {
                std::cout << "rank1PoissonUpdate: Sweep " << sweep
                          << ", leaf " << j
                          << ", solve1D required " << numit
                          << " iterations to reach "
                          << resleaf << std::endl;
            }
        }

        // Check for stagnation
        T h1norm = std::sqrt(computeR1H1InnerProduct(u, nabla, u));
        T h1diff = computeR1H1Difference(u, nabla, uold);
        T stag;
        if (h1diff==0.) stag = 0.;
        else            stag = h1diff/h1norm;
        if (params.verbose) {
            std::cout << "rank1PoissonUpdate: Sweep " << sweep
                      << ", stagnation = " << stag << std::endl;
        }
        if (stag<=params.stag) return sweep;
    }

    std::cerr << "rank1PoissonUpdate: Reached max sweeps " << params.max_sweep
              << std::endl;
    return params.max_sweep;
}






void
R1H1Normalize(SepCoeff& u, LOp_Lapl1D& A)
{
    assert(u.rank()==1);

    for (unsigned j=1; j<=u.dim(); ++j) {
        T scale  = std::sqrt(computeEnergyProduct(u(1, j), A, u(1, j)));
        u(1, j) *= 1./scale;
    }
}


unsigned
greedyPoissonSolve(      Engine                                      *ep,
                         Sepop&                                       A,
                         HTCoeffTree&                                 u,
                         lawa::SeparableRHSD<T, Basis>&               f,
                         IndexSetVec&                                 Lambda,
                         lawa::AdaptiveGreedyParams&                  params)
{
    assert(f.dim() ==A.dim());
    assert(f.dim() ==Lambda.size());
    assert(f.dim() ==(unsigned) u.dim());
    assert(A.type()==lawa::laplace);

    IndexSetVec start   = Lambda;

    // Rank 1 solve
    DiffReaction nabla(A(1, 1), 1., 1.);
    AlsRhs fals(f);
    SepCoeff v(1, f.dim());
    setCoefficientsJ0(v, u.basis());

    params.r1Als.max_sweep    = 100;
    unsigned numit = rank1PoissonUpdate2(nabla, v, f, fals, Lambda,
                                         u.basis(), u.map(), 0, params.r1Als);
    // Fix index set
    params.r1Als.adaptiveLeaf.tol      = 1e-08;
    params.r1Als.adaptiveLeaf.maxit    = 1;
    params.r1Als.adaptiveLeaf.cg_maxit = 300;

    if (params.verbose) {
        std::cout << "greedyPoissonSolve: Update 0 required " << numit
                  << " sweeps\n";
    }

    // Scale core
    R1L2Normalize(v);
    set(u, v);

    SepCoeff Fcp(f.rank(), f.dim());
    genCoefficients(Fcp, f, Lambda);
    HTCoeffTree ftree(u.dim(), u.basis(), u.map());
    ftree.tree().set_tree(u.tree());
    set(ftree, Fcp, Lambda);
    std::vector<GeMat> fred = reduce_rhs(u, ftree);
    std::vector<GeMat> Ared = reduce(A, u, Lambda, Lambda);
    std::vector<GeMat> x0(fred.size()-1);
    IVector ranks(x0.size());
    extract_core(u.tree(), x0, ranks);
    x0 = optTTcoreLaplace(ep, Ared, fred, x0, ranks, params.dmrgTTcore);
    insert_core(u.tree(), x0, ranks);

    // Greedy iterations
    saveCoeffVector1D(v(1, 1), u.basis(), "sol_0.dat");
    for (unsigned k=1; k<=params.maxit; ++k) {
        // Initial guess
        HTCoeffTree vtree(u.dim(), u.basis(), u.map());
        vtree.tree().set_tree(u.tree());
        setCoefficientsJ0(v, u.basis());
        set(vtree, v);

        // ALS residual
        AlsOperatorRhs Au(A, u, vtree);
        Au.setRows(start);
        Au.setCols(Lambda);
        AlsResidual rals(fals, Au);

        // Rank 1 update
        IndexSetVec current = Lambda;
        numit = rank1PoissonUpdate(nabla, v, f, rals, current,
                                     u.basis(), u.map(), k, params.r1Als);
        std::string name = "sol_";
        name += std::to_string(k);
        name += ".dat";
        saveCoeffVector1D(v(1, 1), u.basis(), name.c_str());
        Lambda = unify(Lambda, current);
        if (params.verbose) {
            std::cout << "greedyPoissonSolve: Update " << k <<" required "
                      << numit
                      << " sweeps\n";
        }

        // Update basis
        R1L2Normalize(v);
        set(vtree, v);
        u.tree() = u.tree()+vtree.tree();
        u.tree().orthogonalize();

        // Optimize core
        genCoefficients(Fcp, f, Lambda);
        set(ftree, Fcp, Lambda);
        fred = reduce_rhs(u, ftree);
        Ared = reduce(A, u, Lambda, Lambda);
        x0.clear();
        x0.resize(fred.size()-1);
        ranks.resize(x0.size());
        extract_core(u.tree(), x0, ranks);
        x0 = optTTcoreLaplace(ep, Ared, fred, x0, ranks, params.dmrgTTcore);
        insert_core(u.tree(), x0, ranks);
    }

    std::cerr << "greedyPoissonSolve: Reached max iterations " << params.maxit
              << std::endl;

    return params.maxit;
}


class
Wrapper
{
public:
    Wrapper(const lawa::RHSWithPeaks1D<T, Basis>& _f):f(_f){};

    Coeff1D
    operator()(const IndexSet& Lambda) const
    {
        Coeff1D ret;
        for (const auto& it : Lambda) {
            ret[it] = f(it);
        }

        return ret;
    }

private:
    lawa::RHSWithPeaks1D<T, Basis> f;
};


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

    GeMat peaks(9, 2);
    peaks(1, 1) = 0.1;
    peaks(2, 1) = 0.2;
    peaks(3, 1) = 0.3;
    peaks(4, 1) = 0.4;
    peaks(5, 1) = 0.5;
    peaks(6, 1) = 0.6;
    peaks(7, 1) = 0.7;
    peaks(8, 1) = 0.8;
    peaks(9, 1) = 0.9;
    peaks(1, 2) = -1.;
    peaks(2, 2) = 1.;
    peaks(3, 2) = -1.;
    peaks(4, 2) = 1.;
    peaks(5, 2) = -1.;
    peaks(6, 2) = 1.;
    peaks(7, 2) = -1.;
    peaks(8, 2) = 1.;
    peaks(9, 2) = -1.;

    lawa::RHSWithPeaks1D<T, Basis> diracs(basis, onef, peaks, 4, true, false);

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
    //lawa::NoPreconditioner<T, Index1D>          pnull;
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
    par.tol          = 1e-02;
    par.maxit        = 50;
    par.gamma        = 1e-01;
    par.cg_maxit     = 50;
    par.cg_verbose   = false;
    par.verbose      = true;
    par.alpha        = 0.5;
    par.bulk_verbose = false;
    lawa::Rank1AdaptiveAlsParams pam;
    pam.max_sweep    = 10;
    pam.stag         = 1e-03;
    pam.verbose      = true;
    pam.adaptiveLeaf = par;

    lawa::AdaptiveGreedyParams agparams;
    agparams.r1Als   = pam;
    agparams.verbose = true;
    agparams.maxit   = 4;

    lawa::Rank1UP_Params                    p1;
    lawa::OptTTCoreParams                   p2;
    lawa::GreedyALSParams                   p3;

    /* Start MATLAB session */
    Engine *ep;
    if (!(ep = engOpen("matlab -nojvm"))) {
        std::cerr << "\nCan't start MATLAB engine\n" << std::endl;
        exit(1);
    }

//    std::cout << "adaptive greedy required "
//              << adaptiveGreedy(ep, A, S, dp, u, indexsetvec, Fint, agparams)
//              << std::endl;

    setCoefficientsJ0(coeffs, u.basis());
//    Wrapper wrap(diracs);
//
//    DiffReaction diffusion(lapl, 1e-03, 1.);
//    T omg;
//    (void) solve1D(diffusion, dp, coeffs(1, 1), wrap, omg, indexset, basis,
//                   map, 1, 1, 1, par);
//        plot(basis, coeffs(1, 1), pnull, zero, zero, 0., 1.02, 1./100.,
//             "test.dat");
//    exit(1);

    auto rhss = computeR1RhsIntegrals(coeffs, coeffs, 2);
    rhss(1, 2) = 1.;
    DiffReaction nablaR(lapl, 1., 1.);
    nablaR.setRows(indexset);
    //Coeff1D w = coeffs(1, 1);
    //AlsRhs ff(Fint);
    //ff.setCoeffs(rhss);
    //ff.setDim(1);
    //(void) solve1D(nablaR, dp, w, ff, indexset, basis, par);

    auto start  = std::chrono::system_clock::now();
    (void) greedyPoissonSolve(ep, A, u, Fint, indexsetvec, agparams);
    auto end     = std::chrono::system_clock::now();
    htucker::DimensionIndex idx(1);
    idx[0] = 1;
    coeffs = extract(u, indexsetvec[0], idx);
    //    AlsRhs fals(Fint);
//    auto copy = indexsetvec;
//    (void) rank1PoissonUpdate(nablaR, dp, coeffs, Fint, fals,
//                              indexsetvec, basis, pam);
//    R1L2Normalize(coeffs);
//    set(u, coeffs);
//
//    SepCoeff Fcp(Fint.rank(), Fint.dim());
//    genCoefficients(Fcp, Fint, indexsetvec);
//    set(f, Fcp, indexsetvec);
//    std::vector<GeMat> fred = reduce_rhs(u, f);
//    std::vector<GeMat> Ared = reduce(A, u, indexsetvec, indexsetvec);
//    std::vector<GeMat> x0(fred.size()-1);
//    IVector ranks(x0.size());
//    extract_core(u.tree(), x0, ranks);
//    lawa::OptTTCoreParams dmrgParams;
//    x0 = optTTcoreLaplace(ep, Ared, fred, x0, ranks, dmrgParams);
//    insert_core(u.tree(), x0, ranks);
//
//    setCoefficientsJ0(coeffs, u.basis());
//    HTCoeffTree w(dim, sp, basis, map);
//    set(w, coeffs);
//
//    AlsOperatorRhs Au(A, u, w);
//    Au.setRows(copy);
//    Au.setCols(indexsetvec);
//
//    AlsResidual r(fals, Au);
//    (void) rank1PoissonUpdate(nablaR, dp, coeffs, Fint, r,
//                              copy, basis, pam);
//    R1L2Normalize(coeffs);
    //R1H1Normalize(coeffs, lapl);
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
    auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(end-start);
//    std::cout << "AGALS took " << it << " iterations to reach relative residual "
//            << residual << std::endl;
    std::cout << "It took " << elapsed.count() << " secs\n";
    for (unsigned k=1; k<=coeffs.rank(); ++k) {
        std::string name = "plot_";
        name            += std::to_string(k);
        name            += ".dat";
        plot(basis, coeffs(k, 1), pnull, zero, zero, 0., 1.02, 1./100.,
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

    engClose(ep);

    return 0;
}
