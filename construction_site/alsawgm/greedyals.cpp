#include <iostream>
#include <iomanip>
#include <htucker/htucker.h>
#include <lawa/lawa.h>
#include <vector>
#include <cmath>
#include <cstdlib>

#include <fstream>
#include <utility>
#include <tuple>

#include <cholmod.h>

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
zero(double)
{
    return 0.;
}

lawa::NoPreconditioner<T, Index1D>          pnull;


double
one(double)
{
    return 1.;
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


GeMat
computeR1Norms(const SepCoeff& u, LOp_Lapl1D& A)
{
    assert(u.rank()==1);

    unsigned d = u.dim();
    GeMat ret(2, d);
    for (unsigned i=1; i<=d; ++i) {
        ret(1, i) = u(1, i)*u(1, i);
        ret(2, i) = computeEnergyProduct(u(1, i), A, u(1, i));
    }

    return ret;
}

GeMat
computeR1Norms(const SepCoeff& u, LOp_Lapl1D& A, const SepCoeff& v,
               const unsigned j)
{
    assert(u.rank()==1);
    assert(v.rank()==1);
    assert(u.dim()==v.dim());
    assert(j>=1 && j<=u.dim());

    unsigned d = u.dim();
    GeMat ret(2, d);
    for (unsigned i=1; i<=d; ++i) {
        if (i!=j) {
            ret(1, i) = u(1, i)*v(1, i);
            ret(2, i) = computeEnergyProduct(u(1, i), A, v(1, i));
        }
    }

    return ret;
}

void
updateR1Norms(GeMat& nrms, const Coeff1D& u, LOp_Lapl1D& A,
              const unsigned j)
{
    assert(j>=1 && j<=(unsigned) nrms.numCols());

    nrms(1, j) = u*u;
    nrms(2, j) = computeEnergyProduct(u, A, u);
}

void
updateR1Norms(GeMat& nrms, const Coeff1D& u, LOp_Lapl1D& A,
              const Coeff1D& v, const unsigned j)
{
    assert(j>=1 && j<=(unsigned) nrms.numCols());

    nrms(1, j) = u*v;
    nrms(2, j) = computeEnergyProduct(u, A, v);
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


template <typename PrecT, typename Residual>
unsigned
solve1D(      DiffReaction&                 A,
              PrecT&                        Prec,
              Coeff1D&                      u,
              Residual&                     f,
              T&                            rel,
              IndexSet&                     Lambda,
        const Basis&                        basis,
        const lawa::AdaptiveLeafParams&     params)
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
    T residual = r.norm(2.);
    rel        = residual/b.norm(2.);

    if (params.verbose) {
        std::cout << "solve1D: Iteration 0, r = " << rel << std::endl;
    }
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
        if (t==1) {
            int r = AU.numCols()/2;
            flens::blas::mm(cxxblas::NoTrans, cxxblas::Trans, 1.,
                            AU(_, _(r+1, 2*r)),
                            Pj(_, _(r+1, 2*r)), 0., tmp);
            Coeff1D ret1 = convert(tmp, u, Lambda, dim);
            return ret1;
        }

        if (t==2) {
            int r = AU.numCols()/2;
            flens::blas::mm(cxxblas::NoTrans, cxxblas::Trans, 1.,
                            AU(_, _(1, r)),
                            Pj(_, _(1, r)), 0., tmp);
            Coeff1D ret2 = convert(tmp, u, Lambda, dim);
            return ret2;
        }

        flens::blas::mm(cxxblas::NoTrans, cxxblas::Trans, 1.,
                        AU, Pj, 0., tmp);
        Coeff1D ret = convert(tmp, u, Lambda, dim);

        return ret;
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


class
Penalty
{
public:
    Penalty(const SepCoeff&   _wn,
            const LOp_Lapl1D& _A):
        wn(_wn),
        A(_A)
    {
        assert(wn.rank()==1);
    }

    void
    setDim(const unsigned _dim)
    {
        assert(_dim>=1 && dim<=wn.dim());

        dim      = _dim;
        upToDate = false;
    }

    void
    computeConstants(const GeMat& r1NormsWWn)
    {
        assert(dim>=1 && dim<=wn.dim());
        assert(r1NormsWWn.numRows()==2);
        assert((unsigned) r1NormsWWn.numCols()==wn.dim());

        c1    = computeAlpha(r1NormsWWn, dim);
        c2    = computeBeta(r1NormsWWn, dim);
        std::cout << "c1    = " << c1 << std::endl;
        std::cout << "c2    = " << c2 << std::endl;
        std::cout << "c1/c2 = " << c1/c2 << std::endl;

        upToDate = true;
    }

    void
    setSigma(const T _sigma)
    {
        assert(_sigma>=0.);

        sigma = _sigma;
    }

    T
    getSigma()
    {
        return sigma;
    }

    Coeff1D
    operator()(const IndexSet& Lambda)
    {
        assert(dim>=1 && dim<=wn.dim());

        if (!upToDate) {
            std::cerr << "Penalty: Warning! Constants not up to date! Exiting.\n";
            exit(1);
        }

        Coeff1D ret  = -sigma*c1*apply(A, wn(1, dim), Lambda);
        ret         += -sigma*c2*wn(1, dim);

        return ret;
    }

private:
    SepCoeff   wn;
    LOp_Lapl1D A;
    T          dim;
    T          c1;
    T          c2;
    T          sigma;
    bool       upToDate = false;
};


class
AlsPenalty
{
public:
    AlsPenalty(const AlsResidual& _r,
               const Penalty&     _pen):
        r(_r),
        pen(_pen)
    {}

    void
    setCoeffs(const GeMat& ints)
    {
        r.setCoeffs(ints);
    }

    void
    setDim(const unsigned j)
    {
        assert(j>=1);

        r.setDim(j);
        pen.setDim(j);
    }

    void
    updateR1RhsIntegrals(const SepCoeff& fcp, const Coeff1D& v,
                         const unsigned j)
    {
        assert(j>=1 && j<=fcp.dim());

        r.updateR1RhsIntegrals(fcp, v, j);
    }

    void
    computeConstants(const GeMat& r1NormsWWn)
    {
        pen.computeConstants(r1NormsWWn);
    }

    void
    setSigma(const T _sigma)
    {
        assert(_sigma>=0.);

        pen.setSigma(_sigma);
    }

    T
    getSigma()
    {
        return pen.getSigma();
    }

    Coeff1D
    operator()(const IndexSet& Lambda, const unsigned t = 0)
    {
        Coeff1D ret = r(Lambda, t);
        return ret-pen(Lambda);
    }

private:
    AlsResidual r;
    Penalty     pen;
};


template <typename Residual, typename MapType>
unsigned
rank1PoissonUpdateFixed(      DiffReaction&                  A,
                        const std::vector<SyMat>&            Astiff,
                              SepCoeff&                      u,
                              lawa::SeparableRHSD<T, Basis>& f,
                              Residual&                      r,
                              IndexSetVec&                   Lambda,
                              MapType&                       map,
                        const lawa::Rank1AdaptiveAlsParams&  params)
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

                updateR1Norms(r1Norms, u(1, l), nabla, l);
                alpha = computeAlpha(r1Norms, j);
                beta  = computeBeta(r1Norms, j);
                A.setAlpha(alpha);
                A.setBeta(beta);
                Prec.setAlpha(alpha);
                Prec.setBeta(beta);
                r.updateR1RhsIntegrals(fcp, u(1, l), l);
                r.setDim(j);
            }

            if (params.verbose) {
                std::cout << "rank1PoissonFixed: Sweep " << sweep
                          << ", leaf " << j
                          << ", epsilon = "
                          << alpha/beta << std::endl;
            }

            // Approximate factor
//            T resleaf;
//            A.setRows(Lambda[j-1]);
//            Coeff1D bh = r(Lambda[j-1]);
//            Coeff1D rh = bh - A(u(1, j));
//            Prec(rh);
//            resleaf    = rh.norm(2.);
//            T tol      = resleaf*params.adaptiveLeaf.gamma;
//            unsigned numit = wavPcg(A, Prec, u(1, j), bh, resleaf,
//                                    tol,
//                                    params.adaptiveLeaf.cg_maxit,
//                                    params.adaptiveLeaf.cg_verbose);
            unsigned numit  = 1;
            T resleaf       = 0.;
            Coeff1D bh      = r(Lambda[j-1]);
            DenseVector x   = convert(bh, j, map);
            SyMat Ah        = Astiff[j-1];
            flens::blas::scal(A.getAlpha(), Ah);
            for (int k=1; k<=Ah.numRows(); ++k) {
                Ah(k, k) += A.getBeta();
            }
            (void) flens::lapack::posv(Ah, x);
            u(1, j) = convert(x, Lambda[j-1], j, map);

           if (params.verbose) {
                std::cout << "rank1PoissonUpdateFixed: Sweep " << sweep
                          << ", leaf " << j
                          << ", wavPcg required " << numit
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
            std::cout << "rank1PoissonUpdateFixed: Sweep " << sweep
                      << ", stagnation = " << stag << std::endl;
        }
        if (stag<=params.stag) return sweep;
    }

    std::cerr << "rank1PoissonUpdateFixed: Reached max sweeps " << params.max_sweep
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


template <typename Residual>
unsigned
rank1PoissonUpdate(      DiffReaction&                  A,
                         SepCoeff&                      u,
                         lawa::SeparableRHSD<T, Basis>& f,
                         Residual&                      r,
                         IndexSetVec&                   Lambda,
                   const Basis&                         basis,
                   const lawa::Rank1AdaptiveAlsParams&  params)
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
                r.updateR1RhsIntegrals(fcp, u(1, l), l);
                r.setDim(j);
            }

            if (params.verbose) {
                std::cout << "rank1PoissonUpdate: Sweep " << sweep
                          << ", leaf " << j
                          << ", epsilon = "
                          << alpha/beta << std::endl;
            }

            // Approximate factor
            T resleaf;
            unsigned numit = solve1D(A, Prec, u(1, j), r, resleaf, Lambda[j-1],
                                     basis,
                                     params.adaptiveLeaf);
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


template <typename Residual>
unsigned
rank1PoissonUpdate2(     DiffReaction&                  A,
                         SepCoeff&                      u,
                         lawa::SeparableRHSD<T, Basis>& f,
                         Residual&                      r,
                         IndexSetVec&                   Lambda,
                   const Basis&                         basis,
                   const SepCoeff&                      wn,
                   const lawa::Rank1AdaptiveAlsParams&  params)
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
    auto  nabla    = A.getDiff();
    GeMat r1Norms  = computeR1Norms(u, nabla, 1);
    GeMat wwnNorms = computeR1Norms(u, nabla, wn, 1);
    GeMat rhsInts  = computeR1RhsIntegrals(fcp, u, 1);
    T     alpha    = computeAlpha(r1Norms, 1);
    T     beta     = computeBeta(r1Norms, 1);
    T     sigma    = r.getSigma();

    r.setDim(1);
    r.computeConstants(wwnNorms);
    alpha *= 1.+sigma;
    beta  *= 1.+sigma;
    A.setAlpha(alpha);
    A.setBeta(beta);
    PrecDiffReaction Prec(alpha, beta);
    r.setCoeffs(rhsInts);

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
                updateR1Norms(wwnNorms, u(1, l), nabla, wn(1, l), l);
                alpha = computeAlpha(r1Norms, j);
                beta  = computeBeta(r1Norms, j);
                r.setDim(j);

    r.computeConstants(wwnNorms);
    alpha *= 1.+sigma;
    beta  *= 1.+sigma;

                A.setAlpha(alpha);
                A.setBeta(beta);
                Prec.setAlpha(alpha);
                Prec.setBeta(beta);
                r.updateR1RhsIntegrals(fcp, u(1, l), l);
            }

            if (params.verbose) {
                std::cout << "rank1PoissonUpdate: Sweep " << sweep
                          << ", leaf " << j
                          << ", epsilon = "
                          << alpha/beta << std::endl;
            }

            // Approximate factor
            T resleaf;
            unsigned numit = solve1D(A, Prec, u(1, j), r, resleaf, Lambda[j-1],
                                     basis,
                                     params.adaptiveLeaf);
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
greedyPoissonSolve(      Engine                         *ep,
                         lawa::Sepdiagscal<Basis>&      S,
                         Sepop&                         A,
                         HTCoeffTree&                   u,
                         lawa::SeparableRHSD<T, Basis>& f,
                         IndexSetVec&                   Lambda,
                         lawa::AdaptiveGreedyParams&    params)
{
    assert(f.dim() ==A.dim());
    assert(f.dim() ==Lambda.size());
    assert(f.dim() ==(unsigned) u.dim());
    assert(A.type()==lawa::laplace);

    IndexSetVec sweep = Lambda;
    IndexSetVec total = Lambda;
    IndexSetVec start = Lambda;
    T relative = 1.;
    if (params.verbose) {
        std::cout << "greedyPoissonSolve: Initial,  relative residual = "
                  << relative << std::endl;
    }
    if (relative<=params.tol) {
        return 0;
    }

    // Rank 1 solve
    DiffReaction nabla(A(1, 1), 1., 1.);
    AlsRhs fals(f);
    SepCoeff v(1, f.dim());
    setCoefficientsJ0(v, u.basis());

    unsigned numit = rank1PoissonUpdate(nabla, v, f, fals, Lambda,
                                        u.basis(),
                                        params.r1Als);
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
    HTCoeffTree fScaled(u.dim(), u.basis(), u.map());
    HTCoeffTree rtree(u.dim(), u.basis(), u.map());
    ftree.tree().set_tree(u.tree());
    set(ftree, Fcp, Lambda);
    std::vector<GeMat> fred = reduce_rhs(u, ftree);
    std::vector<GeMat> Ared = reduce(A, u, Lambda, Lambda);
    std::vector<GeMat> x0(fred.size()-1);
    IVector ranks(x0.size());
    extract_core(u.tree(), x0, ranks);
    x0 = optTTcoreLaplace(ep, Ared, fred, x0, ranks, params.dmrgTTcore);
    insert_core(u.tree(), x0, ranks);

    // Residual
    T residual = 1.;
    T trunc  = std::min(1e-02, residual*1e-01);
    sweep    = Lambda;
    total    = Lambda;
    rtree    = szoneres(A, u, ftree, Fcp, f, Lambda, sweep, total);
    rtree    = applyScaleTT(S, rtree, total, trunc);
    fScaled  = applyScaleTT(S, ftree, total, trunc);
    T nrmb   = nrm2(fScaled);
    residual = nrm2(rtree);
    relative = residual/nrmb;

    if (params.verbose) {
        std::cout << "greedyPoissonSolve: Update 0, relative residual = "
                  << relative << std::endl;
    }

    // Check residual
    if (relative<=params.tol) {
        return 0;
    }

    // Greedy iterations
    saveCoeffVector1D(v(1, 1), u.basis(), "sol_0.dat");
//    IndexSetVec Lfirst = Lambda;
//    SepCoeff wn        = v;

    // Assemble fixed
//    std::vector<SyMat>  Astiff;
//    std::cout << "Assembling...\n" << std::endl;
//    for (unsigned j=1; j<=A.dim(); ++j) {
//        Astiff.push_back(assemble_projected_laplace(A, u, Lambda[j-1], j));
//    }
//    std::cout << "Done!\n" << std::endl;

    for (unsigned k=1; k<=params.maxit; ++k) {
        // Initial guess
        HTCoeffTree vtree(u.dim(), u.basis(), u.map());
        vtree.tree().set_tree(u.tree());
        setCoefficientsJ0(v, u.basis());
        //set(vtree, wn);
        set(vtree, v, Lambda);

        // ALS residual
        AlsOperatorRhs Au(A, u, vtree);
        Au.setRows(Lambda);
        Au.setCols(Lambda);
        AlsResidual rals(fals, Au);

        // Compute wn
//        params.r1Als.stag = 1e-02;
//        numit = rank1PoissonUpdateFixed(nabla, Astiff, v, f, rals, Lambda,
//                                        u.map(),
//                                        params.r1Als);
//        Penalty pen(wn, nabla.getDiff());
//
//        v = wn;
//        // Restrict
//        for (unsigned j=1; j<=v.dim(); ++j) {
//            P(start[j-1], v(1, j));
//        }
//
//        // ALS residual + penalty
//        set(vtree, v);
//        AlsOperatorRhs Au2(A, u, vtree);
//        Au2.setRows(start);
//        Au2.setCols(Lambda);
//        AlsResidual rals2(fals, Au2);
//        AlsPenalty  penals(rals2, pen);
//        penals.setSigma(100.);
//
//        // Rank 1 update
        IndexSetVec current = start;
        //numit = rank1PoissonUpdate2(nabla, v, f, penals, current,
        numit = rank1PoissonUpdate(nabla, v, f, rals, current,
                                   u.basis(),
                                   params.r1Als);
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

        std::string name = "sol_";
        name += std::to_string(k);
        name += ".dat";
        saveCoeffVector1D(v(1, 1), u.basis(), name.c_str());

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

        // Residual
        trunc    = std::min(1e-02, residual*1e-01);
        rtree    = eval(A, u, total, Lambda);
        rtree    = applyScaleTT(S, rtree, total, trunc);
        rtree.tree() = fScaled.tree()-rtree.tree();
        residual = nrm2(rtree);
        relative = residual/nrmb;

        if (params.verbose) {
            std::cout << "greedyPoissonSolve: Update " << k
                      << ", relative residual = "
                      << relative << std::endl;
        }

        // Check residual
        if (relative<=params.tol) {
            return k;
        }
    }

    std::cerr << "greedyPoissonSolve: Reached max iterations " << params.maxit
              << std::endl;

    return params.maxit;
}


unsigned
greedyPoissonSolveFixed( Engine                         *ep,
                         lawa::Sepdiagscal<Basis>&      S,
                         Sepop&                         A,
                         HTCoeffTree&                   u,
                         lawa::SeparableRHSD<T, Basis>& f,
                         IndexSetVec&                   Lambda,
                         lawa::AdaptiveGreedyParams&    params)
{
    assert(f.dim() ==A.dim());
    assert(f.dim() ==Lambda.size());
    assert(f.dim() ==(unsigned) u.dim());
    assert(A.type()==lawa::laplace);

    T relative = 1.;
    T residual = 1.;
    if (params.verbose) {
        std::cout << "greedyPoissonSolveFixed: Initial,  relative residual = "
                  << relative << std::endl;
    }
    if (relative<=params.tol) {
        return 0;
    }

    // Assemble fixed
    std::vector<SyMat>  Astiff;
    for (unsigned j=1; j<=A.dim(); ++j) {
        Astiff.push_back(assemble_projected_laplace(A, u, Lambda[j-1], j));
    }

    // Rank 1 solve
    DiffReaction nabla(A(1, 1), 1., 1.);
    AlsRhs fals(f);
    SepCoeff v(1, f.dim());
    setCoefficientsJ0(v, u.basis());

    unsigned numit = rank1PoissonUpdateFixed(nabla, Astiff, v, f, fals, Lambda,
                                             u.map(),
                                             params.r1Als);
    if (params.verbose) {
        std::cout << "greedyPoissonSolveFixed: Update 0 required " << numit
                  << " sweeps\n";
    }

    // Scale core
    R1L2Normalize(v);
    set(u, v);

    SepCoeff Fcp(f.rank(), f.dim());
    genCoefficients(Fcp, f, Lambda);
    HTCoeffTree ftree(u.dim(), u.basis(), u.map());
    HTCoeffTree fScaled(u.dim(), u.basis(), u.map());
    HTCoeffTree rtree(u.dim(), u.basis(), u.map());
    ftree.tree().set_tree(u.tree());
    set(ftree, Fcp, Lambda);
    std::vector<GeMat> fred = reduce_rhs(u, ftree);
    std::vector<GeMat> Ared = reduce(A, u, Lambda, Lambda);
    std::vector<GeMat> x0(fred.size()-1);
    IVector ranks(x0.size());
    extract_core(u.tree(), x0, ranks);
    x0 = optTTcoreLaplace(ep, Ared, fred, x0, ranks, params.dmrgTTcore);
    insert_core(u.tree(), x0, ranks);

    // Residual
    T trunc      = 1e-08;
    rtree        = eval(A, u, Lambda, Lambda);
    rtree        = applyScaleTT(S, rtree, Lambda, trunc);
    fScaled      = applyScaleTT(S, ftree, Lambda, trunc);
    rtree.tree() = fScaled.tree()-rtree.tree();
    T nrmb       = nrm2(fScaled);
    residual     = nrm2(rtree);
    relative     = residual/nrmb;

    if (params.verbose) {
        std::cout << "greedyPoissonSolveFixed: Update 0, relative residual = "
                  << relative << std::endl;
    }

    // Check residual
    if (residual<=params.tol) {
        return 0;
    }

    // Greedy iterations
    for (unsigned k=1; k<=params.maxit; ++k) {
        // Initial guess
        HTCoeffTree vtree(u.dim(), u.basis(), u.map());
        vtree.tree().set_tree(u.tree());
        set(vtree, v);

        // ALS residual
        AlsOperatorRhs Au(A, u, vtree);
        Au.setRows(Lambda);
        Au.setCols(Lambda);
        AlsResidual rals(fals, Au);

        // Compute update
        numit = rank1PoissonUpdateFixed(nabla, Astiff, v, f, rals, Lambda,
                                        u.map(),
                                        params.r1Als);

        if (params.verbose) {
            std::cout << "greedyPoissonSolveFixed: Update " << k <<" required "
                      << numit
                      << " sweeps\n";
        }

        // Update basis
        R1L2Normalize(v);
        set(vtree, v);
        u.tree() = u.tree()+vtree.tree();
        u.tree().orthogonalize();

        // Optimize core
        fred = reduce_rhs(u, ftree);
        Ared = reduce(A, u, Lambda, Lambda);
        x0.clear();
        x0.resize(fred.size()-1);
        ranks.resize(x0.size());
        extract_core(u.tree(), x0, ranks);
        x0 = optTTcoreLaplace(ep, Ared, fred, x0, ranks, params.dmrgTTcore);
        insert_core(u.tree(), x0, ranks);

        // Residual
        rtree        = eval(A, u, Lambda, Lambda);
        rtree        = applyScaleTT(S, rtree, Lambda, trunc);
        rtree.tree() = fScaled.tree()-rtree.tree();
        residual     = nrm2(rtree);
        relative     = residual/nrmb;

        if (params.verbose) {
            std::cout << "greedyPoissonSolve: Update " << k
                      << ", relative residual = "
                      << relative << std::endl;
        }

        // Check residual
        if (residual<=params.tol) {
            return k;
        }
    }

    std::cerr << "greedyPoissonSolve: Reached max iterations " << params.maxit
              << std::endl;

    return params.maxit;
}


void
awgm(Engine                         *ep,
     lawa::Sepdiagscal<Basis>&      S,
     Sepop&                         A,
     HTCoeffTree&                   u,
     lawa::SeparableRHSD<T, Basis>& f,
     IndexSetVec&                   Lambda,
     lawa::AdaptiveGreedyParams&    params,
     T                              alpha,
     T                              gamma,
     unsigned                       maxit)
{
    assert(f.dim() ==A.dim());
    assert(f.dim() ==Lambda.size());
    assert(f.dim() ==(unsigned) u.dim());
    assert(A.type()==lawa::laplace);

    T rel      = 1.;
    T residual = 1.;
    std::cout << "awgm: Initial,  relative residual = " << rel << std::endl;

    IndexSetVec sweep = Lambda;
    IndexSetVec total = Lambda;
    for (unsigned k=1; k<=maxit; ++k) {
        // Galerkin solve
        T tol = residual*gamma;
        params.tol = tol;
        unsigned numit = greedyPoissonSolveFixed(ep, S, A, u, f,
                                                 Lambda, params);
        std::cout << "awgm: Iteration " << k << ", greedy required "
                  << numit << " updates to reach " << tol << std::endl;

        // Estimate residual
        T trunc = std::min(1e-02, residual*1e-01);
        SepCoeff Fcp(f.rank(), f.dim());
        genCoefficients(Fcp, f, total);
        HTCoeffTree ftree(u.dim(), u.basis(), u.map());
        HTCoeffTree rtree(u.dim(), u.basis(), u.map());
        ftree.tree().set_tree(u.tree());
        rtree    = szoneres(A, u, ftree, Fcp, f, Lambda, sweep, total);
        rtree    = applyScaleTT(S, rtree, total, trunc);
        ftree    = applyScaleTT(S, ftree, total, trunc);
        T nrmb   = nrm2(ftree);
        residual = nrm2(rtree);
        rel      = residual/nrmb;

        std::cout << "awgm: Iteration " << k << ", relative residual = "
                  << rel << std::endl;

        // Compute bulk
        if (k<maxit) {
            sweep = bulk(alpha, residual, rtree, Lambda, sweep);
        }
    }
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


void
writeCoordMatrix(      std::ofstream&                         file,
                 const std::vector<std::tuple<int, int, T> >  A,
                 const int                                    nrow,
                 const int                                    ncol,
                 const int                                    nnz,
                 const int                                    stype)
{
    assert((unsigned) nnz==A.size());
    assert(file.is_open());

    file << nrow << " " << ncol << " " << nnz << " " << stype << std::endl;

    for (const auto& it : A) {
        file << std::get<0>(it) << " " << std::get<1>(it) << " "
             << std::get<2>(it) << std::endl;
    }
}


template <typename MapType>
void
assemble_laplace_coord_toFile(      std::ofstream&  file,
                                    LOp_Lapl1D&     a,
                                    MapType&        map,
                              const IndexSet&       Lambda,
                              const unsigned        dim)
{
    assert(dim>=1 && dim<=map.dim());

    int nrow   = maxintindhash(Lambda, dim, map);
    int ncol   = nrow;
    int nnz    = 0;
    int stype  = -1;

    std::vector<std::tuple<int, int, T> > storage;
    for (const auto& mu1 : Lambda) {
        int i = map(mu1, dim);

        for (const auto& mu2: Lambda) {
            int j = map(mu2, dim);
            if (j>i) continue;

            T val = a(mu1, mu2);
            if (val==0.) continue;

            std::tuple<int, int, T> entry(i, j, val);
            storage.push_back(entry);
            ++nnz;
        }
    }

    writeCoordMatrix(file, storage, nrow, ncol, nnz, stype);
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
    Basis                       basis(4);
    basis.template              enforceBoundaryCondition<lawa::DirichletBC>();
    lawa::Mapwavind<Index1D>    map(dim);

    IndexSet                    indexset;
    getFullIndexSet(basis, indexset,  lev);
    IndexSetVec     indexsetvec(dim);
    for (int l=0; (unsigned)l<indexsetvec.size(); ++l) {
        indexsetvec[l] = indexset;
    }

    for (int l=0; (unsigned)l<indexsetvec.size(); ++l) {
        std::cout << "Indexset " << l+1 << " : "
                  << indexsetvec[l].size() << std::endl;
    }

    double sp = 1.;
    HTCoeffTree                     f(dim, sp, basis, map);
    HTCoeffTree                     u(dim, sp, basis, map);
    SepCoeff                        coeffs(rank, dim);

    DenseVector            sings;
    Function               onef(one, sings);
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

    // Assemble
//    SyMat Astiff = assemble_laplace(lapl, map, indexset, 1);
//    for (int k=1; k<=Astiff.numRows(); ++k) {
//        Astiff(k, k) += 1e-04;
//    }
//
//    DenseVector b(Astiff.numRows());
//    for (int i=1; i<=b.length(); ++i) {
//        b(i) = 1.;
//    }
//
//    // Solve dense
//    auto t0  = std::chrono::system_clock::now();
//    int info = flens::lapack::posv(Astiff, b);
//    auto t1  = std::chrono::system_clock::now();
//    std::cout << "info = " << info << std::endl;
//    auto dt  = std::chrono::duration_cast<
//               std::chrono::seconds>(t1-t0);
//    std::cout << "time = " << dt.count() << std::endl;
//
//    // Solve sparse
//    std::ofstream file("testMatrix.txt");
//    assemble_laplace_coord_toFile(file, lapl, map, indexset, 1);
//
//    cholmod_dense *x, *bs, *r;
//    cholmod_sparse *Asparse, *I, *C;
//    cholmod_factor *L;
//    cholmod_common c;
//    cholmod_start(&c);
//    T one_ [2] = {1, 0}, m1 [2] = {-1, 0};
//
//    FILE *foobar = fopen("testMatrix.txt", "r");
//    Asparse = cholmod_read_sparse(foobar, &c);
//    bs = cholmod_ones(Asparse->nrow, 1, Asparse->xtype, &c);
//
//    T alpha[2] = {1., 0.};
//    T beta [2] = {1., 0.};
//    I = cholmod_speye(Asparse->nrow, Asparse->nrow, Asparse->xtype, &c);
//    C = cholmod_add(Asparse, I, alpha, beta, 1, 0, &c);
//    std::cout << "Stype => " << Asparse->stype << std::endl;
//    std::cout << "Stype => " << I->stype << std::endl;
//    std::cout << "Stype => " << C->stype << std::endl;
//
//    t0  = std::chrono::system_clock::now();
//    L = cholmod_analyze(Asparse, &c);
//    cholmod_factorize(Asparse, L, &c);
//    x = cholmod_solve(CHOLMOD_A, L, bs, &c);
//    t1  = std::chrono::system_clock::now();
//    auto dt2  = std::chrono::duration_cast<
//                std::chrono::seconds>(t1-t0);
//
//    r = cholmod_copy_dense(bs, &c);
//    cholmod_sdmult(Asparse, 0, m1, one_, x, r, &c);
//    std::cout << "Residual = " << cholmod_norm_dense(r, 0, &c) << std::endl;
//
//    std::cout << "time = " << dt2.count() << std::endl;
//    //std::cout << "scaling = " << dt.count()/dt2.count() << std::endl;
//
//    cholmod_free_factor(&L, &c);
//    cholmod_free_sparse(&Asparse, &c);
//    cholmod_free_sparse(&C, &c);
//    cholmod_free_dense(&r, &c);
//    cholmod_free_dense(&x, &c);
//    cholmod_free_dense(&bs, &c);
//    cholmod_finish(&c);
//    exit(1);

    /* Test greedy solver */
    lawa::DiagonalLevelPreconditioner1D<T>      p;
    lawa::DiagonalPreconditioner
    <lawa::DiagonalLevelPreconditioner1D<T>, T> dp(p);

    T delta = 1e-01;
    T eta   = 1e-01;
    lawa::Sepdiagscal<Basis> S(u.dim(), u.basis(), map);
    setScaling(S, delta);
    S.set_nu(eta);

    lawa::AdaptiveLeafParams par;
    par.tol          = 1e-02;
    par.maxit        = 50;
    par.gamma        = 1e-01;
    par.cg_maxit     = 50;
    par.cg_verbose   = false;
    par.verbose      = false;
    par.alpha        = 0.5;
    par.bulk_verbose = false;
    lawa::Rank1AdaptiveAlsParams pam;
    pam.max_sweep    = 100;
    pam.stag         = 1e-02;
    pam.verbose      = false;
    pam.adaptiveLeaf = par;

    lawa::AdaptiveGreedyParams agparams;
    agparams.r1Als   = pam;
    agparams.verbose = true;
    agparams.maxit   = 8;
    agparams.tol     = 1e-07;

    lawa::Rank1UP_Params                    p1;
    lawa::OptTTCoreParams                   p2;
    lawa::GreedyALSParams                   p3;

    /* Start MATLAB session */
    Engine *ep;
    if (!(ep = engOpen("matlab -nojvm"))) {
        std::cerr << "\nCan't start MATLAB engine\n" << std::endl;
        exit(1);
    }

    T balpha = 0.5;
    T bgamma = 1e-01;
    auto start  = std::chrono::system_clock::now();
    (void) greedyPoissonSolve(ep, S, A, u, Fint, indexsetvec, agparams);
//    awgm(ep, S, A, u, Fint, indexsetvec, agparams,
//         balpha, bgamma, 100);
    auto end     = std::chrono::system_clock::now();

    htucker::DimensionIndex idx(1);
    idx[0] = 1;
    coeffs = extract(u, indexsetvec[0], idx);

    auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(end-start);
    std::cout << "It took " << elapsed.count() << " secs\n";

    saveCoeffVector1D(coeffs(1, 1), u.basis(), "support_d1.dat");
    for (unsigned k=1; k<=coeffs.rank(); ++k) {
        std::string name = "plot_";
        name            += std::to_string(k);
        name            += ".dat";
        plot(basis, coeffs(k, 1), pnull, zero, zero, 0., 1.02, 1./100.,
             name.c_str());
    }

    engClose(ep);

    return 0;
}
