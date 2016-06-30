/*
  This file is part of LAWA - Library for Adaptive Wavelet Applications.
  Copyright (C) 2008-2011  Mario Rometsch, Alexander Stippler.

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 */

#include <cassert>
#include <cmath>
#include <list>

#include <lawa/aux/arrow.h>
#include <lawa/math/lawa_math.h>
#include <lawa/constructions/realline/dual/bspline.h>
#include <extensions/flens/lapack_flens.h>

namespace lawa {

template <typename T>
MRA<T,Dual,Interval,Primbs>::MRA(FLENS_DEFAULT_INDEXTYPE _d, FLENS_DEFAULT_INDEXTYPE _d_, FLENS_DEFAULT_INDEXTYPE j)
    : d(_d), d_(_d_), mu(d&1),
      min_j0(iceil<FLENS_DEFAULT_INDEXTYPE>(log(d+2*d_-3)/log(2))+1),
      j0((j==-1) ? min_j0 : j), phi_R(d,d_), phi_(),
      _bc(2,0), _j(j0)
{
    assert(d>1);
    assert(_j>=min_j0);

    _calcM0_();
}

template <typename T>
MRA<T,Dual,Interval,Primbs>::~MRA()
{
}

//--- cardinalities of whole, left, inner, right index sets. -------------------

template <typename T>
FLENS_DEFAULT_INDEXTYPE
MRA<T,Dual,Interval,Primbs>::cardI_(FLENS_DEFAULT_INDEXTYPE j) const
{
    assert(j>=min_j0);
    return pow2i<T>(j) + d - 1 - (_bc(0)+_bc(1));
}

template <typename T>
FLENS_DEFAULT_INDEXTYPE
MRA<T,Dual,Interval,Primbs>::cardI_L(FLENS_DEFAULT_INDEXTYPE /*j*/) const
{
    return d + d_ - 2 -_bc(0);
}

template <typename T>
FLENS_DEFAULT_INDEXTYPE
MRA<T,Dual,Interval,Primbs>::cardI_I(FLENS_DEFAULT_INDEXTYPE j) const
{
    assert(j>=min_j0);
    return pow2i<T>(j) + d - 1 - 2*(d+d_-2);
}

template <typename T>
FLENS_DEFAULT_INDEXTYPE
MRA<T,Dual,Interval,Primbs>::cardI_R(FLENS_DEFAULT_INDEXTYPE /*j*/) const
{
    return d + d_ - 2 - _bc(1);
}

//--- ranges of whole, left, inner, right index sets. --------------------------

template <typename T>
flens::Range<FLENS_DEFAULT_INDEXTYPE>
MRA<T,Dual,Interval,Primbs>::rangeI_(FLENS_DEFAULT_INDEXTYPE j) const
{
    assert(j>=min_j0);
    return flens::Range<FLENS_DEFAULT_INDEXTYPE>(1 + _bc(0), pow2i<T>(j) + d - 1 - _bc(1));
}

template <typename T>
flens::Range<FLENS_DEFAULT_INDEXTYPE>
MRA<T,Dual,Interval,Primbs>::rangeI_L(FLENS_DEFAULT_INDEXTYPE /*j*/) const
{
    return flens::Range<FLENS_DEFAULT_INDEXTYPE>(1 + _bc(0), d + d_ - 2);
}

template <typename T>
flens::Range<FLENS_DEFAULT_INDEXTYPE>
MRA<T,Dual,Interval,Primbs>::rangeI_I(FLENS_DEFAULT_INDEXTYPE j) const
{
    assert(j>=min_j0);
    return flens::Range<FLENS_DEFAULT_INDEXTYPE>(d+d_-1, pow2i<T>(j)+d-1-(d+d_-2));
}

template <typename T>
flens::Range<FLENS_DEFAULT_INDEXTYPE>
MRA<T,Dual,Interval,Primbs>::rangeI_R(FLENS_DEFAULT_INDEXTYPE j) const
{
    assert(j>=min_j0);
    return flens::Range<FLENS_DEFAULT_INDEXTYPE>(pow2i<T>(j)+d-1-(d+d_-2)+1, pow2i<T>(j)+d-1-_bc(1));
}

template <typename T>
FLENS_DEFAULT_INDEXTYPE
MRA<T,Dual,Interval,Primbs>::level() const
{
    return _j;
}

template <typename T>
void
MRA<T,Dual,Interval,Primbs>::setLevel(FLENS_DEFAULT_INDEXTYPE j) const
{
    if (j!=_j) {
        assert(j>=min_j0);
        _j = j;
        M0_.setLevel(_j);
    }
}

template <typename T>
template <BoundaryCondition BC>
void
MRA<T,Dual,Interval,Primbs>::enforceBoundaryCondition()
{
    assert(BC==DirichletBC);

    _bc(0) = DirichletBC;
    _bc(1) = DirichletBC;

    _twoScaleDual_2(d, d_);

    flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > A(d+d_-2,d+d_-2);
    for (FLENS_DEFAULT_INDEXTYPE n=1; n<=d+d_-2; ++n) {
        for (FLENS_DEFAULT_INDEXTYPE m=1; m<=d+d_-2; ++m) {
            A(n,m) = MDD(m,n);
        }
    }

    for (FLENS_DEFAULT_INDEXTYPE n=1; n<=d+d_-2; ++n) {
        A(n,n) = A(n,n)-1;
    }
    flens::DenseVector<flens::Array<T> > O(d+d_-2), x, tmp;
    A(d+d_-2,d+d_-2) = 1.;
    O(d+d_-2) = 1.;

    x = inv(A)*O;
    _orthonormalize(MP,MDD);
    tmp = x;
    x = transpose(inv(TT))*tmp;
    x.engine().changeIndexBase(1);
    flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > TR(d+d_-2,d+d_-2);
    for (FLENS_DEFAULT_INDEXTYPE n=1; n<=d+d_-2; ++n) {
        TR(n,n) = 1;
        if (n>1) {
            TR(n,1) = -x(n) / x(1);
        }
    }

    flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > TRD(3*d_+2*d-5,3*d_+2*d-5), MDDD_H, Tmp, InvTR;
    TRD.diag(0) = 1.;

//    TRD(_(1,d+d_-2),_(1,d+d_-2)) = transpose(inv(TR));
    InvTR = inv(TR);
    for (FLENS_DEFAULT_INDEXTYPE r=1; r<=d+d_-2; ++r) {
        for (FLENS_DEFAULT_INDEXTYPE c=1; c<=d+d_-2; ++c) {
            TRD(r,c) = InvTR(c,r);
        }
    }

    flens::blas::mm(cxxblas::NoTrans,cxxblas::NoTrans,1.,TRD,MDDD,0.,Tmp);
    flens::blas::mm(cxxblas::NoTrans,cxxblas::Trans,1.,Tmp,TR,0.,MDDD_H);

    T factor = Const<T>::R_SQRT2;

    flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > Mj0_(pow2i<T>(min_j0+1)+d-3,
                                            pow2i<T>(min_j0)+d-3, 
                                            2,2),
                                       Mj0_Right, Bugfix;
    Bugfix = MDDD_H(_(2,MDDD.numRows()), _(2,MDDD.numCols()));
    flens::blas::scal(factor,Bugfix);
    Mj0_(_(2,MDDD.numRows()), _(2,MDDD.numCols())) = Bugfix;
    FLENS_DEFAULT_INDEXTYPE row = d+d_-1;
    for (FLENS_DEFAULT_INDEXTYPE c=MDDD.lastCol()+1; c<=Mj0_.lastCol()-MDDD.numCols()+1; ++c,row+=2) {
        Mj0_(_(row,row+phi_R.a_.length()-1),c) = factor*phi_R.a_;
    }
    arrow(Bugfix,Mj0_Right);
    Mj0_(_(Mj0_.lastRow()-Mj0_Right.numRows()+1,Mj0_.lastRow()),
         _(Mj0_.lastCol()-Mj0_Right.numCols()+1,Mj0_.lastCol())) = Mj0_Right;
    M0_ = flens::RefinementMatrix<T,Interval,Primbs>(Mj0_Right.numCols(),Mj0_Right.numCols(), 
                                              Mj0_, min_j0, min_j0);
    M0_.setLevel(_j);
}

template <typename T>
void
MRA<T,Dual,Interval,Primbs>::_orthonormalize(const flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > &L3,
                                             const flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > &LL3)
{
    FLENS_DEFAULT_INDEXTYPE n = L3.numRows();
    FLENS_DEFAULT_INDEXTYPE m = L3.numCols();

    flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > L2(n-m,m), LL2(n-m,m), L22, C(m,m);

    for (FLENS_DEFAULT_INDEXTYPE j=1; j<=n-m; ++j) {
        for (FLENS_DEFAULT_INDEXTYPE k=1; k<=m; ++k) {
            L2(j,k) = L3(j+m,k);
            LL2(j,k) = LL3(j+m,k);
        }
    }

    // L22 = transpose(L2);
    L22.engine().resize(L2.cols(),L2.rows());
    for (FLENS_DEFAULT_INDEXTYPE r=L2.firstRow(); r<=L2.lastRow(); ++r) {
        for (FLENS_DEFAULT_INDEXTYPE c=L2.firstCol(); c<=L2.lastCol(); ++c) {
            L22(c,r) = L2(r,c);
        }
    }
    
    flens::blas::mm(cxxblas::NoTrans,cxxblas::NoTrans,0.5,L22,LL2,0.,C);
    flens::DenseVector<flens::Array<T> > c(powii(m,2));

    for (FLENS_DEFAULT_INDEXTYPE j=1; j<=m; ++j) {
        for (FLENS_DEFAULT_INDEXTYPE k=1; k<=m; ++k) {
            c((j-1)*m+k) = C(j,k);
        }
    }
    flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > B(powii(m,2),powii(m,2));

    for (FLENS_DEFAULT_INDEXTYPE j=1; j<=m; ++j) {
        for (FLENS_DEFAULT_INDEXTYPE k=1; k<=m; ++k) {
            for (FLENS_DEFAULT_INDEXTYPE l=1; l<=m; ++l) {
                for (FLENS_DEFAULT_INDEXTYPE s=1; s<=m; ++s) {
                    B((j-1)*m+k,(l-1)*m+s) = -0.5*L3(l,j)*LL3(s,k);
                }
            }
        }
    }

    for (FLENS_DEFAULT_INDEXTYPE j=1; j<=powii(m,2); ++j) {
        B(j,j) += 1;
    }

    flens::DenseVector<flens::Array<T> > t(powii(m,2)), cc, dd;
    bool singular;
    B.engine().changeIndexBase(0,0);
    c.engine().changeIndexBase(0);
    qrf(B,cc,dd,singular);
    qrsolv(B,cc,dd,c);
    t = c;
    t.engine().changeIndexBase(1);
    TT.engine().resize(m,m);

    for (FLENS_DEFAULT_INDEXTYPE j=1; j<=m; ++j) {
        for (FLENS_DEFAULT_INDEXTYPE k=1; k<=m; ++k) {
            TT(j,k) = t((j-1)*m+k);
        }
    }
}

//------------------------------------------------------------------------------

template <typename T>
T
_dividedDifferences(flens::DenseVector<flens::Array<T> > &a, T x)
{
    FLENS_DEFAULT_INDEXTYPE m = a.length();
    flens::DenseVector<flens::Array<T> > b(m);

    for (FLENS_DEFAULT_INDEXTYPE j=1; j<=m; ++j) {
        b(j) = std::pow(std::max(a(j)-x,T(0)),m-2);
    }

    for (FLENS_DEFAULT_INDEXTYPE j=1; j<=m-1; ++j) {
        for (FLENS_DEFAULT_INDEXTYPE k=1; k<=m-j; ++k) {
            if (a(k)==a(k+j)) {
                b(k) = 0;
            } else {
                b(k) = (b(k+1)-b(k)) / (a(k+j)-a(k));
            }
        }
    }
    return b(1);
}

template <typename T>
void
MRA<T,Dual,Interval,Primbs>::_Umrechnung1(FLENS_DEFAULT_INDEXTYPE d)
{
    flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > B(d-1,d-1), C(d-1,d-1);
    FLENS_DEFAULT_INDEXTYPE m = d+1;
    flens::DenseVector<flens::Array<T> > b(m), c(m);
    for (FLENS_DEFAULT_INDEXTYPE j=1; j<=m; ++j) {
        c(j) = j-m;
    }
    for (FLENS_DEFAULT_INDEXTYPE l=1; l<=d-1; ++l) {
        for (FLENS_DEFAULT_INDEXTYPE j=1; j<=m-1; ++j) {
            b(j) = b(j+1); c(j) = c(j)+1;
        }
        ++b(m); ++c(m);
        for (FLENS_DEFAULT_INDEXTYPE i=1; i<=d-1; ++i) {
            B(i,l) = l*_dividedDifferences(b,T(i-1));
            C(i,l) = d*_dividedDifferences(c,T(i-1));
        }
    }
    flens::blas::mm(cxxblas::NoTrans,cxxblas::NoTrans,1.,inv(C),B,0.,W);
}

template <typename T>
void
MRA<T,Dual,Interval,Primbs>::_Umrechnung2(FLENS_DEFAULT_INDEXTYPE d)
{
    flens::DenseVector<flens::Array<T> > a(2*d+1);
    flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > Tmp, InvW, B(2*(d-1),2*(d-1));
    Mj.engine().resize(2*(d-1),d-1);

    for (FLENS_DEFAULT_INDEXTYPE j=1; j<=d+1; ++j) {
        a(2*d+2-j) = pow2i<T>(1-d)*binomial(d,j-1);
    }

    for (FLENS_DEFAULT_INDEXTYPE j=1; j<=d-1; ++j) {
        for (FLENS_DEFAULT_INDEXTYPE k=1; k<=2*j; ++k) {
            Mj(k,j) = a(2*d+1-2*j+k);
        }
    }

    _Umrechnung1(d);
    InvW = inv(W);

    for (FLENS_DEFAULT_INDEXTYPE j=1; j<=d-1; ++j) {
        for (FLENS_DEFAULT_INDEXTYPE k=1; k<=d-1; ++k) {
            B(j,k) = InvW(j,k);
        }
    }

    for (FLENS_DEFAULT_INDEXTYPE j=d; j<=2*(d-1); ++j) {
        B(j,j) = 1;
    }
    flens::blas::mm(cxxblas::NoTrans,cxxblas::NoTrans,1.,B,Mj,0.,Tmp);
    flens::blas::mm(cxxblas::NoTrans,cxxblas::NoTrans,1.,Tmp,W,0.,Mp);
}

template <typename T>
void
MRA<T,Dual,Interval,Primbs>::_twoScalePrimal(FLENS_DEFAULT_INDEXTYPE d, FLENS_DEFAULT_INDEXTYPE d_)
{
    MP.engine().resize(3*d_+2*d-5, d+d_-2);
    _Umrechnung2(d);

    for (FLENS_DEFAULT_INDEXTYPE j=1; j<=2*d-2; ++j) {
        for (FLENS_DEFAULT_INDEXTYPE k=1; k<=d-1; ++k) {
            MP(j,k) = Mp(j,k);
        }
    }

    flens::DenseVector<flens::Array<T> > b = BSpline<T,Primal,R,CDF>(d).a;
    b.engine().changeIndexBase(1);
    for (FLENS_DEFAULT_INDEXTYPE k=d; k<=d+d_-2; ++k) {
        for (FLENS_DEFAULT_INDEXTYPE j=2*k-d; j<=2*k; ++j) {
            MP(j,k) = b(j-2*k+d+1);
        }
    }
}

template <typename T>
void
MRA<T,Dual,Interval,Primbs>::_twoScaleDual_1(FLENS_DEFAULT_INDEXTYPE d, FLENS_DEFAULT_INDEXTYPE d_)
{
    M1.engine().resize(3*d_+2*d-5,d+d_-2) || M1.engine().fill();
    MD.engine().resize(3*d_+2*d-5,d+d_-2) || MD.engine().fill();
    flens::DenseVector<flens::Array<T> > b = BSpline<T,Primal,R,CDF>(d).a;
    flens::DenseVector<flens::Array<T> > a = BSpline<T,Dual,R,CDF>(d,d_).a_;
    b.engine().changeIndexBase(1); a.engine().changeIndexBase(1);

    for (FLENS_DEFAULT_INDEXTYPE l=d/2+d_-1; l<=d/2+3*d_+d-5; ++l) {
        FLENS_DEFAULT_INDEXTYPE gu = std::max((FLENS_DEFAULT_INDEXTYPE) 0,iceil<FLENS_DEFAULT_INDEXTYPE>((d/2+d_-3-l)/2.));
        FLENS_DEFAULT_INDEXTYPE go = std::min(2*d_+d-4, (d+d/2+3*d_-5-l)/2);
        for (FLENS_DEFAULT_INDEXTYPE n=gu; n<=go; ++n) {
            M1(l+iceil<FLENS_DEFAULT_INDEXTYPE>(d/2.),d+d_-2) += a(l+2*n-d/2-d_+4);
        }
        for (FLENS_DEFAULT_INDEXTYPE k=1; k<=d_-1; ++k) {
            for (FLENS_DEFAULT_INDEXTYPE n=gu; n<=go; ++n) {
                M1(l+iceil<FLENS_DEFAULT_INDEXTYPE>(d/2.),d+d_-2-k) += powii(n,k)*a(l+2*n-d/2-d_+4);
            }
        }
    }

    flens::DenseVector<flens::Array<T> > m(d_);
    for (FLENS_DEFAULT_INDEXTYPE n=-1; n<=d+2*d_-3; ++n) {
        m(1) += a(n+2);
    }
    for (FLENS_DEFAULT_INDEXTYPE j=2; j<=d_; ++j) {
        for (FLENS_DEFAULT_INDEXTYPE n=-1; n<=d+2*d_-3; ++n) {
            m(j) += a(n+2)*powii(n,j-1);
        }
    }
    m *= .5;

    for (FLENS_DEFAULT_INDEXTYPE k=0; k<=d_-1; ++k) {
        for (FLENS_DEFAULT_INDEXTYPE r=0; r<=k; ++r) {
            M1(d+d_-2-r,d+d_-2-k) = pow2i<T>(-k)*binomial(k,r)*m(k-r+1);
        }
    }

    flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > D3(d_,d_);
    D3(1,1) = 1;
    for (FLENS_DEFAULT_INDEXTYPE k=2; k<=d_; ++k) {
        for (FLENS_DEFAULT_INDEXTYPE i=1; i<=k; ++i) {
            for (FLENS_DEFAULT_INDEXTYPE s=0; s<=i-1; ++s) {
                D3(k,i) += powii(-1,(i-1)-s)*binomial(i-1,s)*powii(s,k-1);
            }
        }
    }

    flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > J(d_,d_);
    for (FLENS_DEFAULT_INDEXTYPE j=1; j<=d_; ++j) {
        J(d_+1-j,j) = 1;
    }

    flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > D4, Tmp;
    flens::blas::mm(cxxblas::NoTrans,cxxblas::Trans,1.,J,D3,0.,Tmp);
    flens::blas::mm(cxxblas::NoTrans,cxxblas::NoTrans,1.,Tmp,J,0.,D4);

    flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > m1(d_,d_);
    for (FLENS_DEFAULT_INDEXTYPE j=1; j<=d_; ++j) {
        for (FLENS_DEFAULT_INDEXTYPE k=1; k<=d_; ++k) {
            m1(j,k) = M1(d-2+j,d-2+k);
        }
    }

    flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > m2(2*d_+d-3,d_);
    for (FLENS_DEFAULT_INDEXTYPE j=1; j<=2*d_+d-3; ++j) {
        for (FLENS_DEFAULT_INDEXTYPE k=1; k<=d_; ++k) {
            m2(j,k) = M1(d+d_-2+j,d-2+k);
        }
    }
    flens::blas::mm(cxxblas::NoTrans,cxxblas::NoTrans,1.,D4,m1,0.,Tmp);

    flens::blas::mm(cxxblas::NoTrans,cxxblas::NoTrans,1.,Tmp,inv(D4),0.,m1);
    Tmp = m2;
    flens::blas::mm(cxxblas::NoTrans,cxxblas::NoTrans,1.,Tmp,inv(D4),0.,m2);

    for (FLENS_DEFAULT_INDEXTYPE j=1; j<=d_; ++j) {
        for (FLENS_DEFAULT_INDEXTYPE k=1; k<=d_; ++k) {
            MD(d-2+j,d-2+k) = m1(j,k);
        }
    }
    for (FLENS_DEFAULT_INDEXTYPE j=1; j<=2*d_+d-3; ++j) {
        for (FLENS_DEFAULT_INDEXTYPE k=1; k<=d_; ++k) {
            MD(d+d_-2+j,d-2+k) = m2(j,k);
        }
    }
}

template <typename T>
void
MRA<T,Dual,Interval,Primbs>::_twoScaleDual_2(FLENS_DEFAULT_INDEXTYPE d, FLENS_DEFAULT_INDEXTYPE d_)
{
    FLENS_DEFAULT_INDEXTYPE t = std::max(d+d_-2,4*d-6);
    flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > TL(t,t);
    for (FLENS_DEFAULT_INDEXTYPE j=1; j<=t; ++j) {
        TL(j,j) = 1.;
    }
    _twoScalePrimal(d,d_);
    _twoScaleDual_1(d,d_);
    MDD.engine().resize(3*d_+2*d-5,d+d_-2) || MDD.engine().fill();

    for (FLENS_DEFAULT_INDEXTYPE j=1; j<=3*d_+2*d-5; ++j) {
        for (FLENS_DEFAULT_INDEXTYPE k=1; k<=d+d_-2; ++k) {
            MDD(j,k) = MD(j,k);
        }
    }

    for (FLENS_DEFAULT_INDEXTYPE k=d-2; k>=1; --k) {
        flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > B(d-1+k,d-1+k), C;
        _orthonormalize(MP,MDD);
        for (FLENS_DEFAULT_INDEXTYPE j=1; j<=d+d_-2; ++j) {
            for (FLENS_DEFAULT_INDEXTYPE s=1; s<=d+d_-2; ++s) {
                TL(j,s) = TT(j,s);
            }
        }
        for (FLENS_DEFAULT_INDEXTYPE j=1; j<=d-1+k; ++j) {
            for (FLENS_DEFAULT_INDEXTYPE s=1; s<=d-1+k; ++s) {
                for (FLENS_DEFAULT_INDEXTYPE l=1; l<=2*k+2*d-2; ++l) {
                    B(j,s) += MP(l,j)*TL(l,k+s);
                }
                B(j,s) *= .5;
            }
        }

        C = inv(B);

        MDD(k,k) = 0;
        flens::DenseVector<flens::Array<T> > x(d+k-1), y(d+k-1);
        for (FLENS_DEFAULT_INDEXTYPE j=1; j<=d+k-1; ++j) {
            x(j) = MDD(k,k)*.5*MP(k,j);
        }
        y(k) = 1; x *= -1; x += y; y = C*x;
        for (FLENS_DEFAULT_INDEXTYPE j=1; j<=d+k-1; ++j) {
            MDD(k+j,k) = y(j);
        }
    }
    _orthonormalize(MP,MDD);
    flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > InvT, D(2*d_+d-3,2*d_+d-3), Tmp, TmpTT;
    InvT = inv(TT);

    flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > K(3*d_+2*d-5,3*d_+2*d-5);
    K.diag(0) = 1;
    K(TT.rows(),TT.cols()) = TT;
    flens::blas::mm(cxxblas::NoTrans,cxxblas::NoTrans,1.,K,MDD,0.,Tmp);
    flens::blas::mm(cxxblas::NoTrans,cxxblas::NoTrans,1.,Tmp,InvT,0.,MDDD);
}


template <typename T>
void
MRA<T,Dual,Interval,Primbs>::_calcM0_()
{
    _twoScaleDual_2(d,d_);

    flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > Mj0_(pow2i<T>(min_j0+1)+d-1,
                                            pow2i<T>(min_j0)+d-1),
                                       Mj0_Right, Bugfix;
    T factor = Const<T>::R_SQRT2;
    Bugfix = MDDD;
    flens::blas::scal(factor,Bugfix);
    Mj0_(MDDD.rows(), MDDD.cols()) = Bugfix;
    FLENS_DEFAULT_INDEXTYPE row = d+d_-1;
    for (FLENS_DEFAULT_INDEXTYPE c=MDDD.lastCol()+1; c<=Mj0_.lastCol()-MDDD.numCols(); ++c,row+=2) {
        Mj0_(_(row,row+phi_R.a_.length()-1),c) = factor*phi_R.a_;
    }
    arrow(Bugfix,Mj0_Right);
    Mj0_(_(Mj0_.lastRow()-Mj0_Right.numRows()+1,Mj0_.lastRow()),
         _(Mj0_.lastCol()-Mj0_Right.numCols()+1,Mj0_.lastCol())) = Mj0_Right;

    M0_ = flens::RefinementMatrix<T,Interval,Primbs>(Mj0_Right.numCols(),Mj0_Right.numCols(), 
                                              Mj0_, min_j0, min_j0);
    setLevel(_j);
}

} // namespace lawa

