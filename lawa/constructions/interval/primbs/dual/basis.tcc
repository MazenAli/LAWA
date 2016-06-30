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

namespace lawa {

template <typename T>
Basis<T,Dual,Interval,Primbs>::Basis(FLENS_DEFAULT_INDEXTYPE _d, FLENS_DEFAULT_INDEXTYPE _d_, FLENS_DEFAULT_INDEXTYPE j)
    : mra(_d, j), mra_(_d, _d_, j),
      d(_d), d_(_d_), mu(d&1),
      min_j0(mra_.min_j0), j0(mra_.j0), _bc(2,0), _j(j0), psi_(*this)
{
    _calcM1_();
}

template <typename T>
FLENS_DEFAULT_INDEXTYPE
Basis<T,Dual,Interval,Primbs>::level() const
{
    return _j;
}

template <typename T>
void
Basis<T,Dual,Interval,Primbs>::setLevel(FLENS_DEFAULT_INDEXTYPE j) const
{
    if (j!=_j) {
        assert(j>=min_j0);
        _j = j;
        mra.setLevel(_j);
        mra_.setLevel(_j);
        M1_.setLevel(_j);
    }
}

template <typename T>
template <BoundaryCondition BC>
void
Basis<T,Dual,Interval,Primbs>::enforceBoundaryCondition()
{
    typedef flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > 
            FullColMatrix;

    BSpline<T,Primal,R,CDF> phiR(d);

    if ((_bc(0)==0) && (_bc(1)==0)) {
        setLevel(min_j0);
        _bc(0) = _bc(1) = 1;
        mra.template enforceBoundaryCondition<BC>();
        mra_.template enforceBoundaryCondition<BC>();

        FLENS_DEFAULT_INDEXTYPE j = min_j0;
        FullColMatrix H(pow2i<T>(j+1)-d+1,
                        pow2i<T>(j+1)-d+1), 
                      Tmp, Tmp2, Tmp3, Tmp4, Mj1, Mjj1;
        H.diag(0) = 1.;

        flens::DenseVector<flens::Array<T> > x(pow2i<T>(j+1)-d+1), tmp;
        for (FLENS_DEFAULT_INDEXTYPE i=1; i<=d+1; ++i) {
            x(i) = phiR.a(phiR.a.firstIndex()-1+i);
        }

        for (FLENS_DEFAULT_INDEXTYPE i=1; i<=d; ++i) {
            FullColMatrix A(pow2i<T>(j+1)-d+1,
                             pow2i<T>(j+1)-d+1);
            A.diag(0) = 1.;
            if (i-2*ifloor(i/2.)==1) {
                for (FLENS_DEFAULT_INDEXTYPE k=iceil<FLENS_DEFAULT_INDEXTYPE>(1/2.-(i+1)/4.); k<=ifloor(pow2i<T>(j)-d/2.-(i+1)/4.); ++k) {
                    A((i+1)/2.+2*k,(i+1)/2.+2*k+1) = -x(ifloor(i/2.)+1) / x(ifloor(i/2.)+2);
                }
                tmp = x;
                x = A*tmp;
                flens::blas::mm(cxxblas::NoTrans,cxxblas::NoTrans,1.,A,H,0.,Tmp);
                H = Tmp;
            } else if (i-2*ifloor(i/2.)==0) {
                for (FLENS_DEFAULT_INDEXTYPE k=iceil<FLENS_DEFAULT_INDEXTYPE>(1/2.-i/4.); k<=ifloor(pow2i<T>(j)-d/2.-i/4.); ++k) {
                    A(pow2i<T>(j+1)-d+2-i/2.-2*k,pow2i<T>(j+1)-d+1-i/2.-2*k) = -x(d-i/2.+2)/x(d-i/2.+1);
                }
                tmp = x;
                x = A*tmp;
                flens::blas::mm(cxxblas::NoTrans,cxxblas::NoTrans,1.,A,H,0.,Tmp);
                H = Tmp;
            }
        }
        T b = 0;
        for (FLENS_DEFAULT_INDEXTYPE i=1; i<=pow2i<T>(j+1)-d+1; ++i) {
            b += x(i);
        }

        FullColMatrix Hj(pow2i<T>(j+1)+d-3,
                          pow2i<T>(j+1)+d-3), InvHj;
        Hj.diag(0) = 1.;
        for (FLENS_DEFAULT_INDEXTYPE n=1; n<=pow2i<T>(j+1)-d+1; ++n) {
            for (FLENS_DEFAULT_INDEXTYPE m=1; m<=pow2i<T>(j+1)-d+1; ++m) {
                Hj(d-2+n,d-2+m) = H(n,m);
            }
        }
        InvHj = inv(Hj);

        FullColMatrix F(pow2i<T>(j+1)-d+1,pow2i<T>(j)-d+1);
        for (FLENS_DEFAULT_INDEXTYPE k=1; k<=pow2i<T>(j)-d+1; ++k) {
            F(iceil<FLENS_DEFAULT_INDEXTYPE>(d/2.)+2*k-2,k) = 1;
        }

        FullColMatrix Fj(pow2i<T>(j+1)+d-3,pow2i<T>(j));
        for (FLENS_DEFAULT_INDEXTYPE k=1; k<=ifloor(d/2.); ++k) {
            Fj(d-2+k,k) =1;
            Fj(pow2i<T>(j+1)-k,pow2i<T>(j)+1-k) = 1;
        }
        for (FLENS_DEFAULT_INDEXTYPE n=1; n<=pow2i<T>(j+1)-d+1; ++n) {
            for (FLENS_DEFAULT_INDEXTYPE m=1; m<=pow2i<T>(j)-d+1; ++m) {
                Fj(d-2+n,iceil<FLENS_DEFAULT_INDEXTYPE>(d/2.)-1+m) = F(n,m);
            }
        }

        FullColMatrix Pj(pow2i<T>(j+1)+d-3,
                          pow2i<T>(j+1)+d-3), 
                       Identity(pow2i<T>(j+1)+d-3,
                                pow2i<T>(j+1)+d-3);
        Pj.diag(0) = 1.;
        for (FLENS_DEFAULT_INDEXTYPE n=1; n<=2*(d-2)+1; ++n) {
            for (FLENS_DEFAULT_INDEXTYPE m=1; m<=d-2; ++m) {
                Pj(n,m) = mra_.MP(n+1,m+1);
                Pj(pow2i<T>(j+1)+d-2-n,pow2i<T>(j+1)+d-2-m) = mra_.MP(n+1,m+1);
            }
        }

        flens::blas::mm(cxxblas::NoTrans,cxxblas::NoTrans,1.,Pj,InvHj,0.,Tmp);
        flens::blas::mm(cxxblas::NoTrans,cxxblas::NoTrans,1.,Tmp,Fj,0.,Tmp2);

        Identity.diag(0) = 1.;
        FullColMatrix DM0, DM0_;
        densify(cxxblas::NoTrans, mra.M0, DM0);
        densify(cxxblas::NoTrans, mra_.M0_, DM0_);

        T factor = powii(-1,d+1)*1./b;
        flens::blas::mm(cxxblas::NoTrans,cxxblas::Trans,1.,DM0,DM0_,0.,Tmp3);
        Tmp4 = Identity;
        flens::blas::axpy(cxxblas::NoTrans,-1.,Tmp3,Tmp4);
        flens::blas::mm(cxxblas::NoTrans,cxxblas::NoTrans,factor,Tmp4,Tmp2,0.,Mj1);

        flens::blas::mm(cxxblas::Trans,cxxblas::NoTrans,1.,Fj,Hj,0.,Tmp);
        flens::blas::mm(cxxblas::NoTrans,cxxblas::NoTrans,1.,Tmp,inv(Pj),0.,Tmp2);

        Mjj1.engine().resize(Tmp2.cols(),Tmp2.rows());
        for (FLENS_DEFAULT_INDEXTYPE row=Mjj1.firstRow(); row<=Mjj1.lastRow(); ++row) {
            for (FLENS_DEFAULT_INDEXTYPE col=Mjj1.firstCol(); col<=Mjj1.lastCol(); ++col) {
                Mjj1(row,col) = 2./factor*Tmp2(col,row); 
            }
        }

        if ((d/2.)>ifloor(d/2.)) {
            FullColMatrix Nj1(pow2i<T>(j+1)+d-3,pow2i<T>(j-1)),
            Nj1tr(pow2i<T>(j+1)+d-3,pow2i<T>(j-1)),
            Njj1(pow2i<T>(j+1)+d-3,pow2i<T>(j-1)),
            Njj1tr(pow2i<T>(j+1)+d-3,pow2i<T>(j-1));
            for (FLENS_DEFAULT_INDEXTYPE n=1; n<=pow2i<T>(j+1)+d-3; ++n) {
                for (FLENS_DEFAULT_INDEXTYPE m=1; m<=pow2i<T>(j-1); ++m) {
                    Nj1(n,m) = Mj1(n,m);
                    Nj1tr(pow2i<T>(j+1)+d-2-n,pow2i<T>(j-1)+1-m) = Nj1(n,m);
                    Njj1(n,m) = Mjj1(n,m);
                    Njj1tr(pow2i<T>(j+1)+d-2-n,pow2i<T>(j-1)+1-m) = Njj1(n,m);
                }
            }

            FullColMatrix Hj1(pow2i<T>(j+1)+d-3,pow2i<T>(j)),
            Hjj1(pow2i<T>(j+1)+d-3,pow2i<T>(j));
            for (FLENS_DEFAULT_INDEXTYPE n=1; n<=pow2i<T>(j+1)+d-3; ++n) {
                for (FLENS_DEFAULT_INDEXTYPE m=1; m<=pow2i<T>(j-1); ++m) {
                    Hj1(n,m) = Nj1(n,m);
                    Hjj1(n,m) = Njj1(n,m);
                }
                for (FLENS_DEFAULT_INDEXTYPE m=pow2i<T>(j-1)+1; m<=pow2i<T>(j); ++m) {
                    Hj1(n,m) = Nj1tr(n,m-pow2i<T>(j-1));
                    Hjj1(n,m) = Njj1tr(n,m-pow2i<T>(j-1));
                }
            }

            FullColMatrix Kj, InvKj;
            flens::blas::mm(cxxblas::Trans,cxxblas::NoTrans,1.,Hjj1,Hj1,0.,Kj);
            InvKj = inv(Kj);
            flens::blas::mm(cxxblas::NoTrans,cxxblas::NoTrans,2.,Hj1,InvKj,0.,Mj1);
            Mjj1 = Hjj1;
        }

        FLENS_DEFAULT_INDEXTYPE numPrimalCols = (d+d_)/2-1,
            numDualCols = (d+1)/2-1;

        if (d==2) {
            numPrimalCols += 1;
            numDualCols += 1;
        }

        flens::blas::scal(Const<T>::R_SQRT2*.5, Mjj1);
        M1_ = flens::RefinementMatrix<T,Interval,Primbs>(numDualCols, numDualCols, Mjj1, min_j0, min_j0);
        setLevel(_j);
    }
}

template <typename T>
const BasisFunction<T,Dual,Interval,Primbs> &
Basis<T,Dual,Interval,Primbs>::generator(XType xtype) const
{
    if (xtype==XBSpline) {
        return mra.phi_;
    } else {
        return psi_;
    }
}

// cardinalities of whole, left, inner, right index sets (primal).
template <typename T>
FLENS_DEFAULT_INDEXTYPE
Basis<T,Dual,Interval,Primbs>::cardJ_(FLENS_DEFAULT_INDEXTYPE j) const
{
    assert(0);
    assert(j>=min_j0);
    return pow2i<T>(j);
}

template <typename T>
FLENS_DEFAULT_INDEXTYPE
Basis<T,Dual,Interval,Primbs>::cardJ_L(FLENS_DEFAULT_INDEXTYPE j) const
{
    assert(0);
    assert(j>=min_j0);
    return d + d_ - 2;
}

template <typename T>
FLENS_DEFAULT_INDEXTYPE
Basis<T,Dual,Interval,Primbs>::cardJ_I(FLENS_DEFAULT_INDEXTYPE j) const
{
    assert(0);
    assert(j>=min_j0);
    return pow2i<T>(j) - 2*(d + d_ - 2);
}

template <typename T>
FLENS_DEFAULT_INDEXTYPE
Basis<T,Dual,Interval,Primbs>::cardJ_R(FLENS_DEFAULT_INDEXTYPE j) const
{
    assert(0);
    assert(j>=min_j0);
    return d + d_ - 2;
}

// ranges of whole, left, inner, right index sets (primal).
template <typename T>
const flens::Range<FLENS_DEFAULT_INDEXTYPE>
Basis<T,Dual,Interval,Primbs>::rangeJ_(FLENS_DEFAULT_INDEXTYPE j) const
{
    assert(0);
    assert(j>=min_j0);
    return flens::_(1,pow2i<T>(j));
}

template <typename T>
const flens::Range<FLENS_DEFAULT_INDEXTYPE>
Basis<T,Dual,Interval,Primbs>::rangeJ_L(FLENS_DEFAULT_INDEXTYPE j) const
{
    assert(0);
    assert(j>=min_j0);
    return flens::_(1,d+d_-2);
}

template <typename T>
const flens::Range<FLENS_DEFAULT_INDEXTYPE>
Basis<T,Dual,Interval,Primbs>::rangeJ_I(FLENS_DEFAULT_INDEXTYPE j) const
{
    assert(0);
    assert(j>=min_j0);
    return flens::_(d+d_-1,pow2i<T>(j)-(d+d_-3));
    
}

template <typename T>
const flens::Range<FLENS_DEFAULT_INDEXTYPE>
Basis<T,Dual,Interval,Primbs>::rangeJ_R(FLENS_DEFAULT_INDEXTYPE j) const
{
    assert(0);
    assert(j>=min_j0);
    return flens::_(pow2i<T>(j)-(d+d_-2),pow2i<T>(j));

}

template <typename T>
void
Basis<T,Dual,Interval,Primbs>::_calcM1_()
{
    typedef flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > FullColMatrix;
    BSpline<T,Primal,R,CDF> phiR(d);

    setLevel(min_j0);
    FullColMatrix H(pow2i<T>(min_j0+1)-d+1,pow2i<T>(min_j0+1)-d+1), 
                   Tmp, Tmp2, Tmp3, Tmp4, Mj1, Mjj1;
    H.diag(0) = 1.;

    flens::DenseVector<flens::Array<T> > x(pow2i<T>(min_j0+1)-d+1), tmp;
    for (FLENS_DEFAULT_INDEXTYPE i=1; i<=d+1; ++i) {
        x(i) = phiR.a(phiR.a.firstIndex()-1+i);
    }
    for (FLENS_DEFAULT_INDEXTYPE i=1; i<=d; ++i) {
        FullColMatrix A(pow2i<T>(min_j0+1)-d+1,pow2i<T>(min_j0+1)-d+1);
        A.diag(0) = 1.;
        if (i-2*ifloor(i/2.)==1) {
            for (FLENS_DEFAULT_INDEXTYPE k=iceil<FLENS_DEFAULT_INDEXTYPE>(1/2.-(i+1)/4.); k<=ifloor(pow2i<T>(min_j0)-d/2.-(i+1)/4.); ++k) {
                A((i+1)/2.+2*k,(i+1)/2.+2*k+1) = -x(ifloor(i/2.)+1) / x(ifloor(i/2.)+2);
            }
            tmp = x;
            x = A*tmp;
            flens::blas::mm(cxxblas::NoTrans,cxxblas::NoTrans,1.,A,H,0.,Tmp);
            H = Tmp;
        } else if (i-2*ifloor(i/2.)==0) {
            for (FLENS_DEFAULT_INDEXTYPE k=iceil<FLENS_DEFAULT_INDEXTYPE>(1/2.-i/4.); k<=ifloor(pow2i<T>(min_j0)-d/2.-i/4.); ++k) {
                A(pow2i<T>(min_j0+1)-d+2-i/2.-2*k,pow2i<T>(min_j0+1)-d+1-i/2.-2*k) = -x(d-i/2.+2)/x(d-i/2.+1);
            }
            tmp = x;
            x = A*tmp;
            flens::blas::mm(cxxblas::NoTrans,cxxblas::NoTrans,1.,A,H,0.,Tmp);
            H = Tmp;
        }
    }

    T b = 0;
    for (FLENS_DEFAULT_INDEXTYPE i=1; i<=pow2i<T>(min_j0+1)-d+1; ++i) {
        b += x(i);
    }

    FullColMatrix Hj(pow2i<T>(min_j0+1)+d-1,pow2i<T>(min_j0+1)+d-1), InvHj;
    Hj.diag(0) = 1.;
    for (FLENS_DEFAULT_INDEXTYPE n=1; n<=pow2i<T>(min_j0+1)-d+1; ++n) {
        for (FLENS_DEFAULT_INDEXTYPE m=1; m<=pow2i<T>(min_j0+1)-d+1; ++m) {
            Hj(d-1+n,d-1+m) = H(n,m);
        }
    }
    InvHj = inv(Hj);

    FullColMatrix F(pow2i<T>(min_j0+1)-d+1,pow2i<T>(min_j0)-d+1);
    for (FLENS_DEFAULT_INDEXTYPE k=1; k<=pow2i<T>(min_j0)-d+1; ++k) {
        F(iceil<FLENS_DEFAULT_INDEXTYPE>(d/2.)+2*k-2,k) = 1;
    }

    FullColMatrix Fj(pow2i<T>(min_j0+1)+d-1,pow2i<T>(min_j0));
    for (FLENS_DEFAULT_INDEXTYPE k=1; k<=iceil<FLENS_DEFAULT_INDEXTYPE>((d-1)/2.); ++k) {
        Fj(d-1+k,k) =1;
        Fj(pow2i<T>(min_j0+1)+1-k,pow2i<T>(min_j0)+1-k) = 1;
    }
    for (FLENS_DEFAULT_INDEXTYPE n=1; n<=pow2i<T>(min_j0+1)-d+1; ++n) {
        for (FLENS_DEFAULT_INDEXTYPE m=1; m<=pow2i<T>(min_j0)-d+1; ++m) {
            Fj(d-1+n,ifloor((d-1)/2.+m)) = F(n,m);
        }
    }

    FullColMatrix Pj(pow2i<T>(min_j0+1)+d-1,pow2i<T>(min_j0+1)+d-1),
                   Identity(pow2i<T>(min_j0+1)+d-1,pow2i<T>(min_j0+1)+d-1);
    Pj.diag(0) = 1.;
    for (FLENS_DEFAULT_INDEXTYPE n=1; n<=2*(d-1); ++n) {
        for (FLENS_DEFAULT_INDEXTYPE m=1; m<=d-1; ++m) {
            Pj(n,m) = mra_.MP(n,m);
            Pj(pow2i<T>(min_j0+1)+d-n,pow2i<T>(min_j0+1)+d-m) = mra_.MP(n,m);
        }
    }

    flens::blas::mm(cxxblas::NoTrans,cxxblas::NoTrans,1.,Pj,InvHj,0.,Tmp);
    flens::blas::mm(cxxblas::NoTrans,cxxblas::NoTrans,1.,Tmp,Fj,0.,Tmp2);
    Identity.diag(0) = 1.;

    FullColMatrix DM0, DM0_;
    densify(cxxblas::NoTrans, mra.M0, DM0);
    densify(cxxblas::NoTrans, mra_.M0_, DM0_);

    T factor = powii(-1,d+1)*2./b;
    flens::blas::mm(cxxblas::NoTrans,cxxblas::Trans,1.,DM0,DM0_,0.,Tmp3);
    Tmp4 = Identity;
    flens::blas::axpy(cxxblas::NoTrans,-1.,Tmp3,Tmp4);
    flens::blas::mm(cxxblas::NoTrans,cxxblas::NoTrans,factor,Tmp4,Tmp2,0.,Mj1);

    flens::blas::mm(cxxblas::Trans,cxxblas::NoTrans,1.,Fj,Hj,0.,Tmp);
    flens::blas::mm(cxxblas::NoTrans,cxxblas::NoTrans,1.,Tmp,inv(Pj),0.,Tmp2);

//    flens::blas::axpy(cxxblas::Trans,2./factor,Tmp2,Mjj1);
    Mjj1.engine().resize(Tmp2.cols(),Tmp2.rows());
    for (FLENS_DEFAULT_INDEXTYPE row=Mjj1.firstRow(); row<=Mjj1.lastRow(); ++row) {
        for (FLENS_DEFAULT_INDEXTYPE col=Mjj1.firstCol(); col<=Mjj1.lastCol(); ++col) {
            Mjj1(row,col) = 2./factor*Tmp2(col,row); 
        }
    }

    if ((d/2.)>ifloor(d/2.)) {
        FullColMatrix Nj1(pow2i<T>(min_j0+1)+d-1,pow2i<T>(min_j0-1)),
                      Nj1tr(pow2i<T>(min_j0+1)+d-1,pow2i<T>(min_j0-1)),
                      Njj1(pow2i<T>(min_j0+1)+d-1,pow2i<T>(min_j0-1)),
                      Njj1tr(pow2i<T>(min_j0+1)+d-1,pow2i<T>(min_j0-1));
        for (FLENS_DEFAULT_INDEXTYPE n=1; n<=pow2i<T>(min_j0+1)+d-1; ++n) {
            for (FLENS_DEFAULT_INDEXTYPE m=1; m<=pow2i<T>(min_j0-1); ++m) {
                Nj1(n,m) = Mj1(n,m);
                Nj1tr(pow2i<T>(min_j0+1)+d-n,pow2i<T>(min_j0-1)+1-m) = Nj1(n,m);
                Njj1(n,m) = Mjj1(n,m);
                Njj1tr(pow2i<T>(min_j0+1)+d-n,pow2i<T>(min_j0-1)+1-m) = Njj1(n,m);
            }
        }

        FullColMatrix Hj1(pow2i<T>(min_j0+1)+d-1,pow2i<T>(min_j0)),
                      Hjj1(pow2i<T>(min_j0+1)+d-1,pow2i<T>(min_j0));
        for (FLENS_DEFAULT_INDEXTYPE n=1; n<=pow2i<T>(min_j0+1)+d-1; ++n) {
            for (FLENS_DEFAULT_INDEXTYPE m=1; m<=pow2i<T>(min_j0-1); ++m) {
                Hj1(n,m) = Nj1(n,m);
                Hjj1(n,m) = Njj1(n,m);
            }
            for (FLENS_DEFAULT_INDEXTYPE m=pow2i<T>(min_j0-1)+1; m<=pow2i<T>(min_j0); ++m) {
                Hj1(n,m) = Nj1tr(n,m-pow2i<T>(min_j0-1));
                Hjj1(n,m) = Njj1tr(n,m-pow2i<T>(min_j0-1));
            }
        }

        FullColMatrix Kj, InvKj;
        flens::blas::mm(cxxblas::Trans,cxxblas::NoTrans,1.,Hjj1,Hj1,0.,Kj);
        InvKj = inv(Kj);
        flens::blas::mm(cxxblas::NoTrans,cxxblas::NoTrans,2.,Hj1,InvKj,0.,Mj1);
        Mjj1 = Hjj1;
    }

    FLENS_DEFAULT_INDEXTYPE numPrimalCols = d-1 + (d_-d)/2,
        numDualCols = (d+1)/2-1;

    if (d==2) {
        numPrimalCols += 1;
        numDualCols += 1;
    }

    flens::blas::scal(Const<T>::R_SQRT2*.5, Mjj1);
    M1_ = flens::RefinementMatrix<T,Interval,Primbs>(numDualCols,numDualCols,Mjj1,min_j0,min_j0);
    setLevel(_j);
}


} // namespace lawa

