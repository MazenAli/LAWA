/*
  LAWA - Library for Adaptive Wavelet Applications.
  Copyright (C) 2008,2009  Mario Rometsch, Alexander Stippler.

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

#include <complex>
#include <numeric>


namespace lawa {

//--- Traits for return type.---------------------------------------------------

template <>
struct NormTraits<l1, FLENS_DEFAULT_INDEXTYPE>
{
    typedef FLENS_DEFAULT_INDEXTYPE Type;
};

template <typename T>
struct NormTraits<l1, flens::complex<T> >
{
    typedef T Type;
};

template <>
struct NormTraits<l2, FLENS_DEFAULT_INDEXTYPE>
{
    typedef double Type;
};

template <typename T>
struct NormTraits<l2, flens::complex<T> >
{
    typedef T Type;
};

template <>
struct NormTraits<lInfinity, FLENS_DEFAULT_INDEXTYPE>
{
    typedef FLENS_DEFAULT_INDEXTYPE Type;
};

template <typename T>
struct NormTraits<lInfinity, flens::complex<T> >
{
    typedef T Type;
};

//--- VectorNorm: internal class doing actual calculation ----------------------

template <NormType N, typename T>
struct VectorNorm
{
};

template <typename T>
struct VectorNorm<l1, T>
{
    static const typename NormTraits<l1,T>::Type
    apply(const flens::DenseVector<flens::Array<T> > &x)
    {
        return asum(x);
    }

    static const typename NormTraits<l1,T>::Type
    apply(const flens::Array<T> &x)
    {
        return asum(x.length(), x.data(), 1);
    }
};

template <typename T>
struct VectorNorm<l2, T>
{
    static const typename NormTraits<l2,T>::Type
    apply(const flens::DenseVector<flens::Array<T> > &x)
    {
        T norm=0.;
        cxxblas::nrm2<FLENS_DEFAULT_INDEXTYPE>((FLENS_DEFAULT_INDEXTYPE)x.length(),x.engine().data(),(FLENS_DEFAULT_INDEXTYPE)1,norm);
        return norm;
    }

    static const typename NormTraits<l2,T>::Type
    apply(const flens::Array<T> &x)
    {
        return nrm2(x.length(), x.data(), 1);
    }
};

template <typename T>
struct VectorNorm<lInfinity, T>
{
    static const typename NormTraits<lInfinity,T>::Type
    apply(const flens::DenseVector<flens::Array<T> > &x)
    {
        FLENS_DEFAULT_INDEXTYPE i=0;
        cxxblas::iamax<FLENS_DEFAULT_INDEXTYPE>((FLENS_DEFAULT_INDEXTYPE)x.length(),x.engine().data(),(FLENS_DEFAULT_INDEXTYPE)1,i);
        i+=1;       //Flens vectors are index one-based
        return fabs(x(i));
    }

    static const typename NormTraits<lInfinity,T>::Type
    apply(const flens::Array<T> &x)
    {
        return absolute(x(amax(x.length(), x.data(), 1)));
    }
};

template <NormType N, typename X>
typename NormTraits<N,typename X::ElementType>::Type
norm(const flens::DenseVector<X> &x)
{
    return VectorNorm<N,typename X::ElementType>::apply(x);
}

template <NormType N, typename T>
typename NormTraits<N,T>::Type
norm(const flens::Array<T> &x)
{
    return VectorNorm<N, T>::apply(x);
}

} // namespace flens
