/*
  LAWA - Library for Adaptive Wavelet Applications.
  Copyright (C) 2008-2014  Sebastian Kestler, Mario Rometsch, Kristina Steih, 
  Alexander Stippler.

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

#ifndef  LAWA_METHODS_ADAPTIVE_ALGORITHMS_MULTITREEOPERATIONS_H
#define  LAWA_METHODS_ADAPTIVE_ALGORITHMS_MULTITREEOPERATIONS_H 1


#include <iostream>
#include <cstring>
#include <lawa/constructions/constructions.h>
#include <lawa/methods/adaptive/datastructures/datastructures.h>

namespace lawa {

// ----------------- extendMultiTree 1D-3D ------------------- //


template <typename T, typename Basis>
void
extendMultiTree(const Basis &basis, const Coefficients<Lexicographical,T,Index1D>  &v,
                Coefficients<Lexicographical,T,Index1D>  &C_v, const char* residualType,
                bool sparsetree=false);

template <typename T, typename Basis>
void
extendMultiTree(const Basis &basis, const Coefficients<Lexicographical,T,Index2D>  &v,
                Coefficients<Lexicographical,T,Index2D>  &C_v, const char* residualType,
                bool IsMW=false, bool sparsetree=false);

template <typename T, typename Basis>
void
extendMultiTree(const Basis &basis, const Coefficients<Lexicographical,T,Index2D>  &v,
                Coefficients<Lexicographical,T,Index2D>  &C_v, IndexSet<Index2D>& Cdiff_v,
                const char* residualType, bool IsMW=false, bool sparsetree=false);

template <typename T, typename Basis>
void
extendMultiTree(const Basis &basis, const Coefficients<Lexicographical,T,Index3D>  &v,
                Coefficients<Lexicographical,T,Index3D>  &C_v, const char* residualType,
                bool IsMW=false, bool sparsetree=false);

// ----------------- extendMultiTreeAtBoundary 2D,3D ------------------- //


template <typename T, typename Basis>
void
extendMultiTreeAtBoundary(const Basis &basis, const Coefficients<Lexicographical,T,Index2D>  &v,
                          Coefficients<Lexicographical,T,Index2D>  &C_v, int J,
                          bool sparsetree=false);

template <typename T, typename Basis>
void
extendMultiTreeAtBoundary(const Basis &basis, const Coefficients<Lexicographical,T,Index3D>  &v,
                          Coefficients<Lexicographical,T,Index3D>  &C_v, int J,
                          bool sparsetree=false);


/*
 * To partially specialize for periodic basis, we use SFINAE
 */

// ----------------- completeMultiTree 1D ------------------- //

// ORIGINAL: Non-Periodic version
template <typename T, typename Basis>
typename cxxblas::RestrictTo<SFINAE_Wrapper<!IsPeriodic<Basis>::value, T>::value, void>::Type
completeMultiTree(const Basis &basis, const Index1D &index1d,
                  Coefficients<Lexicographical,T,Index1D>  &v, bool sparsetree=false);

// Periodic version
template <typename T, typename Basis>
typename cxxblas::RestrictTo<SFINAE_Wrapper<IsPeriodic<Basis>::value, T>::value, void>::Type
completeMultiTree(const Basis &basis, const Index1D &index1d,
                  Coefficients<Lexicographical,T,Index1D>  &v, bool sparsetree=false);

// ---- Non-Periodic + returning added indizes
template <typename T, typename Basis>
typename cxxblas::RestrictTo<SFINAE_Wrapper<!IsPeriodic<Basis>::value, T>::value, void>::Type
completeMultiTree(const Basis &basis, const Index1D &index1d,
                  Coefficients<Lexicographical,T,Index1D>  &v,
                  IndexSet<Index2D>& diff_v, bool sparsetree=false);

// ---- Periodic + returning added indizes
template <typename T, typename Basis>
typename cxxblas::RestrictTo<SFINAE_Wrapper<IsPeriodic<Basis>::value, T>::value, void>::Type
completeMultiTree(const Basis &basis, const Index1D &index1d,
                  Coefficients<Lexicographical,T,Index1D>  &v,
                  IndexSet<Index2D>& diff_v, bool sparsetree=false);


// ----------------- completeMultiTree 2D ------------------- //

// ORIGINAL: Non-Periodic version
template <typename T, typename Basis>
typename cxxblas::RestrictTo<SFINAE_Wrapper<!IsPeriodic<typename Basis::FirstBasisType>::value
					and !IsPeriodic<typename Basis::SecondBasisType>::value, T>::value, void>::Type
completeMultiTree(const Basis &basis, const Index2D &index2d,
                  Coefficients<Lexicographical,T,Index2D>  &v,
                  int coordDirec=0, bool sparsetree=false, bool isAlreadyMultiTree=true);

// Periodic version
template <typename T, typename Basis>
typename cxxblas::RestrictTo<SFINAE_Wrapper<IsPeriodic<typename Basis::FirstBasisType>::value
					and IsPeriodic<typename Basis::SecondBasisType>::value, T>::value, void>::Type
completeMultiTree(const Basis &basis, const Index2D &index2d,
                  Coefficients<Lexicographical,T,Index2D>  &v,
                  int coordDirec=0, bool sparsetree=false, bool isAlreadyMultiTree=true);

// Periodic-NonPeriodic version
template <typename T, typename Basis>
typename cxxblas::RestrictTo<SFINAE_Wrapper<IsPeriodic<typename Basis::FirstBasisType>::value
					and !IsPeriodic<typename Basis::SecondBasisType>::value, T>::value, void>::Type
completeMultiTree(const Basis &basis, const Index2D &index2d,
                  Coefficients<Lexicographical,T,Index2D>  &v,
                  int coordDirec=0, bool sparsetree=false, bool isAlreadyMultiTree=true);

// NonPeriodic-Periodic version
template <typename T, typename Basis>
typename cxxblas::RestrictTo<SFINAE_Wrapper<!IsPeriodic<typename Basis::FirstBasisType>::value
					and IsPeriodic<typename Basis::SecondBasisType>::value, T>::value, void>::Type
completeMultiTree(const Basis &basis, const Index2D &index2d,
                  Coefficients<Lexicographical,T,Index2D>  &v,
                  int coordDirec=0, bool sparsetree=false, bool isAlreadyMultiTree=true);

// ---- Non-Periodic + returning added indizes
template <typename T, typename Basis>
typename cxxblas::RestrictTo<SFINAE_Wrapper<!IsPeriodic<typename Basis::FirstBasisType>::value
					and !IsPeriodic<typename Basis::SecondBasisType>::value, T>::value, void>::Type
completeMultiTree(const Basis &basis, const Index2D &index2d,
                  Coefficients<Lexicographical,T,Index2D>  &v,
                  IndexSet<Index2D>& diff_v,
                  int coordDirec=0, bool sparsetree=false);

// ---- Periodic + returning added indizes
template <typename T, typename Basis>
typename cxxblas::RestrictTo<SFINAE_Wrapper<IsPeriodic<typename Basis::FirstBasisType>::value
					and IsPeriodic<typename Basis::SecondBasisType>::value, T>::value, void>::Type
completeMultiTree(const Basis &basis, const Index2D &index2d,
                  Coefficients<Lexicographical,T,Index2D>  &v,
                  IndexSet<Index2D>& diff_v,
                  int coordDirec=0, bool sparsetree=false);

// ---- Periodic-NonPeriodic + returning added indizes
template <typename T, typename Basis>
typename cxxblas::RestrictTo<SFINAE_Wrapper<IsPeriodic<typename Basis::FirstBasisType>::value
					and !IsPeriodic<typename Basis::SecondBasisType>::value, T>::value, void>::Type
completeMultiTree(const Basis &basis, const Index2D &index2d,
                  Coefficients<Lexicographical,T,Index2D>  &v,
                  IndexSet<Index2D>& diff_v,
                  int coordDirec=0, bool sparsetree=false);

// ---- NonPeriodic-Periodic + returning added indizes
template <typename T, typename Basis>
typename cxxblas::RestrictTo<SFINAE_Wrapper<!IsPeriodic<typename Basis::FirstBasisType>::value
					and IsPeriodic<typename Basis::SecondBasisType>::value, T>::value, void>::Type
completeMultiTree(const Basis &basis, const Index2D &index2d,
                  Coefficients<Lexicographical,T,Index2D>  &v,
                  IndexSet<Index2D>& diff_v,
                  int coordDirec=0, bool sparsetree=false);

// ----------------- completeMultiTree 3D ------------------- //

// For L2-orth. multiwavelets only!!
/*
template <typename T, typename Basis>
void
completeMultiTree(const Basis &basis, const Index3D &index3d,
                  Coefficients<Lexicographical,T,Index3D>  &v);
*/
template <typename T, typename Basis>
void
completeMultiTree(const Basis &basis, const Index3D &index3d,
                  Coefficients<Lexicographical,T,Index3D>  &v,
                  int coordDirec=0, bool sparsetree=false);



template <typename T, typename Basis>
void
getSparseGridVector(const Basis &basis, Coefficients<Lexicographical,T,Index2D> &v, int j, T gamma);

template <typename T, typename Basis>
void
getSparseGridVector(const Basis &basis, Coefficients<Lexicographical,T,Index3D> &v, int j, T gamma);

// ----------------- extendMultiTree ------------------- //

// ORIGINAL: Non-Periodic version
template <typename Index, typename Basis>
typename cxxblas::RestrictTo<SFINAE_Wrapper<!IsPeriodic<typename Basis::FirstBasisType>::value
					and !IsPeriodic<typename Basis::SecondBasisType>::value, Basis>::value, void>::Type
extendMultiTree(const Basis &basis, const Index &index2d, IndexSet<Index> &Lambda);

// Periodic version
template <typename Index, typename Basis>
typename cxxblas::RestrictTo<SFINAE_Wrapper<IsPeriodic<typename Basis::FirstBasisType>::value
					and IsPeriodic<typename Basis::SecondBasisType>::value, Basis>::value, void>::Type
extendMultiTree(const Basis &basis, const Index &index2d, IndexSet<Index> &Lambda);

// Periodic-NonPeriodic version
template <typename Index, typename Basis>
typename cxxblas::RestrictTo<SFINAE_Wrapper<IsPeriodic<typename Basis::FirstBasisType>::value
					and !IsPeriodic<typename Basis::SecondBasisType>::value, Basis>::value, void>::Type
extendMultiTree(const Basis &basis, const Index &index2d, IndexSet<Index> &Lambda);

// Non-Periodic-Periodic version
template <typename Index, typename Basis>
typename cxxblas::RestrictTo<SFINAE_Wrapper<!IsPeriodic<typename Basis::FirstBasisType>::value
					and IsPeriodic<typename Basis::SecondBasisType>::value, Basis>::value, void>::Type
extendMultiTree(const Basis &basis, const Index &index2d, IndexSet<Index> &Lambda);

// ----------------- extendMultiTree2 ------------------- //

// ORIGINAL: Non-Periodic version
template <typename Index, typename Basis>
typename cxxblas::RestrictTo<SFINAE_Wrapper<!IsPeriodic<typename Basis::FirstBasisType>::value
					and !IsPeriodic<typename Basis::SecondBasisType>::value, Basis>::value, void>::Type
extendMultiTree2(const Basis &basis, const Index &index2d, const int offset, IndexSet<Index> &Lambda);

// Periodic version
template <typename Index, typename Basis>
typename cxxblas::RestrictTo<SFINAE_Wrapper<IsPeriodic<typename Basis::FirstBasisType>::value
					and IsPeriodic<typename Basis::SecondBasisType>::value, Basis>::value, void>::Type
extendMultiTree2(const Basis &basis, const Index &index2d, const int offset, IndexSet<Index> &Lambda);

// Periodic-NonPeriodic version
template <typename Index, typename Basis>
typename cxxblas::RestrictTo<SFINAE_Wrapper<IsPeriodic<typename Basis::FirstBasisType>::value
					and !IsPeriodic<typename Basis::SecondBasisType>::value, Basis>::value, void>::Type
extendMultiTree2(const Basis &basis, const Index &index2d, const int offset, IndexSet<Index> &Lambda);

// NonPeriodic-Periodic version
template <typename Index, typename Basis>
typename cxxblas::RestrictTo<SFINAE_Wrapper<!IsPeriodic<typename Basis::FirstBasisType>::value
					and IsPeriodic<typename Basis::SecondBasisType>::value, Basis>::value, void>::Type
extendMultiTree2(const Basis &basis, const Index &index2d, const int offset, IndexSet<Index> &Lambda);

// ----------------- getCounterpart / getStableExpansion ------------------- //

// To a given indexset in basis_origin, find a corresponding indexset in basis_target.
// This is realized by taking the cone consisting of neighbours of all bfs in indexset_origin
// We have to make sure that this is a MT!!
template <typename T, typename Basis2D_Origin, typename Basis2D_Target>
void
getCounterpart(const Basis2D_Origin& basis_origin, const Basis2D_Target& basis_target,
		IndexSet<Index2D>& indexset_origin, Coefficients<Lexicographical,T,Index2D>& coeffs_target);

// To a given indexset in basis_origin, find a corresponding indexset in basis_target,
// so that the system is "stable" (or at least A^T A is approximated well enough).
// This is realized by taking the cone consisting of neighbours of all bfs in indexset_origin
// plus the HigherWavelet-Neighbours.

// Version A) Include all combinations of HigherWavelet-Neighbours
// (i.e. with levels (lambda_t + 1, lambda_x + 1))
template <typename T, typename Basis2D_Origin, typename Basis2D_Target>
void
getStableExpansion(const Basis2D_Origin& basis_origin, const Basis2D_Target& basis_target,
				   IndexSet<Index2D>& indexset_origin, Coefficients<Lexicographical,T,Index2D>& coeffs_target);

template <typename T, typename Basis2D_Origin, typename Basis2D_Target>
void
getStableExpansion(const Basis2D_Origin& basis_origin, const Basis2D_Target& basis_target,
			Coefficients<Lexicographical,T,Index2D>& coeffs_origin, Coefficients<Lexicographical,T,Index2D>& coeffs_target);

// Version B) Include only HigherWavelet-Neighbours in one coordinate direction
// (i.e. with levels (lambda_t + 1, lambda_x), (lambda_t, lambda_x + 1))
template <typename T, typename Basis2D_Origin, typename Basis2D_Target>
void
getStableExpansion_woMixedHW(const Basis2D_Origin& basis_origin, const Basis2D_Target& basis_target,
				   IndexSet<Index2D>& indexset_origin, Coefficients<Lexicographical,T,Index2D>& coeffs_target);

template <typename T, typename Basis2D_Origin, typename Basis2D_Target>
void
getStableExpansion_woMixedHW(const Basis2D_Origin& basis_origin, const Basis2D_Target& basis_target,
			Coefficients<Lexicographical,T,Index2D>& coeffs_origin, Coefficients<Lexicographical,T,Index2D>& coeffs_target);

// Version C) Include only HigherWavelet-Neighbours in time
// (i.e. with levels (lambda_t + 1, lambda_x))
template <typename T, typename Basis2D_Origin, typename Basis2D_Target>
void
getStableExpansion_onlyTemporalHW(const Basis2D_Origin& basis_origin, const Basis2D_Target& basis_target,
				   IndexSet<Index2D>& indexset_origin, Coefficients<Lexicographical,T,Index2D>& coeffs_target);

template <typename T, typename Basis2D_Origin, typename Basis2D_Target>
void
getStableExpansion_onlyTemporalHW(const Basis2D_Origin& basis_origin, const Basis2D_Target& basis_target,
			Coefficients<Lexicographical,T,Index2D>& coeffs_origin, Coefficients<Lexicographical,T,Index2D>& coeffs_target);

} // namespace lawa

#include <lawa/methods/adaptive/algorithms/multitreeoperations.tcc>

#endif // LAWA_METHODS_ADAPTIVE_ALGORITHMS_MULTITREEOPERATIONS_H
