namespace lawa {

template <typename PrimalBasis>
LocalRefinement<PrimalBasis>::LocalRefinement(const PrimalBasis &_basis)
 : basis(_basis), refinementbasis(basis.refinementbasis)
{

}

// Non-Periodic version
template <typename PrimalBasis>
template <typename T_>
typename cxxblas::RestrictTo<SFINAE_Wrapper<!IsPeriodic<PrimalBasis>::value,T_>::value, void>::Type
LocalRefinement<PrimalBasis>::reconstruct(const Coefficients<Lexicographical,T_,Index1D> &u, FLENS_DEFAULT_INDEXTYPE j_scaling,
                                           Coefficients<Lexicographical,T_,Index1D> &u_loc_single) const
{
    TreeCoefficients1D<T> u_tree(255,basis.j0);
    fromCoefficientsToTreeCoefficients(u, u_tree);
    FLENS_DEFAULT_INDEXTYPE j_bspline = j_scaling;
    FLENS_DEFAULT_INDEXTYPE j_wavelet = j_scaling;

    CoefficientsByLevel<T> u_bspline;
    if ((PrimalBasis::Cons==Multi && basis.d>1) || PrimalBasis::Domain==Periodic) {
        this->reconstructOnlyMultiScaling(u_tree.bylevel[0], j_scaling, u_bspline, j_bspline);
        u_tree.bylevel[0] = u_bspline;
    }

    FLENS_DEFAULT_INDEXTYPE imax = u_tree.getMaxTreeLevel();
    for (FLENS_DEFAULT_INDEXTYPE i=0; i<imax; ++i) {
        FLENS_DEFAULT_INDEXTYPE  j_refinement = j_bspline + i;
        CoefficientsByLevel<T> help;
        help = u_tree[i];
        for (const_coeffbylevel_it it=help.map.begin(); it!=help.map.end(); ++it) {

            FLENS_DEFAULT_INDEXTYPE k_refinement = (*it).first;
            FLENS_DEFAULT_INDEXTYPE test_j_wavelet = 0;
            FLENS_DEFAULT_INDEXTYPE k_first = (FLENS_DEFAULT_INDEXTYPE) 0, k_last = (FLENS_DEFAULT_INDEXTYPE) 0;
            refinementbasis.getWaveletNeighborsForBSpline(j_refinement,k_refinement, basis, test_j_wavelet, k_first, k_last);

            assert(test_j_wavelet==j_wavelet+i);
            bool has_neighbor=false;

            for (FLENS_DEFAULT_INDEXTYPE k=k_first; k<=k_last; ++k) {
                if (   u_tree[i+1].map.find(k)!=u_tree[i+1].map.end()) {
                    has_neighbor = true;
                    break;
                }
            }

            if (!has_neighbor) {
                u_tree[i].map.erase((*it).first);
                u_loc_single[Index1D(j_refinement,k_refinement,XBSpline)] = (*it).second;
            }
        }

        CoefficientsByLevel<T> u_loc_single_jP1;
        FLENS_DEFAULT_INDEXTYPE test_j_refinement = 0;
        this->reconstruct(u_tree[i], j_bspline+i, u_tree[i+1], j_wavelet+i, u_loc_single_jP1, test_j_refinement);
        assert(test_j_refinement==j_refinement+1);
        u_tree[i+1] = u_loc_single_jP1;
    }
    for (const_coeffbylevel_it it=u_tree[imax].map.begin(); it!=u_tree[imax].map.end(); ++it) {
        u_loc_single[Index1D(j_bspline+imax,(*it).first,XBSpline)] = (*it).second;
    }
    return;
}

// Periodic version
template <typename PrimalBasis>
template <typename T_>
typename cxxblas::RestrictTo<SFINAE_Wrapper<IsPeriodic<PrimalBasis>::value,T_>::value, void>::Type
LocalRefinement<PrimalBasis>::reconstruct(const Coefficients<Lexicographical,T_,Index1D> &u, FLENS_DEFAULT_INDEXTYPE j_scaling,
                                           Coefficients<Lexicographical,T_,Index1D> &u_loc_single) const
{
    TreeCoefficients1D<T> u_tree(255,basis.j0);
    fromCoefficientsToTreeCoefficients(u, u_tree);
    FLENS_DEFAULT_INDEXTYPE j_bspline = j_scaling;
    FLENS_DEFAULT_INDEXTYPE j_wavelet = j_scaling;

    CoefficientsByLevel<T> u_bspline;
    if ((PrimalBasis::Cons==Multi && basis.d>1) || PrimalBasis::Domain==Periodic) {
        this->reconstructOnlyMultiScaling(u_tree.bylevel[0], j_scaling, u_bspline, j_bspline);
        u_tree.bylevel[0] = u_bspline;
    }

    FLENS_DEFAULT_INDEXTYPE imax = u_tree.getMaxTreeLevel();
    for (FLENS_DEFAULT_INDEXTYPE i=0; i<imax; ++i) {
        FLENS_DEFAULT_INDEXTYPE  j_refinement = j_bspline + i;
        CoefficientsByLevel<T> help;
        help = u_tree[i];
        for (const_coeffbylevel_it it=help.map.begin(); it!=help.map.end(); ++it) {

            FLENS_DEFAULT_INDEXTYPE k_refinement = (*it).first;
            FLENS_DEFAULT_INDEXTYPE test_j_wavelet = 0;
            FLENS_DEFAULT_INDEXTYPE k_first = (FLENS_DEFAULT_INDEXTYPE) 0, k_last = (FLENS_DEFAULT_INDEXTYPE) 0;
            refinementbasis.getWaveletNeighborsForBSpline(j_refinement,k_refinement, basis, test_j_wavelet, k_first, k_last);

            assert(test_j_wavelet==j_wavelet+i);
            bool has_neighbor=false;

            if(k_first < k_last){
                for (FLENS_DEFAULT_INDEXTYPE k=k_first; k<=k_last; ++k) {
                    if (   u_tree[i+1].map.find(k)!=u_tree[i+1].map.end()) {
                        has_neighbor = true;
                        break;
                    }
                }
            }
            else{
                for (FLENS_DEFAULT_INDEXTYPE k=k_first; k<=(FLENS_DEFAULT_INDEXTYPE) refinementbasis.rangeJ(test_j_wavelet).lastIndex(); ++k) {
                    if (   u_tree[i+1].map.find(k)!=u_tree[i+1].map.end()) {
                        has_neighbor = true;
                        break;
                    }
                }
                if(!has_neighbor){
                	for(FLENS_DEFAULT_INDEXTYPE k = refinementbasis.rangeJ(test_j_wavelet).firstIndex(); k <= k_last; ++k){
                        if (   u_tree[i+1].map.find(k)!=u_tree[i+1].map.end()) {
                            has_neighbor = true;
                            break;
                        }
                	}
                }
            }
            if (!has_neighbor) {
                u_tree[i].map.erase((*it).first);
                u_loc_single[Index1D(j_refinement,k_refinement,XBSpline)] = (*it).second;
            }
        }

        CoefficientsByLevel<T> u_loc_single_jP1;
        FLENS_DEFAULT_INDEXTYPE test_j_refinement = 0;
        this->reconstruct(u_tree[i], j_bspline+i, u_tree[i+1], j_wavelet+i, u_loc_single_jP1, test_j_refinement);
        assert(test_j_refinement==j_refinement+1);
        u_tree[i+1] = u_loc_single_jP1;
    }
    for (const_coeffbylevel_it it=u_tree[imax].map.begin(); it!=u_tree[imax].map.end(); ++it) {
        u_loc_single[Index1D(j_bspline+imax,(*it).first,XBSpline)] = (*it).second;
    }
    return;
}

template <typename PrimalBasis>
void
LocalRefinement<PrimalBasis>::reconstruct(const CoefficientsByLevel<T> &u_bspline, FLENS_DEFAULT_INDEXTYPE j_bspline,
                                           const CoefficientsByLevel<T> &u_wavelet, FLENS_DEFAULT_INDEXTYPE j_wavelet,
                                           CoefficientsByLevel<T> &u_loc_single, FLENS_DEFAULT_INDEXTYPE &j_refinement) const
{
    FLENS_DEFAULT_INDEXTYPE j1_refinement = refinementbasis.mra.phi.getRefinementLevel(j_bspline);
    // pre-initialization need if we do not enter the following loop.
    for (typename CoefficientsByLevel<T>::const_it it=u_bspline.map.begin(); it!=u_bspline.map.end(); ++it) {
        this->reconstructBSpline(j_bspline, (*it).first, (*it).second, u_loc_single, j1_refinement);
    }
    FLENS_DEFAULT_INDEXTYPE j2_refinement = basis.psi.getRefinementLevel(j_wavelet);
    // pre-initialization need if we do not enter the following loop.
    for (typename CoefficientsByLevel<T>::const_it it=u_wavelet.map.begin(); it!=u_wavelet.map.end(); ++it) {
        this->reconstructWavelet(j_wavelet, (*it).first, (*it).second, u_loc_single, j2_refinement);
    }
    //assert(j1_refinement==j2_refinement);
    if(j1_refinement!=j2_refinement) {
        std::cerr << "LocalRefinement<PrimalBasis> ERROR: j_bspline = " << j_bspline << ", j_wavelet = " << j_wavelet << std::endl;
        std::cerr << "                                j1_refinement = " << j1_refinement << ", j2_refinement = " << j2_refinement << std::endl;
    }
    j_refinement = j2_refinement;
}

// Non-Periodic Version
template <typename PrimalBasis>
template <typename T_>
typename cxxblas::RestrictTo<SFINAE_Wrapper<!IsPeriodic<PrimalBasis>::value,T_>::value, void>::Type
LocalRefinement<PrimalBasis>::reconstructOnlyMultiScaling
                               (const CoefficientsByLevel<T_> &u_scaling, FLENS_DEFAULT_INDEXTYPE j,
                                CoefficientsByLevel<T_> &u_loc_single, FLENS_DEFAULT_INDEXTYPE &j_refinement) const
{
    DenseVectorLD *refCoeffs;
    FLENS_DEFAULT_INDEXTYPE k_refinement_first = (FLENS_DEFAULT_INDEXTYPE) 0;
    FLENS_DEFAULT_INDEXTYPE split = 100.L;
    FLENS_DEFAULT_INDEXTYPE k_refinement_restart = 100.L;

    for (const_coeffbylevel_it it=u_scaling.map.begin(); it!=u_scaling.map.end(); ++it) {
        refCoeffs = basis.mra.phi.getRefinement(j,(*it).first,j_refinement,k_refinement_first, split, k_refinement_restart);
        for (FLENS_DEFAULT_INDEXTYPE i=(*refCoeffs).firstIndex(); i<= (FLENS_DEFAULT_INDEXTYPE)(*refCoeffs).lastIndex(); ++i) {
            u_loc_single.map[k_refinement_first+i] += (*refCoeffs).operator()(i) * (*it).second;
        }
    }
}

// Periodic Version
template <typename PrimalBasis>
template <typename T_>
typename cxxblas::RestrictTo<SFINAE_Wrapper<IsPeriodic<PrimalBasis>::value,T_>::value, void>::Type
LocalRefinement<PrimalBasis>::reconstructOnlyMultiScaling
                               (const CoefficientsByLevel<T_> &u_scaling, FLENS_DEFAULT_INDEXTYPE j,
                                CoefficientsByLevel<T_> &u_loc_single, FLENS_DEFAULT_INDEXTYPE &j_refinement) const
{
    DenseVectorLD *refCoeffs;
    FLENS_DEFAULT_INDEXTYPE k_refinement_first = (FLENS_DEFAULT_INDEXTYPE) 0;
    FLENS_DEFAULT_INDEXTYPE split = 100.L;
    FLENS_DEFAULT_INDEXTYPE k_refinement_restart = 100.L;

    for (const_coeffbylevel_it it=u_scaling.map.begin(); it!=u_scaling.map.end(); ++it) {
        refCoeffs = basis.mra.phi.getRefinement(j,(*it).first,j_refinement,k_refinement_first, split, k_refinement_restart);
    	// First part of coefficients (= all coefficients, if basis is not periodic
        for (FLENS_DEFAULT_INDEXTYPE i=(*refCoeffs).firstIndex(); i<= std::min((*refCoeffs).firstIndex()+split-1, (FLENS_DEFAULT_INDEXTYPE)(*refCoeffs).lastIndex()); ++i) {
            u_loc_single.map[k_refinement_first+i] += (*refCoeffs).operator()(i) * (*it).second;
        }
    	// Second part of coefficients: i is still index in coefficient vector
        for (FLENS_DEFAULT_INDEXTYPE i= (*refCoeffs).firstIndex()+split; i <= (*refCoeffs).lastIndex(); ++i) {
            u_loc_single.map[k_refinement_restart+i-((*refCoeffs).firstIndex()+split)] += (*refCoeffs).operator()(i) * (*it).second;
    	}
    }
}


template <typename PrimalBasis>
void
LocalRefinement<PrimalBasis>::reconstructBSpline(FLENS_DEFAULT_INDEXTYPE j, FLENS_DEFAULT_INDEXTYPE k, T coeff,
                                                  CoefficientsByLevel<T> &u_loc_single,
                                                  FLENS_DEFAULT_INDEXTYPE &j_refinement) const
{
    DenseVectorLD *refCoeffs;
    j_refinement = 0;
    FLENS_DEFAULT_INDEXTYPE k_refinement_first = (FLENS_DEFAULT_INDEXTYPE) 0;
    FLENS_DEFAULT_INDEXTYPE split = 100.L;
    FLENS_DEFAULT_INDEXTYPE k_refinement_restart = 100.L;
    refCoeffs = refinementbasis.mra.phi.getRefinement(j,k,j_refinement,k_refinement_first, split, k_refinement_restart);

    for (FLENS_DEFAULT_INDEXTYPE i=(*refCoeffs).firstIndex(); i<= (FLENS_DEFAULT_INDEXTYPE)(*refCoeffs).lastIndex(); ++i) {
        u_loc_single.map[k_refinement_first+i] += (*refCoeffs).operator()(i) * coeff;
    }
}

// Non-Periodic Version
template <typename PrimalBasis>
template <typename T_>
typename cxxblas::RestrictTo<SFINAE_Wrapper<!IsPeriodic<PrimalBasis>::value,T_>::value, void>::Type
LocalRefinement<PrimalBasis>::reconstructWavelet(FLENS_DEFAULT_INDEXTYPE j, FLENS_DEFAULT_INDEXTYPE k, T_ coeff,
                                                  CoefficientsByLevel<T_> &u_loc_single,
                                                  FLENS_DEFAULT_INDEXTYPE &j_refinement) const
{
    DenseVectorLD *refCoeffs;
    j_refinement = 0;
    FLENS_DEFAULT_INDEXTYPE k_refinement_first = (FLENS_DEFAULT_INDEXTYPE) 0;
    FLENS_DEFAULT_INDEXTYPE split = 100.L;
    FLENS_DEFAULT_INDEXTYPE k_refinement_restart = 100.L;
    refCoeffs = basis.psi.getRefinement(j,k,j_refinement,k_refinement_first, split, k_refinement_restart);

    for (FLENS_DEFAULT_INDEXTYPE i=(*refCoeffs).firstIndex(); i<= (FLENS_DEFAULT_INDEXTYPE)(*refCoeffs).lastIndex(); ++i) {
        u_loc_single.map[k_refinement_first+i] += (*refCoeffs).operator()(i) * coeff;
    }
}

//Periodic Version
template <typename PrimalBasis>
template <typename T_>
typename cxxblas::RestrictTo<SFINAE_Wrapper<IsPeriodic<PrimalBasis>::value,T_>::value, void>::Type
LocalRefinement<PrimalBasis>::reconstructWavelet(FLENS_DEFAULT_INDEXTYPE j, FLENS_DEFAULT_INDEXTYPE k, T_ coeff,
                                                  CoefficientsByLevel<T_> &u_loc_single,
                                                  FLENS_DEFAULT_INDEXTYPE &j_refinement) const
{
    DenseVectorLD *refCoeffs;
    j_refinement = 0;
    FLENS_DEFAULT_INDEXTYPE k_refinement_first = (FLENS_DEFAULT_INDEXTYPE) 0;
    FLENS_DEFAULT_INDEXTYPE split = 100.L;
    FLENS_DEFAULT_INDEXTYPE k_refinement_restart = 100.L;
    refCoeffs = basis.psi.getRefinement(j,k,j_refinement,k_refinement_first, split, k_refinement_restart);

	// First part of coefficients (= all coefficients, if basis is not periodic
    for (FLENS_DEFAULT_INDEXTYPE i=(*refCoeffs).firstIndex(); i<= std::min((*refCoeffs).firstIndex()+split-1, (FLENS_DEFAULT_INDEXTYPE)(*refCoeffs).lastIndex()); ++i) {
        u_loc_single.map[k_refinement_first+i] += (*refCoeffs).operator()(i) * coeff;
    }
	// Second part of coefficients: i is still index in coefficient vector
    for (FLENS_DEFAULT_INDEXTYPE i= (*refCoeffs).firstIndex()+split; i <= (*refCoeffs).lastIndex(); ++i) {
        u_loc_single.map[k_refinement_restart+i-((*refCoeffs).firstIndex()+split)] += (*refCoeffs).operator()(i) * coeff;
	}
}

template <typename PrimalBasis>
void
LocalRefinement<PrimalBasis>::decompose_(const CoefficientsByLevel<T>  &u_loc_single,
                                          CoefficientsByLevel<T>  &u_bspline, FLENS_DEFAULT_INDEXTYPE j_bspline,
                                          CoefficientsByLevel<T>  &u_wavelet, FLENS_DEFAULT_INDEXTYPE j_wavelet) const
{
    if (u_loc_single.map.size()==0) return;
    for (coeffbylevel_it it=u_bspline.map.begin(); it!=u_bspline.map.end(); ++it) {
        T coeff = this->decompose_BSpline(u_loc_single, j_bspline, (*it).first);
        (*it).second += coeff;
    }
    for (coeffbylevel_it it=u_wavelet.map.begin(); it!=u_wavelet.map.end(); ++it) {
        T coeff = this->decompose_Wavelet(u_loc_single, j_wavelet, (*it).first);
        (*it).second += coeff;
    }
}

template <typename PrimalBasis>
void
LocalRefinement<PrimalBasis>::decompose_OnlyMultiScaling(const CoefficientsByLevel<T>  &u_loc_single,
                                                          CoefficientsByLevel<T>  &u_scaling, FLENS_DEFAULT_INDEXTYPE j_scaling)
                                                          const
{
    if (u_loc_single.map.size()==0) return;
    for (coeffbylevel_it it=u_scaling.map.begin(); it!=u_scaling.map.end(); ++it) {
        T coeff = this->decompose_Scaling(u_loc_single, j_scaling, (*it).first);
        (*it).second += coeff;
    }
}


// Non-Periodic Version
template <typename PrimalBasis>
template <typename T_>
typename cxxblas::RestrictTo<SFINAE_Wrapper<!IsPeriodic<PrimalBasis>::value,T_>::value, T_>::Type
LocalRefinement<PrimalBasis>::decompose_Scaling(const CoefficientsByLevel<T_> &u_loc_single,
                                                 FLENS_DEFAULT_INDEXTYPE j, FLENS_DEFAULT_INDEXTYPE k) const
{
    const_coeffbylevel_it u_loc_single_end = u_loc_single.map.end();
    const_coeffbylevel_it u_loc_single_ptr;
    DenseVectorLD *refCoeffs;
    FLENS_DEFAULT_INDEXTYPE refinement_j = 0;
    FLENS_DEFAULT_INDEXTYPE refinement_k_first = (FLENS_DEFAULT_INDEXTYPE) 0;
    FLENS_DEFAULT_INDEXTYPE split = 100L;
    FLENS_DEFAULT_INDEXTYPE refinement_k_restart = (FLENS_DEFAULT_INDEXTYPE) 1;
    T val = 0.;
    refCoeffs = basis.mra.phi.getRefinement(j,k,refinement_j,refinement_k_first, split, refinement_k_restart);

    for (FLENS_DEFAULT_INDEXTYPE i=(*refCoeffs).firstIndex(); i<=(*refCoeffs).lastIndex(); ++i) {
        u_loc_single_ptr=u_loc_single.map.find(refinement_k_first+i);
        if (u_loc_single_ptr!=u_loc_single_end) {
            val += (*refCoeffs).operator()(i) * (*u_loc_single_ptr).second;
        }
    }

    return val;
}


// Periodic Version
template <typename PrimalBasis>
template <typename T_>
typename cxxblas::RestrictTo<SFINAE_Wrapper<IsPeriodic<PrimalBasis>::value,T_>::value, T_>::Type
LocalRefinement<PrimalBasis>::decompose_Scaling(const CoefficientsByLevel<T_> &u_loc_single,
                                                 FLENS_DEFAULT_INDEXTYPE j, FLENS_DEFAULT_INDEXTYPE k) const
{
    const_coeffbylevel_it u_loc_single_end = u_loc_single.map.end();
    const_coeffbylevel_it u_loc_single_ptr;
    DenseVectorLD *refCoeffs;
    FLENS_DEFAULT_INDEXTYPE refinement_j = 0;
    FLENS_DEFAULT_INDEXTYPE refinement_k_first = (FLENS_DEFAULT_INDEXTYPE) 0;
    FLENS_DEFAULT_INDEXTYPE split = 100L;
    FLENS_DEFAULT_INDEXTYPE refinement_k_restart = (FLENS_DEFAULT_INDEXTYPE) 1;
    T val = 0.;
    refCoeffs = basis.mra.phi.getRefinement(j,k,refinement_j,refinement_k_first, split, refinement_k_restart);

	// First part of coefficients
	for (FLENS_DEFAULT_INDEXTYPE i=(*refCoeffs).firstIndex(); i<= std::min( ((*refCoeffs).firstIndex()+split-1) , (FLENS_DEFAULT_INDEXTYPE)(*refCoeffs).lastIndex()); ++i) {
    	u_loc_single_ptr=u_loc_single.map.find(refinement_k_first+i);
        if (u_loc_single_ptr!=u_loc_single_end) {
            val += (*refCoeffs).operator()(i) * (*u_loc_single_ptr).second;
        }
	}
	// Second part of coefficients: i is still index in coefficient vector
	for (FLENS_DEFAULT_INDEXTYPE i= (*refCoeffs).firstIndex()+split; i <= (*refCoeffs).lastIndex(); ++i) {
		FLENS_DEFAULT_INDEXTYPE index = refinement_k_restart+i-((*refCoeffs).firstIndex()+split);
        u_loc_single_ptr=u_loc_single.map.find(index);
        if (u_loc_single_ptr!=u_loc_single_end) {
            val += (*refCoeffs).operator()(i) * (*u_loc_single_ptr).second;
        }
	}
    return val;
}

template <typename PrimalBasis>
typename PrimalBasis::T
LocalRefinement<PrimalBasis>::decompose_BSpline(const CoefficientsByLevel<T> &u_loc_single,
                                                 FLENS_DEFAULT_INDEXTYPE j, FLENS_DEFAULT_INDEXTYPE k) const
{
    const_coeffbylevel_it u_loc_single_end = u_loc_single.map.end();
    const_coeffbylevel_it u_loc_single_ptr;
    DenseVectorLD *refCoeffs;
    FLENS_DEFAULT_INDEXTYPE refinement_j = 0;
    FLENS_DEFAULT_INDEXTYPE refinement_k_first = (FLENS_DEFAULT_INDEXTYPE) 0;
    FLENS_DEFAULT_INDEXTYPE split = 100L;
    FLENS_DEFAULT_INDEXTYPE refinement_k_restart = (FLENS_DEFAULT_INDEXTYPE) 0;
    T val = 0.;
    refCoeffs = refinementbasis.mra.phi.getRefinement(j,k,refinement_j,refinement_k_first, split, refinement_k_restart);
    for (FLENS_DEFAULT_INDEXTYPE i=(*refCoeffs).firstIndex(); i<=(*refCoeffs).lastIndex(); ++i) {
        u_loc_single_ptr=u_loc_single.map.find(refinement_k_first+i);
        if (u_loc_single_ptr!=u_loc_single_end) {
            val += (*refCoeffs).operator()(i) * (*u_loc_single_ptr).second;
        }
    }
    return val;
}

// Non-Periodic Version
template <typename PrimalBasis>
template <typename T_>
typename cxxblas::RestrictTo<SFINAE_Wrapper<!IsPeriodic<PrimalBasis>::value,T_>::value, T_>::Type
LocalRefinement<PrimalBasis>::decompose_Wavelet(const CoefficientsByLevel<T_> &u_loc_single,
                                                 FLENS_DEFAULT_INDEXTYPE j, FLENS_DEFAULT_INDEXTYPE k) const
{
    const_coeffbylevel_it u_loc_single_end = u_loc_single.map.end();
    const_coeffbylevel_it u_loc_single_ptr;
    DenseVectorLD *refCoeffs;
    FLENS_DEFAULT_INDEXTYPE refinement_j = 0;
    FLENS_DEFAULT_INDEXTYPE refinement_k_first = (FLENS_DEFAULT_INDEXTYPE) 0;
    FLENS_DEFAULT_INDEXTYPE split = 100L;
    FLENS_DEFAULT_INDEXTYPE refinement_k_restart = (FLENS_DEFAULT_INDEXTYPE) 0;
    T val = 0.;
    refCoeffs = basis.psi.getRefinement(j,k,refinement_j,refinement_k_first,split,refinement_k_restart);

    for (FLENS_DEFAULT_INDEXTYPE i=(*refCoeffs).firstIndex(); i<=(*refCoeffs).lastIndex(); ++i) {
        u_loc_single_ptr=u_loc_single.map.find(refinement_k_first+i);
        if (u_loc_single_ptr!=u_loc_single_end) {
            val += (*refCoeffs).operator()(i) * (*u_loc_single_ptr).second;
        }
    }

    return val;
}


// Periodic Version
template <typename PrimalBasis>
template <typename T_>
typename cxxblas::RestrictTo<SFINAE_Wrapper<IsPeriodic<PrimalBasis>::value,T_>::value, T_>::Type
LocalRefinement<PrimalBasis>::decompose_Wavelet(const CoefficientsByLevel<T_> &u_loc_single,
                                                 FLENS_DEFAULT_INDEXTYPE j, FLENS_DEFAULT_INDEXTYPE k) const
{
    const_coeffbylevel_it u_loc_single_end = u_loc_single.map.end();
    const_coeffbylevel_it u_loc_single_ptr;
    DenseVectorLD *refCoeffs;
    FLENS_DEFAULT_INDEXTYPE refinement_j = 0;
    FLENS_DEFAULT_INDEXTYPE refinement_k_first = (FLENS_DEFAULT_INDEXTYPE) 0;
    FLENS_DEFAULT_INDEXTYPE split = 100L;
    FLENS_DEFAULT_INDEXTYPE refinement_k_restart = (FLENS_DEFAULT_INDEXTYPE) 0;
    T val = 0.;
    refCoeffs = basis.psi.getRefinement(j,k,refinement_j,refinement_k_first,split,refinement_k_restart);

    /*for (FLENS_DEFAULT_INDEXTYPE i=(*refCoeffs).firstIndex(); i<=(*refCoeffs).lastIndex(); ++i) {
        u_loc_single_ptr=u_loc_single.map.find(refinement_k_first+i);
        if (u_loc_single_ptr!=u_loc_single_end) {
            val += (*refCoeffs).operator()(i) * (*u_loc_single_ptr).second;
        }
    }*/

	// First part of coefficients
	for (FLENS_DEFAULT_INDEXTYPE i=(*refCoeffs).firstIndex(); i<= std::min( ((*refCoeffs).firstIndex()+split-1) , (FLENS_DEFAULT_INDEXTYPE)(*refCoeffs).lastIndex()); ++i) {
    	u_loc_single_ptr=u_loc_single.map.find(refinement_k_first+i);
        if (u_loc_single_ptr!=u_loc_single_end) {
            val += (*refCoeffs).operator()(i) * (*u_loc_single_ptr).second;
        }
	}
	// Second part of coefficients: i is still index in coefficient vector
	for (FLENS_DEFAULT_INDEXTYPE i= (*refCoeffs).firstIndex()+split; i <= (*refCoeffs).lastIndex(); ++i) {
		FLENS_DEFAULT_INDEXTYPE index = refinement_k_restart+i-((*refCoeffs).firstIndex()+split);
        u_loc_single_ptr=u_loc_single.map.find(index);
        if (u_loc_single_ptr!=u_loc_single_end) {
            val += (*refCoeffs).operator()(i) * (*u_loc_single_ptr).second;
        }
	}
    return val;
}

}   // namespace lawa
