namespace lawa {

template <typename T, typename Basis2D>
void
getSparseGridIndexSet(const Basis2D &basis, IndexSet<Index2D> &Lambda, FLENS_DEFAULT_INDEXTYPE j, FLENS_DEFAULT_INDEXTYPE deltaL, T gamma)
{
    FLENS_DEFAULT_INDEXTYPE j0_1 = basis.first.j0;
    FLENS_DEFAULT_INDEXTYPE j0_2 = basis.second.j0;
    for (FLENS_DEFAULT_INDEXTYPE k1=basis.first.mra.rangeI(j0_1).firstIndex(); k1<=basis.first.mra.rangeI(j0_1).lastIndex(); ++k1) {
        for (FLENS_DEFAULT_INDEXTYPE k2=basis.second.mra.rangeI(j0_2).firstIndex(); k2<=basis.second.mra.rangeI(j0_2).lastIndex(); ++k2) {
            Index1D row(j0_1,k1,XBSpline);
            Index1D col(j0_2,k2,XBSpline);
            Lambda.insert(Index2D(row,col));
        }
        for (FLENS_DEFAULT_INDEXTYPE i2=1; i2<=j; ++i2) {
            FLENS_DEFAULT_INDEXTYPE j2=j0_2+i2-1;
            for (FLENS_DEFAULT_INDEXTYPE k2=basis.second.rangeJ(j2).firstIndex(); k2<=basis.second.rangeJ(j2).lastIndex(); ++k2) {
                Index1D row(j0_1,k1,XBSpline);
                Index1D col(j2,k2,XWavelet);
                Lambda.insert(Index2D(row,col));
            }
        }
    }
    for (FLENS_DEFAULT_INDEXTYPE k2=basis.second.mra.rangeI(j0_2).firstIndex(); k2<=basis.second.mra.rangeI(j0_2).lastIndex(); ++k2) {
        for (FLENS_DEFAULT_INDEXTYPE i1=1; i1<=j+deltaL; ++i1) {
            FLENS_DEFAULT_INDEXTYPE j1=j0_1+i1-1;
            for (FLENS_DEFAULT_INDEXTYPE k1=basis.first.rangeJ(j1).firstIndex(); k1<=basis.first.rangeJ(j1).lastIndex(); ++k1) {
                Index1D row(j1,k1,XWavelet);
                Index1D col(j0_2,k2,XBSpline);
                Lambda.insert(Index2D(row,col));
            }
        }
    }
    for (FLENS_DEFAULT_INDEXTYPE i1=1; i1<=j+deltaL; ++i1) {
        FLENS_DEFAULT_INDEXTYPE j1=j0_1+i1-1;
        for (FLENS_DEFAULT_INDEXTYPE i2=1; i1+i2<=j; ++i2) {
            if (T(i1-deltaL+i2)-gamma*std::max(i1-deltaL,i2)>(1-gamma)*j) continue;
            FLENS_DEFAULT_INDEXTYPE j2=j0_2+i2-1;
            for (FLENS_DEFAULT_INDEXTYPE k1=basis.first.rangeJ(j1).firstIndex(); k1<=basis.first.rangeJ(j1).lastIndex(); ++k1) {
                for (FLENS_DEFAULT_INDEXTYPE k2=basis.second.rangeJ(j2).firstIndex(); k2<=basis.second.rangeJ(j2).lastIndex(); ++k2) {
                    Index1D row(j1,k1,XWavelet);
                    Index1D col(j2,k2,XWavelet);
                    Lambda.insert(Index2D(row,col));
                }
            }
        }
    }
    return;
}

template <typename Basis1D>
void
getFullIndexSet(const Basis1D &basis, IndexSet<Index1D> &Lambda, FLENS_DEFAULT_INDEXTYPE J)
{
	FLENS_DEFAULT_INDEXTYPE j0 = basis.j0;

	// Scaling Functions
	for(FLENS_DEFAULT_INDEXTYPE k = basis.mra.rangeI(j0).firstIndex(); k <= basis.mra.rangeI(j0).lastIndex(); ++k){
		Lambda.insert(Index1D(j0, k, XBSpline));
	}

	// Wavelets
	for(FLENS_DEFAULT_INDEXTYPE j = j0; j < J; ++j){
		for(FLENS_DEFAULT_INDEXTYPE k = basis.rangeJ(j).firstIndex(); k <= basis.rangeJ(j).lastIndex(); ++k){
			Lambda.insert(Index1D(j, k, XWavelet));
		}
	}
}



template <typename Basis2D>
void
getFullIndexSet(const Basis2D &basis, IndexSet<Index2D> &Lambda, FLENS_DEFAULT_INDEXTYPE J1, FLENS_DEFAULT_INDEXTYPE J2, FLENS_DEFAULT_INDEXTYPE deltaL)
{
    FLENS_DEFAULT_INDEXTYPE j0_1 = basis.first.j0;
    FLENS_DEFAULT_INDEXTYPE j0_2 = basis.second.j0;
    for (FLENS_DEFAULT_INDEXTYPE k1=basis.first.mra.rangeI(j0_1).firstIndex(); k1<=basis.first.mra.rangeI(j0_1).lastIndex(); ++k1) {
    	// Scaling x Scaling
        for (FLENS_DEFAULT_INDEXTYPE k2=basis.second.mra.rangeI(j0_2).firstIndex(); k2<=basis.second.mra.rangeI(j0_2).lastIndex(); ++k2) {
            Index1D row(j0_1,k1,XBSpline);
            Index1D col(j0_2,k2,XBSpline);
            Lambda.insert(Index2D(row,col));
        }
    	// Scaling x Wavelet
        for (FLENS_DEFAULT_INDEXTYPE j2= j0_2; j2 <= J2; ++j2) {
            for (FLENS_DEFAULT_INDEXTYPE k2=basis.second.rangeJ(j2).firstIndex(); k2<=basis.second.rangeJ(j2).lastIndex(); ++k2) {
                Index1D row(j0_1,k1,XBSpline);
                Index1D col(j2,k2,XWavelet);
                Lambda.insert(Index2D(row,col));
            }
        }
    }
    for (FLENS_DEFAULT_INDEXTYPE k2=basis.second.mra.rangeI(j0_2).firstIndex(); k2<=basis.second.mra.rangeI(j0_2).lastIndex(); ++k2) {
    	// Wavelet x Scaling
        for (FLENS_DEFAULT_INDEXTYPE j1=j0_1; j1<= J1+deltaL; ++j1) {
            for (FLENS_DEFAULT_INDEXTYPE k1=basis.first.rangeJ(j1).firstIndex(); k1<=basis.first.rangeJ(j1).lastIndex(); ++k1) {
                Index1D row(j1,k1,XWavelet);
                Index1D col(j0_2,k2,XBSpline);
                Lambda.insert(Index2D(row,col));
            }
        }
    }
	// Wavelet x Wavelet
    for (FLENS_DEFAULT_INDEXTYPE j1=j0_1; j1<= J1+deltaL; ++j1) {
        for (FLENS_DEFAULT_INDEXTYPE j2= j0_2; j2 <= J2; ++j2) {
            for (FLENS_DEFAULT_INDEXTYPE k1=basis.first.rangeJ(j1).firstIndex(); k1<=basis.first.rangeJ(j1).lastIndex(); ++k1) {
                for (FLENS_DEFAULT_INDEXTYPE k2=basis.second.rangeJ(j2).firstIndex(); k2<=basis.second.rangeJ(j2).lastIndex(); ++k2) {
                    Index1D row(j1,k1,XWavelet);
                    Index1D col(j2,k2,XWavelet);
                    Lambda.insert(Index2D(row,col));
                }
            }
        }
    }
    return;
}

template <typename Basis2D>
void
getScalingFctIndexSet(const Basis2D &basis, IndexSet<Index2D> &Lambda, FLENS_DEFAULT_INDEXTYPE J1, FLENS_DEFAULT_INDEXTYPE J2)
{
    for (FLENS_DEFAULT_INDEXTYPE k1=basis.first.mra.rangeI(J1).firstIndex(); k1<=basis.first.mra.rangeI(J1).lastIndex(); ++k1) {
    	// Scaling x Scaling
        for (FLENS_DEFAULT_INDEXTYPE k2=basis.second.mra.rangeI(J2).firstIndex(); k2<=basis.second.mra.rangeI(J2).lastIndex(); ++k2) {
            Index1D row(J1,k1,XBSpline);
            Index1D col(J2,k2,XBSpline);
            Lambda.insert(Index2D(row,col));
        }
    }
}

} // namespace lawa
