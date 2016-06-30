namespace lawa {

template <typename T, typename Basis2D, typename S1_x, typename S1_y, typename S2_x, typename S2_y>
TensorSparseGrid2D<T, Basis2D, S1_x, S1_y, S2_x, S2_y>::TensorSparseGrid2D
                                                        (const Basis2D &basis,
                                                         const S1_x &s1_x, const S1_y &s1_y,
                                                         const S2_x &s2_x, const S2_y &s2_y,
                                                         FLENS_DEFAULT_INDEXTYPE I, T gamma)
: _basis(basis), _s1_x(s1_x), _s1_y(s1_y), _s2_x(s2_x), _s2_y(s2_y),
  _j0_x(basis.first.j0), _j0_y(basis.second.j0), _I(I), _gamma(gamma),
  _blockassembler1d(basis.first), _dim(0)
{
    FLENS_DEFAULT_INDEXTYPE fi=1, li=1;

    for (FLENS_DEFAULT_INDEXTYPE ix=0; ix<=_I; ++ix) {
        for (FLENS_DEFAULT_INDEXTYPE iy=0; iy<=_I; ++iy) {

            //if (ix+iy>_I) break;
            if (ix+iy> (1-_gamma)*_I + _gamma*std::max(ix,iy)) break;

            FLENS_DEFAULT_INDEXTYPE *pair = new FLENS_DEFAULT_INDEXTYPE[6];
            pair[0] = ix; pair[1] = iy;

            if (ix==0) {
                pair[4] = _basis.first.mra.cardI(_j0_x);
                if (iy==0) {
                    li = fi+_basis.first.mra.cardI(_j0_x)*_basis.second.mra.cardI(_j0_y)-1;
                    pair[5] = _basis.second.mra.cardI(_j0_y);
                }
                else {
                    li = fi+_basis.first.mra.cardI(_j0_x)*_basis.second.cardJ(_j0_y+iy-1)-1;
                    pair[5] = _basis.second.cardJ(_j0_y+iy-1);
                }
            }
            else {
                pair[4] = _basis.first.cardJ(_j0_x+ix-1);
                if (iy==0) {
                    li = fi+_basis.first.cardJ(_j0_x+ix-1)*_basis.second.mra.cardI(_j0_y)-1;
                    pair[5] = _basis.second.mra.cardI(_j0_y);
                }
                else {
                    li = fi+_basis.first.cardJ(_j0_x+ix-1)*_basis.second.cardJ(_j0_y+iy-1)-1;
                    pair[5] = _basis.second.cardJ(_j0_y+iy-1);
                }
            }
            pair[2] = fi; pair[3] = li;
            fi = li+1;
            _dim += pair[4] * pair[5];
            _sg_blocks.push_back(pair);
        }
    }

    FLENS_DEFAULT_INDEXTYPE levelpair_pos = 0;
    for (FLENS_DEFAULT_INDEXTYPE level1=0; level1<=_I; ++level1) {
        for (FLENS_DEFAULT_INDEXTYPE level2=0; level2<=_I; ++level2) {
            _levelpair_map[std::pair<FLENS_DEFAULT_INDEXTYPE,FLENS_DEFAULT_INDEXTYPE>(level1,level2)] = levelpair_pos;
            ++levelpair_pos;
        }
    }
    this->assembleMatrices();
}

template<typename T, typename Basis2D, typename S1_x, typename S1_y, typename S2_x, typename S2_y>
FLENS_DEFAULT_INDEXTYPE
TensorSparseGrid2D<T,Basis2D,S1_x,S1_y,S2_x,S2_y>::getDimension() const
{
    return _dim;
}

template<typename T, typename Basis2D, typename S1_x, typename S1_y, typename S2_x, typename S2_y>
FLENS_DEFAULT_INDEXTYPE
TensorSparseGrid2D<T,Basis2D,S1_x,S1_y,S2_x,S2_y>::numCols() const
{
    return _dim;
}

template<typename T, typename Basis2D, typename S1_x, typename S1_y, typename S2_x, typename S2_y>
FLENS_DEFAULT_INDEXTYPE
TensorSparseGrid2D<T,Basis2D,S1_x,S1_y,S2_x,S2_y>::numRows() const
{
    return _dim;
}

template<typename T, typename Basis2D, typename S1_x, typename S1_y, typename S2_x, typename S2_y>
IndexSet<Index2D>
TensorSparseGrid2D<T,Basis2D,S1_x,S1_y,S2_x,S2_y>::getIndexSet() const
{
    IndexSet<Index2D> sparsegridindexset;
    for (FLENS_DEFAULT_INDEXTYPE block=0; block<(FLENS_DEFAULT_INDEXTYPE)_sg_blocks.size(); ++block) {
        FLENS_DEFAULT_INDEXTYPE i1 = _sg_blocks[block][0]-1;
        FLENS_DEFAULT_INDEXTYPE i2 = _sg_blocks[block][1]-1;

        if (i1<0) {
            for (FLENS_DEFAULT_INDEXTYPE k1=_basis.first.mra.rangeI(_j0_x).firstIndex(); k1<=_basis.first.mra.rangeI(_j0_x).lastIndex(); ++k1) {
                if (i2<0) {
                    for (FLENS_DEFAULT_INDEXTYPE k2=_basis.second.mra.rangeI(_j0_y).firstIndex(); k2<=_basis.second.mra.rangeI(_j0_y).lastIndex(); ++k2) {
                        Index2D index2d(Index1D(_j0_x,k1,XBSpline),Index1D(_j0_y,k2,XBSpline));
                        sparsegridindexset.insert(index2d);
                    }
                }
                else {
                    for (FLENS_DEFAULT_INDEXTYPE k2=_basis.second.rangeJ(_j0_y+i2).firstIndex(); k2<=_basis.second.rangeJ(_j0_y+i2).lastIndex(); ++k2) {
                        Index2D index2d(Index1D(_j0_x,k1,XBSpline),Index1D(_j0_y+i2,k2,XWavelet));
                        sparsegridindexset.insert(index2d);
                    }
                }
            }
        }
        else {
            for (FLENS_DEFAULT_INDEXTYPE k1=_basis.first.rangeJ(_j0_x+i1).firstIndex(); k1<=_basis.first.rangeJ(_j0_x+i1).lastIndex(); ++k1) {
                if (i2<0) {
                    for (FLENS_DEFAULT_INDEXTYPE k2=_basis.second.mra.rangeI(_j0_y).firstIndex(); k2<=_basis.second.mra.rangeI(_j0_y).lastIndex(); ++k2) {
                        Index2D index2d(Index1D(_j0_x+i1,k1,XWavelet),Index1D(_j0_y,k2,XBSpline));
                        sparsegridindexset.insert(index2d);
                    }
                }
                else {
                    for (FLENS_DEFAULT_INDEXTYPE k2=_basis.second.rangeJ(_j0_y+i2).firstIndex(); k2<=_basis.second.rangeJ(_j0_y+i2).lastIndex(); ++k2) {
                        Index2D index2d(Index1D(_j0_x+i1,k1,XWavelet),Index1D(_j0_y+i2,k2,XWavelet));
                        sparsegridindexset.insert(index2d);
                    }
                }
            }
        }
    }
    return sparsegridindexset;
}

template <typename T, typename Basis2D, typename S1_x, typename S1_y, typename S2_x, typename S2_y>
void
TensorSparseGrid2D<T, Basis2D, S1_x, S1_y, S2_x, S2_y>::toCoefficients
                                   (const DenseVectorT &vec,
                                    Coefficients<Lexicographical,T,Index2D> &sparsegridcoefficients)
{
    FLENS_DEFAULT_INDEXTYPE pos=1;
    for (FLENS_DEFAULT_INDEXTYPE block=0; block<(FLENS_DEFAULT_INDEXTYPE)_sg_blocks.size(); ++block) {
        FLENS_DEFAULT_INDEXTYPE ix = _sg_blocks[block][0]-1;
        FLENS_DEFAULT_INDEXTYPE iy = _sg_blocks[block][1]-1;

        if (ix<0) {
            for (FLENS_DEFAULT_INDEXTYPE k1=_basis.first.mra.rangeI(_j0_x).firstIndex(); k1<=_basis.first.mra.rangeI(_j0_x).lastIndex(); ++k1) {
                if (iy<0) {
                    for (FLENS_DEFAULT_INDEXTYPE k2=_basis.second.mra.rangeI(_j0_y).firstIndex(); k2<=_basis.second.mra.rangeI(_j0_y).lastIndex(); ++k2) {
                        Index2D index2d(Index1D(_j0_x,k1,XBSpline),Index1D(_j0_y,k2,XBSpline));
                        sparsegridcoefficients[index2d] = vec(pos);
                        ++pos;
                    }
                }
                else {
                    for (FLENS_DEFAULT_INDEXTYPE k2=_basis.second.rangeJ(_j0_y+iy).firstIndex(); k2<=_basis.second.rangeJ(_j0_y+iy).lastIndex(); ++k2) {
                        Index2D index2d(Index1D(_j0_x,k1,XBSpline),Index1D(_j0_y+iy,k2,XWavelet));
                        sparsegridcoefficients[index2d] = vec(pos);
                        ++pos;
                    }
                }
            }
        }
        else {
            for (FLENS_DEFAULT_INDEXTYPE k1=_basis.first.rangeJ(_j0_x+ix).firstIndex(); k1<=_basis.first.rangeJ(_j0_x+ix).lastIndex(); ++k1) {
                if (iy<0) {
                    for (FLENS_DEFAULT_INDEXTYPE k2=_basis.second.mra.rangeI(_j0_y).firstIndex(); k2<=_basis.second.mra.rangeI(_j0_y).lastIndex(); ++k2) {
                        Index2D index2d(Index1D(_j0_x+ix,k1,XWavelet),Index1D(_j0_y,k2,XBSpline));
                        sparsegridcoefficients[index2d] = vec(pos);
                        ++pos;
                    }
                }
                else {
                    for (FLENS_DEFAULT_INDEXTYPE k2=_basis.second.rangeJ(_j0_y+iy).firstIndex(); k2<=_basis.second.rangeJ(_j0_y+iy).lastIndex(); ++k2) {
                        Index2D index2d(Index1D(_j0_x+ix,k1,XWavelet),Index1D(_j0_y+iy,k2,XWavelet));
                        sparsegridcoefficients[index2d] = vec(pos);
                        ++pos;
                    }
                }
            }
        }
    }
}

template <typename T, typename Basis2D, typename S1_x, typename S1_y, typename S2_x, typename S2_y>
T
TensorSparseGrid2D<T, Basis2D, S1_x, S1_y, S2_x, S2_y>::evaluate(const DenseVectorT &u, T x, T y,
                                                                 FLENS_DEFAULT_INDEXTYPE deriv_x, FLENS_DEFAULT_INDEXTYPE deriv_y)
const
{
    FLENS_DEFAULT_INDEXTYPE pos=1;
    T ret=0.;
    for (FLENS_DEFAULT_INDEXTYPE block=0; block<(FLENS_DEFAULT_INDEXTYPE)_sg_blocks.size(); ++block) {
        FLENS_DEFAULT_INDEXTYPE ix=_sg_blocks[block][0]-1, iy = _sg_blocks[block][1]-1;
        //FLENS_DEFAULT_INDEXTYPE dim_x = _sg_blocks[block][4], dim_y = _sg_blocks[block][5];

        if (ix<0) {
            if (iy<0) {    //SF*SF
                for(FLENS_DEFAULT_INDEXTYPE kx = _basis.first.mra.rangeI(_j0_x).firstIndex(); kx <= _basis.first.mra.rangeI(_j0_x).lastIndex(); ++kx){
                    for(FLENS_DEFAULT_INDEXTYPE ky = _basis.second.mra.rangeI(_j0_y).firstIndex(); ky <= _basis.second.mra.rangeI(_j0_y).lastIndex(); ++ky){
                        ret += u(pos) * _basis.first.generator(XBSpline)(x,_j0_x,kx,deriv_x)
                                      * _basis.second.generator(XBSpline)(y,_j0_y,ky,deriv_y);
                        ++pos;
                    }
                }
            }
            else {    // SF * W
                for(FLENS_DEFAULT_INDEXTYPE kx = _basis.first.mra.rangeI(_j0_x).firstIndex(); kx <= _basis.first.mra.rangeI(_j0_x).lastIndex(); ++kx){
                    for(FLENS_DEFAULT_INDEXTYPE ky = _basis.second.rangeJ(_j0_y+iy).firstIndex(); ky <= _basis.second.rangeJ(_j0_y+iy).lastIndex(); ++ky){
                        ret += u(pos) * _basis.first.generator(XBSpline)(x,_j0_x,kx,deriv_x)
                                      * _basis.second.generator(XWavelet)(y,_j0_y+iy,ky,deriv_y);
                        ++pos;
                    }
                }
            }
        }

        else {
            if (iy<0) { // W * SF
                for(FLENS_DEFAULT_INDEXTYPE kx = _basis.first.rangeJ(_j0_x+ix).firstIndex(); kx <= _basis.first.rangeJ(_j0_x+ix).lastIndex(); ++kx){
                    for(FLENS_DEFAULT_INDEXTYPE ky = _basis.second.mra.rangeI(_j0_y).firstIndex(); ky <= _basis.second.mra.rangeI(_j0_y).lastIndex(); ++ky){
                        ret += u(pos) * _basis.first.generator(XWavelet)(x,_j0_x+ix,kx,deriv_x)
                                      * _basis.second.generator(XBSpline)(y,_j0_y,ky,deriv_y);
                        ++pos;
                    }
                }
            }
            else {    // W * W
                for(FLENS_DEFAULT_INDEXTYPE kx = _basis.first.rangeJ(_j0_x+ix).firstIndex(); kx <= _basis.first.rangeJ(_j0_x+ix).lastIndex(); ++kx){
                    for(FLENS_DEFAULT_INDEXTYPE ky = _basis.second.rangeJ(_j0_y+iy).firstIndex(); ky <= _basis.second.rangeJ(_j0_y+iy).lastIndex(); ++ky){
                        ret += u(pos) * _basis.first.generator(XWavelet)(x,_j0_x+ix,kx,deriv_x)
                                      * _basis.second.generator(XWavelet)(y,_j0_y+iy,ky,deriv_y);
                        ++pos;
                    }
                }
            }
        }
    }
    return ret;
}

template <typename T, typename Basis2D, typename S1_x, typename S1_y, typename S2_x, typename S2_y>
flens::DiagonalMatrix<T>
TensorSparseGrid2D<T, Basis2D, S1_x, S1_y, S2_x, S2_y>::assembleDiagonalMatrixPreconditioner()
{
    std::vector<DenseVectorT> diags_s1_x, diags_s1_y, diags_s2_x, diags_s2_y;

    for (FLENS_DEFAULT_INDEXTYPE l=0; l<=_I; ++l) {
        DenseVectorT diag_s1_x = _blockassembler1d.assembleStiffnessMatrixBlockDiagonal(_s1_x, l-1, 0.);
        DenseVectorT diag_s1_y = _blockassembler1d.assembleStiffnessMatrixBlockDiagonal(_s1_y, l-1, 0.);
        DenseVectorT diag_s2_x = _blockassembler1d.assembleStiffnessMatrixBlockDiagonal(_s2_x, l-1, 0.);
        DenseVectorT diag_s2_y = _blockassembler1d.assembleStiffnessMatrixBlockDiagonal(_s2_y, l-1, 0.);
        diags_s1_x.push_back(diag_s1_x);
        diags_s1_y.push_back(diag_s1_y);
        diags_s2_x.push_back(diag_s2_x);
        diags_s2_y.push_back(diag_s2_y);
    }

    flens::DenseVector<flens::Array<T> > prec(_dim);

    FLENS_DEFAULT_INDEXTYPE pos=1;
    for (FLENS_DEFAULT_INDEXTYPE block=0; block<(FLENS_DEFAULT_INDEXTYPE)_sg_blocks.size(); ++block) {
        FLENS_DEFAULT_INDEXTYPE ix=_sg_blocks[block][0], iy = _sg_blocks[block][1];
        DenseVectorT diag_s1_x = diags_s1_x[ix];
        DenseVectorT diag_s1_y = diags_s1_y[iy];
        DenseVectorT diag_s2_x = diags_s2_x[ix];
        DenseVectorT diag_s2_y = diags_s2_y[iy];

        FLENS_DEFAULT_INDEXTYPE dim_x = _sg_blocks[block][4], dim_y = _sg_blocks[block][5];

        for (FLENS_DEFAULT_INDEXTYPE i=1; i<=dim_x; ++i) {
            for (FLENS_DEFAULT_INDEXTYPE j=1; j<=dim_y; ++j) {
                prec(pos) = 1./sqrt(fabs(diag_s1_x(i)*diag_s1_y(j) + diag_s2_x(i)*diag_s2_y(j)));
                ++pos;
            }
        }
    }

    flens::DiagonalMatrix<T> diagprec(prec);

    return diagprec;
}

template <typename T, typename Basis2D, typename S1_x, typename S1_y, typename S2_x, typename S2_y>
template <typename RHSIntegral_x, typename RHSIntegral_y>
flens::DenseVector<flens::Array<T> >
TensorSparseGrid2D<T, Basis2D, S1_x, S1_y, S2_x, S2_y>::assembleRHS
                                                       (const RHSIntegral_x &rhs_integral_x,
                                                        const RHSIntegral_y &rhs_integral_y)
{
    std::vector<DenseVectorT> rhss_x, rhss_y;

    for (FLENS_DEFAULT_INDEXTYPE l=0; l<=_I; ++l) {
        DenseVectorT rhs_x = _blockassembler1d.assembleRHSBlock(rhs_integral_x, l-1);
        DenseVectorT rhs_y = _blockassembler1d.assembleRHSBlock(rhs_integral_y, l-1);
        rhss_x.push_back(rhs_x);
        rhss_y.push_back(rhs_y);
    }

    flens::DenseVector<flens::Array<T> > f(_dim);

    FLENS_DEFAULT_INDEXTYPE pos=1;
    for (FLENS_DEFAULT_INDEXTYPE block=0; block<(FLENS_DEFAULT_INDEXTYPE)_sg_blocks.size(); ++block) {
        FLENS_DEFAULT_INDEXTYPE ix=_sg_blocks[block][0], iy = _sg_blocks[block][1];
        DenseVectorT rhs_x = rhss_x[ix];
        DenseVectorT rhs_y = rhss_y[iy];

        FLENS_DEFAULT_INDEXTYPE dim_x = _sg_blocks[block][4], dim_y = _sg_blocks[block][5];

        for (FLENS_DEFAULT_INDEXTYPE i=1; i<=dim_x; ++i) {
            for (FLENS_DEFAULT_INDEXTYPE j=1; j<=dim_y; ++j) {
                f(pos) = rhs_x(i)*rhs_y(j);
                ++pos;
            }
        }
    }

    return f;
}

template<typename T, typename Basis2D, typename S1_x, typename S1_y, typename S2_x, typename S2_y>
flens::DenseVector<flens::Array<T> >
TensorSparseGrid2D<T,Basis2D,S1_x,S1_y,S2_x,S2_y>::operator*(const DenseVectorT &x) const
{
    using flens::_;

    std::vector<DenseMatrixT> y;
    for (FLENS_DEFAULT_INDEXTYPE row_block=0; row_block<(FLENS_DEFAULT_INDEXTYPE)_sg_blocks.size(); ++row_block) {
        FLENS_DEFAULT_INDEXTYPE dim_x = _sg_blocks[row_block][4], dim_y = _sg_blocks[row_block][5];
        DenseMatrixT Yi(dim_y,dim_x);
        y.push_back(Yi);
    }


    for (FLENS_DEFAULT_INDEXTYPE col_block=0; col_block<(FLENS_DEFAULT_INDEXTYPE)_sg_blocks.size(); ++col_block) {
        FLENS_DEFAULT_INDEXTYPE dim_x = _sg_blocks[col_block][4], dim_y = _sg_blocks[col_block][5];
        DenseMatrixT Xj(dim_y,dim_x);
        FLENS_DEFAULT_INDEXTYPE fi = _sg_blocks[col_block][2], li=fi+dim_y-1;
        for (FLENS_DEFAULT_INDEXTYPE l=1; l<=dim_x; ++l) {
            Xj(_,l) = x(_(fi,li));
            fi = li+1;
            li = fi+dim_y-1;
        }

        for (FLENS_DEFAULT_INDEXTYPE row_block=0; row_block<(FLENS_DEFAULT_INDEXTYPE)_sg_blocks.size(); ++row_block) {
            //std::cout << "Multiplying (" << _sg_blocks[row_block][0] << ", " << _sg_blocks[row_block][1]
            //          << ") * (" << _sg_blocks[col_block][0] << ", " << _sg_blocks[col_block][1] << ")" << std::endl;
            y[row_block] += this->block_multiplication(row_block,col_block,Xj);
        }
    }

    DenseVectorT ret(_dim);
    FLENS_DEFAULT_INDEXTYPE pos=1;
    for (FLENS_DEFAULT_INDEXTYPE block=0; block<(FLENS_DEFAULT_INDEXTYPE)_sg_blocks.size(); ++block) {
        FLENS_DEFAULT_INDEXTYPE dim_x = _sg_blocks[block][4], dim_y = _sg_blocks[block][5];
        for (FLENS_DEFAULT_INDEXTYPE i=1; i<=dim_x; ++i) {
            for (FLENS_DEFAULT_INDEXTYPE j=1; j<=dim_y; ++j) {
                ret(pos) = y[block](j,i);
                ++pos;
            }
        }
    }
    return ret;
}

template <typename T, typename Basis2D, typename S1_x, typename S1_y, typename S2_x, typename S2_y>
flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >
TensorSparseGrid2D<T, Basis2D, S1_x, S1_y, S2_x, S2_y>::block_multiplication(FLENS_DEFAULT_INDEXTYPE row_block,
                                                                             FLENS_DEFAULT_INDEXTYPE col_block,
                                                                             const DenseMatrixT &Xj)
                                                                             const
{
    using flens::_;

    FLENS_DEFAULT_INDEXTYPE dim_Lambda_x =     _sg_blocks[row_block][4], dim_Lambda_y =     _sg_blocks[row_block][5];
    FLENS_DEFAULT_INDEXTYPE dim_Lambda_prime_x=_sg_blocks[col_block][4], dim_Lambda_prime_y=_sg_blocks[col_block][5];


    std::pair<FLENS_DEFAULT_INDEXTYPE,FLENS_DEFAULT_INDEXTYPE> lp_x(_sg_blocks[row_block][0],_sg_blocks[col_block][0]);
    std::pair<FLENS_DEFAULT_INDEXTYPE,FLENS_DEFAULT_INDEXTYPE> lp_y(_sg_blocks[row_block][1],_sg_blocks[col_block][1]);

    LevelPairMap::const_iterator it;
    it = _levelpair_map.find(lp_x);
    FLENS_DEFAULT_INDEXTYPE pos_x  = (*it).second;
    it = _levelpair_map.find(lp_y);
    FLENS_DEFAULT_INDEXTYPE pos_y  = (*it).second;

    FLENS_DEFAULT_INDEXTYPE NumRows_x = _matrixblocks_s1_x[pos_x].numRows();
    FLENS_DEFAULT_INDEXTYPE NumCols_x = _matrixblocks_s1_x[pos_x].numCols();
    FLENS_DEFAULT_INDEXTYPE NumRows_y = _matrixblocks_s1_y[pos_y].numRows();
    FLENS_DEFAULT_INDEXTYPE NumCols_y = _matrixblocks_s1_y[pos_y].numCols();

    assert(NumRows_x == dim_Lambda_x);
    assert(NumRows_y == dim_Lambda_y);
    assert(NumCols_x == dim_Lambda_prime_x);
    assert(NumCols_y == dim_Lambda_prime_y);

    assert(Xj.numCols() == dim_Lambda_prime_x);
    assert(Xj.numRows() == dim_Lambda_prime_y);


    DenseMatrixT Yi(dim_Lambda_y,dim_Lambda_x);

    if (dim_Lambda_y*dim_Lambda_prime_x <= dim_Lambda_prime_y*dim_Lambda_x) {
        //std::cout << "->: " << dim_Lambda_y << " x " << dim_Lambda_prime_x << std::endl;
        DenseMatrixT S_x_otimes_I_y_X(dim_Lambda_y,dim_Lambda_prime_x);
        for (FLENS_DEFAULT_INDEXTYPE j=1; j<=dim_Lambda_prime_x; ++j) {
            DenseVectorT tmp = _matrixblocks_s1_y[pos_y]*Xj(_,j);
            S_x_otimes_I_y_X(_,j) = tmp;
        }

        for (FLENS_DEFAULT_INDEXTYPE i=1; i<=dim_Lambda_y; ++i) {
            DenseVectorT v=S_x_otimes_I_y_X(i,_);
            DenseVectorT tmp = _matrixblocks_s1_x[pos_x]*v;
            Yi(i,_) = tmp;
        }

        for (FLENS_DEFAULT_INDEXTYPE j=1; j<=dim_Lambda_prime_x; ++j) {
            DenseVectorT tmp = _matrixblocks_s2_y[pos_y]*Xj(_,j);
            S_x_otimes_I_y_X(_,j) = tmp;
        }

        for (FLENS_DEFAULT_INDEXTYPE i=1; i<=dim_Lambda_y; ++i) {
            DenseVectorT v=S_x_otimes_I_y_X(i,_);
            DenseVectorT tmp = _matrixblocks_s2_x[pos_x]*v;
            Yi(i,_) += tmp;
        }
    }

    else {
        //std::cout << "->: " << dim_Lambda_prime_y << " x " << dim_Lambda_x << std::endl;
        DenseMatrixT I_x_otimes_S_y_X(dim_Lambda_prime_y,dim_Lambda_x);
        for (FLENS_DEFAULT_INDEXTYPE i=1; i<=dim_Lambda_prime_y; ++i) {
            DenseVectorT v=Xj(i,_);
            DenseVectorT tmp = _matrixblocks_s1_x[pos_x]*v;
            I_x_otimes_S_y_X(i,_) = tmp;
        }
        for (FLENS_DEFAULT_INDEXTYPE j=1; j<=dim_Lambda_x; ++j) {
            DenseVectorT tmp = _matrixblocks_s1_y[pos_y]*I_x_otimes_S_y_X(_,j);
            Yi(_,j) = tmp;
        }

        for (FLENS_DEFAULT_INDEXTYPE i=1; i<=dim_Lambda_prime_y; ++i) {
            DenseVectorT v=Xj(i,_);
            DenseVectorT tmp = _matrixblocks_s2_x[pos_x]*v;
            I_x_otimes_S_y_X(i,_) = tmp;
        }
        for (FLENS_DEFAULT_INDEXTYPE j=1; j<=dim_Lambda_x; ++j) {
            DenseVectorT tmp = _matrixblocks_s2_y[pos_y]*I_x_otimes_S_y_X(_,j);
            Yi(_,j) += tmp;
        }
    }

    return Yi;
}

template <typename T, typename Basis2D, typename S1_x, typename S1_y, typename S2_x, typename S2_y>
void
TensorSparseGrid2D<T, Basis2D, S1_x, S1_y, S2_x, S2_y>::assembleMatrices()
{
    for (LevelPairMap::const_iterator it=_levelpair_map.begin(); it!=_levelpair_map.end(); ++it) {
        FLENS_DEFAULT_INDEXTYPE i1 = (*it).first.first;
        FLENS_DEFAULT_INDEXTYPE i2 = (*it).first.second;
        //FLENS_DEFAULT_INDEXTYPE pos = (*it).second;
        //std::cout << "(" << i1 << ", " << i2 << "): " << pos << std::endl;
        SparseMatrixT A_s1_x = _blockassembler1d.assembleStiffnessMatrixBlock(_s1_x,i1-1,i2-1,0.);
        SparseMatrixT A_s1_y = _blockassembler1d.assembleStiffnessMatrixBlock(_s1_y,i1-1,i2-1,0.);
        SparseMatrixT A_s2_x = _blockassembler1d.assembleStiffnessMatrixBlock(_s2_x,i1-1,i2-1,0.);
        SparseMatrixT A_s2_y = _blockassembler1d.assembleStiffnessMatrixBlock(_s2_y,i1-1,i2-1,0.);
        _matrixblocks_s1_x.push_back(A_s1_x);
        _matrixblocks_s1_y.push_back(A_s1_y);
        _matrixblocks_s2_x.push_back(A_s2_x);
        _matrixblocks_s2_y.push_back(A_s2_y);
    }
}




}   // namespace lawa
