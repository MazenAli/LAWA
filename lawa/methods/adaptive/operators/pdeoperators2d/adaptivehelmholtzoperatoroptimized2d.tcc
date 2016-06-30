namespace lawa {


template <typename T, DomainType Domain>
AdaptiveHelmholtzOperatorOptimized2D<T,Orthogonal,Domain,Multi,Orthogonal,Domain,Multi>
::AdaptiveHelmholtzOperatorOptimized2D(const Basis2D &_basis, T _c)
: basis(_basis), c(_c),
  cA(0.), CA(0.), kappa(0.),
  compression_1d_x(basis.first), compression_1d_y(basis.second), compression(basis),
  laplace_data1d(basis.first),
  P_data()
{
    T cA_x=0., CA_x=0., cA_y=0., CA_y = 0.;
    if (Domain==R) {
        //todo: change to c/2.
        AdaptiveHelmholtzOperatorOptimized1D<T,Orthogonal,Domain,Multi> helmholtzOp1D_x(basis.first,c);
        cA_x = helmholtzOp1D_x.cA;
        CA_x = helmholtzOp1D_x.CA;
        AdaptiveHelmholtzOperatorOptimized1D<T,Orthogonal,Domain,Multi> helmholtzOp1D_y(basis.second,c);
        cA_y = helmholtzOp1D_y.cA;
        CA_y = helmholtzOp1D_y.CA;
    }
    else if (Domain==Interval) {
        AdaptiveHelmholtzOperatorOptimized1D<T,Orthogonal,Domain,Multi> helmholtzOp1D_x(basis.first,0.);
        cA_x = helmholtzOp1D_x.cA;
        CA_x = helmholtzOp1D_x.CA;
        AdaptiveHelmholtzOperatorOptimized1D<T,Orthogonal,Domain,Multi> helmholtzOp1D_y(basis.second,0.);
        cA_y = helmholtzOp1D_y.cA;
        CA_y = helmholtzOp1D_y.CA;
    }
    else {
        assert(0);
    }
/*

    std::cout << "AdaptiveHelmholtzOperatorOptimized2D<T,Orthogonal,Domain1,Multi,Orthogonal,Domain2,Multi>: "
              << "j0_x=" << basis.first.j0 << ", j0_y=" << basis.second.j0 << std::endl;
    if (basis.first.d==2 && basis.second.d==2 && c==1.) {

        if      (basis.first.j0==0)  {    cA_x = 0.14;  CA_x = 2.9;    }
        else if (basis.first.j0==-1) {    cA_x = 0.14;  CA_x = 2.9;    }
        else if (basis.first.j0==-2) {    cA_x = 0.18;  CA_x = 2.9;    }
        else if (basis.first.j0==-3) {    cA_x = 0.18;  CA_x = 2.9;    }
        else if (basis.first.j0==-4) {    cA_x = 0.18;  CA_x = 2.9;    }
        else if (basis.first.j0==-5) {    cA_x = 0.18;  CA_x = 2.9;    }
        else if (basis.first.j0==-6) {    cA_x = 0.18;  CA_x = 2.9;    }
        else if (basis.first.j0==-7) {    cA_x = 0.18;  CA_x = 2.9;    }
        else if (basis.first.j0==-8) {    cA_x = 0.18;  CA_x = 2.9;    }
        else if (basis.first.j0==-9) {    cA_x = 0.18;  CA_x = 2.9;    }
        else if (basis.first.j0==-10) {    cA_x = 0.18;  CA_x = 2.9;    }
        else assert(0);

        if      (basis.second.j0==0)  {    cA_y = 0.14;  CA_y = 2.9;    }
        else if (basis.second.j0==-1) {    cA_y = 0.14;  CA_y = 2.9;    }
        else if (basis.second.j0==-2) {    cA_y = 0.18;  CA_y = 2.9;    }
        else if (basis.second.j0==-3) {    cA_y = 0.18;  CA_y = 2.9;    }
        else if (basis.second.j0==-4) {    cA_y = 0.18;  CA_y = 2.9;    }
        else if (basis.second.j0==-5) {    cA_y = 0.18;  CA_y = 2.9;    }
        else if (basis.second.j0==-6) {    cA_y = 0.18;  CA_y = 2.9;    }
        else if (basis.second.j0==-7) {    cA_x = 0.18;  CA_x = 2.9;    }
        else if (basis.second.j0==-8) {    cA_x = 0.18;  CA_x = 2.9;    }
        else if (basis.second.j0==-9) {    cA_x = 0.18;  CA_x = 2.9;    }
        else if (basis.second.j0==-10) {    cA_x = 0.18;  CA_x = 2.9;    }
        else assert(0);

    }
*/
    cA = std::min(cA_x,cA_y);
    CA = std::max(CA_x,CA_y);
    kappa = CA/cA;
}

template <typename T, DomainType Domain>
T
AdaptiveHelmholtzOperatorOptimized2D<T,Orthogonal,Domain,Multi,Orthogonal,Domain,Multi>
::operator()(const Index2D &row_index, const Index2D &col_index)
{
    T id_x = 0., dd_y = 0.;
    if ((row_index).index1.k==(col_index).index1.k && (row_index).index1.j==(col_index).index1.j &&
        (row_index).index1.xtype==(col_index).index1.xtype) {
        id_x = 1.;
        dd_y = laplace_data1d(row_index.index2,col_index.index2);
    }
    T id_y = 0., dd_x = 0.;
    if ((row_index).index2.k==(col_index).index2.k && (row_index).index2.j==(col_index).index2.j &&
        (row_index).index2.xtype==(col_index).index2.xtype) {
        id_y = 1.;
        dd_x = laplace_data1d(row_index.index1,col_index.index1);
    }

    T val = (dd_x*id_y + id_x*dd_y + c*id_x*id_y);
    if (fabs(val)>0) {
        return this->prec(row_index) * val * this->prec(col_index);
    }
    else {
        return 0.;
    }
}

template <typename T, DomainType Domain>
T
AdaptiveHelmholtzOperatorOptimized2D<T,Orthogonal,Domain,Multi,Orthogonal,Domain,Multi>
::prec(const Index2D &index)
{
    T prec = 1.;
    //const_coeff2d_it it_P_end   = P_data.end();
    //const_coeff2d_it it_index   = P_data.find(index);
    //if (it_index != it_P_end) {
    //    prec *= (*it_index).second;
    //}
    //else {
        T prec_dd_x = laplace_data1d(index.index1,index.index1);
        T prec_dd_y = laplace_data1d(index.index2,index.index2);
        T tmp = 1./std::sqrt(fabs(prec_dd_x + prec_dd_y + c ));
    //    P_data[index] = tmp;
        prec *= tmp;
    //}
    return prec;
}

template <typename T, DomainType Domain>
Coefficients<Lexicographical,T,Index2D>
AdaptiveHelmholtzOperatorOptimized2D<T,Orthogonal,Domain,Multi,Orthogonal,Domain,Multi>
::mv(const IndexSet<Index2D> &LambdaRow, const Coefficients<Lexicographical,T,Index2D> &v)
{
    Coefficients<Lexicographical,T,Index2D> ret;

    IndexSet<Index1D> LambdaRow_x, LambdaRow_y;
    split(LambdaRow, LambdaRow_x, LambdaRow_y);

    std::map<Index1D,IndexSet<Index1D>,lt<Lexicographical,Index1D> > sparsitypatterns_x,
                                                                     sparsitypatterns_y;

    compression_1d_x.setParameters(LambdaRow_x);
    compression_1d_y.setParameters(LambdaRow_y);


    for (const_coeff2d_it col=v.begin(); col!=v.end(); ++col) {
        Index2D col_index = (*col).first;
        T prec_val_col_index = this->prec(col_index) * (*col).second;

        if (LambdaRow.count(col_index)>0) {
            ret[col_index] += c * prec_val_col_index;
        }


        IndexSet<Index1D> LambdaRowSparse_x, LambdaRowSparse_y;
        if (sparsitypatterns_x.count(col_index.index1) == 0) {
            LambdaRowSparse_x = this->compression_1d_x.SparsityPattern(col_index.index1, LambdaRow_x);
            sparsitypatterns_x[col_index.index1] = LambdaRowSparse_x;
        }
        else {
            LambdaRowSparse_x = sparsitypatterns_x[col_index.index1];
        }

//        if (col_index.index1.xtype==XBSpline) {
//            std::cout << "LambdaRow_x for " << col_index.index1 << " : " << LambdaRow_x<< std::endl;
//            std::cout << "LambdaRowSparse for " << col_index.index1 << " : " << LambdaRowSparse_x<< std::endl;
//        }

        if (sparsitypatterns_y.count(col_index.index2) == 0) {
            LambdaRowSparse_y = this->compression_1d_y.SparsityPattern(col_index.index2, LambdaRow_y);
            sparsitypatterns_y[col_index.index2] = LambdaRowSparse_y;
        }
        else {
            LambdaRowSparse_y = sparsitypatterns_y[col_index.index2];
        }

//        LambdaRowSparse_x = lambdaTilde1d_PDE(col_index.index1, basis.first, 10, basis.first.j0, basis.first.j0+10, false);
        for (const_set1d_it row_x=LambdaRowSparse_x.begin(); row_x!=LambdaRowSparse_x.end(); ++row_x) {
            Index2D tmp_row_index(*row_x,col_index.index2);
            if (LambdaRow.count(tmp_row_index)>0) {

                T val_x = laplace_data1d(*row_x,col_index.index1);
                if (fabs(val_x)>0.) {
                    ret[tmp_row_index] += val_x * prec_val_col_index;
                }

            }
        }
//        LambdaRowSparse_y = lambdaTilde1d_PDE(col_index.index2, basis.second, 10, basis.second.j0, basis.second.j0+10, false);
        for (const_set1d_it row_y=LambdaRowSparse_y.begin(); row_y!=LambdaRowSparse_y.end(); ++row_y) {
            Index2D tmp_row_index(col_index.index1, *row_y);
            if (LambdaRow.count(tmp_row_index)>0) {

                T val_y = laplace_data1d(*row_y,col_index.index2);
                if (fabs(val_y)>0.) {
                    ret[tmp_row_index] += val_y * prec_val_col_index;
                }

            }
        }
    }
    for (coeff2d_it it=ret.begin(); it!=ret.end(); ++it) {
        (*it).second *=  this->prec((*it).first);
    }
    return ret;
}

template <typename T, DomainType Domain>
void
AdaptiveHelmholtzOperatorOptimized2D<T,Orthogonal,Domain,Multi,Orthogonal,Domain,Multi>
::toFlensSparseMatrix(const IndexSet<Index2D>& LambdaRow, const IndexSet<Index2D>& LambdaCol,
                      SparseMatrixT &A_flens, FLENS_DEFAULT_INDEXTYPE J)
{
    std::cerr << "  -> toFlensSparseMatrix called with J= " << J << std::endl;
    IndexSet<Index1D> LambdaRow_x, LambdaRow_y;
    split(LambdaRow, LambdaRow_x, LambdaRow_y);
    std::map<Index1D,IndexSet<Index1D>,lt<Lexicographical,Index1D> > sparsitypatterns_x,
                                                                     sparsitypatterns_y;

    compression_1d_x.setParameters(LambdaRow_x);
    compression_1d_y.setParameters(LambdaRow_y);

    const_set2d_it LambdaRow_end=LambdaRow.end();


    std::cerr << "AdaptiveHelmholtzOperatorOptimized2D<T,Orthogonal,Domain1,Multi,Orthogonal,Domain2,Multi>::toFlensSparseMatrix called." << std::endl;

    std::map<Index2D,FLENS_DEFAULT_INDEXTYPE,lt<Lexicographical,Index2D> > row_indices;
    FLENS_DEFAULT_INDEXTYPE row_count = 1;
    for (const_set2d_it row=LambdaRow.begin(); row!=LambdaRow.end(); ++row, ++row_count) {
        row_indices[(*row)] = row_count;
    }

    FLENS_DEFAULT_INDEXTYPE col_count = 1;
    for (const_set2d_it col=LambdaCol.begin(); col!=LambdaCol.end(); ++col, ++col_count) {
        Index2D col_index = *col;

        T prec_col_index = this->prec(col_index);
        if (LambdaRow.count(col_index)>0) {
            A_flens(col_count,col_count) += c * prec_col_index*prec_col_index;
        }

        IndexSet<Index1D> LambdaRowSparse_x, LambdaRowSparse_y;
        if (sparsitypatterns_x.count((*col).index1) == 0) {
            LambdaRowSparse_x = this->compression_1d_x.SparsityPattern((*col).index1, LambdaRow_x, J);
            sparsitypatterns_x[(*col).index1] = LambdaRowSparse_x;
        }
        else {
            LambdaRowSparse_x = sparsitypatterns_x[(*col).index1];
        }


        if (sparsitypatterns_y.count((*col).index2) == 0) {
            LambdaRowSparse_y = this->compression_1d_y.SparsityPattern((*col).index2, LambdaRow_y, J);
            sparsitypatterns_y[(*col).index2] = LambdaRowSparse_y;
        }
        else {
            LambdaRowSparse_y = sparsitypatterns_y[(*col).index2];
        }

        for (const_set1d_it row_x=LambdaRowSparse_x.begin(); row_x!=LambdaRowSparse_x.end(); ++row_x) {
            Index2D tmp_row_index(*row_x,(*col).index2);
            if (LambdaRow.count(tmp_row_index)>0) {
                T val_x = laplace_data1d(*row_x,(*col).index1);
                if (fabs(val_x)>0.) {
                    A_flens(row_indices[tmp_row_index],col_count) += this->prec(tmp_row_index)*val_x*prec_col_index;
                }
            }
        }
        for (const_set1d_it row_y=LambdaRowSparse_y.begin(); row_y!=LambdaRowSparse_y.end(); ++row_y) {
            Index2D tmp_row_index((*col).index1, *row_y);
            if (LambdaRow.count(tmp_row_index)>0) {
                 T val_y = laplace_data1d(*row_y,(*col).index2);
                if (fabs(val_y)>0.) {
                    A_flens(row_indices[tmp_row_index],col_count) += this->prec(tmp_row_index)*val_y*prec_col_index;
                }
            }
        }
    }
    A_flens.finalize();
}

template <typename T, DomainType Domain>
void
AdaptiveHelmholtzOperatorOptimized2D<T,Orthogonal,Domain,Multi,Orthogonal,Domain,Multi>
::toFlensSparseMatrix(const IndexSet<Index2D>& LambdaRow, const IndexSet<Index2D>& LambdaCol,
                      SparseMatrixT &A_flens, T eps)
{
    FLENS_DEFAULT_INDEXTYPE d = basis.first.d;
    if (basis.first.d!=basis.second.d) {
        std::cerr << "AdaptiveHelmholtzOperatorOptimized2D::toFlensSparseMatrix not implemented "
                  << "for different polynomial orders " << basis.first.d << ", " << basis.second.d
                  << "." << std::endl;
        exit(1);
    }
    FLENS_DEFAULT_INDEXTYPE J=0;        //compression
    J = (FLENS_DEFAULT_INDEXTYPE)std::ceil(-1./(d-1.5)*log(eps/CA)/log(2.));
    std::cerr << "   -> toFlensSparseMatrix: Estimated compression level for "
              << "tol = " << eps << " : " << J << std::endl;
    this->toFlensSparseMatrix(LambdaRow,LambdaCol,A_flens,J);
}

template <typename T, DomainType Domain>
Coefficients<Lexicographical,T,Index2D>
AdaptiveHelmholtzOperatorOptimized2D<T,Orthogonal,Domain,Multi,Orthogonal,Domain,Multi>
::apply(const Coefficients<Lexicographical,T,Index2D> &v, FLENS_DEFAULT_INDEXTYPE k, FLENS_DEFAULT_INDEXTYPE J, cxxblas::Transpose /*trans*/)
{
    Coefficients<Lexicographical,T,Index2D> ret;
    if (v.size() == 0) return ret;

    Coefficients<AbsoluteValue,T,Index2D> temp;
    temp = v;

    FLENS_DEFAULT_INDEXTYPE s = 0, count = 0;
    for (const_abs_coeff2d_it it = temp.begin(); (it != temp.end()) && (s<=k); ++it) {
        ret[(*it).second] += c * (*it).first * this->prec((*it).second);
        Index1D col_index_x = (*it).second.index1;
        Index1D col_index_y = (*it).second.index2;

        T prec_col_index = this->prec(Index2D(col_index_x, col_index_y));

        IndexSet<Index1D> Lambda_x, Lambda_y;
        FLENS_DEFAULT_INDEXTYPE maxlevel_x, maxlevel_y;

        J==-1000 ? maxlevel_x=col_index_x.j+(k-s)+1 : maxlevel_x=J;
        maxlevel_x =36;
        maxlevel_y =36;
        Lambda_x=lambdaTilde1d_PDE(col_index_x, basis.first, (k-s), basis.first.j0,
                                   maxlevel_x,false);

        J==-1000 ? maxlevel_y=col_index_y.j+(k-s)+1 : maxlevel_y=J;
        Lambda_y=lambdaTilde1d_PDE(col_index_y, basis.second, (k-s), basis.second.j0,
                                   maxlevel_y,false);

        for (const_set1d_it row_x = Lambda_x.begin(); row_x != Lambda_x.end(); ++row_x) {
            Index2D row_index(*row_x, col_index_y);
            ret[row_index] += this->laplace_data1d(*row_x, col_index_x) * prec_col_index * (*it).first;
        }

        for (const_set1d_it row_y = Lambda_y.begin(); row_y != Lambda_y.end(); ++row_y) {
            Index2D row_index(col_index_x, *row_y);
            ret[row_index] += this->laplace_data1d(*row_y, col_index_y) * prec_col_index * (*it).first;
        }
        ++count;
        s = (FLENS_DEFAULT_INDEXTYPE)(log(T(count))/log(T(2))) + 1;
    }
    for (const_coeff2d_it it=ret.begin(); it!=ret.end(); ++it) {
        ret[(*it).first] *= this->prec((*it).first);
    }

    return ret;
}

template <typename T, DomainType Domain>
void
AdaptiveHelmholtzOperatorOptimized2D<T,Orthogonal,Domain,Multi,Orthogonal,Domain,Multi>
::apply(const Coefficients<Lexicographical,T,Index2D> &v, T eps,
        Coefficients<Lexicographical,T,Index2D> &ret, cxxblas::Transpose /*trans*/)
{
//    Coefficients<Lexicographical,T,Index2D> ret;
/*
    Coefficients<AbsoluteValue,T,Index2D> v_abs;
    v_abs = v;
    FLENS_DEFAULT_INDEXTYPE k = this->findK(v_abs, eps);
    ret = this->apply(v, k);
    return ret;
*/
    Coefficients<Bucket,T,Index2D> v_bucket;
    T tol = 0.5*eps/CA;
    v_bucket.bucketsort(v,tol);
    long double squared_v_norm = (long double)std::pow(v.norm(2.),(T)2.);
    long double squared_v_bucket_norm = 0.;
    T delta=0.;
    FLENS_DEFAULT_INDEXTYPE l=0;
    FLENS_DEFAULT_INDEXTYPE support_size_all_buckets=0;
    for (FLENS_DEFAULT_INDEXTYPE i=0; i<(FLENS_DEFAULT_INDEXTYPE)v_bucket.buckets.size(); ++i) {
        squared_v_bucket_norm += v_bucket.bucket_ell2norms[i]*v_bucket.bucket_ell2norms[i];
        T squared_delta = fabs(squared_v_norm - squared_v_bucket_norm);
        support_size_all_buckets += v_bucket.buckets[i].size();
        delta = std::sqrt(squared_delta);
        l = i+1;
        if (squared_delta<tol*tol) {
            break;
        }
    }

    if (delta>tol) delta = eps/2.;
    //std::cerr << "APPLY: squared_v_norm=" << squared_v_norm << ", squared_v_bucket_norm=" << squared_v_bucket_norm << std::endl;

    for (FLENS_DEFAULT_INDEXTYPE i=0; i<l; ++i) {
        Coefficients<Lexicographical,T,Index2D> w_p;
        v_bucket.addBucketToCoefficients(w_p,i);
        if (w_p.size()==0) continue;
        T numerator = w_p.norm(2.) * support_size_all_buckets;
        T denominator = w_p.size() * (eps-delta) / CA;
        //std::cout << "Bucket " << i << ": size=" << w_p.size() << ", (eps-delta) " << fabs(eps-delta) << std::endl;
        //FLENS_DEFAULT_INDEXTYPE jp = (FLENS_DEFAULT_INDEXTYPE)std::max(std::log(numerator/denominator) / std::log(2.) / (basis.d-1.5), 0.);
        FLENS_DEFAULT_INDEXTYPE jp = (FLENS_DEFAULT_INDEXTYPE)std::max(std::log(numerator/denominator) / std::log(2.), (T)0.);
        //std::cout << "Bucket " << i << ": #wp= " << w_p.size() << ", jp=" << jp << std::endl;
        for (const_coeff2d_it it=w_p.begin(); it!=w_p.end(); ++it) {
            Index1D col_index_x = (*it).first.index1;
            Index1D col_index_y = (*it).first.index2;
            T prec_col_index = this->prec(Index2D(col_index_x, col_index_y));

            ret[(*it).first] += c * prec_col_index * (*it).second;

            IndexSet<Index1D> Lambda_x, Lambda_y;

            FLENS_DEFAULT_INDEXTYPE maxlevel_x = std::min(col_index_x.j+jp,(FLENS_DEFAULT_INDEXTYPE) 60);
            FLENS_DEFAULT_INDEXTYPE maxlevel_y = std::min(col_index_y.j+jp,(FLENS_DEFAULT_INDEXTYPE) 60);
            Lambda_x=lambdaTilde1d_PDE(col_index_x, basis.first, jp, basis.first.j0,
                                       maxlevel_x,false);
            Lambda_y=lambdaTilde1d_PDE(col_index_y, basis.second,jp, basis.second.j0,
                                       maxlevel_y,false);
            for (const_set1d_it row_x = Lambda_x.begin(); row_x != Lambda_x.end(); ++row_x) {
                Index2D row_index(*row_x, col_index_y);
                ret[row_index] += this->laplace_data1d(*row_x, col_index_x) * prec_col_index * (*it).second;
            }
            for (const_set1d_it row_y = Lambda_y.begin(); row_y != Lambda_y.end(); ++row_y) {
                Index2D row_index(col_index_x, *row_y);
                ret[row_index] += this->laplace_data1d(*row_y, col_index_y) * prec_col_index * (*it).second;
            }
        }
    }
    for (coeff2d_it it=ret.begin(); it!=ret.end(); ++it) {
        (*it).second *= this->prec((*it).first);
    }

    //return ret;
}

template <typename T, DomainType Domain>
void
AdaptiveHelmholtzOperatorOptimized2D<T,Orthogonal,Domain,Multi,Orthogonal,Domain,Multi>
::apply(const Coefficients<Lexicographical,T,Index2D> &v, T eps,
        const IndexSet<Index2D> &Lambda, Coefficients<Lexicographical,T,Index2D> &ret, cxxblas::Transpose /*trans*/)
{

/*
    Coefficients<AbsoluteValue,T,Index2D> v_abs;
    v_abs = v;
    FLENS_DEFAULT_INDEXTYPE k = this->findK(v_abs, eps);
    ret = this->apply(v, k);
    return ret;
*/

    Coefficients<Bucket,T,Index2D> v_bucket;
    T tol = 0.5*eps/CA;
    v_bucket.bucketsort(v,tol);
    //std::cerr << "APPLY: NumOfBuckets=" << v_bucket.buckets.size() << std::endl;
    long double squared_v_norm = (long double)std::pow(v.norm(2.),(T)2.);
    long double squared_v_bucket_norm = 0.;
    T delta=0.;
    FLENS_DEFAULT_INDEXTYPE l=0;
    FLENS_DEFAULT_INDEXTYPE support_size_all_buckets=0;
    for (FLENS_DEFAULT_INDEXTYPE i=0; i<(FLENS_DEFAULT_INDEXTYPE)v_bucket.buckets.size(); ++i) {
        squared_v_bucket_norm += v_bucket.bucket_ell2norms[i]*v_bucket.bucket_ell2norms[i];
        T squared_delta = fabs(squared_v_norm - squared_v_bucket_norm);
        support_size_all_buckets += v_bucket.buckets[i].size();
        delta = std::sqrt(squared_delta);
        l = i+1;
        if (squared_delta<tol*tol) {
            break;
        }
    }
    //std::cerr << "APPLY: squared_v_norm=" << squared_v_norm << ", squared_v_bucket_norm=" << squared_v_bucket_norm << std::endl;

    for (FLENS_DEFAULT_INDEXTYPE i=0; i<l; ++i) {
        Coefficients<Lexicographical,T,Index2D> w_p;
        v_bucket.addBucketToCoefficients(w_p,i);
        if (w_p.size()==0) continue;
        T numerator = w_p.norm(2.) * support_size_all_buckets;
        T denominator = w_p.size() * (eps-delta) / CA;
        //std::cout << "Bucket " << i << ": size=" << w_p.size() << ", (eps-delta) " << fabs(eps-delta) << std::endl;
        //FLENS_DEFAULT_INDEXTYPE jp = (FLENS_DEFAULT_INDEXTYPE)std::max(std::log(numerator/denominator) / std::log(2.) / (basis.d-1.5), 0.);
        FLENS_DEFAULT_INDEXTYPE jp = (FLENS_DEFAULT_INDEXTYPE)std::max(std::log(numerator/denominator) / std::log(2.), (T)0.)+2;
        //std::cout << "Bucket " << i << ": #wp= " << w_p.size() << ", jp=" << jp << std::endl;
        for (const_coeff2d_it it=w_p.begin(); it!=w_p.end(); ++it) {
            Index1D col_index_x = (*it).first.index1;
            Index1D col_index_y = (*it).first.index2;
            T prec_col_index = this->prec(Index2D(col_index_x, col_index_y));

            ret[(*it).first] += c * prec_col_index * (*it).second;

            IndexSet<Index1D> Lambda_x, Lambda_y;

            FLENS_DEFAULT_INDEXTYPE maxlevel_x = std::min(col_index_x.j+jp,(FLENS_DEFAULT_INDEXTYPE) 60);
            FLENS_DEFAULT_INDEXTYPE maxlevel_y = std::min(col_index_y.j+jp,(FLENS_DEFAULT_INDEXTYPE) 60);
            Lambda_x=lambdaTilde1d_PDE(col_index_x, basis.first, jp, basis.first.j0,
                                       maxlevel_x,false);

            Lambda_y=lambdaTilde1d_PDE(col_index_y, basis.second,jp, basis.second.j0,
                                       maxlevel_y,false);

            for (const_set1d_it row_x = Lambda_x.begin(); row_x != Lambda_x.end(); ++row_x) {
                Index2D row_index(*row_x, col_index_y);
                if (Lambda.count(row_index)>0) {
                    ret[row_index] += this->laplace_data1d(*row_x, col_index_x) * prec_col_index * (*it).second;
                }
            }

            for (const_set1d_it row_y = Lambda_y.begin(); row_y != Lambda_y.end(); ++row_y) {
                Index2D row_index(col_index_x, *row_y);
                if (Lambda.count(row_index)>0) {
                    ret[row_index] += this->laplace_data1d(*row_y, col_index_y) * prec_col_index * (*it).second;
                }

            }
        }
    }
    for (coeff2d_it it=ret.begin(); it!=ret.end(); ++it) {
        (*it).second *= this->prec((*it).first);
    }

    return;// ret;
}

template <typename T, DomainType Domain>
FLENS_DEFAULT_INDEXTYPE
AdaptiveHelmholtzOperatorOptimized2D<T,Orthogonal,Domain,Multi,Orthogonal,Domain,Multi>
::findK(const Coefficients<AbsoluteValue,T,Index2D> &v, T eps)
{
    FLENS_DEFAULT_INDEXTYPE d=basis.first.d;
    if (v.size() == 0) return 1;
    T s=d-1.5;    //s = gamma-1, gamma the smoothness index of the wavelet basis

    T tau = 1.0 / (s + 0.5);
    // here the constant in (7.27) (KU-Wavelet) is estimated with 10
    FLENS_DEFAULT_INDEXTYPE k_eps = static_cast<FLENS_DEFAULT_INDEXTYPE>(10*log(std::pow(eps, -1.0/s)*std::pow(v.wtauNorm(tau), 1.0/s)) / log(2.0));
    flens::DenseVector<flens::Array<T> > normsec = v.norm_sections();
    T ErrorEstimateFactor = 1.;
    //std::cout << "eps = " << eps << ", k_eps = " << k_eps << std::endl;

    for (FLENS_DEFAULT_INDEXTYPE k=1; k<=k_eps; ++k) {
        //std::cout << "At k = " << setw(3) << k;

        T R_k = 0.0;
        for (FLENS_DEFAULT_INDEXTYPE i=k; i<=normsec.lastIndex()-1; ++i) {
            R_k += normsec(i+1);
        }
        R_k *= 2.8;//parameters.CA;
        //std::cout << ", R_k = " << setw(10) << R_k;
        R_k += std::pow(2.,-k*s) * normsec(1);
        //std::cout << ", R_k = " << setw(10) << R_k;

        for (FLENS_DEFAULT_INDEXTYPE l=0; l<=k-1; ++l) {
            if (k-l<=normsec.lastIndex()-1) {
                //R_k += std::pow(l,-1.01)*std::pow(2.,-l*s) * normsec(k-l+1);
                R_k += std::pow(2.,-l*s) * normsec(k-l+1);
            }
        }
        //std::cout << ", R_k = " << setw(10) << R_k;
        R_k *= ErrorEstimateFactor;
        //std::cout << ", R_k = " << setw(10) << R_k << ", eps = " << setw(10) << eps << endl;

        if (R_k<=eps) {
            std::cout << "   findK ==> k = " << k << " for eps =  " << eps << std::endl;
            FLENS_DEFAULT_INDEXTYPE maxlevel=100;
            if (d==2)         {    maxlevel=100; }
            else if (d==3)    {   maxlevel=100; }    //for non-singular examples, also lower values are possible
            return std::min(std::max(k,(FLENS_DEFAULT_INDEXTYPE) 1),maxlevel);
        }
    }
    return std::min(k_eps,(FLENS_DEFAULT_INDEXTYPE) 60);    //higher level differences result in translation indices that cannot be stored in FLENS_DEFAULT_INDEXTYPE.
}

template <typename T, DomainType Domain>
void
AdaptiveHelmholtzOperatorOptimized2D<T,Orthogonal,Domain,Multi,Orthogonal,Domain,Multi>
::clear()
{

}



template <typename T, DomainType Domain1, DomainType Domain2>
AdaptiveHelmholtzOperatorOptimized2D<T,Primal,Domain1,SparseMulti,Primal,Domain2,SparseMulti>
::AdaptiveHelmholtzOperatorOptimized2D(const Basis2D &_basis, T _c, T _diffusion_y)
: basis(_basis), c(_c), diffusion_y(_diffusion_y), thresh(0.),
  cA(0.), CA(0.), kappa(0.),
  compression_1d_x(basis.first), compression_1d_y(basis.second), compression(basis),
  laplace_data1d_x(basis.first), identity_data1d_x(basis.first),
  laplace_data1d_y(basis.second), identity_data1d_y(basis.second),
  P_data()
{
    T cA_x=0., CA_x=0., cA_y=0., CA_y = 0.;
    std::cout << "AdaptiveHelmholtzOperatorOptimized2D<T,Primal,Domain1,SparseMulti,"
              << "Primal,Domain2,SparseMulti>: j0_x=" << basis.first.j0 << ", j0_y="
              << basis.second.j0 << std::endl;
    if (basis.first.d==4 && basis.second.d==4 && c==1.) {

        if      (basis.first.j0==0)  {    cA_x = 0.14;  CA_x = 2.9;    }
        else if (basis.first.j0== 1) {    cA_x = 0.14;  CA_x = 2.9;    }
        else if (basis.first.j0==-1) {    cA_x = 0.14;  CA_x = 2.9;    }
        else if (basis.first.j0==-2) {    cA_x = 0.19;  CA_x = 2.9;    }
        else if (basis.first.j0==-3) {    cA_x = 0.19;  CA_x = 2.9;    }
        else if (basis.first.j0==-4) {    cA_x = 0.19;  CA_x = 2.9;    }
        else assert(0);

        if      (basis.second.j0==0)  {    cA_y = 0.14;  CA_y = 2.9;    }
        else if (basis.second.j0== 1) {    cA_y = 0.14;  CA_y = 2.9;    }
        else if (basis.second.j0==-1) {    cA_y = 0.14;  CA_y = 2.9;    }
        else if (basis.second.j0==-2) {    cA_y = 0.19;  CA_y = 2.9;    }
        else if (basis.second.j0==-3) {    cA_y = 0.19;  CA_y = 2.9;    }
        else if (basis.second.j0==-4) {    cA_y = 0.19;  CA_y = 2.9;    }
        else assert(0);

    }
    cA = std::min(cA_x*cA_x,cA_y*cA_y);
    CA = std::max(CA_x*CA_x,CA_y*CA_y);
    kappa = CA/cA;
}

template <typename T, DomainType Domain1, DomainType Domain2>
T
AdaptiveHelmholtzOperatorOptimized2D<T,Primal,Domain1,SparseMulti,Primal,Domain2,SparseMulti>
::operator()(const Index2D &row_index, const Index2D &col_index)
{
    T id_x = 0., dd_y = 0.;
    if (fabs((row_index).index1.j-(col_index).index1.j)<=1) {
        id_x = identity_data1d_x(row_index.index1,col_index.index1);
        if (fabs(id_x)>1e-14) {
            dd_y = laplace_data1d_y(row_index.index2,col_index.index2);
        }
    }
    T id_y = 0., dd_x = 0.;
    if (fabs((row_index).index2.j-(col_index).index2.j)<=1) {
        id_y = identity_data1d_y(row_index.index2,col_index.index2);
        if (fabs(id_y)>1e-14) {
            dd_x = laplace_data1d_x(row_index.index1,col_index.index1);
        }
    }

    T val = (dd_x*id_y + diffusion_y*id_x*dd_y + c*id_x*id_y);
    val *= this->prec(row_index) * this->prec(col_index);
    if (fabs(val)>thresh) {
        return val;
    }
    else {
        return 0.;
    }
}

template <typename T, DomainType Domain1, DomainType Domain2>
T
AdaptiveHelmholtzOperatorOptimized2D<T,Primal,Domain1,SparseMulti,Primal,Domain2,SparseMulti>
::prec(const Index2D &index)
{
    T prec = 1.;
    const_coeff2d_it it_P_end   = P_data.end();
    const_coeff2d_it it_index   = P_data.find(index);
    if (it_index != it_P_end) {
        prec *= (*it_index).second;
    }
    else {
        T prec_id_x = identity_data1d_x(index.index1,index.index1);
        T prec_id_y = identity_data1d_y(index.index2,index.index2);
        T prec_dd_x = laplace_data1d_x(index.index1,index.index1);
        T prec_dd_y = laplace_data1d_y(index.index2,index.index2);
        T tmp = 1./std::sqrt(fabs(prec_dd_x*prec_id_y + diffusion_y*prec_id_x*prec_dd_y + c*prec_id_x*prec_id_y ));
        P_data[index] = tmp;
        prec *= tmp;
    }
    return prec;
}

template <typename T, DomainType Domain1, DomainType Domain2>
Coefficients<Lexicographical,T,Index2D>
AdaptiveHelmholtzOperatorOptimized2D<T,Primal,Domain1,SparseMulti,Primal,Domain2,SparseMulti>
::mv(const IndexSet<Index2D> &LambdaRow, const Coefficients<Lexicographical,T,Index2D> &v)
{
    Coefficients<Lexicographical,T,Index2D> ret;
    this->apply(v, 0., LambdaRow, ret);
    return ret;
}

template <typename T, DomainType Domain1, DomainType Domain2>
void
AdaptiveHelmholtzOperatorOptimized2D<T,Primal,Domain1,SparseMulti,Primal,Domain2,SparseMulti>
::toFlensSparseMatrix(const IndexSet<Index2D>& LambdaRow, const IndexSet<Index2D>& LambdaCol,
                      SparseMatrixT &A_flens, FLENS_DEFAULT_INDEXTYPE J)
{
    std::cerr << "  -> toFlensSparseMatrix called with J= " << J << std::endl;

    FLENS_DEFAULT_INDEXTYPE maxSizeSparsityPattern=0;

    std::map<Index1D,IndexSet<Index1D>,lt<Lexicographical,Index1D> > sparsitypatterns_x,
                                                                             sparsitypatterns_y;
    std::map<Index1D,Coefficients<Lexicographical,T,Index1D>,lt<Lexicographical,Index1D> >
        sparsitypatterns_identity_x, sparsitypatterns_laplace_x,
        sparsitypatterns_identity_y, sparsitypatterns_laplace_y;

    IndexSet<Index1D> LambdaRow_x, LambdaRow_y;
    split(LambdaRow, LambdaRow_x, LambdaRow_y);

    compression_1d_x.setParameters(LambdaRow_x);
    compression_1d_y.setParameters(LambdaRow_y);


    std::map<Index2D,FLENS_DEFAULT_INDEXTYPE,lt<Lexicographical,Index2D> > row_indices;



    FLENS_DEFAULT_INDEXTYPE row_count = 1;
    for (const_set2d_it row=LambdaRow.begin(); row!=LambdaRow.end(); ++row, ++row_count) {
        row_indices[(*row)] = row_count;
    }
    std::map<Index2D,FLENS_DEFAULT_INDEXTYPE,lt<Lexicographical,Index2D> >::const_iterator row_indices_end;
    row_indices_end = row_indices.end();
    const_set2d_it LambdaRow_end=LambdaRow.end();

    FLENS_DEFAULT_INDEXTYPE col_count = 1;
    for (const_set2d_it col=LambdaCol.begin(); col!=LambdaCol.end(); ++col, ++col_count) {
        Index2D col_index = *col;
        T prec_col_index = this->prec(col_index);

        IndexSet<Index1D> LambdaRowSparse_x, LambdaRowSparse_y;
        Coefficients<Lexicographical,T,Index1D> LambdaRowSparseIdentity_x, LambdaRowSparseLaplace_x,
                                                LambdaRowSparseIdentity_y, LambdaRowSparseLaplace_y;
        if (sparsitypatterns_x.count((*col).index1) == 0) {
            LambdaRowSparse_x = this->compression_1d_x.SparsityPattern((*col).index1, LambdaRow_x, 1);
            sparsitypatterns_x[(*col).index1] = LambdaRowSparse_x;
            for (const_set1d_it row_x=LambdaRowSparse_x.begin(); row_x!=LambdaRowSparse_x.end(); ++row_x) {
                LambdaRowSparseIdentity_x[*row_x] = identity_data1d_x(*row_x,col_index.index1);
                LambdaRowSparseLaplace_x[*row_x]  =  laplace_data1d_x(*row_x,col_index.index1);
            }
            sparsitypatterns_identity_x[(*col).index1] = LambdaRowSparseIdentity_x;
            sparsitypatterns_laplace_x[(*col).index1]  = LambdaRowSparseLaplace_x;
        }
        else {
            LambdaRowSparse_x         = sparsitypatterns_x[(*col).index1];
            LambdaRowSparseIdentity_x = sparsitypatterns_identity_x[(*col).index1];
            LambdaRowSparseLaplace_x  = sparsitypatterns_laplace_x[(*col).index1];
        }

        if (sparsitypatterns_y.count((*col).index2) == 0) {
            LambdaRowSparse_y = this->compression_1d_y.SparsityPattern((*col).index2, LambdaRow_y, 1);
            sparsitypatterns_y[(*col).index2] = LambdaRowSparse_y;
            for (const_set1d_it row_y=LambdaRowSparse_y.begin(); row_y!=LambdaRowSparse_y.end(); ++row_y) {
                LambdaRowSparseIdentity_y[*row_y] = identity_data1d_y(*row_y,col_index.index2);
                LambdaRowSparseLaplace_y[*row_y]  =  laplace_data1d_y(*row_y,col_index.index2);
            }
            sparsitypatterns_identity_y[(*col).index2] = LambdaRowSparseIdentity_y;
            sparsitypatterns_laplace_y[(*col).index2]  = LambdaRowSparseLaplace_y;
        }
        else {
            LambdaRowSparse_y = sparsitypatterns_y[(*col).index2];
            LambdaRowSparseIdentity_y = sparsitypatterns_identity_y[(*col).index2];
            LambdaRowSparseLaplace_y  = sparsitypatterns_laplace_y[(*col).index2];
        }

        maxSizeSparsityPattern=std::max(maxSizeSparsityPattern,(FLENS_DEFAULT_INDEXTYPE)(LambdaRowSparse_x.size()*LambdaRowSparse_y.size()));

        for (const_set1d_it row_x=LambdaRowSparse_x.begin(); row_x!=LambdaRowSparse_x.end(); ++row_x) {
            T id_x = 0., dd_x = 0.;
            id_x = LambdaRowSparseIdentity_x[*row_x];
            dd_x = LambdaRowSparseLaplace_x[*row_x];
            if (fabs(id_x)<1e-13 && fabs(dd_x)<1e-13) continue;
            for (const_set1d_it row_y=LambdaRowSparse_y.begin(); row_y!=LambdaRowSparse_y.end(); ++row_y) {
                Index2D row_index(*row_x,*row_y);
                std::map<Index2D,FLENS_DEFAULT_INDEXTYPE,lt<Lexicographical,Index2D> >::const_iterator row_count_ptr;
                row_count_ptr = row_indices.find(row_index);
                if (row_count_ptr!=row_indices_end) {
                    T id_y = 0., dd_y = 0.;
                    id_y = LambdaRowSparseIdentity_y[*row_y];
                    dd_y = LambdaRowSparseLaplace_y[*row_y];

                    T val = (dd_x*id_y + diffusion_y*id_x*dd_y + c*id_x*id_y);

                    T prec_row_index = this->prec(row_index);
                    T prec_val = prec_row_index* val * prec_col_index;
                    if (fabs(prec_val)>thresh) {
                        A_flens((*row_count_ptr).second,col_count) = prec_val;
                    }
                }
            }
        }
    }
    std::cout << "Size of LambdaRow: " << LambdaRow.size()
              << ", size of LambdaRowSparse: " << maxSizeSparsityPattern << std::endl;

    A_flens.finalize();
    std::cerr << "   Number of non-zeros: " << A_flens.numNonZeros() << std::endl;
}

template <typename T, DomainType Domain1, DomainType Domain2>
void
AdaptiveHelmholtzOperatorOptimized2D<T,Primal,Domain1,SparseMulti,Primal,Domain2,SparseMulti>
::toFlensSparseMatrix(const IndexSet<Index2D>& LambdaRow, const IndexSet<Index2D>& LambdaCol,
                      SparseMatrixT &A_flens, T /*eps*/)
{
    this->toFlensSparseMatrix(LambdaRow,LambdaCol,A_flens,1);
}

template <typename T, DomainType Domain1, DomainType Domain2>
Coefficients<Lexicographical,T,Index2D>
AdaptiveHelmholtzOperatorOptimized2D<T,Primal,Domain1,SparseMulti,Primal,Domain2,SparseMulti>
::apply(const Coefficients<Lexicographical,T,Index2D> &v, FLENS_DEFAULT_INDEXTYPE /*k*/, FLENS_DEFAULT_INDEXTYPE /*J*/,
        cxxblas::Transpose /*trans*/)
{
    Coefficients<Lexicographical,T,Index2D> ret;
    this->apply(v,0.,ret);
    return ret;
}

template <typename T, DomainType Domain1, DomainType Domain2>
void
AdaptiveHelmholtzOperatorOptimized2D<T,Primal,Domain1,SparseMulti,Primal,Domain2,SparseMulti>
::apply(const Coefficients<Lexicographical,T,Index2D> &v, T /*eps*/,
        Coefficients<Lexicographical,T,Index2D> &ret, cxxblas::Transpose /*trans*/)
{
    if (v.size()==0) return;

    Index1D_Coefficients1D_Hash y_v;
    Coefficients<Lexicographical,T,Index2D> I_S_v(SIZEHASHINDEX2D);
    Index1D_Coefficients1D_Hash x_I_S_v;

    Timer time;
    time.start();

    for (const_coeff2d_it col=v.begin(); col!=v.end(); ++col) {
        Index1D_Coefficients1D_Hash_it p_y_v=y_v.find((*col).first.index2);
        if (p_y_v!=y_v.end()) {
            (*p_y_v).second.operator[]((*col).first.index1) = this->prec((*col).first) * (*col).second;
        }
        else {
            Coefficients<Lexicographical,T,Index1D> coeff_x;
            coeff_x[(*col).first.index1] = this->prec((*col).first) * (*col).second;
            y_v[(*col).first.index2] = coeff_x;
        }
    }

    for (const_Index1D_Coefficients1D_Hash_it it=y_v.begin(); it!=y_v.end(); ++it) {
        Index1D col_index_y = (*it).first;

        IndexSet<Index1D> LambdaRowSparse_y;
        LambdaRowSparse_y = lambdaTilde1d_PDE(col_index_y, basis.second, 1, std::max(col_index_y.j-1,basis.second.j0),
                                                          col_index_y.j+1, false);
        for (const_set1d_it row_y=LambdaRowSparse_y.begin(); row_y!=LambdaRowSparse_y.end(); ++row_y) {
            T dd_y = laplace_data1d_y(*row_y,col_index_y);
            T id_y = identity_data1d_y(*row_y,col_index_y);
            T val2_y = (diffusion_y*dd_y + 0.5*c*id_y);

            for (const_coeff1d_it coeff_x_it=(*it).second.begin(); coeff_x_it!=(*it).second.end(); ++coeff_x_it) {
                Index2D row_index((*coeff_x_it).first,*row_y);
                T val2 = val2_y * (*coeff_x_it).second;
                if (fabs(val2)>0) I_S_v[row_index] += val2;
            }
        }
    }
    for (const_coeff2d_it col=I_S_v.begin(); col!=I_S_v.end(); ++col) {
        Index1D_Coefficients1D_Hash_it p_x_I_S_v=x_I_S_v.find((*col).first.index1);
        if (p_x_I_S_v!=x_I_S_v.end()) {
            (*p_x_I_S_v).second.operator[]((*col).first.index2) = (*col).second;
        }
        else {
            Coefficients<Lexicographical,T,Index1D> coeff_y;
            coeff_y[(*col).first.index2] = (*col).second;
            x_I_S_v[(*col).first.index1] = coeff_y;
        }
    }
    for (const_Index1D_Coefficients1D_Hash_it it=x_I_S_v.begin(); it!=x_I_S_v.end(); ++it) {
        Index1D col_index_x = (*it).first;
        IndexSet<Index1D> LambdaRowSparse_x;
        LambdaRowSparse_x = lambdaTilde1d_PDE(col_index_x, basis.first, 1, std::max(col_index_x.j-1,basis.first.j0),
                                              col_index_x.j+1, false);
        for (const_set1d_it row_x=LambdaRowSparse_x.begin(); row_x!=LambdaRowSparse_x.end(); ++row_x) {
            T id_x = identity_data1d_x(*row_x,col_index_x);
            for (const_coeff1d_it coeff_y_it=(*it).second.begin(); coeff_y_it!=(*it).second.end(); ++coeff_y_it) {
                Index2D row_index(*row_x,(*coeff_y_it).first);
                T val = id_x * (*coeff_y_it).second;
                if (fabs(val)>0.) ret[row_index] += val;
            }
        }
    }

    for (coeff2d_it it=I_S_v.begin(); it!=I_S_v.end(); ++it) {
                (*it).second = 0.;
    }
    for (Index1D_Coefficients1D_Hash_it it=x_I_S_v.begin(); it!=x_I_S_v.end(); ++it) {
        for (coeff1d_it coeff_it=(*it).second.begin(); coeff_it!=(*it).second.end(); ++coeff_it) {
            (*coeff_it).second = 0.;
        }
    }


    for (const_Index1D_Coefficients1D_Hash_it it=y_v.begin(); it!=y_v.end(); ++it) {
        Index1D col_index_y = (*it).first;

        IndexSet<Index1D> LambdaRowSparse_y;
        LambdaRowSparse_y = lambdaTilde1d_PDE(col_index_y, basis.second, 1, std::max(col_index_y.j-1,basis.second.j0),
                                                          col_index_y.j+1, false);
        for (const_set1d_it row_y=LambdaRowSparse_y.begin(); row_y!=LambdaRowSparse_y.end(); ++row_y) {
            T id_y = identity_data1d_y(*row_y,col_index_y);
            T val1_y =  id_y;

            for (const_coeff1d_it coeff_x_it=(*it).second.begin(); coeff_x_it!=(*it).second.end(); ++coeff_x_it) {
                Index2D row_index((*coeff_x_it).first,*row_y);
                T val1 = val1_y * (*coeff_x_it).second;
                if (fabs(val1)>0) I_S_v[row_index] += val1;
            }
        }
    }
    for (const_coeff2d_it col=I_S_v.begin(); col!=I_S_v.end(); ++col) {
        Index1D_Coefficients1D_Hash_it p_x_I_S_v=x_I_S_v.find((*col).first.index1);
        if (p_x_I_S_v!=x_I_S_v.end()) {
            (*p_x_I_S_v).second.operator[]((*col).first.index2) = (*col).second;
        }
        else {
            Coefficients<Lexicographical,T,Index1D> coeff_y;
            coeff_y[(*col).first.index2] = (*col).second;
            x_I_S_v[(*col).first.index1] = coeff_y;
        }
    }

    for (const_Index1D_Coefficients1D_Hash_it it=x_I_S_v.begin(); it!=x_I_S_v.end(); ++it) {

        Index1D col_index_x = (*it).first;
        IndexSet<Index1D> LambdaRowSparse_x;
        LambdaRowSparse_x = lambdaTilde1d_PDE(col_index_x, basis.first, 1, std::max(col_index_x.j-1,basis.first.j0),
                                              col_index_x.j+1, false);
        for (const_set1d_it row_x=LambdaRowSparse_x.begin(); row_x!=LambdaRowSparse_x.end(); ++row_x) {
            T id_x = identity_data1d_x(*row_x,col_index_x);
            T dd_x = laplace_data1d_x(*row_x,col_index_x);
            T val_x = dd_x + 0.5*c*id_x;
            for (const_coeff1d_it coeff_y_it=(*it).second.begin(); coeff_y_it!=(*it).second.end(); ++coeff_y_it) {
                Index2D row_index(*row_x,(*coeff_y_it).first);
                T val = val_x * (*coeff_y_it).second;
                if (fabs(val)>0.) ret[row_index] += val;
            }
        }
    }


    for (coeff2d_it it=ret.begin(); it!=ret.end(); ++it) {

        (*it).second *=  this->prec((*it).first);
    }

    time.stop();
}

template <typename T, DomainType Domain1, DomainType Domain2>
void
AdaptiveHelmholtzOperatorOptimized2D<T,Primal,Domain1,SparseMulti,Primal,Domain2,SparseMulti>
::apply(const Coefficients<Lexicographical,T,Index2D> &v, T /*eps*/,
        const IndexSet<Index2D> &Lambda, Coefficients<Lexicographical,T,Index2D> &ret,
        cxxblas::Transpose /*trans*/)
{
    if (v.size()==0) return;


    IndexSet<Index1D> Lambda_x(SIZEHASHINDEX1D), Lambda_y(SIZEHASHINDEX1D);
    split(Lambda, Lambda_x, Lambda_y);
    const_set1d_it Lambda_x_end = Lambda_x.end();
    const_set1d_it Lambda_y_end = Lambda_y.end();
    const_set2d_it Lambda_end   = Lambda.end();

    //Index1D_Coefficients1D_Hash y_v(SIZEHASHINDEX1D);
    Index1D_Coefficients1D_Hash y_v;
    Coefficients<Lexicographical,T,Index2D> I_S_v(SIZEHASHINDEX2D);
    //Index1D_Coefficients1D_Hash x_I_S_v(SIZEHASHINDEX1D);
    Index1D_Coefficients1D_Hash x_I_S_v;

    Timer time;
    time.start();

    for (const_coeff2d_it col=v.begin(); col!=v.end(); ++col) {
        Index1D_Coefficients1D_Hash_it p_y_v=y_v.find((*col).first.index2);
        if (p_y_v!=y_v.end()) {
            (*p_y_v).second.operator[]((*col).first.index1) = this->prec((*col).first) * (*col).second;
        }
        else {
            Coefficients<Lexicographical,T,Index1D> coeff_x;
            y_v[(*col).first.index2] = coeff_x;
            p_y_v=y_v.find((*col).first.index2);
            //(*p_y_v).second.resize(500);
            (*p_y_v).second.operator[]((*col).first.index1) = this->prec((*col).first) * (*col).second;
            //coeff_x[(*col).first.index1] = this->prec((*col).first) * (*col).second;
            //y_v[(*col).first.index2] = coeff_x;
        }
    }

    for (const_Index1D_Coefficients1D_Hash_it it=y_v.begin(); it!=y_v.end(); ++it) {
        Index1D col_index_y = (*it).first;

        IndexSet<Index1D> LambdaRowSparse_y;
        LambdaRowSparse_y = lambdaTilde1d_PDE(col_index_y, basis.second, 1, std::max(col_index_y.j-1,basis.second.j0),
                                                          col_index_y.j+1, false);
        for (const_set1d_it row_y=LambdaRowSparse_y.begin(); row_y!=LambdaRowSparse_y.end(); ++row_y) {
            if (Lambda_y.find(*row_y)==Lambda_y_end) continue;
            T dd_y = laplace_data1d_y(*row_y,col_index_y);
            T id_y = identity_data1d_y(*row_y,col_index_y);
            T val2_y = (diffusion_y*dd_y + 0.5*c*id_y);

            for (const_coeff1d_it coeff_x_it=(*it).second.begin(); coeff_x_it!=(*it).second.end(); ++coeff_x_it) {
                Index2D row_index((*coeff_x_it).first,*row_y);
                T val2 = val2_y * (*coeff_x_it).second;
                if (fabs(val2)>0) I_S_v[row_index] += val2;
            }
        }
    }
    for (const_coeff2d_it col=I_S_v.begin(); col!=I_S_v.end(); ++col) {
        Index1D_Coefficients1D_Hash_it p_x_I_S_v=x_I_S_v.find((*col).first.index1);
        if (p_x_I_S_v!=x_I_S_v.end()) {
            (*p_x_I_S_v).second.operator[]((*col).first.index2) = (*col).second;
        }
        else {
            Coefficients<Lexicographical,T,Index1D> coeff_y;
            x_I_S_v[(*col).first.index1] = coeff_y;
            p_x_I_S_v=x_I_S_v.find((*col).first.index1);
            //(*p_x_I_S_v).second.resize(500);
            (*p_x_I_S_v).second.operator[]((*col).first.index2) = (*col).second;
            //coeff_y[(*col).first.index2] = (*col).second;
            //x_I_S_v[(*col).first.index1] = coeff_y;
        }
    }
    for (const_Index1D_Coefficients1D_Hash_it it=x_I_S_v.begin(); it!=x_I_S_v.end(); ++it) {
        Index1D col_index_x = (*it).first;
        IndexSet<Index1D> LambdaRowSparse_x;
        LambdaRowSparse_x = lambdaTilde1d_PDE(col_index_x, basis.first, 1, std::max(col_index_x.j-1,basis.first.j0),
                                              col_index_x.j+1, false);
        for (const_set1d_it row_x=LambdaRowSparse_x.begin(); row_x!=LambdaRowSparse_x.end(); ++row_x) {
            if (Lambda_x.find(*row_x)==Lambda_x_end) continue;
            T id_x = identity_data1d_x(*row_x,col_index_x);
            for (const_coeff1d_it coeff_y_it=(*it).second.begin(); coeff_y_it!=(*it).second.end(); ++coeff_y_it) {
                Index2D row_index(*row_x,(*coeff_y_it).first);
                if (Lambda.find(row_index)==Lambda_end) continue;
                T val = id_x * (*coeff_y_it).second;
                if (fabs(val)>0.) ret[row_index] += val;
            }
        }
    }

    for (coeff2d_it it=I_S_v.begin(); it!=I_S_v.end(); ++it) {
                (*it).second = 0.;
    }
    for (Index1D_Coefficients1D_Hash_it it=x_I_S_v.begin(); it!=x_I_S_v.end(); ++it) {
        for (coeff1d_it coeff_it=(*it).second.begin(); coeff_it!=(*it).second.end(); ++coeff_it) {
            (*coeff_it).second = 0.;
        }
    }


    for (const_Index1D_Coefficients1D_Hash_it it=y_v.begin(); it!=y_v.end(); ++it) {
        Index1D col_index_y = (*it).first;

        IndexSet<Index1D> LambdaRowSparse_y;
        LambdaRowSparse_y = lambdaTilde1d_PDE(col_index_y, basis.second, 1, std::max(col_index_y.j-1,basis.second.j0),
                                                          col_index_y.j+1, false);
        for (const_set1d_it row_y=LambdaRowSparse_y.begin(); row_y!=LambdaRowSparse_y.end(); ++row_y) {
            if (Lambda_y.find(*row_y)==Lambda_y_end) continue;
            T id_y = identity_data1d_y(*row_y,col_index_y);
            T val1_y =  id_y;

            for (const_coeff1d_it coeff_x_it=(*it).second.begin(); coeff_x_it!=(*it).second.end(); ++coeff_x_it) {
                Index2D row_index((*coeff_x_it).first,*row_y);
                T val1 = val1_y * (*coeff_x_it).second;
                if (fabs(val1)>0) I_S_v[row_index] += val1;
            }
        }
    }
    for (const_coeff2d_it col=I_S_v.begin(); col!=I_S_v.end(); ++col) {
        Index1D_Coefficients1D_Hash_it p_x_I_S_v=x_I_S_v.find((*col).first.index1);
        if (p_x_I_S_v!=x_I_S_v.end()) {
            (*p_x_I_S_v).second.operator[]((*col).first.index2) = (*col).second;
        }
        else {
            Coefficients<Lexicographical,T,Index1D> coeff_y;
            x_I_S_v[(*col).first.index1] = coeff_y;
            p_x_I_S_v=x_I_S_v.find((*col).first.index1);
            //(*p_x_I_S_v).second.resize(500);
            (*p_x_I_S_v).second.operator[]((*col).first.index2) = (*col).second;
            //coeff_y[(*col).first.index2] = (*col).second;
            //x_I_S_v[(*col).first.index1] = coeff_y;
        }
    }

    for (const_Index1D_Coefficients1D_Hash_it it=x_I_S_v.begin(); it!=x_I_S_v.end(); ++it) {

        Index1D col_index_x = (*it).first;
        IndexSet<Index1D> LambdaRowSparse_x;
        LambdaRowSparse_x = lambdaTilde1d_PDE(col_index_x, basis.first, 1, std::max(col_index_x.j-1,basis.first.j0),
                                              col_index_x.j+1, false);
        for (const_set1d_it row_x=LambdaRowSparse_x.begin(); row_x!=LambdaRowSparse_x.end(); ++row_x) {
            if (Lambda_x.find(*row_x)==Lambda_x_end) continue;
            T id_x = identity_data1d_x(*row_x,col_index_x);
            T dd_x = laplace_data1d_x(*row_x,col_index_x);
            T val_x = dd_x + 0.5*c*id_x;
            for (const_coeff1d_it coeff_y_it=(*it).second.begin(); coeff_y_it!=(*it).second.end(); ++coeff_y_it) {
                Index2D row_index(*row_x,(*coeff_y_it).first);
                if (Lambda.find(row_index)==Lambda_end) continue;
                T val = val_x * (*coeff_y_it).second;
                if (fabs(val)>0.) ret[row_index] += val;
            }
        }
    }
    for (coeff2d_it it=ret.begin(); it!=ret.end(); ++it) {

        (*it).second *=  this->prec((*it).first);
    }
    time.stop();
}

template <typename T, DomainType Domain1, DomainType Domain2>
void
AdaptiveHelmholtzOperatorOptimized2D<T,Primal,Domain1,SparseMulti,Primal,Domain2,SparseMulti>
::clear()
{
    Coefficients<Lexicographical,T,Index2D> tmp;
    P_data = tmp;
}

}   // namespace lawa
