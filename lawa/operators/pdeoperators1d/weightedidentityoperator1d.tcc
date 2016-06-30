namespace lawa{

template <typename T, typename Basis, QuadratureType Quad>
WeightedIdentityOperator1D<T, Basis, Quad>::WeightedIdentityOperator1D(const Basis& _basis,
                                                                       Function<T> weightFct,
                                                                       FLENS_DEFAULT_INDEXTYPE order,
                                                                       const T left, const T right)
    : basis(_basis), W(weightFct), integral(W, _basis, _basis, left, right)
{
    integral.quadrature.setOrder(order);
}

template <typename T, typename Basis, QuadratureType Quad>
T
WeightedIdentityOperator1D<T, Basis, Quad>::operator()(XType xtype1, FLENS_DEFAULT_INDEXTYPE j1, FLENS_DEFAULT_INDEXTYPE k1,
                                                       XType xtype2, FLENS_DEFAULT_INDEXTYPE j2, FLENS_DEFAULT_INDEXTYPE k2) const
{   
    return integral(j1, k1, xtype1, 0, j2, k2, xtype2, 0);
}

template <typename T, typename Basis, QuadratureType Quad>
T
WeightedIdentityOperator1D<T, Basis, Quad>::operator()(const Index1D &row_index,
                                                       const Index1D &col_index) const
{
    return this->operator()(row_index.xtype, row_index.j, row_index.k,
                            col_index.xtype, col_index.j, col_index.k);
}


} // namespace lawa

