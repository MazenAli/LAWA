namespace lawa{
    
template<MethodType Method, typename FirstBasis, typename SecondBasis>
TensorBasis2D<Method, FirstBasis, SecondBasis>::TensorBasis2D(const FirstBasis &_basis1,
                                                              const SecondBasis &_basis2)
    : first(_basis1), second(_basis2), d(std::max(_basis1.d, _basis2.d))
{}

} // namespace lawa

