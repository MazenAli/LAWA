namespace lawa {
    
template <typename T>
T 
HomogeneousRHS<T>::operator()(T /*time*/, XType /*xtype*/, FLENS_DEFAULT_INDEXTYPE /*j*/, FLENS_DEFAULT_INDEXTYPE /*k*/) const{
    return 0;
}

} // namespace lawa

