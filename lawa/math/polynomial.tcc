namespace lawa {

template <typename T>
Polynomial<T>::Polynomial(FLENS_DEFAULT_INDEXTYPE n)
    : _coefficients(flens::DenseVector<flens::Array<T> >(flens::_(0,n)))
{
}

template <typename T>
Polynomial<T>::Polynomial(const flens::DenseVector<flens::Array<T> > &coefficients)
    : _coefficients(coefficients)
{
}


template <typename T>
const T &
Polynomial<T>::operator()(FLENS_DEFAULT_INDEXTYPE n) const
{
    assert(0<=n);
    assert(n<=this->degree());

    return _coefficients(n);
}

template <typename T>
T &
Polynomial<T>::operator()(FLENS_DEFAULT_INDEXTYPE n)
{
    assert(0<=n);
    assert(n<=this->degree());

    return _coefficients(n);
}


template <typename T>
Polynomial<T> &
Polynomial<T>::operator+=(const Polynomial<T> &rhs)
{
    if (this->degree()>=rhs.degree()) {
        for (FLENS_DEFAULT_INDEXTYPE i=0; i<=rhs.degree(); ++i) {
            _coefficients(i) += rhs._coefficients(i);
        }
    } else {
        flens::DenseVector<flens::Array<T> > tmp = rhs._coefficients;
        for (FLENS_DEFAULT_INDEXTYPE i=0; i<=this->degree(); ++i) {
            tmp(i) += _coefficients(i);
        }
        this->_coefficients = tmp;
    }
    return *this;
}

template <typename T>
FLENS_DEFAULT_INDEXTYPE
Polynomial<T>::degree() const
{
    return _coefficients.lastIndex();
}

//------------------------------------------------------------------------------

template <typename T>
Polynomial<T>
operator*(const Polynomial<T> &lhs, const Polynomial<T> &rhs)
{
    FLENS_DEFAULT_INDEXTYPE degree = lhs.degree() + rhs.degree();

    Polynomial<T> res(degree);
    for (FLENS_DEFAULT_INDEXTYPE i=0; i<=lhs.degree(); i++) {
        for (FLENS_DEFAULT_INDEXTYPE j=0; j<=rhs.degree(); j++) {
            res(i+j) += lhs(i)*rhs(j);
        }
    }
    return res;
}

template <typename S, typename T>
Polynomial<T>
operator*(const S &lhs, const Polynomial<T> &rhs)
{
    Polynomial<T> res(rhs);
    for (FLENS_DEFAULT_INDEXTYPE i=0; i<=rhs.degree(); ++i) {
        res(i) *= lhs;
    }
    return res;
}

template <typename T>
Polynomial<T>
pow(const Polynomial<T> &p, FLENS_DEFAULT_INDEXTYPE n)
{
    if (!n) {
        Polynomial<T> res(0);
        res(0) = 1.;
        return res;
    }
    Polynomial<T> res(p);
    for (FLENS_DEFAULT_INDEXTYPE k=1; k<n; ++k) {
        res = res*p;
    }
    return res;
}

template <typename T>
std::ostream &
operator<<(std::ostream &out, const Polynomial<T> &rhs)
{
    out << "(";
    for (FLENS_DEFAULT_INDEXTYPE i=0; i<=rhs.degree(); ++i) {
        out << rhs(i) << " ";
    }
    out << ")";
    return out;
}

} // namespace lawa

