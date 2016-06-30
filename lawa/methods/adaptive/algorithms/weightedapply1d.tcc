namespace lawa {

template <typename T, typename Basis1D, typename Parameters, typename MA>
SYM_WEIGHTED_APPLY_1D<T,Basis1D,Parameters,MA>::SYM_WEIGHTED_APPLY_1D(const Parameters &_parameters,
                                                                      const Basis1D &_basis, MA &_A)
: parameters(_parameters), basis(_basis), A(_A)
{
}

template <typename T, typename Basis1D, typename Parameters, typename MA>
Coefficients<Lexicographical,T,Index1D>
SYM_WEIGHTED_APPLY_1D<T,Basis1D,Parameters,MA>::operator()
                                                (const Coefficients<Lexicographical,T,Index1D> &v,
                                                 FLENS_DEFAULT_INDEXTYPE k)
{
    //FLENS_DEFAULT_INDEXTYPE d =basis.d;
    //FLENS_DEFAULT_INDEXTYPE d_=basis.d_;
    Coefficients<Lexicographical,T,Index1D> ret;
    if (v.size() == 0) return ret;

    Coefficients<AbsoluteValue,T,Index1D> temp;
    temp = v;
    FLENS_DEFAULT_INDEXTYPE s = 0, count = 0;
    for (const_coeff_abs_it it = temp.begin(); (it != temp.end()) && (s<=k); ++it) {
        IndexSet<Index1D> Lambda_v;
        FLENS_DEFAULT_INDEXTYPE s_tilde_level    = k-s;
        //FLENS_DEFAULT_INDEXTYPE s_tilde_singsupp = (FLENS_DEFAULT_INDEXTYPE)((d-1.5)*s_tilde_level/(1.+d_))+1;
        Lambda_v=lambdaTilde1d_WeightedPDE((*it).second, basis, s_tilde_level,
                                           basis.j0, (*it).second.j+(k-s)+1);
        for (const_set_it mu = Lambda_v.begin(); mu != Lambda_v.end(); ++mu) {
            ret[*mu] += A(*mu, (*it).second) * (*it).first;
        }
        ++count;
        s = (FLENS_DEFAULT_INDEXTYPE)(log(T(count))/log(T(2))) + 1;
    }

    return ret;
}

template <typename T, typename Basis1D, typename Parameters, typename MA>
Coefficients<Lexicographical,T,Index1D>
SYM_WEIGHTED_APPLY_1D<T,Basis1D,Parameters,MA>::operator()
                                                (const Coefficients<Lexicographical,T,Index1D> &v,
                                                 T eps)
{
    Coefficients<AbsoluteValue,T,Index1D> v_abs;
    v_abs = v;
    FLENS_DEFAULT_INDEXTYPE k = this->findK(v_abs, eps);
    Coefficients<Lexicographical,T,Index1D> ret;
    ret = this->operator()(v, k);
    return ret;
}

template <typename T, typename Basis1D, typename Parameters, typename MA>
FLENS_DEFAULT_INDEXTYPE
SYM_WEIGHTED_APPLY_1D<T,Basis1D,Parameters,MA>::findK(const Coefficients<AbsoluteValue,T,Index1D> &v, T eps) {
    FLENS_DEFAULT_INDEXTYPE d=basis.d;
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
        R_k *= parameters.CA;
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
            std::cout << "   findK ==> k = " << k << ", k_eps = " << k_eps << std::endl;
            FLENS_DEFAULT_INDEXTYPE maxlevel=22;
            if (d==2)         {    maxlevel=25; }
            else if (d==3)    {   maxlevel=19; }    //for non-singular examples, also lower values are possible
            return std::min(k,maxlevel);
        }
    }
    return std::min(k_eps,(FLENS_DEFAULT_INDEXTYPE) 60);    //higher level differences result in translation indices that cannot be stored in FLENS_DEFAULT_INDEXTYPE.
}

}   // namespace lawa
