namespace lawa {

template <typename SingularIntegral>
flens::GeMatrix<flens::FullStorage<long double,cxxblas::ColMajor> >
SingularQuadrature<SingularIntegral>::_legendreweights;

template <typename SingularIntegral>
flens::GeMatrix<flens::FullStorage<long double,cxxblas::ColMajor> >
SingularQuadrature<SingularIntegral>::_legendreknots;

template <typename SingularIntegral>
flens::DenseVector<flens::Array<FLENS_DEFAULT_INDEXTYPE> >
SingularQuadrature<SingularIntegral>::_hp_legendrenumofpoints;

template <typename SingularIntegral>
flens::GeMatrix<flens::FullStorage<long double,cxxblas::ColMajor> >
SingularQuadrature<SingularIntegral>::_hp_legendreweights;

template <typename SingularIntegral>
flens::GeMatrix<flens::FullStorage<long double,cxxblas::ColMajor> >
SingularQuadrature<SingularIntegral>::_hp_legendreknots;

template <typename SingularIntegral>
FLENS_DEFAULT_INDEXTYPE
SingularQuadrature<SingularIntegral>::_precalculated_order = 30;

template <typename SingularIntegral>
FLENS_DEFAULT_INDEXTYPE
SingularQuadrature<SingularIntegral>::_precalculated_n = 10;

template <typename SingularIntegral>
double
SingularQuadrature<SingularIntegral>::_precalculated_sigma = 0.2;

template <typename SingularIntegral>
double
SingularQuadrature<SingularIntegral>::_precalculated_mu = 0.5;

template <typename SingularIntegral>
SingularQuadrature<SingularIntegral>::SingularQuadrature(const SingularIntegral &_singularintegral)
: singularintegral(_singularintegral), _order_eta(8), _order(6), _n(6), _sigma(_precalculated_sigma),
  _mu(_precalculated_mu), _omega(0.01)
{
    _order_eta = (singularintegral.first.d+singularintegral.second.d)/2+2;
    std::cerr << "_order_eta = " << _order_eta << std::endl;
    _legendre(_precalculated_order);
    _hp_composite_legendre(_precalculated_n, _precalculated_sigma, _precalculated_mu);
}

template <typename SingularIntegral>
void
SingularQuadrature<SingularIntegral>::setLegendreOrder(FLENS_DEFAULT_INDEXTYPE order_eta)
{
    if (order_eta>=_precalculated_order) {
        _legendre(order_eta);
        _precalculated_order = order_eta;
    }
    _order_eta = order_eta;
}

template <typename SingularIntegral>
void
SingularQuadrature<SingularIntegral>::setParameters(FLENS_DEFAULT_INDEXTYPE order, FLENS_DEFAULT_INDEXTYPE n, double sigma,
                                                    double mu, double omega)
{
    if (order>=_precalculated_order) {
        _legendre(order);
        _precalculated_order = order;
    }
    _order = order;

    if (n > _precalculated_n || sigma != _precalculated_sigma || mu != _precalculated_mu) {
        _hp_composite_legendre(n, sigma, mu);
        _precalculated_n = n; _precalculated_sigma = sigma; _precalculated_mu = mu;
    }
    _n = n;
    _sigma = sigma;
    _mu = mu;
    _omega = omega;
}

template <typename SingularIntegral>
long double
SingularQuadrature<SingularIntegral>::operator()(long double a1, long double b1,
                                                 long double a2, long double b2, long double eps)
{
    //std::cerr << "[" << a1 << ", " << b1 << "], [" << a2 << ", " << b2 << "]" << std::endl;
    if (fabs(a1-a2)<1e-15 && fabs(b1-b2)<1e-15) {
        //std::cerr << "   _integrate_singular_diagonal [" << a1 << ", " << b1 << "], [" << a2 << ", " << b2 << "]" << std::endl;
        return _integrate_singular_diagonal(a1, b1, a2, b2, eps);
    }
    else if (fabs(a1-b2)<1e-15) {   // Corner singularity in (a1,b2)
        //std::cerr << "   _integrate_singular_corner for [" << a1 << ", " << b1 << "], [" << a2 << ", " << b2 << "]" << std::endl;
        long double h1 = b1-a1, h2 = b2-a2;
        if (h1 < h2) {
            //std::cerr << "     First case" << std::endl;
            return   _integrate_singular_corner(a1, b1, b2-h1, b2, (long double)1e-10)
                   + _integrate_nonsingular(a1, b1, a2, b2-h1, (long double)1e-10);
        }
        else if (h1 > h2) {
            //std::cerr << "     Second case" << std::endl;
            return   _integrate_singular_corner(a1, a1+h2, a2, b2, (long double)1e-10)
                   + _integrate_nonsingular(a1+h2, b1, a2, b2, (long double)1e-10);
        }
        else {  // Quadratic domain
            //std::cerr << "      quadratic domain" << std::endl;
            return _integrate_singular_corner(a1, b1, a2, b2, eps);
        }
    }
    else if (fabs(a2-b1)<1e-15) {   // Corner singularity in (b1,a2)
        //std::cerr << "   _integrate_singular_corner for [" << a1 << ", " << b1 << "], [" << a2 << ", " << b2 << "]" << std::endl;
        long double h1 = b1-a1, h2 = b2-a2;
        if (h1 < h2) {
            return   _integrate_singular_corner(a1, b1, a2, a2+h1, (long double)1e-10)
                   + _integrate_nonsingular(a1, b1, a2+h1, b2, (long double)1e-10);
        }
        else if (h1 > h2) {
            return   _integrate_singular_corner(b1-h2, b1, a2, b2, (long double)1e-10)
                   + _integrate_nonsingular(a1, b1-h2, a2, b2, (long double)1e-10);
        }
        else {  // Quadratic domain
            return _integrate_singular_corner(a1, b1, a2, b2, eps);
        }

    }
    else {
        //std::cerr << "   _integrate_nonsingular for [" << a1 << ", " << b1 << "], [" << a2 << ", " << b2 << "]" << std::endl;
        return _integrate_nonsingular(a1, b1, a2, b2, eps);
    }

}

template <typename SingularIntegral>
long double
SingularQuadrature<SingularIntegral>::_integrate_singular_diagonal(long double a1, long double b1,
                                                                   long double a2, long double b2,
                                                                   long double eps)
{
    long double ret = 0.L;
    FLENS_DEFAULT_INDEXTYPE hp_legendre_N = _hp_legendrenumofpoints(_n);
    long double h = (b1-a1);
    for (FLENS_DEFAULT_INDEXTYPE i=2; i<=hp_legendre_N; ++i) {      // skipping the first point = skip the first interval
        long double xi = _hp_legendreknots(_n,i);
        long double kernel_val_xi_Times_One_M_xi = singularintegral.kernel(h*xi) * (1-xi);
        long double kernel_val_Minus_xi_Times_One_M_xi = singularintegral.kernel(-h*xi)  * (1-xi);
        for (FLENS_DEFAULT_INDEXTYPE j=1; j<=_order_eta; ++j) {
            long double eta = _legendreknots(_order_eta,j);
            long double productweight = _hp_legendreweights(_n,i) * _legendreweights(_order_eta,j);

            long double val1 =   singularintegral.p1(h*(xi+eta*(1-xi))+a1)
                               * singularintegral.p2(h*(eta*(1-xi))+a2) * kernel_val_xi_Times_One_M_xi;
            ret += productweight * val1;
            long double val2 =   singularintegral.p1(h*(eta*(1-xi))+a1)
                               * singularintegral.p2(h*(xi+eta*(1-xi))+a2) * kernel_val_Minus_xi_Times_One_M_xi;
            ret += productweight * val2;
        }
    }
    return h*h*ret;
}

template <typename SingularIntegral>
long double
SingularQuadrature<SingularIntegral>::_integrate_singular_corner(long double a1, long double b1,
                                                                 long double a2, long double b2,
                                                                 long double eps)
{
    long double ret = 0.L;
    long double h = (b1-a1);
    long double _a1, _a2, _b1, _b2, _h1, _h2;
    long double factor = 1.L;
    if (fabs(a1-b2)<1e-15) {
        //std::cerr << "        corner (a1,b2)" << std::endl;
        _a1 = a1, _b1 = b1, _a2 = b2, _b2 = a2; _h1 = h, _h2 = -h;
    }
    else {
        _a1 = b1, _b1 = a1, _a2 = a2, _b2 = b2; _h1 = -h, _h2 = h;  factor = -1.L;
    }
    FLENS_DEFAULT_INDEXTYPE hp_legendre_N = _hp_legendrenumofpoints(_n);
    for (FLENS_DEFAULT_INDEXTYPE i=1; i<=hp_legendre_N; ++i) {
        long double xi = _hp_legendreknots(_n,i);
        long double kernel_val_xi_Times_xi = singularintegral.kernel(factor*h*xi) * xi;
        long double kernel_val_One_Plus_xi_Times_One_Minus_xi = singularintegral.kernel(factor*h*(1+xi)) * (1-xi);
        for (FLENS_DEFAULT_INDEXTYPE j=1; j<=_order_eta; ++j) {
            long double eta = _legendreknots(_order_eta,j);
            long double productweight = _hp_legendreweights(_n,i) * _legendreweights(_order_eta,j);

            long double val1 =   singularintegral.p1(_h1*(xi*eta)+_a1)
                               * singularintegral.p2(_h2*(xi*(1-eta))+_a2) * kernel_val_xi_Times_xi;
            ret += productweight * val1;
            long double val2 =  singularintegral.p1(_h1*(xi+eta*(1-xi))+_a1)
                              * singularintegral.p2(_h2*(1+eta*(xi-1))+_a2) * kernel_val_One_Plus_xi_Times_One_Minus_xi;
            ret += productweight * val2;
        }
    }
    //std::cerr << "   Num of points: " << _order_eta * hp_legendre_N << std::endl;
    return h*h*ret;
}

template <typename SingularIntegral>
long double
SingularQuadrature<SingularIntegral>::_integrate_nonsingular(long double a1, long double b1,
                                                             long double a2, long double b2,
                                                             long double eps)
{
    long double ret = 0.L;
    if ( (fabs(fabs(b1-a1)-fabs(b2-a2))<1e-15) ) {
        //std::cerr << "      quadratic domain." << std::endl;
        long double _a1, _a2, _b1, _b2, _h1, _h2, _h, c;
        long double h = (b1-a1);
        if (a1>b2) {
            //std::cerr << "      Special formula 1 called." << std::endl;
            _a1 = a1, _b1 = b1, _a2 = b2, _b2 = a2; _h1 = h, _h2 = -h, _h = h, c = (a1-b2)/h;
        }
        else if (a2>b1){
            //std::cerr << "      Special formula 2 called." << std::endl;
            _a1 = b1, _b1 = a1, _a2 = a2, _b2 = b2; _h1 = -h, _h2 = h, _h = -h, c = (a2-b1)/h;
        }
        else {
            std::cerr << "SingularQuadrature<SingularIntegral>::_integrate_nonsingular: "
                      << "unknown configuration." << std::endl;
            exit(1);
        }
        for (FLENS_DEFAULT_INDEXTYPE i=1; i<=_order; ++i) {
            long double xi = _legendreknots(_order,i);
            long double kernel_val_xi_Times_xi = singularintegral.kernel(_h*(xi+c)) * xi;
            if (xi<1e-14) kernel_val_xi_Times_xi = 0.L;
            long double kernel_val_One_Plus_xi_Times_One_Minus_xi = singularintegral.kernel(_h*(1+xi+c)) * (1-xi);
            for (FLENS_DEFAULT_INDEXTYPE j=1; j<=_order_eta; ++j) {
                long double eta = _legendreknots(_order_eta,j);
                long double productweight = _legendreweights(_order,i) * _legendreweights(_order_eta,j);

                long double val1 =   singularintegral.p1(_h1*(xi*eta)+_a1)
                                   * singularintegral.p2(_h2*(xi*(1-eta))+_a2) * kernel_val_xi_Times_xi;
                ret += productweight * val1;
                long double val2 =  singularintegral.p1(_h1*(xi+eta*(1-xi))+_a1)
                                  * singularintegral.p2(_h2*(1+eta*(xi-1))+_a2) * kernel_val_One_Plus_xi_Times_One_Minus_xi;
                ret += productweight * val2;
            }
        }
//        std::cerr << "   Num of points: " << 2 * _order_eta * _order << std::endl;
        return h*h*ret;
    }
    else {
        long double h1 = (b1-a1), h2 = (b2-a2);
        //std::cerr << "      non-quadratic domain." << std::endl;
        for (FLENS_DEFAULT_INDEXTYPE i=1; i<=_order; ++i) {
            long double xi = _legendreknots(_order,i);
            for (FLENS_DEFAULT_INDEXTYPE j=1; j<=_order; ++j) {
                long double eta = _legendreknots(_order,j);
                long double productweight = _legendreweights(_order,i) * _legendreweights(_order,j);
                long double val =   singularintegral.p1(h1*xi+a1) * singularintegral.p2(h2*eta+a2)
                                  * singularintegral.kernel(h1*xi+a1-h2*eta-a2);
                ret += productweight * val;
            }
        }
        return h1*h2*ret;
    }
}

template <typename SingularIntegral>
void
SingularQuadrature<SingularIntegral>::_legendre(FLENS_DEFAULT_INDEXTYPE order)
{
    //std::cerr << "_legendre called for order = " << order << std::endl;
    long double eps = Const<long double>::EQUALITY_EPS;
    flens::GeMatrix<flens::FullStorage<long double,cxxblas::ColMajor> > initiallegendreweights;
    flens::GeMatrix<flens::FullStorage<long double,cxxblas::ColMajor> > initiallegendreknots;
    initiallegendreweights.engine().resize(order,order);
    initiallegendreknots.engine().resize(order,order);
    _legendreknots.engine().resize(order, order);
    _legendreweights.engine().resize(order, order);

    long double x1 = -1.L, x2 =  1.L;
    for (FLENS_DEFAULT_INDEXTYPE k=1; k<=order; ++k) {
        FLENS_DEFAULT_INDEXTYPE     m = (k+1)/2;
        long double xm = 0.5L * (x2+x1), xl = 0.5L * (x2-x1);
        for (FLENS_DEFAULT_INDEXTYPE i=1; i<=m; ++i) {
            long double z = cos(M_PI*(i-0.25L)/(k+0.5L)), z1, pp;
            do {
                long double p1 = 1.0L, p2 = 2.0L;
                for (FLENS_DEFAULT_INDEXTYPE j=1; j<=k; ++j) {
                    long double p3 = p2;
                    p2 = p1;
                    p1 = ((2.0L*j-1.0L)*z*p2-(j-1.0L)*p3)/j;
                }
                pp = k * (z*p1-p2)/(z*z-1.0L);
                z1 = z;
                z = z1-p1/pp;
            } while (fabs(z-z1) > eps);
            initiallegendreknots(k,i)     = xm - xl*z;
            initiallegendreknots(k,k+1-i) = xm + xl*z;
            initiallegendreweights(k,i)     = 2.0L*xl/((1.0L-z*z)*pp*pp);
            initiallegendreweights(k,k+1-i) = initiallegendreweights(k,i);
        }
    }

    for (FLENS_DEFAULT_INDEXTYPE k=1; k<=order; ++k) {
        for (FLENS_DEFAULT_INDEXTYPE l=1; l<=order; ++l) {
            _legendreknots(k,l)   = 0.5L*initiallegendreknots(k,l)+0.5L;
            _legendreweights(k,l) = 0.5L*initiallegendreweights(k,l);
        }
    }
}

template <typename SingularIntegral>
void
SingularQuadrature<SingularIntegral>::_hp_composite_legendre(FLENS_DEFAULT_INDEXTYPE max_n, double sigma, double mu)
{
    using flens::_;

    FLENS_DEFAULT_INDEXTYPE N=0;
    for (FLENS_DEFAULT_INDEXTYPE j=0; j<=max_n; ++j) {
        N += std::max((FLENS_DEFAULT_INDEXTYPE)(mu*j)+1,(FLENS_DEFAULT_INDEXTYPE) 1);
    }
    _hp_legendrenumofpoints.engine().resize(max_n);
    _hp_legendreknots.engine().resize(max_n,N);
    _hp_legendreweights.engine().resize(max_n,N);

    FLENS_DEFAULT_INDEXTYPE max_q = std::max((FLENS_DEFAULT_INDEXTYPE)(mu*max_n)+1,(FLENS_DEFAULT_INDEXTYPE) 1);
    if (max_q>=_precalculated_order) {
        _legendre(max_q);
        _precalculated_order = max_q;
    }

    for (FLENS_DEFAULT_INDEXTYPE n=1; n<=max_n; ++n) {
        N = 0;
        DenseVector q(_(0,n));
        for (FLENS_DEFAULT_INDEXTYPE j=0; j<=n; ++j) {
            FLENS_DEFAULT_INDEXTYPE qj = std::max((FLENS_DEFAULT_INDEXTYPE)(mu*j)+1,(FLENS_DEFAULT_INDEXTYPE) 1);
            q(j) = qj;
            N += qj;
        }
        _hp_legendrenumofpoints(n) = N;
        //std::cerr << "q = " << q << std::endl;
        //std::cerr << "_hp_N = " << _hp_legendrenumofpoints(n) << std::endl;
        FLENS_DEFAULT_INDEXTYPE i=1;
        for (FLENS_DEFAULT_INDEXTYPE j=0; j<=n; ++j) {
            FLENS_DEFAULT_INDEXTYPE qj = q(j);
            long double a = std::pow((long double)sigma,n+1-j), b = std::pow((long double)sigma,n-j);
            if (j==0) {
                a = 0.L; b = std::pow((long double)sigma,n);
            }
            long double h = b-a;
            for (FLENS_DEFAULT_INDEXTYPE l=1; l<=qj; ++l) {
                _hp_legendreknots(n,i)  =  h*_legendreknots(qj,l) + a;
                _hp_legendreweights(n,i) = h*_legendreweights(qj,l);
                ++i;
            }
        }
    }
    return;
}

}   // namespace lawa


/*
 *
 *     for (FLENS_DEFAULT_INDEXTYPE i=1; i<=_order; ++i) {
        long double xi = _legendreknots(_order,i);
        long double kernel_val_xi_Times_xi = singularintegral.kernel(h*xi) * xi;
        if (xi<1e-14) kernel_val_xi_Times_xi = 0.L;
        long double kernel_val_One_Plus_xi_Times_One_Minus_xi = singularintegral.kernel(h*(1+xi)) * (1-xi);
        for (FLENS_DEFAULT_INDEXTYPE j=1; j<=_order_eta; ++j) {
            long double eta = _legendreknots(_order_eta,j);
            long double productweight = _legendreweights(_order,i) * _legendreweights(_order_eta,j);

            long double val1 =   singularintegral.p1(_h1*(xi*eta)+_a1)
                               * singularintegral.p2(_h2*(xi*(1-eta))+_a2) * kernel_val_xi_Times_xi;
            ret += productweight * val1;
            long double val2 =  singularintegral.p1(_h1*(xi+eta*(1-xi))+_a1)
                              * singularintegral.p2(_h2*(1+eta*(xi-1))+_a2) * kernel_val_One_Plus_xi_Times_One_Minus_xi;
            ret += productweight * val2;
        }
    }
    std::cerr << "   Num of points: " << _order_eta * _order << std::endl;
    return h*h*ret;
 *
 */
