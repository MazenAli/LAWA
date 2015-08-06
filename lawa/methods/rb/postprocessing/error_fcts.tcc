namespace lawa {

template<typename T, typename Basis2D, typename Prec>
T
L2_H1_error(Basis2D& basis2d, const Coefficients<Lexicographical,T,Index2D> &u, Prec& P,
            T (*u_ref)(T,T), T (*dx_u_ref)(T,T), T a_t, T b_t, int n_t, T a_x, T b_x, int n_x)
{
  typedef typename Coefficients<Lexicographical,T,Index2D>::const_iterator coeff_it;
  
  T L2_H1_error = 0.;
  
  T h_t = (b_t - a_t)/(n_t-1);
  T h_x = (b_x - a_x)/(n_x-1);
  for (int i = 0; i < n_t; ++i) {
      T H1_error = 0.;
      T t = a_t + i*h_t;
      for (int j = 0; j < n_x; ++j) {
          T x = a_x + j*h_x;
          
          T u_exact= u_ref(t,x);
          T dx_u_exact= dx_u_ref(t,x);
          
          T u_appr = 0.0;
          T dx_u_appr = 0.0;
          for (coeff_it it = u.begin(); it != u.end(); ++it) {
              XType xtype_t = (*it).first.index1.xtype;
              XType xtype_x = (*it).first.index2.xtype;
              int j_t = (*it).first.index1.j, k_t = (*it).first.index1.k;
              int j_x = (*it).first.index2.j, k_x = (*it).first.index2.k;

              T coeff = (*it).second;
              T prec = P((*it).first);
              T u_t_appr = basis2d.first.generator(xtype_t)(t,j_t,k_t,0);
              
              u_appr += prec * coeff * u_t_appr * basis2d.second.generator(xtype_x)(x,j_x,k_x,0);
              dx_u_appr += prec * coeff * u_t_appr * basis2d.second.generator(xtype_x)(x,j_x,k_x,1);
          }
          
          T factor_x = 1.;
          if((j == 0)||(j == n_x)){
            factor_x = 0.5;
          }
          
          H1_error += factor_x * ((u_appr - u_exact)*(u_appr - u_exact) + (dx_u_appr - dx_u_exact)*(dx_u_appr - dx_u_exact));
          //std::cout << " x = " << x << ": H1-error = " << H1_error << "  (err = " << (u_appr - u_exact) 
          //     << ", dx_err = " << dx_u_appr - dx_u_exact << ")"<<std::endl;
      }
      H1_error /= n_x;
      
      T factor_t = 1.;
      if((i == 0)||(i == n_t)){
        factor_t = 0.5;
      }
      L2_H1_error += factor_t * H1_error;
      //std::cout << "t = " << t << ": L2_H1_error = " << L2_H1_error << std::endl;
  } 
  
  L2_H1_error /= n_t;
  
  return sqrt(L2_H1_error);
}

} // namespace lawa
