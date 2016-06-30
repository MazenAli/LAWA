# include <fstream>
# include <string>

namespace webbur
{
  void binary_vector_next ( FLENS_DEFAULT_INDEXTYPE n, FLENS_DEFAULT_INDEXTYPE bvec[] );

  void chebyshev1_compute ( FLENS_DEFAULT_INDEXTYPE order, double x[], double w[] );
  void chebyshev1_compute_np ( FLENS_DEFAULT_INDEXTYPE order, FLENS_DEFAULT_INDEXTYPE np, double p[], double x[], double w[] );
  void chebyshev1_compute_points ( FLENS_DEFAULT_INDEXTYPE order, double x[] );
  void chebyshev1_compute_points_np ( FLENS_DEFAULT_INDEXTYPE order, FLENS_DEFAULT_INDEXTYPE np, double p[], double x[] );
  void chebyshev1_compute_weights ( FLENS_DEFAULT_INDEXTYPE order, double w[] );
  void chebyshev1_compute_weights_np ( FLENS_DEFAULT_INDEXTYPE order, FLENS_DEFAULT_INDEXTYPE np, double p[], double w[] );
  double chebyshev1_integral ( FLENS_DEFAULT_INDEXTYPE expon );

  void chebyshev2_compute ( FLENS_DEFAULT_INDEXTYPE order, double x[], double w[] );
  void chebyshev2_compute_np ( FLENS_DEFAULT_INDEXTYPE order, FLENS_DEFAULT_INDEXTYPE np, double p[], double x[], double w[] );
  void chebyshev2_compute_points ( FLENS_DEFAULT_INDEXTYPE order, double x[] );
  void chebyshev2_compute_points_np ( FLENS_DEFAULT_INDEXTYPE order, FLENS_DEFAULT_INDEXTYPE np, double p[], double x[] );
  void chebyshev2_compute_weights ( FLENS_DEFAULT_INDEXTYPE order, double w[] );
  void chebyshev2_compute_weights_np ( FLENS_DEFAULT_INDEXTYPE order, FLENS_DEFAULT_INDEXTYPE np, double p[], double w[] );
  double chebyshev2_integral ( FLENS_DEFAULT_INDEXTYPE expon );

  void clenshaw_curtis_compute ( FLENS_DEFAULT_INDEXTYPE order, double x[], double w[] );
  void clenshaw_curtis_compute_np ( FLENS_DEFAULT_INDEXTYPE order, FLENS_DEFAULT_INDEXTYPE np, double p[], double x[], double w[] );
  void clenshaw_curtis_compute_points ( FLENS_DEFAULT_INDEXTYPE order, double x[] );
  void clenshaw_curtis_compute_points_np ( FLENS_DEFAULT_INDEXTYPE order, FLENS_DEFAULT_INDEXTYPE np, double p[], double x[] );
  void clenshaw_curtis_compute_weights ( FLENS_DEFAULT_INDEXTYPE order, double w[] );
  void clenshaw_curtis_compute_weights_np ( FLENS_DEFAULT_INDEXTYPE order, FLENS_DEFAULT_INDEXTYPE np, double p[], double w[] );

  void comp_next ( FLENS_DEFAULT_INDEXTYPE n, FLENS_DEFAULT_INDEXTYPE k, FLENS_DEFAULT_INDEXTYPE a[], bool *more, FLENS_DEFAULT_INDEXTYPE *h, FLENS_DEFAULT_INDEXTYPE *t );

  double cpu_time ( );

  void fejer2_compute ( FLENS_DEFAULT_INDEXTYPE order, double x[], double w[] );
  void fejer2_compute_np ( FLENS_DEFAULT_INDEXTYPE order, FLENS_DEFAULT_INDEXTYPE np, double p[], double x[], double w[] );
  void fejer2_compute_points ( FLENS_DEFAULT_INDEXTYPE order, double x[] );
  void fejer2_compute_points_np ( FLENS_DEFAULT_INDEXTYPE order, FLENS_DEFAULT_INDEXTYPE np, double p[], double x[] );
  void fejer2_compute_weights ( FLENS_DEFAULT_INDEXTYPE order, double w[] );
  void fejer2_compute_weights_np ( FLENS_DEFAULT_INDEXTYPE order, FLENS_DEFAULT_INDEXTYPE np, double p[], double w[] );

  void gegenbauer_compute ( FLENS_DEFAULT_INDEXTYPE order, double alpha, double x[], double w[] );
  void gegenbauer_compute_np ( FLENS_DEFAULT_INDEXTYPE order, FLENS_DEFAULT_INDEXTYPE np, double p[], double x[], double w[] );
  void gegenbauer_compute_points ( FLENS_DEFAULT_INDEXTYPE order, double alpha, double x[] );
  void gegenbauer_compute_points_np ( FLENS_DEFAULT_INDEXTYPE order, FLENS_DEFAULT_INDEXTYPE np, double p[], double x[] );
  void gegenbauer_compute_weights ( FLENS_DEFAULT_INDEXTYPE order, double alpha, double w[] );
  void gegenbauer_compute_weights_np ( FLENS_DEFAULT_INDEXTYPE order, FLENS_DEFAULT_INDEXTYPE np, double p[], double w[] );
  double gegenbauer_integral ( FLENS_DEFAULT_INDEXTYPE expon, double alpha );
  void gegenbauer_recur ( double *p2, double *dp2, double *p1, double x, 
    FLENS_DEFAULT_INDEXTYPE order, double alpha, double c[] );
  void gegenbauer_root ( double *x, FLENS_DEFAULT_INDEXTYPE order, double alpha,  double *dp2, 
    double *p1, double c[] );

  void gen_hermite_compute ( FLENS_DEFAULT_INDEXTYPE order, double alpha, double x[], double w[] );
  void gen_hermite_compute_np ( FLENS_DEFAULT_INDEXTYPE order, FLENS_DEFAULT_INDEXTYPE np, double p[], double x[], double w[] );
  void gen_hermite_compute_points ( FLENS_DEFAULT_INDEXTYPE order, double alpha, double x[] );
  void gen_hermite_compute_points_np ( FLENS_DEFAULT_INDEXTYPE order, FLENS_DEFAULT_INDEXTYPE np, double p[], double x[] );
  void gen_hermite_compute_weights ( FLENS_DEFAULT_INDEXTYPE order, double alpha, double w[] );
  void gen_hermite_compute_weights_np ( FLENS_DEFAULT_INDEXTYPE order, FLENS_DEFAULT_INDEXTYPE np, double p[], double w[] );
  double gen_hermite_integral ( FLENS_DEFAULT_INDEXTYPE expon, double alpha );

  void gen_laguerre_compute ( FLENS_DEFAULT_INDEXTYPE order, double alpha, double x[], double w[] );
  void gen_laguerre_compute_np ( FLENS_DEFAULT_INDEXTYPE order, FLENS_DEFAULT_INDEXTYPE np, double p[], double x[], double w[] );
  void gen_laguerre_compute_points ( FLENS_DEFAULT_INDEXTYPE order, double alpha, double x[] );
  void gen_laguerre_compute_points_np ( FLENS_DEFAULT_INDEXTYPE order, FLENS_DEFAULT_INDEXTYPE np, double p[], double x[] );
  void gen_laguerre_compute_weights ( FLENS_DEFAULT_INDEXTYPE order, double alpha, double w[] );
  void gen_laguerre_compute_weights_np ( FLENS_DEFAULT_INDEXTYPE order, FLENS_DEFAULT_INDEXTYPE np, double p[], double w[] );
  double gen_laguerre_integral ( FLENS_DEFAULT_INDEXTYPE expon, double alpha );
  void gen_laguerre_recur ( double *p2, double *dp2, double *p1, double x, 
    FLENS_DEFAULT_INDEXTYPE order, double alpha, double b[], double c[] );
  void gen_laguerre_root ( double *x, FLENS_DEFAULT_INDEXTYPE order, double alpha, double *dp2, 
    double *p1, double b[], double c[] );

  void hermite_compute ( FLENS_DEFAULT_INDEXTYPE order, double x[], double w[] );
  void hermite_compute_np ( FLENS_DEFAULT_INDEXTYPE order, FLENS_DEFAULT_INDEXTYPE np, double p[], double x[], double w[] );
  void hermite_compute_points ( FLENS_DEFAULT_INDEXTYPE order, double x[] );
  void hermite_compute_points_np ( FLENS_DEFAULT_INDEXTYPE order, FLENS_DEFAULT_INDEXTYPE np, double p[], double x[] );
  void hermite_compute_weights ( FLENS_DEFAULT_INDEXTYPE order, double w[] );
  void hermite_compute_weights_np ( FLENS_DEFAULT_INDEXTYPE order, FLENS_DEFAULT_INDEXTYPE np, double p[], double w[] );

  void hermite_genz_keister_lookup ( FLENS_DEFAULT_INDEXTYPE n, double x[], double w[] );
  void hermite_genz_keister_lookup_points ( FLENS_DEFAULT_INDEXTYPE n, double x[] );
  void hermite_genz_keister_lookup_points_np ( FLENS_DEFAULT_INDEXTYPE n, FLENS_DEFAULT_INDEXTYPE np, double p[], 
    double x[] );
  void hermite_genz_keister_lookup_weights ( FLENS_DEFAULT_INDEXTYPE n, double w[] );
  void hermite_genz_keister_lookup_weights_np ( FLENS_DEFAULT_INDEXTYPE n, FLENS_DEFAULT_INDEXTYPE np, double p[], 
    double w[] );

  double hermite_integral ( FLENS_DEFAULT_INDEXTYPE n );
  void hermite_lookup ( FLENS_DEFAULT_INDEXTYPE n, double x[], double w[] );
  void hermite_lookup_points ( FLENS_DEFAULT_INDEXTYPE n, double x[] );
  void hermite_lookup_weights ( FLENS_DEFAULT_INDEXTYPE n, double w[] );
  void hermite_recur ( double *p2, double *dp2, double *p1, double x, 
    FLENS_DEFAULT_INDEXTYPE order );
  void hermite_root ( double *x, FLENS_DEFAULT_INDEXTYPE order, double *dp2, double *p1 );

  FLENS_DEFAULT_INDEXTYPE i4_max ( FLENS_DEFAULT_INDEXTYPE i1, FLENS_DEFAULT_INDEXTYPE i2 );
  FLENS_DEFAULT_INDEXTYPE i4_min ( FLENS_DEFAULT_INDEXTYPE i1, FLENS_DEFAULT_INDEXTYPE i2 );
  FLENS_DEFAULT_INDEXTYPE i4_power ( FLENS_DEFAULT_INDEXTYPE i, FLENS_DEFAULT_INDEXTYPE j );

  void i4mat_write ( std::string output_filename, FLENS_DEFAULT_INDEXTYPE m, FLENS_DEFAULT_INDEXTYPE n, FLENS_DEFAULT_INDEXTYPE table[] );

  FLENS_DEFAULT_INDEXTYPE *i4vec_add_new ( FLENS_DEFAULT_INDEXTYPE n, FLENS_DEFAULT_INDEXTYPE a[], FLENS_DEFAULT_INDEXTYPE b[] );
  bool i4vec_any_lt ( FLENS_DEFAULT_INDEXTYPE n, FLENS_DEFAULT_INDEXTYPE a[], FLENS_DEFAULT_INDEXTYPE b[] );
  void i4vec_print ( FLENS_DEFAULT_INDEXTYPE n, FLENS_DEFAULT_INDEXTYPE a[], std::string title );
  FLENS_DEFAULT_INDEXTYPE i4vec_product ( FLENS_DEFAULT_INDEXTYPE n, FLENS_DEFAULT_INDEXTYPE a[] );
  FLENS_DEFAULT_INDEXTYPE i4vec_sum ( FLENS_DEFAULT_INDEXTYPE n, FLENS_DEFAULT_INDEXTYPE a[] );
  void i4vec_zero ( FLENS_DEFAULT_INDEXTYPE n, FLENS_DEFAULT_INDEXTYPE a[] );
  FLENS_DEFAULT_INDEXTYPE *i4vec_zero_new ( FLENS_DEFAULT_INDEXTYPE n );

  void jacobi_compute ( FLENS_DEFAULT_INDEXTYPE order, double alpha, double beta, double x[], 
    double w[] );
  void jacobi_compute_np ( FLENS_DEFAULT_INDEXTYPE order, FLENS_DEFAULT_INDEXTYPE np, double p[], double x[], double w[] );
  void jacobi_compute_points ( FLENS_DEFAULT_INDEXTYPE order, double alpha, double beta, 
    double x[] );
  void jacobi_compute_points_np ( FLENS_DEFAULT_INDEXTYPE order, FLENS_DEFAULT_INDEXTYPE np, double p[], double x[] );
  void jacobi_compute_weights ( FLENS_DEFAULT_INDEXTYPE order, double alpha, double beta, 
    double w[] );
  void jacobi_compute_weights_np ( FLENS_DEFAULT_INDEXTYPE order, FLENS_DEFAULT_INDEXTYPE np, double p[], double w[] );
  double jacobi_integral ( FLENS_DEFAULT_INDEXTYPE expon, double alpha, double beta );
  void jacobi_recur ( double *p2, double *dp2, double *p1, double x, FLENS_DEFAULT_INDEXTYPE order, 
    double alpha, double beta, double b[], double c[] );
  void jacobi_root ( double *x, FLENS_DEFAULT_INDEXTYPE order, double alpha, double beta, 
    double *dp2, double *p1, double b[], double c[] );

  void laguerre_compute ( FLENS_DEFAULT_INDEXTYPE order, double x[], double w[] );
  void laguerre_compute_np ( FLENS_DEFAULT_INDEXTYPE order, FLENS_DEFAULT_INDEXTYPE np, double p[], double x[], double w[] );
  void laguerre_compute_points ( FLENS_DEFAULT_INDEXTYPE order, double x[] );
  void laguerre_compute_points_np ( FLENS_DEFAULT_INDEXTYPE order, FLENS_DEFAULT_INDEXTYPE np, double p[], double x[] );
  void laguerre_compute_weights ( FLENS_DEFAULT_INDEXTYPE order, double w[] );
  void laguerre_compute_weights_np ( FLENS_DEFAULT_INDEXTYPE order, FLENS_DEFAULT_INDEXTYPE np, double p[], double w[] );
  double laguerre_integral ( FLENS_DEFAULT_INDEXTYPE expon );
  void laguerre_lookup ( FLENS_DEFAULT_INDEXTYPE n, double x[], double w[] );
  void laguerre_lookup_points ( FLENS_DEFAULT_INDEXTYPE n, double x[] );
  void laguerre_lookup_weights ( FLENS_DEFAULT_INDEXTYPE n, double w[] );
  void laguerre_recur ( double *p2, double *dp2, double *p1, double x, 
    FLENS_DEFAULT_INDEXTYPE order, double b[], double c[] );
  void laguerre_root ( double *x, FLENS_DEFAULT_INDEXTYPE order, double *dp2, double *p1, 
    double b[], double c[] );

  void legendre_compute ( FLENS_DEFAULT_INDEXTYPE order, double x[], double w[] );
  void legendre_compute_np ( FLENS_DEFAULT_INDEXTYPE order, FLENS_DEFAULT_INDEXTYPE np, double p[], double x[], double w[] );
  void legendre_compute_points ( FLENS_DEFAULT_INDEXTYPE order, double x[] );
  void legendre_compute_points_np ( FLENS_DEFAULT_INDEXTYPE order, FLENS_DEFAULT_INDEXTYPE np, double p[], double x[] );
  void legendre_compute_weights ( FLENS_DEFAULT_INDEXTYPE order, double w[] );
  void legendre_compute_weights_np ( FLENS_DEFAULT_INDEXTYPE order, FLENS_DEFAULT_INDEXTYPE np, double p[], double w[] );
  double legendre_integral ( FLENS_DEFAULT_INDEXTYPE expon );
  void legendre_lookup ( FLENS_DEFAULT_INDEXTYPE n, double x[], double w[] );
  void legendre_lookup_points ( FLENS_DEFAULT_INDEXTYPE n, double x[] );
  void legendre_lookup_weights ( FLENS_DEFAULT_INDEXTYPE n, double w[] );

  void level_growth_to_order ( FLENS_DEFAULT_INDEXTYPE dim_num, FLENS_DEFAULT_INDEXTYPE level[], FLENS_DEFAULT_INDEXTYPE rule[], FLENS_DEFAULT_INDEXTYPE growth[],
    FLENS_DEFAULT_INDEXTYPE order[] );
  void level_to_order_default ( FLENS_DEFAULT_INDEXTYPE dim_num, FLENS_DEFAULT_INDEXTYPE level[], FLENS_DEFAULT_INDEXTYPE rule[], 
    FLENS_DEFAULT_INDEXTYPE order[] );
  void level_to_order_exponential ( FLENS_DEFAULT_INDEXTYPE dim_num, FLENS_DEFAULT_INDEXTYPE level[], FLENS_DEFAULT_INDEXTYPE rule[], 
    FLENS_DEFAULT_INDEXTYPE order[] );
  void level_to_order_exponential_slow ( FLENS_DEFAULT_INDEXTYPE dim_num, FLENS_DEFAULT_INDEXTYPE level[], FLENS_DEFAULT_INDEXTYPE rule[], 
    FLENS_DEFAULT_INDEXTYPE order[] );
  void level_to_order_linear ( FLENS_DEFAULT_INDEXTYPE dim_num, FLENS_DEFAULT_INDEXTYPE level[], FLENS_DEFAULT_INDEXTYPE rule[], 
    FLENS_DEFAULT_INDEXTYPE order[] );

  void nc_compute ( FLENS_DEFAULT_INDEXTYPE n, double x_min, double x_max, double x[], double w[] );

  void ncc_compute_points ( FLENS_DEFAULT_INDEXTYPE n, double x[] );
  void ncc_compute_weights ( FLENS_DEFAULT_INDEXTYPE n, double w[] );

  void nco_compute_points ( FLENS_DEFAULT_INDEXTYPE n, double x[] );
  void nco_compute_weights ( FLENS_DEFAULT_INDEXTYPE n, double w[] );

  void patterson_lookup ( FLENS_DEFAULT_INDEXTYPE n, double x[], double w[] );
  void patterson_lookup_points ( FLENS_DEFAULT_INDEXTYPE n, double x[] );
  void patterson_lookup_points_np ( FLENS_DEFAULT_INDEXTYPE n, FLENS_DEFAULT_INDEXTYPE np, double p[], double x[] );
  void patterson_lookup_weights ( FLENS_DEFAULT_INDEXTYPE n, double w[] );
  void patterson_lookup_weights_np ( FLENS_DEFAULT_INDEXTYPE n, FLENS_DEFAULT_INDEXTYPE np, double p[], double w[] );

  FLENS_DEFAULT_INDEXTYPE point_radial_tol_unique_count ( FLENS_DEFAULT_INDEXTYPE m, FLENS_DEFAULT_INDEXTYPE n, double a[], double tol, 
    FLENS_DEFAULT_INDEXTYPE *seed );
  FLENS_DEFAULT_INDEXTYPE point_radial_tol_unique_index ( FLENS_DEFAULT_INDEXTYPE m, FLENS_DEFAULT_INDEXTYPE n, double a[], double tol, 
    FLENS_DEFAULT_INDEXTYPE *seed, FLENS_DEFAULT_INDEXTYPE undx[], FLENS_DEFAULT_INDEXTYPE xdnu[] );

  void point_unique_index ( FLENS_DEFAULT_INDEXTYPE m, FLENS_DEFAULT_INDEXTYPE n, double a[], FLENS_DEFAULT_INDEXTYPE unique_num, FLENS_DEFAULT_INDEXTYPE undx[], 
    FLENS_DEFAULT_INDEXTYPE xdnu[] );

  void product_mixed_weight ( FLENS_DEFAULT_INDEXTYPE dim_num, FLENS_DEFAULT_INDEXTYPE order_1d[], FLENS_DEFAULT_INDEXTYPE order_nd, 
    FLENS_DEFAULT_INDEXTYPE rule[], double alpha[], double beta[], double weight_nd[] );

  double r8_abs ( double x );
  double r8_ceiling ( double x );
  double r8_choose ( FLENS_DEFAULT_INDEXTYPE n, FLENS_DEFAULT_INDEXTYPE k );
  double r8_epsilon ( );
  double r8_factorial ( FLENS_DEFAULT_INDEXTYPE n );
  double r8_factorial2 ( FLENS_DEFAULT_INDEXTYPE n );
  double r8_floor ( double x );
  double r8_gamma ( double x );
  double r8_huge ( );
  double r8_hyper_2f1 ( double a, double b, double c, double x );
  double r8_max ( double x, double y );
  double r8_min ( double x, double y );
  double r8_mop ( FLENS_DEFAULT_INDEXTYPE i );
  double r8_psi ( double xx );

  FLENS_DEFAULT_INDEXTYPE r8col_compare ( FLENS_DEFAULT_INDEXTYPE m, FLENS_DEFAULT_INDEXTYPE n, double a[], FLENS_DEFAULT_INDEXTYPE i, FLENS_DEFAULT_INDEXTYPE j );
  void r8col_sort_heap_a ( FLENS_DEFAULT_INDEXTYPE m, FLENS_DEFAULT_INDEXTYPE n, double a[] );
  FLENS_DEFAULT_INDEXTYPE *r8col_sort_heap_index_a ( FLENS_DEFAULT_INDEXTYPE m, FLENS_DEFAULT_INDEXTYPE n, double a[] );
  FLENS_DEFAULT_INDEXTYPE r8col_sorted_unique_count ( FLENS_DEFAULT_INDEXTYPE m, FLENS_DEFAULT_INDEXTYPE n, double a[], double tol );
  void r8col_swap ( FLENS_DEFAULT_INDEXTYPE m, FLENS_DEFAULT_INDEXTYPE n, double a[], FLENS_DEFAULT_INDEXTYPE j1, FLENS_DEFAULT_INDEXTYPE j2 );
  void r8col_tol_undex ( FLENS_DEFAULT_INDEXTYPE x_dim, FLENS_DEFAULT_INDEXTYPE x_num, double x_val[], FLENS_DEFAULT_INDEXTYPE x_unique_num, 
    double tol, FLENS_DEFAULT_INDEXTYPE undx[], FLENS_DEFAULT_INDEXTYPE xdnu[] );
  FLENS_DEFAULT_INDEXTYPE r8col_tol_unique_count ( FLENS_DEFAULT_INDEXTYPE m, FLENS_DEFAULT_INDEXTYPE n, double a[], double tol );
  void r8col_undex ( FLENS_DEFAULT_INDEXTYPE x_dim, FLENS_DEFAULT_INDEXTYPE x_num, double x_val[], FLENS_DEFAULT_INDEXTYPE x_unique_num, 
    double tol, FLENS_DEFAULT_INDEXTYPE undx[], FLENS_DEFAULT_INDEXTYPE xdnu[] );
  void r8col_unique_index ( FLENS_DEFAULT_INDEXTYPE m, FLENS_DEFAULT_INDEXTYPE n, double a[], double tol, 
    FLENS_DEFAULT_INDEXTYPE unique_index[] );

  void r8mat_transpose_print ( FLENS_DEFAULT_INDEXTYPE m, FLENS_DEFAULT_INDEXTYPE n, double a[], std::string title );
  void r8mat_transpose_print_some ( FLENS_DEFAULT_INDEXTYPE m, FLENS_DEFAULT_INDEXTYPE n, double a[], FLENS_DEFAULT_INDEXTYPE ilo, FLENS_DEFAULT_INDEXTYPE jlo, 
    FLENS_DEFAULT_INDEXTYPE ihi, FLENS_DEFAULT_INDEXTYPE jhi, std::string title );
  void r8mat_write ( std::string output_filename, FLENS_DEFAULT_INDEXTYPE m, FLENS_DEFAULT_INDEXTYPE n, double table[] );

  FLENS_DEFAULT_INDEXTYPE r8vec_compare ( FLENS_DEFAULT_INDEXTYPE n, double a[], double b[] );
  void r8vec_copy ( FLENS_DEFAULT_INDEXTYPE n, double a1[], double a2[] );
  double r8vec_diff_norm_li ( FLENS_DEFAULT_INDEXTYPE n, double a[], double b[] );
  void r8vec_direct_product2 ( FLENS_DEFAULT_INDEXTYPE factor_index, FLENS_DEFAULT_INDEXTYPE factor_order, 
    double factor_value[], FLENS_DEFAULT_INDEXTYPE factor_num, FLENS_DEFAULT_INDEXTYPE point_num, double w[] );
  double r8vec_dot_product ( FLENS_DEFAULT_INDEXTYPE n, double a1[], double a2[] );
  double r8vec_i4vec_dot_product ( FLENS_DEFAULT_INDEXTYPE n, double r8vec[], FLENS_DEFAULT_INDEXTYPE i4vec[] );
  double r8vec_min ( FLENS_DEFAULT_INDEXTYPE n, double r8vec[] );
  double r8vec_min_pos ( FLENS_DEFAULT_INDEXTYPE n, double a[] );
  void r8vec_print ( FLENS_DEFAULT_INDEXTYPE n, double a[], std::string title );
  FLENS_DEFAULT_INDEXTYPE *r8vec_sort_heap_index_a ( FLENS_DEFAULT_INDEXTYPE n, double a[] );
  double r8vec_sum ( FLENS_DEFAULT_INDEXTYPE n, double a[] );
  void r8vec_uniform_01 ( FLENS_DEFAULT_INDEXTYPE n, FLENS_DEFAULT_INDEXTYPE *seed, double r[] );
  double *r8vec_uniform_01_new ( FLENS_DEFAULT_INDEXTYPE n, FLENS_DEFAULT_INDEXTYPE *seed );
  void r8vec_zero ( FLENS_DEFAULT_INDEXTYPE n, double a[] );

  void sort_heap_external ( FLENS_DEFAULT_INDEXTYPE n, FLENS_DEFAULT_INDEXTYPE *indx, FLENS_DEFAULT_INDEXTYPE *i, FLENS_DEFAULT_INDEXTYPE *j, FLENS_DEFAULT_INDEXTYPE isgn );

  void timestamp ( );

  void vec_colex_next3 ( FLENS_DEFAULT_INDEXTYPE dim_num, FLENS_DEFAULT_INDEXTYPE base[], FLENS_DEFAULT_INDEXTYPE a[], bool *more );
}

