#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
IntegerMatrix find_k_maxima(NumericMatrix matr,
                            int k,
                            bool is_suppression = false,
                            int suppression_window = 6) {

  NumericMatrix acc = clone(matr);
  int nrow = acc.nrow(), ncol = acc.ncol();
  IntegerMatrix result(k, 2);

  for (int m = 0; m < k; m++) {
    double maxVal = R_NegInf;
    int maxRow = -1, maxCol = -1;

    
    for (int i = 0; i < nrow; i++)
      for (int j = 0; j < ncol; j++)
        if (acc(i, j) > maxVal) {
          maxVal = acc(i, j);
          maxRow = i;
          maxCol = j;
        }

    if (maxRow < 0) { 
      for (int c = m; c < k; c++) {
        result(c, 0) = NA_INTEGER;
        result(c, 1) = NA_INTEGER;
      }
      break;
    }

    result(m, 0) = maxRow + 1;
    result(m, 1) = maxCol + 1;

    if (is_suppression) {
      int rmin = std::max(0, maxRow - suppression_window);
      int rmax = std::min(nrow - 1, maxRow + suppression_window);
      int cmin = std::max(0, maxCol - suppression_window);
      int cmax = std::min(ncol - 1, maxCol + suppression_window);
      for (int i = rmin; i <= rmax; i++)
        for (int j = cmin; j <= cmax; j++)
          acc(i, j) = R_NegInf;
    } else {
      acc(maxRow, maxCol) = R_NegInf;
    }
  }
  return result;
}

// [[Rcpp::export]]
NumericVector compute_hough_error(NumericMatrix accumulator,
                                  NumericVector quant_rho,
                                  NumericVector quant_theta,
                                  NumericMatrix true_lines_rt,
                                  int k = 2,
                                  bool is_suppression = false,
                                  int suppression_window = 6) {

  IntegerMatrix maxima = find_k_maxima(accumulator, k,
                                       is_suppression, suppression_window);

  double sum_rho = 0.0, sum_theta = 0.0;
  int valid_found = 0;

  for (int i = 0; i < maxima.nrow(); i++) {
    if (IntegerVector::is_na(maxima(i, 0))) continue;

    double rho_est   = quant_rho[maxima(i, 0) - 1];
    double theta_est = quant_theta[maxima(i, 1) - 1];

    double best_r = R_PosInf, best_t = R_PosInf, best_dist = R_PosInf;

    for (int j = 0; j < true_lines_rt.nrow(); j++) {
      double dr = std::abs(rho_est - true_lines_rt(j, 0));
      double dt = std::abs(theta_est - true_lines_rt(j, 1));
      double d  = std::sqrt(dr * dr + dt * dt);
      if (d < best_dist) {
        best_dist = d;
        best_r = dr;
        best_t = dt;
      }
    }

    if (R_finite(best_dist)) {
      sum_rho   += best_r;
      sum_theta += best_t;
      valid_found++;
    }
  }

  if (valid_found == 0)
    return NumericVector::create(NA_REAL, NA_REAL);

  return NumericVector::create(
    Named("rho_error")   = sum_rho   / valid_found,
    Named("theta_error") = sum_theta / valid_found
  );
}