#include <Rcpp.h>
#include <vector>
#include <cmath>
#include <limits>
#include <algorithm>

using namespace Rcpp;

namespace {

constexpr double kPi = 3.141592653589793238462643383279502884;
constexpr double kTwoPi = 2.0 * kPi;

inline double wrap_theta(double theta) {
  theta = std::fmod(theta, kTwoPi);
  if (theta < 0.0) {
    theta += kTwoPi;
  }
  return theta;
}

inline void canonicalize_rt(double& rho, double& theta) {
  theta = wrap_theta(theta);
  if (theta > kPi) {
    theta -= kPi;
    rho = -rho;
  }
}

inline NumericMatrix recycle_pair_vectors(const NumericVector& x, const NumericVector& y) {
  int n = std::max(x.size(), y.size());
  if (n == 0) {
    stop("Пустые векторы параметров.");
  }
  if (x.size() != n && x.size() != 1) {
    stop("Длины параметров не согласованы.");
  }
  if (y.size() != n && y.size() != 1) {
    stop("Длины параметров не согласованы.");
  }

  NumericMatrix out(n, 2);
  for (int i = 0; i < n; ++i) {
    out(i, 0) = x[x.size() == 1 ? 0 : i];
    out(i, 1) = y[y.size() == 1 ? 0 : i];
  }
  return out;
}

inline NumericVector pair_err_internal(double tr, double tt, double pr, double pt) {
  double pt_alt = std::fmod(pt + kPi, kTwoPi);
  double dtheta = std::pow(std::atan2(std::sin(tt - pt), std::cos(tt - pt)), 2.0);
  double dtheta_alt = std::pow(std::atan2(std::sin(tt - pt_alt), std::cos(tt - pt_alt)), 2.0);

  double loss_1 = std::pow(tr - pr, 2.0) + dtheta;
  double loss_2 = std::pow(tr + pr, 2.0) + dtheta_alt;

  NumericVector out = NumericVector::create(
    Named("dr") = std::pow(tr - pr, 2.0),
    Named("dtheta") = dtheta,
    Named("total") = loss_1
  );

  if (loss_2 < loss_1) {
    out["dr"] = std::pow(tr + pr, 2.0);
    out["dtheta"] = dtheta_alt;
    out["total"] = loss_2;
  }

  return out;
}

void permutation_search(
  const NumericMatrix& true_lines,
  const NumericMatrix& pred_lines,
  std::vector<int>& current,
  std::vector<int>& used,
  int depth,
  double& best_cost,
  std::vector<int>& best_perm
) {
  int n = true_lines.nrow();
  if (depth == n) {
    double total_cost = 0.0;
    for (int i = 0; i < n; ++i) {
      NumericVector cur = pair_err_internal(
        true_lines(i, 0), true_lines(i, 1),
        pred_lines(current[i], 0), pred_lines(current[i], 1)
      );
      total_cost += as<double>(cur["total"]);
      if (total_cost >= best_cost) {
        return;
      }
    }
    best_cost = total_cost;
    best_perm = current;
    return;
  }

  for (int j = 0; j < n; ++j) {
    if (used[j]) {
      continue;
    }
    used[j] = 1;
    current[depth] = j;
    permutation_search(true_lines, pred_lines, current, used, depth + 1, best_cost, best_perm);
    used[j] = 0;
  }
}

}  // namespace

// [[Rcpp::export]]
List make_accumulator_processed(
  NumericMatrix bound_matrix,
  double rho_step = 1.0,
  double theta_step = 0.01,
  bool weighted = false
) {
  if (rho_step <= 0.0 || theta_step <= 0.0) {
    stop("Шаги rho_step и theta_step должны быть положительными.");
  }

  int n_row = bound_matrix.nrow();
  int n_col = bound_matrix.ncol();

  std::vector<int> rows;
  std::vector<int> cols;
  std::vector<double> weights;
  rows.reserve(n_row * n_col / 10);
  cols.reserve(n_row * n_col / 10);
  weights.reserve(n_row * n_col / 10);

  double active_weight = 0.0;
  for (int i = 0; i < n_row; ++i) {
    for (int j = 0; j < n_col; ++j) {
      double val = bound_matrix(i, j);
      if (!std::isfinite(val) || val <= 0.0) {
        continue;
      }
      rows.push_back(i + 1);
      cols.push_back(j + 1);
      double w = weighted ? val : 1.0;
      weights.push_back(w);
      active_weight += w;
    }
  }

  int active_pixels = static_cast<int>(weights.size());
  int n_theta = static_cast<int>(std::floor(kPi / theta_step)) + 1;
  NumericVector theta(n_theta);
  std::vector<double> cos_theta(n_theta);
  std::vector<double> sin_theta(n_theta);
  for (int j = 0; j < n_theta; ++j) {
    theta[j] = j * theta_step;
    cos_theta[j] = std::cos(theta[j]);
    sin_theta[j] = std::sin(theta[j]);
  }

  int rho_max = static_cast<int>(std::ceil(std::sqrt(
    static_cast<double>(n_row) * static_cast<double>(n_row) +
    static_cast<double>(n_col) * static_cast<double>(n_col)
  )));
  int n_rho = static_cast<int>(std::floor((2.0 * rho_max) / rho_step)) + 1;
  NumericVector rho(n_rho);
  for (int i = 0; i < n_rho; ++i) {
    rho[i] = -rho_max + i * rho_step;
  }

  NumericMatrix accumulator(n_rho, n_theta);
  if (active_pixels == 0) {
    return List::create(
      _["accumulator"] = accumulator,
      _["rho"] = rho,
      _["theta"] = theta,
      _["active_pixels"] = active_pixels,
      _["active_weight"] = active_weight
    );
  }

  double rho_origin = rho[0];
  for (int t = 0; t < n_theta; ++t) {
    for (int idx = 0; idx < active_pixels; ++idx) {
      double rho_value = cols[idx] * cos_theta[t] + rows[idx] * sin_theta[t];
      int rho_idx = static_cast<int>(std::llround((rho_value - rho_origin) / rho_step));
      if (rho_idx >= 0 && rho_idx < n_rho) {
        accumulator(rho_idx, t) += weights[idx];
      }
    }
  }

  return List::create(
    _["accumulator"] = accumulator,
    _["rho"] = rho,
    _["theta"] = theta,
    _["active_pixels"] = active_pixels,
    _["active_weight"] = active_weight
  );
}

// [[Rcpp::export]]
NumericMatrix convert_ab_to_rt_cpp(NumericVector a, NumericVector b) {
  NumericMatrix params = recycle_pair_vectors(a, b);
  int n = params.nrow();
  NumericMatrix out(n, 2);
  colnames(out) = CharacterVector::create("rho", "theta");

  for (int i = 0; i < n; ++i) {
    double cur_a = params(i, 0);
    double cur_b = params(i, 1);
    double theta = std::atan2(1.0, -cur_a);
    double rho = cur_b / std::sqrt(cur_a * cur_a + 1.0);
    canonicalize_rt(rho, theta);
    out(i, 0) = rho;
    out(i, 1) = theta;
  }

  return out;
}

// [[Rcpp::export]]
NumericMatrix convert_rt_to_ab_cpp(NumericVector rho, NumericVector theta, double eps = 1e-12) {
  NumericMatrix params = recycle_pair_vectors(rho, theta);
  int n = params.nrow();
  NumericMatrix out(n, 2);
  colnames(out) = CharacterVector::create("a", "b");

  for (int i = 0; i < n; ++i) {
    double cur_rho = params(i, 0);
    double cur_theta = params(i, 1);
    canonicalize_rt(cur_rho, cur_theta);
    double s = std::sin(cur_theta);
    if (std::fabs(s) <= eps) {
      out(i, 0) = NA_REAL;
      out(i, 1) = NA_REAL;
      continue;
    }
    out(i, 0) = -std::cos(cur_theta) / s;
    out(i, 1) = cur_rho / s;
  }

  return out;
}

// [[Rcpp::export]]
NumericMatrix find_k_max_cpp(
  NumericMatrix acc,
  int k,
  NumericVector qrho,
  NumericVector qtheta,
  bool suppress = false,
  int window = 6
) {
  if (k <= 0) {
    stop("k должно быть положительным.");
  }
  if (acc.nrow() != qrho.size() || acc.ncol() != qtheta.size()) {
    stop("Размеры accumulator, rho и theta не согласованы.");
  }

  NumericMatrix acc_copy = clone(acc);
  NumericMatrix out(k, 2);
  colnames(out) = CharacterVector::create("rho", "theta");

  for (int iter = 0; iter < k; ++iter) {
    double best_val = -std::numeric_limits<double>::infinity();
    int best_row = 0;
    int best_col = 0;

    for (int r = 0; r < acc_copy.nrow(); ++r) {
      for (int c = 0; c < acc_copy.ncol(); ++c) {
        double cur = acc_copy(r, c);
        if (cur > best_val) {
          best_val = cur;
          best_row = r;
          best_col = c;
        }
      }
    }

    out(iter, 0) = qrho[best_row];
    out(iter, 1) = qtheta[best_col];

    if (suppress) {
      int rmin = std::max(0, best_row - window);
      int rmax = std::min(acc_copy.nrow() - 1, best_row + window);
      int cmin = std::max(0, best_col - window);
      int cmax = std::min(acc_copy.ncol() - 1, best_col + window);
      for (int r = rmin; r <= rmax; ++r) {
        for (int c = cmin; c <= cmax; ++c) {
          acc_copy(r, c) = 0.0;
        }
      }
    } else {
      acc_copy(best_row, best_col) = 0.0;
    }
  }

  return out;
}

// [[Rcpp::export]]
NumericVector compute_err_cpp(NumericMatrix true_lines, NumericMatrix pred_lines) {
  if (true_lines.ncol() != 2 || pred_lines.ncol() != 2) {
    stop("Матрицы true_lines и pred_lines должны иметь два столбца: rho и theta.");
  }
  if (true_lines.nrow() != pred_lines.nrow()) {
    stop("true_lines и pred_lines должны содержать одинаковое число прямых.");
  }

  int n = true_lines.nrow();
  NumericVector result = NumericVector::create(
    Named("dr") = 0.0,
    Named("dtheta") = 0.0
  );

  if (n == 1) {
    NumericVector cur = pair_err_internal(
      true_lines(0, 0), true_lines(0, 1),
      pred_lines(0, 0), pred_lines(0, 1)
    );
    result["dr"] = cur["dr"];
    result["dtheta"] = cur["dtheta"];
    return result;
  }

  std::vector<int> current(n, -1);
  std::vector<int> used(n, 0);
  std::vector<int> best_perm(n, -1);
  double best_cost = std::numeric_limits<double>::infinity();
  permutation_search(true_lines, pred_lines, current, used, 0, best_cost, best_perm);

  double dr_total = 0.0;
  double dtheta_total = 0.0;
  for (int i = 0; i < n; ++i) {
    NumericVector cur = pair_err_internal(
      true_lines(i, 0), true_lines(i, 1),
      pred_lines(best_perm[i], 0), pred_lines(best_perm[i], 1)
    );
    dr_total += as<double>(cur["dr"]);
    dtheta_total += as<double>(cur["dtheta"]);
  }

  result["dr"] = dr_total;
  result["dtheta"] = dtheta_total;
  return result;
}

// [[Rcpp::export]]
List ideal_discretization_error_cpp(NumericMatrix true_lines, double rho_step_ht, double theta_step_ht) {
  if (rho_step_ht <= 0.0 || theta_step_ht <= 0.0) {
    stop("Шаги rho_step_ht и theta_step_ht должны быть положительными.");
  }

  NumericMatrix pred(true_lines.nrow(), 2);
  colnames(pred) = CharacterVector::create("rho", "theta");
  for (int i = 0; i < true_lines.nrow(); ++i) {
    double rho = std::round(true_lines(i, 0) / rho_step_ht) * rho_step_ht;
    double theta = std::round(true_lines(i, 1) / theta_step_ht) * theta_step_ht;
    if (theta < 0.0) {
      theta = 0.0;
    }
    if (theta > kPi) {
      theta = kPi;
    }
    pred(i, 0) = rho;
    pred(i, 1) = theta;
  }

  NumericVector err = compute_err_cpp(true_lines, pred);
  double dr = as<double>(err["dr"]);
  double dtheta = as<double>(err["dtheta"]);

  return List::create(
    _["pred"] = pred,
    _["dr"] = dr,
    _["dtheta"] = dtheta,
    _["baseline_dr"] = std::max(dr, std::pow(rho_step_ht / 2.0, 2.0)),
    _["baseline_dtheta"] = std::max(dtheta, std::pow(theta_step_ht / 2.0, 2.0))
  );
}
