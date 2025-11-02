#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List make_accumulator(NumericMatrix matrix,
    Function detector,
    double theta_step = 1.0,
    double rho_step = 1.0) {

    NumericMatrix bound_matrix = detector(matrix);

    int n_row = bound_matrix.nrow();
    int n_col = bound_matrix.ncol();

    std::vector<int> non_zero_row, non_zero_col;

    for (int i = 0; i < n_row; ++i) {
        for (int j = 0; j < n_col; ++j) {
            if (bound_matrix(i, j) > 0) {
                non_zero_row.push_back(i + 1);
                non_zero_col.push_back(j + 1);
            }
        }
    }

    // theta
    int n_theta = std::floor(M_PI / theta_step);
    NumericVector theta(n_theta);
    for (int j = 0; j < n_theta; j++) {
        theta[j] = (j * M_PI) / (n_theta - 1);
    }

    // rho
    int rho_max = std::ceil(std::sqrt(n_row * n_row + n_col * n_col));
    int n_rho = std::floor((2 * rho_max) / rho_step) + 1;
    NumericVector rho(n_rho);
    for (int i = 0; i < n_rho; ++i) {
        rho[i] = -rho_max + i * rho_step;
    }

    // accumulator
    IntegerMatrix accumulator(n_rho, n_theta);

    for (size_t i = 0; i < non_zero_row.size(); ++i) {
        int xi = non_zero_row[i];
        int yi = non_zero_col[i];

        for (int j = 0; j < n_theta; ++j) {
            double rho_val = xi * std::cos(theta[j]) + yi * std::sin(theta[j]);
            int cur_rho = std::round((rho_val - rho[0]) / rho_step);

            if (cur_rho >= 0 && cur_rho < n_rho) {
                accumulator(cur_rho, j) += 1;
            }
        }
    }

    return List::create(
        _["accumulator"] = accumulator,
        _["rho"] = rho,
        _["theta"] = theta
    );
}
