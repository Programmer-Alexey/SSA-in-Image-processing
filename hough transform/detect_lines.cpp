#include <Rcpp.h>
#include <algorithm>
using namespace Rcpp;


void add_line_to_matrix(NumericMatrix& mat, double a, double b) {
    int nrow = mat.nrow();
    int ncol = mat.ncol();

    for (int x = 0; x < ncol; ++x) {
        double y = a * x + b;
        int yi = (int)std::round(y);
        if (yi >= 0 && yi < nrow) {
            mat(yi, x) = 1.0;
        }
    }
}

// [[Rcpp::export]]
NumericMatrix detect_lines(IntegerMatrix accumulator,
    NumericVector quant_rho,
    NumericVector quant_theta,
    int N,
    int ncol = 100,
    int nrow = 100,
    bool nms = false,
    int suppression_window = 6) {

    double eps = 1.0 / std::max(ncol, nrow);
    int n_rho = accumulator.nrow();
    int n_theta = accumulator.ncol();

    NumericMatrix output(nrow, ncol);

    if (!nms) {
        int total = n_rho * n_theta;
        std::vector<int> indices(total);
        std::iota(indices.begin(), indices.end(), 0); // 0,1,...,total-1

        std::sort(indices.begin(), indices.end(), [&](int a, int b) {
            int ra = a % n_rho;
            int ca = a / n_rho;
            int rb = b % n_rho;
            int cb = b / n_rho;
            return accumulator(ra, ca) > accumulator(rb, cb);
            });

        for (int k = 0; k < N && k < total; ++k) {
            int idx = indices[k];
            int row = idx % n_rho;
            int col = idx / n_rho;

            double rho = quant_rho[row];
            double theta = quant_theta[col];

            if (std::fabs(std::sin(theta)) > eps) {
                double a = -std::tan(theta);
                double b = rho / std::cos(theta);
                add_line_to_matrix(output, a, b);
            }
        }

    }
    else {
        IntegerMatrix suppressed = clone(accumulator);

        for (int k = 0; k < N; ++k) {
            int max_val = 0;
            int max_row = -1, max_col = -1;
            for (int r = 0; r < n_rho; ++r) {
                for (int c = 0; c < n_theta; ++c) {
                    if (suppressed(r, c) > max_val) {
                        max_val = suppressed(r, c);
                        max_row = r;
                        max_col = c;
                    }
                }
            }
            if (max_val == 0) break;

            double rho = quant_rho[max_row];
            double theta = quant_theta[max_col];

            if (std::fabs(std::sin(theta)) > eps) {
                double a = -std::tan(theta);
                double b = rho / std::cos(theta);
                add_line_to_matrix(output, a, b);
            }

            int rmin = std::max(0, max_row - suppression_window);
            int rmax = std::min(n_rho - 1, max_row + suppression_window);
            int cmin = std::max(0, max_col - suppression_window);
            int cmax = std::min(n_theta - 1, max_col + suppression_window);

            for (int r = rmin; r <= rmax; ++r) {
                for (int c = cmin; c <= cmax; ++c) {
                    suppressed(r, c) = 0;
                }
            }
        }
    }

    return output;
}