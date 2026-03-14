#include <Rcpp.h>
#include <algorithm>
#include <numeric>
#include <cmath>
using namespace Rcpp;


void add_line_to_matrix(NumericMatrix& mat, double a, double b) {
    int nrow = mat.nrow();
    int ncol = mat.ncol();

    // Bresenham rasterization to match add.line() behavior in R.
    int x0 = 1;
    int x1 = ncol;
    int y0 = (int)std::round(a * x0 + b);
    int y1 = (int)std::round(a * x1 + b);

    int dx = std::abs(x1 - x0);
    int dy = std::abs(y1 - y0);
    int sx = (x0 < x1) ? 1 : -1;
    int sy = (y0 < y1) ? 1 : -1;
    int err = dx - dy;

    while (true) {
        if (x0 >= 1 && x0 <= ncol && y0 >= 1 && y0 <= nrow) {
            mat(y0 - 1, x0 - 1) = 1.0;
        }
        if (x0 == x1 && y0 == y1) break;
        int e2 = 2 * err;
        if (e2 > -dy) {
            err -= dy;
            x0 += sx;
        }
        if (e2 < dx) {
            err += dx;
            y0 += sy;
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
                double a = -1/std::tan(theta);
                double b = rho / std::sin(theta);
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
                double a = -1/std::tan(theta);
                double b = rho / std::sin(theta);
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
