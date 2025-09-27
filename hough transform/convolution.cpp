#include <Rcpp.h>
#include <vector>
#include <cmath>
#include <algorithm>

using namespace Rcpp;

// --------------------
// Convolution
// [[Rcpp::export]]
NumericMatrix convolution(NumericMatrix matr, NumericMatrix kernel) {
    int matr_nrow = matr.nrow(), matr_ncol = matr.ncol();
    int kernel_nrow = kernel.nrow(), kernel_ncol = kernel.ncol();
    int pad_y = kernel_nrow / 2, pad_x = kernel_ncol / 2;

    NumericMatrix padded(matr_nrow + 2 * pad_y, matr_ncol + 2 * pad_x);
    for (int i = 0; i < matr_nrow; i++)
        for (int j = 0; j < matr_ncol; j++)
            padded(i + pad_y, j + pad_x) = matr(i, j);

    NumericMatrix result(matr_nrow, matr_ncol);
    for (int i = 0; i < matr_nrow; i++) {
        for (int j = 0; j < matr_ncol; j++) {
            double s = 0.0;
            for (int ki = 0; ki < kernel_nrow; ki++)
                for (int kj = 0; kj < kernel_ncol; kj++)
                    s += padded(i + ki, j + kj) * kernel(ki, kj);
            result(i, j) = s;
        }
    }
    return result;
}

// --------------------
// Intensity Detector
// [[Rcpp::export]]
NumericMatrix intensity_detector(NumericMatrix matr, double threshold = 0.8) {
    NumericVector vals = as<NumericVector>(matr);
    std::vector<double> v = Rcpp::as<std::vector<double>>(vals);
    std::sort(v.begin(), v.end());
    int idx = std::floor(threshold * (v.size() - 1));
    double q = v[idx];

    NumericMatrix res(matr.nrow(), matr.ncol());
    for (int i = 0; i < matr.nrow(); i++)
        for (int j = 0; j < matr.ncol(); j++)
            res(i, j) = matr(i, j) > q ? matr(i, j) : 0;
    return res;
}

// --------------------
// Gradient Detector
// [[Rcpp::export]]
NumericMatrix intensity_gradient_detector(NumericMatrix matr, double threshold = 2.0) {
    NumericMatrix sobel_x(3, 3);
    sobel_x(0, 0) = -1; sobel_x(0, 1) = 0; sobel_x(0, 2) = 1;
    sobel_x(1, 0) = -2; sobel_x(1, 1) = 0; sobel_x(1, 2) = 2;
    sobel_x(2, 0) = -1; sobel_x(2, 1) = 0; sobel_x(2, 2) = 1;

    NumericMatrix Gx = convolution(matr, sobel_x);
    NumericMatrix Gy = convolution(matr, transpose(sobel_x));

    NumericMatrix res(matr.nrow(), matr.ncol());
    for (int i = 0; i < matr.nrow(); i++)
        for (int j = 0; j < matr.ncol(); j++) {
            double val = std::sqrt(Gx(i, j) * Gx(i, j) + Gy(i, j) * Gy(i, j));
            res(i, j) = val > threshold ? val : 0;
        }
    return res;
}

// --------------------
// Gaussian Kernel
// [[Rcpp::export]]
NumericMatrix gaussian_kernel(int n, double sigma) {
    NumericMatrix kernel(n, n);
    int center = n / 2;
    double factor = 1.0 / (2 * M_PI * sigma * sigma);
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++) {
            double x = i - center;
            double y = j - center;
            kernel(i, j) = std::exp(-(x * x + y * y) / (2 * sigma * sigma)) * factor;
        }
    return kernel;
}

// --------------------
// Laplace Detector
// [[Rcpp::export]]
NumericMatrix laplace_detector(NumericMatrix matr, double sigma = 2.0) {
    int n = (int)(6 * sigma) - 1;
    NumericMatrix gauss = gaussian_kernel(n, sigma);

    NumericMatrix laplace(3, 3);
    laplace(0, 0) = 1; laplace(0, 1) = 4; laplace(0, 2) = 1;
    laplace(1, 0) = 4; laplace(1, 1) = -20; laplace(1, 2) = 4;
    laplace(2, 0) = 1; laplace(2, 1) = 4; laplace(2, 2) = 1;

    NumericMatrix temp = convolution(matr, gauss);
    NumericMatrix res = convolution(temp, laplace);

    int nrow = res.nrow(), ncol = res.ncol();
    NumericMatrix edges(nrow, ncol);
    for (int i = 1; i < nrow - 1; i++)
        for (int j = 1; j < ncol - 1; j++) {
            int neg = 0, pos = 0;
            for (int ii = i - 1; ii <= i + 1; ii++)
                for (int jj = j - 1; jj <= j + 1; jj++) {
                    if (res(ii, jj) < 0) neg++;
                    else if (res(ii, jj) > 0) pos++;
                }
            edges(i, j) = (neg > 0 && pos > 0) ? 1 : 0;
        }
    return edges;
}
