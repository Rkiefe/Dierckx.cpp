#pragma once
#include <vector>
#include <cmath>
#include <stdexcept>
#include <algorithm>

// ---------------------------------------------------------------------------
//  Spline2D — Smoothing / least-squares B-spline surface fitting
//  Port of P. Dierckx's SURFIT (FITPACK) to modern C++.
//
//  Usage:
//      Spline2D spl(x, y, z);      // k = 3 (bicubic), s = 0 (interpolate)
//      Spline2D spl(x, y, z, 1);   // bilinear
//      Spline2D spl(x, y, z, 3, s);// bicubic, explicit smoothing factor s>0
//      double zq = spl(xq, yq);
// ---------------------------------------------------------------------------
class Spline2D {
public:
    // k = kx = ky
    Spline2D(const std::vector<double>& x,
             const std::vector<double>& y,
             const std::vector<double>& z,
             int k = 3,
             double s = 0.0,
             const std::vector<double>& w = {});

    // separate kx, ky
    Spline2D(const std::vector<double>& x,
             const std::vector<double>& y,
             const std::vector<double>& z,
             int kx, int ky,
             double s = 0.0,
             const std::vector<double>& w = {});

    double operator()(double xq, double yq) const;
    double eval(double xq, double yq) const { return (*this)(xq, yq); }

    int    degreeX()      const { return kx_; }
    int    degreeY()      const { return ky_; }
    int    numKnotsX()    const { return nx_; }
    int    numKnotsY()    const { return ny_; }
    int    numCoeffs()    const { return ncof_; }
    double residual()     const { return fp_; }
    double smoothing()    const { return s_; }
    int    errorCode()    const { return ier_; }

    const std::vector<double>& knotsX()       const { return tx_; }
    const std::vector<double>& knotsY()       const { return ty_; }
    const std::vector<double>& coefficients() const { return c_; }
    double xb() const { return xb_; }  double xe() const { return xe_; }
    double yb() const { return yb_; }  double ye() const { return ye_; }

private:
    std::vector<double> x_, y_, z_, w_;
    std::vector<double> tx_, ty_, c_;
    int    m_ = 0, kx_ = 3, ky_ = 3, nx_ = 0, ny_ = 0, ncof_ = 0;
    double xb_ = 0, xe_ = 0, yb_ = 0, ye_ = 0;
    double s_ = 0, fp_ = 0;
    int    ier_ = 0;

    static constexpr double EPS_DEFAULT = 1e-16;

    // ---- Kernels (ports of FITPACK) ----
    static void fpbspl(const std::vector<double>& t, int n, int k,
                       double x, int l, double h[6]);
    static void fpgivs(double piv, double& ww, double& cos, double& sin);
    static void fprota(double cos, double sin, double& a, double& b);
    static void fpdisc(const std::vector<double>& t, int n, int k2,
                       std::vector<double>& b, int nest);
    static void fporde(const std::vector<double>& x, const std::vector<double>& y,
                       int m, int kx, int ky,
                       const std::vector<double>& tx, int nx,
                       const std::vector<double>& ty, int ny,
                       std::vector<int>& nummer, std::vector<int>& index, int nreg);
    static void fpback(const std::vector<double>& a, int lda,
                       const std::vector<double>& z, int n, int k,
                       std::vector<double>& c);
    static double fprati(double& p1, double& f1, double p2, double f2,
                         double& p3, double& f3);
    static void fprank(std::vector<double>& a, int lda,
                       std::vector<double>& f, int n, int m, double tol,
                       std::vector<double>& c, double& sq, int& rank,
                       std::vector<double>& aa, std::vector<double>& ff,
                       std::vector<double>& h);

    // Row-major matrix helper: a(i,j) with 1-based i,j.
    static double& AM(std::vector<double>& a, int ncols, int i, int j) {
        return a[(i-1)*ncols + (j-1)];
    }
    static double  AM(const std::vector<double>& a, int ncols, int i, int j) {
        return a[(i-1)*ncols + (j-1)];
    }

    void fpsurf(int iopt);
    void initAndFit(int iopt);
};