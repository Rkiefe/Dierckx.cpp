#pragma once
#include <vector>
#include <cmath>
#include <stdexcept>
#include <algorithm>
#include <numeric>

// ---------------------------------------------------------------------------
//  Spline1D — Smoothing / least-squares B-spline curve fitting
//
//  Port of P. Dierckx's CURFIT (FITPACK) routines to modern C++.
//
//  Usage:
//      Spline1D spl(x, y);              // default k=3, s = smoothing (auto)
//      Spline1D spl(x, y, 1);           // k=1 linear
//      Spline1D spl(x, y, 3, s);        // cubic, explicit smoothing factor
//      double yq = spl(xq);             // evaluate
// ---------------------------------------------------------------------------
class Spline1D {
public:
    // ---- Construction -----------------------------------------------------

    // Fit with automatic smoothing choice.  w defaults to 1.
    // iopt = 0 (restart).  s < 0 triggers an automatic s search.
    explicit Spline1D(const std::vector<double>& x,
                      const std::vector<double>& y,
                      int k = 3,
                      double s = 0.0,
                      const std::vector<double>& w = {});

    // Convenience: iterators / any container with .size() and operator[]
    template <class XC, class YC>
    Spline1D(const XC& x, const YC& y, int k = 3)
        : Spline1D(toVec(x), toVec(y), k) {}

    // ---- Evaluation -------------------------------------------------------

    double operator()(double xq) const;          // evaluate spline
    double eval(double xq) const { return (*this)(xq); }

    // Derivative (order = 0 returns value)
    double derivative(double xq, int order = 1) const;

    // ---- Introspection ----------------------------------------------------

    int    degree()       const { return k_; }
    int    numKnots()     const { return n_; }
    int    numCoeffs()    const { return nc_; }
    double residual()     const { return fp_; }
    double smoothing()    const { return s_; }
    int    errorCode()    const { return ier_; }

    const std::vector<double>& knots()       const { return t_; }
    const std::vector<double>& coefficients() const { return c_; }
    double xb() const { return xb_; }
    double xe() const { return xe_; }

private:
    // ---- Data -------------------------------------------------------------
    std::vector<double> x_, y_, w_;
    std::vector<double> t_;   // knots
    std::vector<double> c_;   // B-spline coefficients
    int    m_ = 0, k_ = 3, n_ = 0, nc_ = 0;
    double xb_ = 0.0, xe_ = 0.0, s_ = 0.0, fp_ = 0.0;
    int    ier_ = 0;

    // ---- Static helpers ---------------------------------------------------
    template <class C>
    static std::vector<double> toVec(const C& c) {
        return std::vector<double>(c.begin(), c.end());
    }

    // ---- Internal kernels (port of FITPACK) --------------------------------
    static void fpbspl(const std::vector<double>& t, int n, int k,
                       double x, int l, double h[6]);
    static void fpgivs(double piv, double& ww, double& cos, double& sin);
    static void fprota(double cos, double sin, double& a, double& b);
    static void fpback(const std::vector<double>& a, const std::vector<double>& z,
                       int n, int k, std::vector<double>& c);
    static void fpdisc(const std::vector<double>& t, int n, int k2,
                       std::vector<double>& b, int nest);
    static void fpknot(const std::vector<double>& x, int m,
                       std::vector<double>& t, int& n,
                       std::vector<double>& fpint, std::vector<int>& nrdata,
                       int& nrint, int nest, int istart);
    static double fprati(double p1, double f1, double p2, double f2,
                         double p3, double f3);

    void fpcurf(int iopt, int maxit, double tol);
    void fpchec(int& ier) const;

    // Matrix layout: a is (nest x k1) row-major flattened
    static double& A(std::vector<double>& a, int nest, int k1, int i, int j) {
        return a[(i - 1) * k1 + (j - 1)];
    }
    static double A(const std::vector<double>& a, int nest, int k1, int i, int j) {
        return a[(i - 1) * k1 + (j - 1)];
    }

    // // Matrix layout: a is (nest x k1) row-major flattened
    // static double& A(std::vector<double>& a, int /*nest*/, int k1, int i, int j) {
    //     return a[(i - 1) * k1 + (j - 1)];
    // }
    // static double  A(const std::vector<double>& a, int /*nest*/, int k1, int i, int j) {
    //     return a[(i - 1) * k1 + (j - 1)];
    // }
};