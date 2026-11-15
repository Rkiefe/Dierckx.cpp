// Header-Only version of 'Dierckx.cpp'

#pragma once
#include <vector>
#include <cmath>
#include <stdexcept>
#include <algorithm>
#include <numeric>
#include <limits>
#include <cassert>

// ---------------------------------------------------------------------------
//  Spline1D — Smoothing / least-squares B-spline curve fitting
//
//  Port of P. Dierckx's CURFIT (FITPACK) routines to modern C++.
//  Header-only: just #include "spline1d.hpp".
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

    explicit Spline1D(const std::vector<double>& x,
                      const std::vector<double>& y,
                      int k = 3,
                      double s = 0.0,
                      const std::vector<double>& w = {});

    // Convenience: any container with .begin() / .end()
    template <class XC, class YC>
    Spline1D(const XC& x, const YC& y, int k = 3)
        : Spline1D(toVec(x), toVec(y), k) {}

    // ---- Evaluation -------------------------------------------------------

    double operator()(double xq) const;
    double eval(double xq) const { return (*this)(xq); }
    double derivative(double xq, int order = 1) const;

    // ---- Introspection ----------------------------------------------------

    int    degree()       const { return k_; }
    int    numKnots()     const { return n_; }
    int    numCoeffs()    const { return nc_; }
    double residual()     const { return fp_; }
    double smoothing()    const { return s_; }
    int    errorCode()    const { return ier_; }

    const std::vector<double>& knots()        const { return t_; }
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

    // ---- Internal kernels (port of FITPACK) -------------------------------
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
    static double fprati(double &p1, double &f1, double p2, double f2,
                         double &p3, double &f3);

    void fpcurf(int iopt, int maxit, double tol);
    void fpchec(int& ier) const;

    // Matrix layout: a is (nest x k1) row-major flattened
    static double& A(std::vector<double>& a, int nest, int k1, int i, int j) {
        return a[(i - 1) * k1 + (j - 1)];
    }
    static double A(const std::vector<double>& a, int nest, int k1, int i, int j) {
        return a[(i - 1) * k1 + (j - 1)];
    }
};

// ============================================================================
//  Implementation (all inline)
// ============================================================================

// ---- Low-level FITPACK kernels ---------------------------------------------

inline void Spline1D::fpbspl(const std::vector<double>& t, int /*n*/, int k,
                             double x, int l, double h[6])
{
    const double one = 1.0;
    double hh[5];

    h[0] = one;                       // h(1)
    for (int j = 1; j <= k; ++j) {
        for (int i = 1; i <= j; ++i)
            hh[i - 1] = h[i - 1];
        h[0] = 0.0;
        for (int i = 1; i <= j; ++i) {
            int li = l + i;
            int lj = li - j;
            double f = hh[i - 1] / (t[li - 1] - t[lj - 1]);
            h[i - 1]     = h[i - 1] + f * (t[li - 1] - x);
            h[i]         = f * (x - t[lj - 1]);
        }
    }
}

inline void Spline1D::fpgivs(double piv, double& ww, double& cos, double& sin)
{
    const double one = 1.0;
    double store = std::abs(piv);

    if (store == 0.0 && ww == 0.0) { cos = 1.0; sin = 0.0; return; }

    double dd;
    if (store >= ww) dd = store * std::sqrt(one + (ww / piv) * (ww / piv));
    else             dd = ww    * std::sqrt(one + (piv / ww) * (piv / ww));
    cos = ww / dd;
    sin = piv / dd;
    ww  = dd;
}

inline void Spline1D::fprota(double cos, double sin, double& a, double& b)
{
    double stor1 = a, stor2 = b;
    b = cos * stor2 + sin * stor1;
    a = cos * stor1 - sin * stor2;
}

inline void Spline1D::fpback(const std::vector<double>& a, const std::vector<double>& z,
                             int n, int k, std::vector<double>& c)
{
    int k1 = k - 1;
    c[n - 1] = z[n - 1] / A(a, 0, k, n, 1);
    int i = n - 1;
    if (i == 0) return;
    for (int j = 2; j <= n; ++j) {
        double store = z[i - 1];
        int i1 = k1;
        if (j <= k1) i1 = j - 1;
        int mm = i;
        for (int l = 1; l <= i1; ++l) {
            mm = mm + 1;
            store = store - c[mm - 1] * A(a, 0, k, i, l + 1);
        }
        c[i - 1] = store / A(a, 0, k, i, 1);
        i = i - 1;
    }
}

inline void Spline1D::fpdisc(const std::vector<double>& t, int n, int k2,
                             std::vector<double>& b, int nest)
{
    int k1   = k2 - 1;
    int k    = k1 - 1;
    int nk1  = n - k1;
    int nrint = nk1 - k;
    double an  = static_cast<double>(nrint);
    double fac = an / (t[nk1] - t[k1 - 1]);
    double h[12];

    auto T = [&](int idx1) -> double { return t[idx1 - 1]; };

    for (int l = k2; l <= nk1; ++l) {
        int lmk = l - k1;
        for (int j = 1; j <= k1; ++j) {
            int ik = j + k1;
            int lj = l + j;
            int lk = lj - k2;
            h[j - 1]  = T(l) - T(lk);
            h[ik - 1] = T(l) - T(lj);
        }
        int lp = lmk;
        for (int j = 1; j <= k2; ++j) {
            int jk = j;
            double prod = h[j - 1];
            for (int i = 1; i <= k; ++i) {
                jk = jk + 1;
                prod = prod * h[jk - 1] * fac;
            }
            int lk = lp + k1;
            b[(lmk - 1) * k2 + (j - 1)] = (T(lk) - T(lp)) / prod;
            lp = lp + 1;
        }
    }
}

inline void Spline1D::fpknot(const std::vector<double>& x, int /*m*/,
                             std::vector<double>& t, int& n,
                             std::vector<double>& fpint, std::vector<int>& nrdata,
                             int& nrint, int /*nest*/, int istart)
{
    int k = (n - nrint - 1) / 2;

    double fpmax = 0.0;
    int jbegin = istart, number = 0, maxpt = 0, maxbeg = 0;
    for (int j = 1; j <= nrint; ++j) {
        int jpoint = nrdata[j - 1];
        if (fpmax >= fpint[j - 1] || jpoint == 0) { jbegin += jpoint + 1; continue; }
        fpmax  = fpint[j - 1];
        number = j;
        maxpt  = jpoint;
        maxbeg = jbegin;
        jbegin += jpoint + 1;
    }

    int ihalf = maxpt / 2 + 1;
    int nrx   = maxbeg + ihalf;
    int next  = number + 1;
    if (next <= nrint) {
        for (int j = next; j <= nrint; ++j) {
            int jj = next + nrint - j;
            fpint[jj]      = fpint[jj - 1];
            nrdata[jj]     = nrdata[jj - 1];
            int jk = jj + k;
            t[jk] = t[jk - 1];
        }
    }
    nrdata[number - 1] = ihalf - 1;
    nrdata[next - 1]   = maxpt - ihalf;
    double am = static_cast<double>(maxpt);
    double an = static_cast<double>(nrdata[number - 1]);
    fpint[number - 1] = fpmax * an / am;
    an = static_cast<double>(nrdata[next - 1]);
    fpint[next - 1]  = fpmax * an / am;
    int jk = next + k;
    t[jk - 1] = x[nrx - 1];
    ++n;
    ++nrint;
}

inline double Spline1D::fprati(double& p1, double& f1,
                               double  p2, double  f2,
                               double& p3, double& f3)
{
    double p;
    if (p3 > 0.0) {
        double h1 = f1 * (f2 - f3);
        double h2 = f2 * (f3 - f1);
        double h3 = f3 * (f1 - f2);
        p = -(p1 * p2 * h3 + p2 * p3 * h1 + p3 * p1 * h2) /
             (p1 * h1 + p2 * h2 + p3 * h3);
    } else {
        p = (p1 * (f1 - f3) * f2 - p2 * (f2 - f3) * f1) / ((f1 - f2) * f3);
    }

    if (f2 >= 0.0) { p1 = p2; f1 = f2; }
    else           { p3 = p2; f3 = f2; }

    return p;
}

// ---- Data validity check (port of fpchec) ----------------------------------

inline void Spline1D::fpchec(int& ier) const
{
    int m = m_, n = n_, k = k_;
    int k1 = k + 1;
    int k2 = k1 + 1;
    int nk1 = n - k1;
    int nk2 = nk1 + 1;
    ier = 10;

    if (nk1 < k1 || nk1 > m) return;

    {
        int j = n;
        for (int i = 1; i <= k; ++i) {
            if (t_[i - 1] > t_[i]) return;
            if (t_[j - 1] < t_[j - 2]) return;
            --j;
        }
    }

    for (int i = k2; i <= nk2; ++i)
        if (t_[i - 1] <= t_[i - 2]) return;

    if (x_[0] < t_[k1 - 1] || x_[m - 1] > t_[nk2 - 1]) return;

    if (x_[0] >= t_[k2 - 1] || x_[m - 1] <= t_[nk1 - 1]) return;
    {
        int i = 1, l = k2;
        int nk3 = nk1 - 1;
        if (nk3 < 2) { ier = 0; return; }
        for (int j = 2; j <= nk3; ++j) {
            double tj = t_[j - 1];
            ++l;
            double tl = t_[l - 1];
            do {
                ++i;
                if (i >= m) return;
                if (x_[i - 1] > tj) {
                    if (x_[i - 1] >= tl) return;
                    break;
                }
            } while (true);
        }
    }
    ier = 0;
}

// ---- Main fitting driver (port of fpcurf) ----------------------------------

inline void Spline1D::fpcurf(int iopt, int maxit, double tol)
{
    const double one  = 1.0;
    const double con1 = 0.1;
    const double con9 = 0.9;
    const double con4 = 0.04;
    const double half = 0.5;

    int k1 = k_ + 1;
    int k2 = k1 + 1;
    int m  = m_;
    int nest = static_cast<int>(t_.size());

    int n = n_;
    struct Writeback {
        int& dst_n;  int& dst_nc;
        const int& src_n;  int k1;
        ~Writeback() {
            dst_n  = src_n;
            dst_nc = (src_n > k1) ? (src_n - k1) : 0;
        }
    } _wb{ n_, nc_, n, k1 };

    std::vector<double> fpint(nest, 0.0);
    std::vector<double> z(nest, 0.0);
    std::vector<double> a((size_t)nest * k1, 0.0);
    std::vector<double> b((size_t)nest * k2, 0.0);
    std::vector<double> g((size_t)nest * k2, 0.0);
    std::vector<double> q((size_t)m * k1, 0.0);
    std::vector<int>    nrdata(nest, 0);

    double h[7];

    int nmin = 2 * k1;
    int nmax = m + k1, nrint = 0, nk1 = 0, nplus = 0;
    double fp = 0.0, fp0 = 0.0, fpold = 0.0, fpms = 0.0;
    double acc = 0.0;

    if (iopt >= 0) {
        acc = tol * s_;

        if (s_ > 0.0) {
            if (iopt == 0) {
                n = nmin;
                fpold = 0.0;
                nplus = 0;
                nrdata[0] = m - 2;
            } else {
                if (n == nmin) {
                    fpold = 0.0;
                    nplus = 0;
                    nrdata[0] = m - 2;
                } else {
                    fp0   = fpint[n - 1];
                    fpold = fpint[n - 2];
                    nplus = nrdata[n - 1];
                    if (fp0 <= s_) {
                        n = nmin;
                        fpold = 0.0;
                        nplus = 0;
                        nrdata[0] = m - 2;
                    }
                }
            }
        } else {
            n = nmax;
            if (nmax > nest) { ier_ = 1; return; }
            int mk1 = m - k1;
            if (mk1 == 0) {
                n = nmin;
            } else {
                int k3 = k_ / 2;
                int i = k2;
                int j = k3 + 2;
                if (k3 * 2 == k_) {
                    for (int l = 1; l <= mk1; ++l) {
                        t_[i - 1] = (x_[j - 1] + x_[j - 2]) * half;
                        ++i; ++j;
                    }
                } else {
                    for (int l = 1; l <= mk1; ++l) {
                        t_[i - 1] = x_[j - 1];
                        ++i; ++j;
                    }
                }
            }
        }
    }

    for (int iter = 1; iter <= m; ++iter) {
        if (n == nmin) ier_ = -2;
        nrint = n - nmin + 1;

        nk1 = n - k1;
        {
            int i = n;
            for (int j = 1; j <= k1; ++j) {
                t_[j - 1] = xb_;
                t_[i - 1] = xe_;
                --i;
            }
        }

        fp = 0.0;
        for (int i = 1; i <= nk1; ++i) {
            z[i - 1] = 0.0;
            for (int j = 1; j <= k1; ++j)
                A(a, nest, k1, i, j) = 0.0;
        }

        int l = k1;
        for (int it = 1; it <= m; ++it) {
            double xi = x_[it - 1];
            double wi = w_[it - 1];
            double yi = y_[it - 1] * wi;

            while (!(xi < t_[l] || l == nk1)) ++l;

            double hloc[6];
            fpbspl(t_, n, k_, xi, l, hloc);
            for (int i = 1; i <= k1; ++i) {
                q[(size_t)(it - 1) * k1 + (i - 1)] = hloc[i - 1];
                hloc[i - 1] *= wi;
            }

            int j = l - k1;
            for (int i = 1; i <= k1; ++i) {
                ++j;
                double piv = hloc[i - 1];
                if (piv == 0.0) continue;

                double cos, sin;
                fpgivs(piv, A(a, nest, k1, j, 1), cos, sin);
                fprota(cos, sin, yi, z[j - 1]);
                if (i == k1) continue;

                int i2 = 1;
                for (int i1 = i + 1; i1 <= k1; ++i1) {
                    ++i2;
                    fprota(cos, sin, hloc[i1 - 1], A(a, nest, k1, j, i2));
                }
            }
            fp += yi * yi;
        }

        if (ier_ == -2) fp0 = fp;
        fpint[n - 1]   = fp0;
        if (n - 1 >= 1) fpint[n - 2] = fpold;
        nrdata[n - 1]  = nplus;

        c_.assign(nest, 0.0);
        fpback(a, z, nk1, k1, c_);

        if (iopt < 0) { fp_ = fp; return; }

        fpms = fp - s_;
        if (std::abs(fpms) < acc) { fp_ = fp; return; }
        if (fpms < 0.0)           { fp_ = fp; break; }
        if (n == nmax)            { ier_ = -1; fp_ = fp; return; }
        if (n == nest)            { ier_ = 1;  fp_ = fp; return; }

        if (ier_ == 0) {
            int npl1 = nplus * 2;
            double rn = static_cast<double>(nplus);
            if (fpold - fp > acc)
                npl1 = static_cast<int>(rn * fpms / (fpold - fp));
            nplus = std::min(nplus * 2, std::max({ npl1, nplus / 2, 1 }));
        } else {
            nplus = 1;
            ier_ = 0;
        }
        fpold = fp;

        {
            double fpart = 0.0;
            int i = 1;
            l = k2;
            int isnew = 0;
            for (int it = 1; it <= m; ++it) {
                if (x_[it - 1] >= t_[l - 1] && l <= nk1) {
                    isnew = 1;
                    ++l;
                }
                double term = 0.0;
                int l0 = l - k2;
                for (int j = 1; j <= k1; ++j) {
                    ++l0;
                    term += c_[l0 - 1] * q[(size_t)(it - 1) * k1 + (j - 1)];
                }
                term = std::pow(w_[it - 1] * (term - y_[it - 1]), 2);
                fpart += term;
                if (isnew) {
                    double store = term * half;
                    fpint[i - 1] = fpart - store;
                    ++i;
                    fpart = store;
                    isnew = 0;
                }
            }
            fpint[nrint - 1] = fpart;
        }

        for (int li = 1; li <= nplus; ++li) {
            fpknot(x_, m, t_, n, fpint, nrdata, nrint, nest, 1);
            if (n == nmax) {
                int mk1 = m - k1;
                int k3 = k_ / 2;
                int i = k2;
                int j = k3 + 2;
                if (k3 * 2 == k_) {
                    for (int l2 = 1; l2 <= mk1; ++l2) {
                        t_[i - 1] = (x_[j - 1] + x_[j - 2]) * half;
                        ++i; ++j;
                    }
                } else {
                    for (int l2 = 1; l2 <= mk1; ++l2) {
                        t_[i - 1] = x_[j - 1];
                        ++i; ++j;
                    }
                }
                break;
            }
            if (n == nest) break;
        }
        if (n == nest && fpms >= 0.0) { ier_ = 1; fp_ = fp; return; }
    }

    if (ier_ == -2) { fp_ = fp; return; }

    std::vector<double> bmat((size_t)nest * k2, 0.0);
    fpdisc(t_, n, k2, bmat, nest);

    double p1 = 0.0, f1 = fp0 - s_;
    double p3 = -one,  f3 = fpms;
    double p  = 0.0;
    nk1 = n - k1;
    for (int i = 1; i <= nk1; ++i) p += A(a, nest, k1, i, 1);
    p = static_cast<double>(nk1) / p;

    int ich1 = 0, ich3 = 0;
    int n8 = n - nmin;

    for (int iter = 1; iter <= maxit; ++iter) {
        double pinv = one / p;

        c_.assign(nest, 0.0);
        for (int i = 1; i <= nk1; ++i) {
            c_[i - 1] = z[i - 1];
            A(g, nest, k2, i, k2) = 0.0;
            for (int j = 1; j <= k1; ++j)
                A(g, nest, k2, i, j) = A(a, nest, k1, i, j);
        }

        for (int it = 1; it <= n8; ++it) {
            for (int i = 1; i <= k2; ++i)
                h[i - 1] = bmat[(size_t)(it - 1) * k2 + (i - 1)] * pinv;

            double yi = 0.0;
            for (int j = it; j <= nk1; ++j) {
                double piv = h[0];
                double cos, sin;
                fpgivs(piv, A(g, nest, k2, j, 1), cos, sin);
                fprota(cos, sin, yi, c_[j - 1]);
                if (j == nk1) break;

                int i2 = k1;
                if (j > n8) i2 = nk1 - j;
                for (int i = 1; i <= i2; ++i) {
                    int i1 = i + 1;
                    fprota(cos, sin, h[i1 - 1], A(g, nest, k2, j, i1));
                    h[i - 1] = h[i1 - 1];
                }
                h[i2] = 0.0;
            }
        }

        fpback(g, c_, nk1, k2, c_);

        fp = 0.0;
        int l = k2;
        for (int it = 1; it <= m; ++it) {
            if (x_[it - 1] >= t_[l - 1] && l <= nk1) ++l;
            int l0 = l - k2;
            double term = 0.0;
            for (int j = 1; j <= k1; ++j) {
                ++l0;
                term += c_[l0 - 1] * q[(size_t)(it - 1) * k1 + (j - 1)];
            }
            fp += std::pow(w_[it - 1] * (term - y_[it - 1]), 2);
        }

        fpms = fp - s_;
        if (std::abs(fpms) < acc) { fp_ = fp; return; }
        if (iter == maxit)        { ier_ = 3; fp_ = fp; return; }

        double p2 = p, f2 = fpms;

        if (ich3 == 0) {
            if ((f2 - f3) <= acc) {
                p3 = p2; f3 = f2;
                p = p * con4;
                if (p <= p1) p = p1 * con9 + p2 * con1;
                continue;
            }
            if (f2 < 0.0) ich3 = 1;
        }
        if (ich1 == 0) {
            if ((f1 - f2) <= acc) {
                p1 = p2; f1 = f2;
                p = p / con4;
                if (p3 >= 0.0 && p >= p3) p = p2 * con1 + p3 * con9;
                continue;
            }
            if (f2 > 0.0) ich1 = 1;
        }

        if (f2 >= f1 || f2 <= f3) { ier_ = 2; fp_ = fp; return; }

        p = fprati(p1, f1, p2, f2, p3, f3);
    }
    ier_ = 3;
    fp_ = fp;
}

// ---- Public constructor ----------------------------------------------------

inline Spline1D::Spline1D(const std::vector<double>& x,
                          const std::vector<double>& y,
                          int k,
                          double s,
                          const std::vector<double>& w)
    : x_(x), y_(y), k_(k), s_(s)
{
    if (x.size() != y.size())
        throw std::invalid_argument("Spline1D: x and y must have the same size");
    if (k < 1 || k > 5)
        throw std::invalid_argument("Spline1D: degree k must be in [1,5]");

    m_ = static_cast<int>(x.size());
    if (m_ <= k)
        throw std::invalid_argument("Spline1D: need m > k data points");

    if (w.empty()) w_.assign(m_, 1.0);
    else {
        if (static_cast<int>(w.size()) != m_)
            throw std::invalid_argument("Spline1D: w must have size m");
        w_ = w;
    }
    for (double wi : w_)
        if (wi <= 0.0)
            throw std::invalid_argument("Spline1D: weights must be positive");

    for (int i = 1; i < m_; ++i)
        if (!(x_[i - 1] < x_[i]))
            throw std::invalid_argument("Spline1D: x must be strictly increasing");

    xb_ = x_.front();
    xe_ = x_.back();

    int nest;
    if (s_ == 0.0) {
        nest = m_ + k + 1;
    } else {
        nest = std::max(2 * (k + 1) + 2,
                        std::min(m_ + k + 1,
                                 std::max(m_ / 2, 2 * (k + 1) + 2)));
    }
    t_.assign(nest, 0.0);

    if (s_ < 0.0) {
        int  n_save   = n_;
        int  ier_save = ier_;
        std::vector<double> t_save = t_;

        n_  = 2 * (k + 1);
        ier_ = 10;
        fpchec(ier_);
        if (ier_ == 0) {
            fpcurf(-1, 20, 0.001);
            double fp0 = fp_;
            double m_d = static_cast<double>(m_);
            s_ = std::min(fp0, m_d);
        } else {
            s_ = 0.0;
        }

        t_  = t_save;
        n_  = n_save;
        ier_ = ier_save;
    }

    ier_ = 0;
    fpcurf(0, 20, 0.001);
}

// ---- Evaluation ------------------------------------------------------------

inline double Spline1D::operator()(double xq) const
{
    if (c_.empty() || n_ == 0)
        return std::numeric_limits<double>::quiet_NaN();

    int k1 = k_ + 1;
    int l;

    if (xq < t_[0]) {
        l = k1;
    } else if (xq >= t_[n_ - 1]) {
        l = n_ - k1;
    } else {
        int lo = k1, hi = n_ - k1;
        while (lo < hi) {
            int mid = (lo + hi + 1) / 2;
            if (xq < t_[mid - 1]) hi = mid - 1;
            else                  lo = mid;
        }
        l = lo;
    }

    double h[6];
    fpbspl(t_, n_, k_, xq, l, h);

    double yq = 0.0;
    int j = l - k1;
    for (int i = 1; i <= k1; ++i) {
        ++j;
        yq += c_[j - 1] * h[i - 1];
    }
    return yq;
}

// ---- Derivative ------------------------------------------------------------

inline double Spline1D::derivative(double xq, int order) const
{
    if (order <= 0) return (*this)(xq);
    double span = xe_ - xb_;
    double h = std::max(span * 1e-5, 1e-8);
    if (order == 1) {
        return ((*this)(xq + h) - (*this)(xq - h)) / (2.0 * h);
    }
    if (order == 2) {
        return ((*this)(xq + h) - 2.0 * (*this)(xq) + (*this)(xq - h)) / (h * h);
    }
    double f1 = derivative(xq + h, order - 1);
    double f2 = derivative(xq - h, order - 1);
    return (f1 - f2) / (2.0 * h);
}