#pragma once
#include <vector>
#include <cmath>
#include <stdexcept>
#include <algorithm>

// ---------------------------------------------------------------------------
//  Spline2D — Smoothing / least-squares B-spline surface fitting
//  Port of P. Dierckx's SURFIT (FITPACK) to modern C++.
//  Header-only: just #include "spline2d.h".
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

    static constexpr double EPS_DEFAULT = 1e-16;   // requires C++17

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
};

// ============================================================================
//  Implementation (all inline)
// ============================================================================

// ---- Low-level kernels -----------------------------------------------------

inline void Spline2D::fpbspl(const std::vector<double>& t, int /*n*/, int k,
                             double x, int l, double h[6]) {
    const double one = 1.0;
    double hh[5];
    h[0] = one;
    for (int j = 1; j <= k; ++j) {
        for (int i = 1; i <= j; ++i) hh[i-1] = h[i-1];
        h[0] = 0.0;
        for (int i = 1; i <= j; ++i) {
            int li = l + i, lj = li - j;
            double f = hh[i-1] / (t[li-1] - t[lj-1]);
            h[i-1] = h[i-1] + f * (t[li-1] - x);
            h[i]   = f * (x - t[lj-1]);
        }
    }
}

inline void Spline2D::fpgivs(double piv, double& ww, double& cos, double& sin) {
    const double one = 1.0;
    double store = std::abs(piv);
    if (store == 0.0 && ww == 0.0) { cos = 1.0; sin = 0.0; return; }
    double dd;
    if (store >= ww) dd = store * std::sqrt(one + (ww/piv)*(ww/piv));
    else             dd = ww * std::sqrt(one + (piv/ww)*(piv/ww));
    cos = ww / dd;  sin = piv / dd;  ww = dd;
}

inline void Spline2D::fprota(double cos, double sin, double& a, double& b) {
    double s1 = a, s2 = b;
    b = cos*s2 + sin*s1;
    a = cos*s1 - sin*s2;
}

inline void Spline2D::fpdisc(const std::vector<double>& t, int n, int k2,
                             std::vector<double>& b, int /*nest*/) {
    int k1 = k2-1, k = k1-1;
    int nk1 = n - k1, nrint = nk1 - k;
    double an = (double)nrint;
    double fac = an / (t[nk1] - t[k1-1]);
    double h[12];
    auto T = [&](int idx1) -> double { return t[idx1-1]; };
    for (int l = k2; l <= nk1; ++l) {
        int lmk = l - k1;
        for (int j = 1; j <= k1; ++j) {
            int ik = j + k1, lj = l + j, lk = lj - k2;
            h[j-1]  = T(l) - T(lk);
            h[ik-1] = T(l) - T(lj);
        }
        int lp = lmk;
        for (int j = 1; j <= k2; ++j) {
            int jk = j;
            double prod = h[j-1];
            for (int i = 1; i <= k; ++i) { jk++; prod = prod * h[jk-1] * fac; }
            int lk = lp + k1;
            b[(lmk-1)*k2 + (j-1)] = (T(lk) - T(lp)) / prod;
            lp++;
        }
    }
}

inline void Spline2D::fporde(const std::vector<double>& x, const std::vector<double>& y,
                             int m, int kx, int ky,
                             const std::vector<double>& tx, int nx,
                             const std::vector<double>& ty, int ny,
                             std::vector<int>& nummer, std::vector<int>& index, int nreg) {
    int kx1 = kx+1, ky1 = ky+1;
    int nk1x = nx - kx1, nk1y = ny - ky1;
    int nyy = nk1y - ky;
    for (int i = 0; i < nreg; ++i) index[i] = -1;
    for (int im = 0; im < m; ++im) {
        double xi = x[im], yi = y[im];
        int l = kx1, l1 = l+1;
        while (!(xi < tx[l1-1] || l == nk1x)) { l = l1; l1++; }
        int k = ky1, k1 = k+1;
        while (!(yi < ty[k1-1] || k == nk1y)) { k = k1; k1++; }
        int num = (l-kx1)*nyy + k-ky;
        nummer[im]     = index[num-1];
        index[num-1]   = im;
    }
}

inline void Spline2D::fpback(const std::vector<double>& a, int lda,
                             const std::vector<double>& z, int n, int k,
                             std::vector<double>& c) {
    int k1 = k-1;
    c[n-1] = z[n-1] / AM(a, lda, n, 1);
    int i = n-1;
    if (i == 0) return;
    for (int j = 2; j <= n; ++j) {
        double store = z[i-1];
        int i1 = k1;  if (j <= k1) i1 = j-1;
        int mm = i;
        for (int l = 1; l <= i1; ++l) {
            ++mm;
            store -= c[mm-1] * AM(a, lda, i, l+1);
        }
        c[i-1] = store / AM(a, lda, i, 1);
        --i;
    }
}

inline double Spline2D::fprati(double& p1, double& f1, double p2, double f2,
                               double& p3, double& f3) {
    double p;
    if (p3 > 0.0) {
        double h1 = f1*(f2-f3), h2 = f2*(f3-f1), h3 = f3*(f1-f2);
        p = -(p1*p2*h3 + p2*p3*h1 + p3*p1*h2) / (p1*h1 + p2*h2 + p3*h3);
    } else {
        p = (p1*(f1-f3)*f2 - p2*(f2-f3)*f1) / ((f1-f2)*f3);
    }
    if (f2 >= 0.0) { p1 = p2; f1 = f2; }
    else           { p3 = p2; f3 = f2; }
    return p;
}

// Minimum-norm least-squares solve for rank-deficient system.
inline void Spline2D::fprank(std::vector<double>& a, int lda,
                             std::vector<double>& f, int n, int m, double tol,
                             std::vector<double>& c, double& sq, int& rank,
                             std::vector<double>& aa, std::vector<double>& ff,
                             std::vector<double>& h) {
    int m1 = m-1;
    int nl = 0;  sq = 0.0;

    // Phase 1: find rank deficiency, zero out small diagonals.
    for (int i = 1; i <= n; ++i) {
        if (AM(a, lda, i, 1) > tol) continue;
        ++nl;
        if (i == n) continue;
        double yi = f[i-1];
        for (int j = 1; j <= m1; ++j) h[j-1] = AM(a, lda, i, j+1);
        h[m-1] = 0.0;
        for (int ii = i+1; ii <= n; ++ii) {
            int i2 = std::min(n-ii, m1);
            double piv = h[0];
            if (piv == 0.0) {
                if (i2 == 0) break;
                for (int j = 1; j <= i2; ++j) h[j-1] = h[j];
                h[i2] = 0.0;
                continue;
            }
            double cos, sin;
            fpgivs(piv, AM(a, lda, ii, 1), cos, sin);
            fprota(cos, sin, yi, f[ii-1]);
            if (i2 == 0) break;
            for (int j = 1; j <= i2; ++j) {
                int j1 = j+1;
                fprota(cos, sin, h[j1-1], AM(a, lda, ii, j1));
                h[j-1] = h[j1-1];
            }
            h[i2] = 0.0;
        }
        sq += yi * yi;
    }
    rank = n - nl;

    // Phase 2: form reduced upper trapezoidal matrix in aa.
    for (int i = 1; i <= rank; ++i)
        for (int j = 1; j <= m; ++j) AM(aa, m, i, j) = 0.0;
    int ii = 0;
    for (int i = 1; i <= n; ++i) {
        if (AM(a, lda, i, 1) <= tol) continue;
        ++ii;
        ff[ii-1] = f[i-1];
        AM(aa, m, ii, 1) = AM(a, lda, i, 1);
        int jj = ii, kk = 1, j = i;
        int j1 = std::min(j-1, m1);
        for (int k = 1; k <= j1; ++k) {
            --j;
            if (AM(a, lda, j, 1) <= tol) continue;
            ++kk; --jj;
            AM(aa, m, jj, kk) = AM(a, lda, j, k+1);
        }
    }

    // Phase 3: rotate zero-diagonal columns into aa.
    ii = 0;
    for (int i = 1; i <= n; ++i) {
        ++ii;
        if (AM(a, lda, i, 1) > tol) continue;
        --ii;
        if (ii == 0) continue;
        int jj_ = 1, j = i;
        int j1 = std::min(j-1, m1);
        for (int k = 1; k <= j1; ++k) {
            --j;
            if (AM(a, lda, j, 1) <= tol) continue;
            h[jj_-1] = AM(a, lda, j, k+1);
            ++jj_;
        }
        for (int kk = jj_; kk <= m; ++kk) h[kk-1] = 0.0;
        int jj = ii;
        bool exit_inner = false;
        for (int i1 = 1; i1 <= ii; ++i1) {
            int j1_ = std::min(jj-1, m1);
            double piv = h[0];
            int j3 = 0;
            if (piv == 0.0) {
                if (j1_ == 0) { exit_inner = true; break; }
                for (int j2 = 1; j2 <= j1_; ++j2) {
                    j3 = j2+1; h[j2-1] = h[j3-1];
                }
            } else {
                double cos, sin;
                fpgivs(piv, AM(aa, m, jj, 1), cos, sin);
                if (j1_ == 0) { exit_inner = true; break; }
                int kk = jj;
                for (int j2 = 1; j2 <= j1_; ++j2) {
                    j3 = j2+1; --kk;
                    fprota(cos, sin, h[j3-1], AM(aa, m, kk, j3));
                    h[j2-1] = h[j3-1];
                }
            }
            --jj;
            h[j3-1] = 0.0;
        }
        if (exit_inner) continue;
    }

    // Solve (aa)(f1) = ff
    ff[rank-1] = ff[rank-1] / AM(aa, m, rank, 1);
    int i = rank-1;
    if (i > 0) {
        for (int j = 2; j <= rank; ++j) {
            double store = ff[i-1];
            int i1 = std::min(j-1, m1), k = i;
            for (int ii2 = 1; ii2 <= i1; ++ii2) {
                ++k;
                store -= ff[k-1] * AM(aa, m, i, ii2+1);
            }
            ff[i-1] = store / AM(aa, m, i, 1);
            --i;
        }
    }
    // Solve (aa)'(f2) = f1
    ff[0] = ff[0] / AM(aa, m, 1, 1);
    if (rank > 1) {
        for (int j = 2; j <= rank; ++j) {
            double store = ff[j-1];
            int i1 = std::min(j-1, m1), k = j;
            for (int ii2 = 1; ii2 <= i1; ++ii2) {
                --k;
                store -= ff[k-1] * AM(aa, m, k, ii2+1);
            }
            ff[j-1] = store / AM(aa, m, j, 1);
        }
    }
    // Premultiply by transpose of a.
    int k = 0;
    for (int i2 = 1; i2 <= n; ++i2) {
        double store = 0.0;
        if (AM(a, lda, i2, 1) > tol) ++k;
        int j1 = std::min(i2, m), kk = k, ij = i2+1;
        for (int j = 1; j <= j1; ++j) {
            --ij;
            if (AM(a, lda, ij, 1) <= tol) continue;
            store += AM(a, lda, ij, j) * ff[kk-1];
            --kk;
        }
        c[i2-1] = store;
    }
    // Add rank-reduction contribution to sq.
    double stor3 = 0.0;
    for (int i2 = 1; i2 <= n; ++i2) {
        if (AM(a, lda, i2, 1) > tol) continue;
        double store = f[i2-1];
        int i1 = std::min(n-i2, m1);
        for (int j = 1; j <= i1; ++j) {
            store -= c[i2+j-1] * AM(a, lda, i2, j+1);
        }
        double stor1 = AM(a, lda, i2, 1) * c[i2-1];
        stor3 += stor1 * (stor1 - store - store);
    }
    sq += stor3;
}

// ---- Main fitting driver ---------------------------------------------------

inline void Spline2D::fpsurf(int iopt) {
    const double one = 1.0, con1 = 0.1, con9 = 0.9, con4 = 0.04, ten = 10.0;
    const double eps = EPS_DEFAULT, tol = 0.001;
    const int maxit = 20;

    int kx = kx_, ky = ky_;
    int kx1 = kx+1, ky1 = ky+1;
    int nminx = 2*kx1, nminy = 2*ky1;
    int kmax = std::max(kx, ky);
    int km1 = kmax+1, km2 = km1+1;

    // Workspace sizing.
    int sqm  = (int)(std::sqrt((double)m_) + 0.5);
    int nxest = std::max(2*kx1 + 2, 2*kx1 + 2*sqm);
    int nyest = std::max(2*ky1 + 2, 2*ky1 + 2*sqm);

    int nmax = std::max(nxest, nyest);
    int nxk  = nxest - kx1;
    int nyk  = nyest - ky1;
    int ncest = nxk * nyk;
    int nmx  = nxest - nminx + 1;
    int nmy  = nyest - nminy + 1;
    int nrint_max = nmx + nmy;
    int nreg_max  = nmx * nmy;

    int ib1, ib3;
    {
        int a1 = kx*nyk + ky + 1;
        int b1 = ky*nxk + kx + 1;
        if (a1 <= b1) { ib1 = a1; ib3 = kx1*nyk + 1; }
        else          { ib1 = b1; ib3 = ky1*nxk + 1; }
    }
    int lwrk2 = ncest*(std::max(ib1, ib3) + 1) + std::max(ib1, ib3) + 16;

    // Allocate working storage.
    std::vector<double> f(ncest, 0.0), ff(ncest, 0.0);
    std::vector<double> fpint(nrint_max, 0.0), coord(nrint_max, 0.0);
    std::vector<double> a((size_t)ncest * ib1, 0.0);
    std::vector<double> q((size_t)ncest * ib3, 0.0);
    std::vector<double> bx((size_t)nmax * km2, 0.0);
    std::vector<double> by((size_t)nmax * km2, 0.0);
    std::vector<double> spx((size_t)m_ * km1, 0.0);
    std::vector<double> spy((size_t)m_ * km1, 0.0);
    std::vector<double> h(std::max(ib3, 12), 0.0);
    std::vector<int>    index(nreg_max, -1), nummer(m_, -1);
    std::vector<double> wrk(lwrk2, 0.0);

    // -------------------------------------------------------------------
    // Data scale — used to detect numerical-zero residual for s == 0.
    // Without this, s = 0 falls through the |fp - s| <= acc test forever
    // (acc = tol * 0 = 0) and never gets a chance to report ier = -1.
    // -------------------------------------------------------------------
    double fp_scale = 0.0;
    for (int i = 0; i < m_; ++i) {
        double t = w_[i] * z_[i];
        fp_scale += t * t;
    }
    const double fp_zero_thresh = std::max(1e-20 * fp_scale, 1e-30);

    // Local fitting state.
    int nx, ny;
    int iband = 0;
    double fp = 0.0, fp0 = 0.0, fpms = 0.0;
    double acc = tol * s_;
    int ncof = 0, rank = 0;
    int nyy = 0, nreg_cur = 0;

    // Initial knot setup.
    if (iopt < 0) {
        nx = nx_; ny = ny_;
    } else if (iopt == 0 || fp0 <= s_) {
        nx = nminx; ny = nminy;
        ier_ = -2;
    } else {
        nx = nx_; ny = ny_;
    }

    bool go_part2 = false;
    bool done     = false;

    for (int iter = 1; iter <= m_ && !done; ++iter) {
        if (nx > nxest || ny > nyest) { ier_ = 1; break; }
        if (nx == nminx && ny == nminy) ier_ = -2;

        // Boundary knots.
        {
            int l = nx;
            for (int i = 1; i <= kx1; ++i) { tx_[i-1] = xb_; tx_[l-1] = xe_; --l; }
            l = ny;
            for (int i = 1; i <= ky1; ++i) { ty_[i-1] = yb_; ty_[l-1] = ye_; --l; }
        }

        int nxx = nx - 2*kx1 + 1;
        nyy = ny - 2*ky1 + 1;
        int nrint_cur = nxx + nyy;
        nreg_cur  = nxx * nyy;

        // Bandwidth ("no swap" value; we skip FITPACK's coordinate swap).
        iband = kx*(ny-ky1) + ky + 1;

        fporde(x_, y_, m_, kx, ky, tx_, nx, ty_, ny, nummer, index, nreg_cur);

        int nk1x = nx - kx1;
        int nk1y = ny - ky1;
        ncof = nk1x * nk1y;
        if (ncof > ncest) { ier_ = 1; break; }

        // Init a and f.
        for (int i = 1; i <= ncof; ++i) {
            f[i-1] = 0.0;
            for (int j = 1; j <= iband; ++j) AM(a, ib1, i, j) = 0.0;
        }
        fp = 0.0;

        // Main data loop over panels.
        for (int num = 1; num <= nreg_cur; ++num) {
            int num1 = num-1;
            int lx_ = num1 / nyy;
            int l1_ = lx_ + kx1;
            int ly_ = num1 - lx_*nyy;
            int l2_ = ly_ + ky1;
            int jrot = lx_ * nk1y + ly_;

            int in = index[num-1];
            while (in >= 0) {
                double wi = w_[in];
                double zi = z_[in] * wi;

                double hx[6] = {0}, hy[6] = {0};
                fpbspl(tx_, nx, kx, x_[in], l1_, hx);
                fpbspl(ty_, ny, ky, y_[in], l2_, hy);

                for (int i = 1; i <= kx1; ++i) spx[(size_t)in*km1 + (i-1)] = hx[i-1];
                for (int j = 1; j <= ky1; ++j) spy[(size_t)in*km1 + (j-1)] = hy[j-1];

                for (int i = 0; i < iband; ++i) h[i] = 0.0;

                int i1 = 0;
                for (int i = 1; i <= kx1; ++i) {
                    double hxi = hx[i-1];
                    int j1 = i1;
                    for (int j = 1; j <= ky1; ++j) {
                        ++j1;
                        h[j1-1] = hxi * hy[j-1] * wi;
                    }
                    i1 += nk1y;
                }

                int irot = jrot;
                for (int i = 1; i <= iband; ++i) {
                    ++irot;
                    double piv = h[i-1];
                    if (piv == 0.0) continue;
                    double cos, sin;
                    fpgivs(piv, AM(a, ib1, irot, 1), cos, sin);
                    fprota(cos, sin, zi, f[irot-1]);
                    if (i == iband) break;
                    int i2 = 1;
                    for (int j = i+1; j <= iband; ++j) {
                        ++i2;
                        fprota(cos, sin, h[j-1], AM(a, ib1, irot, i2));
                    }
                }
                fp += zi * zi;
                in = nummer[in];
            }
        }

        // Rank check + back-substitution.
        double dmax = 0.0;
        for (int i = 1; i <= ncof; ++i)
            if (AM(a, ib1, i, 1) > dmax) dmax = AM(a, ib1, i, 1);

        double sigma = eps * dmax;
        bool fullrank = true;
        for (int i = 1; i <= ncof; ++i)
            if (AM(a, ib1, i, 1) <= sigma) { fullrank = false; break; }

        std::vector<double> ccur(ncest, 0.0);
        double sq = 0.0;
        if (fullrank) {
            fpback(a, ib1, f, ncof, iband, ccur);
            rank = ncof;
            for (int i = 1; i <= ncof; ++i)
                AM(q, ib3, i, 1) = AM(a, ib1, i, 1) / dmax;
        } else {
            for (int i = 1; i <= ncof; ++i)
                for (int j = 1; j <= iband; ++j)
                    AM(q, ib3, i, j) = AM(a, ib1, i, j);
            std::vector<double> aa((size_t)ncof * iband, 0.0);
            std::vector<double> ffl(ncof, 0.0), hh(iband, 0.0);
            fprank(q, ib3, f, ncof, iband, sigma, ccur, sq, rank, aa, ffl, hh);
            for (int i = 1; i <= ncof; ++i)
                AM(q, ib3, i, 1) = AM(q, ib3, i, 1) / dmax;
            fp += sq;
        }
        for (int i = 0; i < ncof; ++i) c_[i] = ccur[i];
        nx_ = nx; ny_ = ny; ncof_ = ncof; fp_ = fp;

        if (iopt < 0) { done = true; break; }

        // Acceptance check.  For s > 0, |fp - s| <= acc = tol*s.  For
        // s == 0, acc is 0 and we additionally treat fp below the
        // numerical-zero threshold as exact interpolation.
        fpms = fp - s_;
        bool eff_zero = (s_ == 0.0) && (fp <= fp_zero_thresh);
        if (std::abs(fpms) <= acc || eff_zero) {
            if (eff_zero || fp == 0.0) {
                ier_ = -1;
                fp_  = 0.0;
            } else {
                ier_ = (ncof != rank) ? -rank : 0;
                fp_  = fp;
            }
            done = true; break;
        }
        if (fpms < 0.0) { go_part2 = true; break; }
        if (ncof > m_)  { ier_ = 4; done = true; break; }
        ier_ = 0;

        // Compute per-panel residuals and coordinates.
        for (int i = 1; i <= nrint_cur; ++i) { fpint[i-1] = 0.0; coord[i-1] = 0.0; }
        for (int num = 1; num <= nreg_cur; ++num) {
            int num1 = num-1;
            int lx_ = num1 / nyy;
            int l1_ = lx_ + 1;
            int ly_ = num1 - lx_*nyy;
            int l2_ = ly_ + 1 + nxx;
            int jrot = lx_ * nk1y + ly_;

            int in = index[num-1];
            while (in >= 0) {
                double store = 0.0;
                int i1 = jrot;
                for (int i = 1; i <= kx1; ++i) {
                    double hxi = spx[(size_t)in*km1 + (i-1)];
                    int j1 = i1;
                    for (int j = 1; j <= ky1; ++j) {
                        ++j1;
                        store += hxi * spy[(size_t)in*km1 + (j-1)] * c_[j1-1];
                    }
                    i1 += nk1y;
                }
                double t = w_[in] * (z_[in] - store);
                store = t * t;
                fpint[l1_-1] += store;
                coord[l1_-1] += store * x_[in];
                fpint[l2_-1] += store;
                coord[l2_-1] += store * y_[in];
                in = nummer[in];
            }
        }

        // ---------------------------------------------------------------
        //  Knot addition with retry (matches Fortran's "go to 370").
        //  If the candidate panel's knot-position test fails, we zero its
        //  fpint entry and re-pick the max WITHOUT recomputing fpint (the
        //  old `continue` bug re-ran the whole outer iteration and lost
        //  the fpint mutation, so the same bad panel was re-selected).
        // ---------------------------------------------------------------
        bool knot_added = false;
        for (int attempt = 1; attempt <= nrint_cur; ++attempt) {
            int l = 0;
            double fpmax = 0.0;
            int l1_ = 1, l2_ = nrint_cur;
            if (nx == nxest) l1_ = nxx + 1;
            if (ny == nyest) l2_ = nxx;
            if (l1_ > l2_) { ier_ = 1; done = true; break; }

            for (int i = l1_; i <= l2_; ++i) {
                if (fpmax >= fpint[i-1]) continue;
                l = i; fpmax = fpint[i-1];
            }
            if (l == 0) { ier_ = 5; done = true; break; }

            double arg = coord[l-1] / fpint[l-1];

            if (l > nxx) {
                int jxy = l + ky1 - nxx;
                fpint[l-1] = 0.0;
                double fac1 = ty_[jxy-1] - arg, fac2 = arg - ty_[jxy-2];
                if (fac1 > ten*fac2 || fac2 > ten*fac1) continue;
                for (int i = jxy; i <= ny; ++i) ty_[i] = ty_[i-1];
                ty_[jxy-1] = arg;
                ++ny;
                knot_added = true;
                break;
            } else {
                int jxy = l + kx1;
                fpint[l-1] = 0.0;
                double fac1 = tx_[jxy-1] - arg, fac2 = arg - tx_[jxy-2];
                if (fac1 > ten*fac2 || fac2 > ten*fac1) continue;
                for (int i = jxy; i <= nx; ++i) tx_[i] = tx_[i-1];
                tx_[jxy-1] = arg;
                ++nx;
                knot_added = true;
                break;
            }
        }
        if (done) break;
        if (!knot_added) { ier_ = 5; done = true; break; }
    }

    if (ier_ == -2) return;
    if (!go_part2) return;

    // ---- Part 2: smoothing spline ----
    int nk1x = nx - kx1, nk1y = ny - ky1;
    int iband3 = kx1*nk1y;
    int iband4 = iband3 + 1;

    int kx2 = kx1+1, ky2 = ky1+1;
    if (nk1x != kx1) fpdisc(tx_, nx, kx2, bx, nmax);
    if (nk1y != ky1) fpdisc(ty_, ny, ky2, by, nmax);

    double p1 = 0.0, f1 = fp0 - s_;
    double p3 = -one, f3 = fpms;
    double p = 0.0;
    for (int i = 1; i <= ncof; ++i) p += AM(a, ib1, i, 1);
    p = (double)ncof / p;

    int ich1 = 0, ich3 = 0;

    for (int iter = 1; iter <= maxit; ++iter) {
        double pinv = one / p;

        for (int i = 1; i <= ncof; ++i) {
            ff[i-1] = f[i-1];
            for (int j = 1; j <= iband; ++j) AM(q, ib3, i, j) = AM(a, ib1, i, j);
            for (int j = iband+1; j <= iband4; ++j) AM(q, ib3, i, j) = 0.0;
        }

        // Extend with y-polynomial rows.
        if (nk1y != ky1) {
            for (int i = ky2; i <= nk1y; ++i) {
                int ii = i - ky1;
                for (int j = 1; j <= nk1x; ++j) {
                    for (int l = 1; l <= iband; ++l) h[l-1] = 0.0;
                    for (int l = 1; l <= ky2; ++l) h[l-1] = by[(size_t)(ii-1)*km2 + (l-1)] * pinv;
                    double zi = 0.0;
                    int jrot = (j-1)*nk1y + ii;
                    for (int irot = jrot; irot <= ncof; ++irot) {
                        double piv = h[0];
                        int i2 = std::min(iband - 1, ncof - irot);
                        if (piv == 0.0 && i2 <= 0) break;
                        double cos, sin;
                        fpgivs(piv, AM(q, ib3, irot, 1), cos, sin);
                        fprota(cos, sin, zi, ff[irot-1]);
                        if (i2 == 0) break;
                        for (int l = 1; l <= i2; ++l) {
                            int l1 = l+1;
                            fprota(cos, sin, h[l1-1], AM(q, ib3, irot, l1));
                        }
                        for (int l = 1; l <= i2; ++l) h[l-1] = h[l];
                        h[i2] = 0.0;
                    }
                }
            }
        }

        // Extend with x-polynomial rows.
        if (nk1x != kx1) {
            for (int i = kx2; i <= nk1x; ++i) {
                int ii = i - kx1;
                for (int j = 1; j <= nk1y; ++j) {
                    for (int l = 1; l <= iband4; ++l) h[l-1] = 0.0;
                    int j1 = 1;
                    for (int l = 1; l <= kx2; ++l) {
                        h[j1-1] = bx[(size_t)(ii-1)*km2 + (l-1)] * pinv;
                        j1 += nk1y;
                    }
                    double zi = 0.0;
                    int jrot = (i-kx2)*nk1y + j;
                    for (int irot = jrot; irot <= ncof; ++irot) {
                        double piv = h[0];
                        int i2 = std::min(iband3, ncof - irot);
                        if (piv == 0.0 && i2 <= 0) break;
                        double cos, sin;
                        fpgivs(piv, AM(q, ib3, irot, 1), cos, sin);
                        fprota(cos, sin, zi, ff[irot-1]);
                        if (i2 == 0) break;
                        for (int l = 1; l <= i2; ++l) {
                            int l1 = l+1;
                            fprota(cos, sin, h[l1-1], AM(q, ib3, irot, l1));
                        }
                        for (int l = 1; l <= i2; ++l) h[l-1] = h[l];
                        h[i2] = 0.0;
                    }
                }
            }
        }

        // Rank check + solve.
        double dmax = 0.0;
        for (int i = 1; i <= ncof; ++i)
            if (AM(q, ib3, i, 1) > dmax) dmax = AM(q, ib3, i, 1);
        double sigma = eps * dmax;
        bool fullrank = true;
        for (int i = 1; i <= ncof; ++i)
            if (AM(q, ib3, i, 1) <= sigma) { fullrank = false; break; }

        std::vector<double> ccur(ncest, 0.0);
        double sq = 0.0;
        if (fullrank) {
            fpback(q, ib3, ff, ncof, iband4, ccur);
            rank = ncof;
        } else {
            std::vector<double> aa((size_t)ncof * iband4, 0.0);
            std::vector<double> ffl(ncof, 0.0), hh(iband4, 0.0);
            fprank(q, ib3, ff, ncof, iband4, sigma, ccur, sq, rank, aa, ffl, hh);
        }
        for (int i = 1; i <= ncof; ++i) AM(q, ib3, i, 1) /= dmax;

        // Evaluate fp.
        fp = 0.0;
        for (int num = 1; num <= nreg_cur; ++num) {
            int num1 = num-1;
            int lx_ = num1 / nyy;
            int ly_ = num1 - lx_*nyy;
            int jrot = lx_ * nk1y + ly_;
            int in = index[num-1];
            while (in >= 0) {
                double store = 0.0;
                int i1 = jrot;
                for (int i = 1; i <= kx1; ++i) {
                    double hxi = spx[(size_t)in*km1 + (i-1)];
                    int j1 = i1;
                    for (int j = 1; j <= ky1; ++j) {
                        ++j1;
                        store += hxi * spy[(size_t)in*km1 + (j-1)] * ccur[j1-1];
                    }
                    i1 += nk1y;
                }
                double t = w_[in] * (z_[in] - store);
                fp += t * t;
                in = nummer[in];
            }
        }

        for (int i = 0; i < ncof; ++i) c_[i] = ccur[i];
        ncof_ = ncof; fp_ = fp;

        fpms = fp - s_;
        if (std::abs(fpms) <= acc) { ier_ = (ncof != rank) ? -rank : 0; return; }
        if (iter == maxit)        { ier_ = 3; return; }

        double p2 = p, f2 = fpms;
        if (ich3 == 0) {
            if ((f2 - f3) <= acc) {
                p3 = p2; f3 = f2;
                p = p * con4;
                if (p <= p1) p = p1*con9 + p2*con1;
                continue;
            }
            if (f2 < 0.0) ich3 = 1;
        }
        if (ich1 == 0) {
            if ((f1 - f2) <= acc) {
                p1 = p2; f1 = f2;
                p = p / con4;
                if (p3 >= 0.0 && p >= p3) p = p2*con1 + p3*con9;
                continue;
            }
            if (f2 > 0.0) ich1 = 1;
        }
        if (f2 >= f1 || f2 <= f3) { ier_ = 2; return; }
        p = fprati(p1, f1, p2, f2, p3, f3);
    }
    ier_ = 3;
}

// ---- Constructors ----------------------------------------------------------

inline Spline2D::Spline2D(const std::vector<double>& x, const std::vector<double>& y,
                          const std::vector<double>& z, int k,
                          double s, const std::vector<double>& w)
    : Spline2D(x, y, z, k, k, s, w) {}

inline Spline2D::Spline2D(const std::vector<double>& x, const std::vector<double>& y,
                          const std::vector<double>& z, int kx, int ky,
                          double s, const std::vector<double>& w)
    : x_(x), y_(y), z_(z), kx_(kx), ky_(ky), s_(s)
{
    if (x.size() != y.size() || x.size() != z.size())
        throw std::invalid_argument("Spline2D: x, y, z must have the same size");
    if (kx < 1 || kx > 5 || ky < 1 || ky > 5)
        throw std::invalid_argument("Spline2D: degrees must be in [1,5]");
    m_ = (int)x.size();
    if (m_ < (kx+1)*(ky+1))
        throw std::invalid_argument("Spline2D: m must be >= (kx+1)*(ky+1)");

    if (w.empty()) w_.assign(m_, 1.0);
    else {
        if ((int)w.size() != m_)
            throw std::invalid_argument("Spline2D: w must have size m");
        w_ = w;
    }
    for (double wi : w_)
        if (wi <= 0.0) throw std::invalid_argument("Spline2D: weights must be positive");

    xb_ = *std::min_element(x_.begin(), x_.end());
    xe_ = *std::max_element(x_.begin(), x_.end());
    yb_ = *std::min_element(y_.begin(), y_.end());
    ye_ = *std::max_element(y_.begin(), y_.end());
    if (xb_ == xe_ || yb_ == ye_)
        throw std::invalid_argument("Spline2D: x, y must span nonzero ranges");

    int sqm  = (int)(std::sqrt((double)m_) + 0.5);
    int nxest = std::max(2*(kx+1) + 2, 2*(kx+1) + 2*sqm);
    int nyest = std::max(2*(ky+1) + 2, 2*(ky+1) + 2*sqm);
    int nmax  = std::max(nxest, nyest);
    tx_.assign(nmax, 0.0);
    ty_.assign(nmax, 0.0);
    c_.assign((nxest-kx-1)*(nyest-ky-1), 0.0);

    ier_ = 0;
    fpsurf(0);
}

// ---- Evaluation ------------------------------------------------------------

inline double Spline2D::operator()(double xq, double yq) const {
    if (c_.empty() || nx_ == 0 || ny_ == 0) return std::nan("");

    int kx1 = kx_+1, ky1 = ky_+1;
    int nk1x = nx_ - kx1, nk1y = ny_ - ky1;

    auto find_span = [](double q, const std::vector<double>& t, int k1, int nk1) -> int {
        if (q < t[k1-1])   return k1;
        if (q >= t[nk1-1]) return nk1;
        int lo = k1, hi = nk1;
        while (lo < hi) {
            int mid = (lo + hi + 1) / 2;
            if (q < t[mid-1]) hi = mid - 1;
            else              lo = mid;
        }
        return lo;
    };
    int lx = find_span(xq, tx_, kx1, nk1x);
    int ly = find_span(yq, ty_, ky1, nk1y);

    double hx[6] = {0}, hy[6] = {0};
    fpbspl(tx_, nx_, kx_, xq, lx, hx);
    fpbspl(ty_, ny_, ky_, yq, ly, hy);

    double sum = 0.0;
    for (int bi = 0; bi < kx1; ++bi) {
        int i = lx - kx1 + bi + 1;
        for (int bj = 0; bj < ky1; ++bj) {
            int j = ly - ky1 + bj + 1;
            sum += c_[(size_t)(i-1)*nk1y + (j-1)] * hx[bi] * hy[bj];
        }
    }
    return sum;
}