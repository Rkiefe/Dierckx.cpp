#include <iostream>
#include <vector>

// #include "spline1d.hpp"
// #include "spline2d.hpp"

#include "src/dierckx1d.h"
#include "src/dierckx2d.h"


double trial_1d_fun(double xq){
    return pow(xq, 3);
}

double trial_2d_fun(double xq, double yq){
    return pow(xq, 2) + pow(yq, 3);
}

void test_1d_linear(){
    std::cout << "TESTING LINEAR SPLINES" << std::endl;

    std::vector<double> x(100);
    std::vector<double> y(100);

    for(uint i = 0; i<x.size(); i++){
        x[i] = (double) i/x.size();
        y[i] = trial_1d_fun(x[i]);
    }

    // Test linear fit over a range of points
    Spline1D spl(x, y, 1);
    std::vector<double> xq = {0.0, 0.25, 0.333, 0.5, 0.666, 0.75, 1.0};
    for(uint i = 0; i<xq.size(); i++){
        std::cout << "f = " << trial_1d_fun(xq[i]) << " | spl = " << spl(xq[i]) << std::endl;
    }

    return;
}

void test_1d_quadratic(){
    std::cout << "TESTING QUADRATIC SPLINES" << std::endl;

    std::vector<double> x(100);
    std::vector<double> y(100);

    for(uint i = 0; i<x.size(); i++){
        x[i] = (double) i/x.size();
        y[i] = trial_1d_fun(x[i]);
    }

    // Test linear fit over a range of points
    Spline1D spl(x, y, 2);
    std::vector<double> xq = {0.0, 0.25, 0.333, 0.5, 0.666, 0.75, 1.0};
    for(uint i = 0; i<xq.size(); i++){
        std::cout << "f = " << trial_1d_fun(xq[i]) << " | spl = " << spl(xq[i]) << std::endl;
    }

    return;
}

void test_1d_cubic(){
    std::cout << "TESTING CUBIC SPLINES" << std::endl;

    std::vector<double> x(100);
    std::vector<double> y(100);

    for(uint i = 0; i<x.size(); i++){
        x[i] = (double) i/x.size();
        y[i] = trial_1d_fun(x[i]);
    }

    // Test linear fit over a range of points
    Spline1D spl(x, y, 3);
    std::vector<double> xq = {0.0, 0.25, 0.333, 0.5, 0.666, 0.75, 1.0};
    for(uint i = 0; i<xq.size(); i++){
        std::cout << "f = " << trial_1d_fun(xq[i]) << " | spl = " << spl(xq[i]) << std::endl;
    }

    return;
}

void test_2d_spline(){

    const int N = 50;
    std::vector<double> x(N*N);
    std::vector<double> y(N*N);
    std::vector<double> z(N*N);
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            double xi = (double)i / (N - 1);   // 0 .. 1
            double yj = (double)j / (N - 1);   // 0 .. 1
            x[i*N + j] = xi;
            y[i*N + j] = yj;
            z[i*N + j] = trial_2d_fun(xi, yj);
        }
    }

    int kx = 2, ky = 3;
    Spline2D spl(x, y, z, kx, ky);

    std::cout << "  ier      = " << spl.errorCode()    << "\n"
              << "  knots    = " << spl.numKnotsX()
              << " x "        << spl.numKnotsY()   << "\n"
              << "  ncof     = " << spl.numCoeffs()    << "\n"
              << "  residual = " << spl.residual()     << "\n"
              << "  smoothing= " << spl.smoothing()    << "\n";


    std::vector<double> xq = {0.0, 0.25, 0.333, 0.5, 0.666, 0.75, 1.0};
    std::vector<double> yq = {0.0, 0.25, 0.333, 0.5, 0.666, 0.75, 1.0};
    for(uint i = 0; i<xq.size(); i++){
        for(uint j = 0; j<xq.size(); j++){
            std::cout << "f = " << trial_2d_fun(xq[i], yq[j]);
            std::cout << " | spl = " << spl(xq[i], yq[j]) << std::endl;
        }
    }
    return;
}

int main() {

    test_1d_linear();
    test_1d_quadratic();
    test_1d_cubic();
    test_2d_spline();

    return 0;
}