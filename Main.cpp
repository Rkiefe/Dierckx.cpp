#include "spline1d.hpp"
#include <iostream>
#include <vector>

void test_linear(){
    std::cout << "TESTING LINEAR SPLINES" << std::endl;

    std::vector<double> x(100);
    std::vector<double> y(100);

    for(uint i = 0; i<x.size(); i++){
        x[i] = (double) i/x.size();
        y[i] = x[i]*x[i]*x[i];
    }

    // Test linear fit over a range of points
    Spline1D spl(x, y, 1);
    std::vector<double> xq = {0.0, 0.25, 0.333, 0.5, 0.666, 0.75, 1.0};
    for(uint i = 0; i<xq.size(); i++){
        std::cout << "f = " << xq[i]*xq[i]*xq[i] << " | spl = " << spl(xq[i]) << std::endl;
    }

    return;
}

void test_quadratic(){
    std::cout << "TESTING QUADRATIC SPLINES" << std::endl;

    std::vector<double> x(100);
    std::vector<double> y(100);

    for(uint i = 0; i<x.size(); i++){
        x[i] = (double) i/x.size();
        y[i] = x[i]*x[i]*x[i];
    }

    // Test linear fit over a range of points
    Spline1D spl(x, y, 2);
    std::vector<double> xq = {0.0, 0.25, 0.333, 0.5, 0.666, 0.75, 1.0};
    for(uint i = 0; i<xq.size(); i++){
        std::cout << "f = " << xq[i]*xq[i]*xq[i] << " | spl = " << spl(xq[i]) << std::endl;
    }

    return;
}

void test_cubic(){
    std::cout << "TESTING CUBIC SPLINES" << std::endl;

    std::vector<double> x(100);
    std::vector<double> y(100);

    for(uint i = 0; i<x.size(); i++){
        x[i] = (double) i/x.size();
        y[i] = x[i]*x[i]*x[i];
    }

    // Test linear fit over a range of points
    Spline1D spl(x, y, 3);
    std::vector<double> xq = {0.0, 0.25, 0.333, 0.5, 0.666, 0.75, 1.0};
    for(uint i = 0; i<xq.size(); i++){
        std::cout << "f = " << xq[i]*xq[i]*xq[i] << " | spl = " << spl(xq[i]) << std::endl;
    }

    return;
}

int main() {
    // std::vector<double> x = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
    // std::vector<double> y = {0, 0.8, 0.9, 0.1, -0.8, -1.0, -0.6, 0.2, 0.7, 1.0};

    test_linear();
    test_quadratic();
    test_cubic();





    // // Cubic with automatic smoothing
    // Spline1D cub(x, y, 3);
    // std::cout << "cub(2.5) = " << cub(2.5) << "\n";
    // std::cout << "knots = " << cub.numKnots()
    //           << "  residual = " << cub.residual() << "\n";

    // // Cubic with explicit smoothing
    // Spline1D cub2(x, y, 3, /*s=*/1.0);
    // std::cout << "cub2(2.5) = " << cub2(2.5) << "\n";

    // // Evaluate at query point
    // double yq = cub(2.5);
    // std::cout << "yq = " << yq << "\n";

    // std::cout << "cub2 residual = " << cub2.residual() << "\n";   // ≈ 1.0
    // std::cout << "cub2 ier = " << cub2.errorCode() << "\n";       // should be 0

    // for (int i = 0; i < (int)x.size(); ++i)
    //     std::cout << i << ": " << cub(x[i]) - y[i] << "\n";   // all ≈ 0

    // std::cout << "numCoeffs = " << cub.numCoeffs() << "\n";   // 10

    // for (int i = 0; i < cub.numKnots(); ++i)
    //     std::cout << "t[" << i << "] = " << cub.knots()[i] << "\n";

    return 0;
}