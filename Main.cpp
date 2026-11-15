// #include "spline1d.hpp"
#include <iostream>
#include <vector>

#include "dierckx.h"

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

    test_linear();
    test_quadratic();
    test_cubic();

    return 0;
}