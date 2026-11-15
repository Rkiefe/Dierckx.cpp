## A modern C++ port of Fortran Dierckx Spline and Interpolation methods.

### Examples

Interpolate a cubic function with a linear fit
```cpp
std::vector<double> x(100);
std::vector<double> y(100);

for(uint i = 0; i<x.size(); i++){
    x[i] = (double) i/x.size();
    y[i] = pow(x[i], 3);
}

Spline1D spl(x, y, /*k = */ 1); // Create a linear fit | k = 1, linear; k = 2, quadratic; k = 3, cubic
double yq = spl(0.5); // Evaluate the spline on target point
```

Compile with:
- **note:** if you use the header-only version:
    `g++ -o main Main.cpp`

- **note:** if you want separate files: 
    `g++ -std=c++17 -Wall -O2 -o main Main.cpp spline1d.cpp`

