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
`g++ -o main Main.cpp`

