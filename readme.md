## A modern C++ port of Fotran Dierckx Spline and Interpolation methods.

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

// Evaluate the spline over a query range
std::vector<double> xq = {0.0, 0.25, 0.333, 0.5, 0.666, 0.75, 1.0};
for(uint i = 0; i<xq.size(); i++){
    std::cout << "f = " << pow(xq[i], 3) << " | spl = " << spl(xq[i]) << std::endl;
}

```

Compile with:
`g++ -std=c++17 -Wall -O2 -o main Main.cpp spline1d.cpp`

