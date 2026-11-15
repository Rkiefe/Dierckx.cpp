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

Interpolate `f = x^2 + y^3` with a 2D spline
```cpp
// Initialize 3 vectors
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
// Initialize the Spline
int kx = 2, ky = 3;
Spline2D spl(x, y, z, kx, ky);
double zq = spl(0.5, 0.5); // outputs: 0.375, theoretical: 0.375
```