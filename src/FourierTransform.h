#pragma once

#include <complex>
#include <QuSim.h>

enum class FourierTransformLibrary {
    KISS,
    FFTW,
    CUDA,
};

struct FourierTransform1D {
    // factory
    static FourierTransform1D *Create(size_t n, bool inverse, FourierTransformLibrary lib);

    // Fourier Transform definition:
    // X_k = sum_n x_n Exp{-2 pi I k n/N}
    // Inverse Fourier Transform definition:
    // x_n = sum_n X_k Exp{+2 pi I k n/N}
    virtual void Transform(Complex64 const* from, Complex64* to) = 0;
    virtual ~FourierTransform1D() {}
};

struct FourierTransform2D {
    // factory
    static FourierTransform2D* Create(size_t rows, size_t cols, bool inverse, FourierTransformLibrary lib);
    virtual void Transform(Complex const* from, Complex* to) = 0;
    virtual ~FourierTransform2D() {}
};

struct FourierTransform3D {
    // factory
    static FourierTransform3D *Create(size_t nz, size_t ny, size_t nx, bool inverse, FourierTransformLibrary lib);
    virtual void Transform(Complex const *from, Complex*to) = 0;
    virtual ~FourierTransform3D() {}
};
