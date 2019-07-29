#pragma once

#include <complex>

enum class FourierTransformLibrary {
	KISS,
	CUDA,
};

struct FourierTransform {
	// factory
	static FourierTransform *Create(size_t n, bool inverse, FourierTransformLibrary lib);
	virtual void Transform(std::complex<double> const *from, std::complex<double> *to) = 0;
	virtual ~FourierTransform(){ }
};

struct FourierTransform2D {
	// factory
	static FourierTransform2D *Create(size_t rows, size_t cols, bool inverse, FourierTransformLibrary lib);
	virtual void Transform(std::complex<double> const *from, std::complex<double> *to) = 0;
	virtual ~FourierTransform2D() { }
};

struct FourierTransform3D {
	// factory
	static FourierTransform3D *Create(size_t nz, size_t ny, size_t nx, bool inverse, FourierTransformLibrary lib);
	virtual void Transform(std::complex<double> const *from, std::complex<double> *to) = 0;
	virtual ~FourierTransform3D() {}
};
