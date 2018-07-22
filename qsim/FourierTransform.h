#pragma once

#include <complex>

enum class FourierTransformLibrary {
	KISS,
	FFTW,
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
