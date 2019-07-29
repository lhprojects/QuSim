#include "FourierTransform.h"
#include "kissfft.hh"

struct KissFourierTransform : FourierTransform {
	kissfft<double> impl;

	KissFourierTransform(size_t n, bool inverse) : impl(n, inverse)
	{

	}
	void Transform(std::complex<double> const *from, std::complex<double> *to) override
	{
		impl.transform(from, to);
	}

};

FourierTransform * FourierTransform::Create(size_t n, bool inverse, FourierTransformLibrary lib)
{
	//return CreateCUDAFourierTransform(n, inverse);
	if (lib == FourierTransformLibrary::KISS) {
		return new KissFourierTransform(n, inverse);
	} else {
		throw std::runtime_error("not supported fft library");
	}
}

struct KissFourierTransform2D : FourierTransform2D {
	kissfft<double> impl1;
	kissfft<double> impl2;
	std::vector<std::complex<double> > tmp;
	std::vector<std::complex<double> > tmp2;
	std::vector<std::complex<double> > tmp3;
	size_t rows;
	size_t cols;
	KissFourierTransform2D(size_t rows, size_t cols, bool inverse) : impl1(rows, inverse), impl2(cols, inverse)
	{
		this->rows = rows;
		this->cols = cols;
		tmp.resize(rows*cols);
		tmp2.resize(cols);
		tmp3.resize(cols);
	}
	void Transform(std::complex<double> const *from, std::complex<double> *to) override
	{

		std::complex<double> *t = tmp.data();
		std::complex<double> *t2 = tmp2.data();
		std::complex<double> *t3 = tmp3.data();

		for (size_t i = 0; i < cols; ++i) {
			impl1.transform(from + i * rows, t + i * rows);
		}

		for (size_t i = 0; i < rows; ++i) {
			for (size_t j = 0; j < cols; ++j) {
				t2[j] = t[j*rows + i];
			}
			impl2.transform(t2, t3);
			for (size_t j = 0; j < cols; ++j) {
				to[j*rows + i] = t3[j];
			}
		}

	}

	~KissFourierTransform2D()
	{
	}

};


FourierTransform2D * FourierTransform2D::Create(size_t rows, size_t cols, bool inverse, FourierTransformLibrary lib)
{
	if (lib == FourierTransformLibrary::KISS) {
		return new KissFourierTransform2D(rows, cols, inverse);
	} else {
		throw std::runtime_error("not supported fft library");
	}
}

struct KissFourierTransform3D : FourierTransform3D {
	kissfft<double> impl1;
	kissfft<double> impl2;
	std::vector<std::complex<double> > tmp;
	std::vector<std::complex<double> > tmp2;
	std::vector<std::complex<double> > tmp3;
	size_t nz;
	size_t rows;
	size_t cols;

	KissFourierTransform3D(size_t nz, size_t rows, size_t cols, bool inverse) : impl1(rows, inverse), impl2(cols, inverse)
	{
		this->rows = rows;
		this->cols = cols;
		this->nz = nz;
		tmp.resize(rows*cols);
		tmp2.resize(cols);
		tmp3.resize(cols);
	}
	void Transform(std::complex<double> const *from, std::complex<double> *to) override
	{

		std::complex<double> *t = tmp.data();
		std::complex<double> *t2 = tmp2.data();
		std::complex<double> *t3 = tmp3.data();

		for (size_t i = 0; i < cols; ++i) {
			impl1.transform(from + i * rows, t + i * rows);
		}

		for (size_t i = 0; i < rows; ++i) {
			for (size_t j = 0; j < cols; ++j) {
				t2[j] = t[j*rows + i];
			}
			impl2.transform(t2, t3);
			for (size_t j = 0; j < cols; ++j) {
				to[j*rows + i] = t3[j];
			}
		}

	}

	~KissFourierTransform3D()
	{
	}

};



FourierTransform3D * FourierTransform3D::Create(size_t nz, size_t ny, size_t nx, bool inverse, FourierTransformLibrary lib)
{
	return new KissFourierTransform3D(nz, ny, nx, inverse);
}
