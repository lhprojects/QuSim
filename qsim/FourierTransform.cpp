#include "FourierTransform.h"
#include "kissfft.hh"

#if defined(_M_X64) 
#include "../fftw-3.3.5-dll64/fftw3.h"
#include "CUDAFourierTransform.h"
#endif

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

#if defined(_M_X64) 
struct FFTWFourierTransform : FourierTransform {
	fftw_plan plan;
	fftw_complex *in;
	fftw_complex *out;
	size_t fn;
	FFTWFourierTransform(size_t n, bool inverse)
	{
		in = fftw_alloc_complex(n);
		out = fftw_alloc_complex(n);
		plan = fftw_plan_dft_1d((int)n, in, out, inverse ? FFTW_BACKWARD : FFTW_FORWARD, FFTW_ESTIMATE);
		fn = n;
	}
	void Transform(std::complex<double> const *from, std::complex<double> *to) override
	{
		memcpy(in, from, sizeof(std::complex<double>)*fn);
		fftw_execute(plan);
		memcpy(to, out, sizeof(std::complex<double>)*fn);
	}

	~FFTWFourierTransform()
	{
		fftw_destroy_plan(plan);
		fftw_free(in);
		fftw_free(out);
	}
};
#endif

FourierTransform * FourierTransform::Create(size_t n, bool inverse, FourierTransformLibrary lib)
{
	//return CreateCUDAFourierTransform(n, inverse);
	if (lib == FourierTransformLibrary::KISS) {
		return new KissFourierTransform(n, inverse);
#if defined(_M_X64) 
	} else if (lib == FourierTransformLibrary::FFTW) {
		return new FFTWFourierTransform(n, inverse);
	} else if (lib == FourierTransformLibrary::CUDA) {
		return CreateCUDAFourierTransform(n, inverse);
#endif
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


struct FFTWFourierTransform2D : FourierTransform2D {

	fftw_plan plan;
	fftw_complex *in;
	fftw_complex *out;
	size_t fn1;
	size_t fn2;
	FFTWFourierTransform2D(size_t n1, size_t n2, bool inverse)
	{
		in = (fftw_complex*)fftw_malloc(n1*n2*sizeof(fftw_complex));
		out = (fftw_complex*)fftw_malloc(n1*n2 * sizeof(fftw_complex));
		// we are colum major, reverse n1 n2
		plan = fftw_plan_dft_2d((int)n2, (int)n1, in, out, inverse ? FFTW_BACKWARD : FFTW_FORWARD, FFTW_ESTIMATE);
		fn1 = n1;
		fn2 = n2;
	}

	void Transform(std::complex<double> const *from, std::complex<double> *to) override
	{
		//for (size_t i = 0; i < fn2*fn1; ++i) {
		//	(std::complex<double>&)in[i] = from[i];
		//}

		memcpy(in, from, sizeof(std::complex<double>)*fn1*fn2);
		fftw_execute(plan);
		memcpy(to, out, sizeof(std::complex<double>)*fn1*fn2);
		//for (size_t i = 0; i < fn2*fn1; ++i) {
		//	to[i] = (std::complex<double>&)out[i];
		//}

	}

	~FFTWFourierTransform2D()
	{
		fftw_destroy_plan(plan);
		fftw_free(in);
		fftw_free(out);
	}

};

FourierTransform2D * FourierTransform2D::Create(size_t n1, size_t n2, bool inverse, FourierTransformLibrary lib)
{
	if (lib == FourierTransformLibrary::KISS) {
		return new KissFourierTransform2D(n1, n2, inverse);
#if defined(_M_X64) 
	} else if (lib == FourierTransformLibrary::FFTW) {
		return new FFTWFourierTransform2D(n1, n2, inverse);
#endif
	} else {
		throw std::runtime_error("not supported fft library");
	}
}
