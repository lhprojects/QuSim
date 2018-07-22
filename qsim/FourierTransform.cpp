#include "FourierTransform.h"
#include "kissfft.hh"
#if defined(_M_X64) 
#include "../fftw-3.3.5-dll64/fftw3.h"
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
		fftw_free(in);
		fftw_free(out);
	}
};
#endif

FourierTransform * FourierTransform::Create(size_t n, bool inverse, FourierTransformLibrary lib)
{
	if (lib == FourierTransformLibrary::KISS) {
		return new KissFourierTransform(n, inverse);
#if defined(_M_X64) 
	} else if (lib == FourierTransformLibrary::FFTW) {
		return new FFTWFourierTransform(n, inverse);
#endif
	} else {
		throw std::runtime_error("not supported fft library");
	}
}
