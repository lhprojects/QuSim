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
	if (lib == FourierTransformLibrary::KISS) {
		return new KissFourierTransform(n, inverse);
	} else {
		throw std::runtime_error("not supported fft library");
	}
}
