#pragma once

#include "FourierTransform.h"
#include "QuSim.h"

struct FourierTransformOptions {

	FourierTransformLibrary const fLib;

	FourierTransformOptions() : fLib() {
	}

	void Init(OptionsImpl const &opts)
	{
		{
			FourierTransformLibrary lib = FourierTransformLibrary::KISS;
			std::string lib_str;

			if (opts.Get("fft_lib", lib_str)) {
				if (lib_str == "kiss") {
					lib = FourierTransformLibrary::KISS;
				} else if (lib_str == "FFTW") {
					lib = FourierTransformLibrary::FFTW;
				} else if (lib_str == "cuda") {
					lib = FourierTransformLibrary::CUDA;
				} else {
					throw std::runtime_error("unkown lib");
				}
			}
			const_cast<FourierTransformLibrary&>(fLib) = lib;
		}
	}
};
