#pragma once

#include "FourierTransform.h"
#include "QuSim.h"

struct FourierTransformOptions {

	FourierTransformLibrary const fLib;

	FourierTransformOptions() : fLib() {
	}

	void Init(std::map<std::string, std::string> const &opts)
	{
		{
			FourierTransformLibrary lib = FourierTransformLibrary::KISS;
			auto it = opts.find("fft_lib");
			if (it != opts.end()) {
				if (it->second == "kiss") {
					lib = FourierTransformLibrary::KISS;
				} else if (it->second == "FFTW") {
					lib = FourierTransformLibrary::FFTW;
				} else if (it->second == "cuda") {
					lib = FourierTransformLibrary::CUDA;
				}
			}
			const_cast<FourierTransformLibrary&>(fLib) = lib;
		}
	}
};
