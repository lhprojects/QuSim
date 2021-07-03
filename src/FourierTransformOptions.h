#pragma once

#include "FourierTransform.h"
#include "QuSim.h"
#include "Device.h"

struct FourierTransformOptions {

    FourierTransformLibrary const fLib;

    FourierTransformOptions() : fLib() {
    }

    void Init(OptionsImpl const &opts, DeviceType devType)
    {
        FourierTransformLibrary lib = FourierTransformLibrary::KISS;

        if (devType == DeviceType::GPU_CUDA) {
            std::string lib_str;
            if (opts.Get("fft_lib", lib_str)) {
                if (lib_str != "CUDA") {
                    throw std::invalid_argument("invalid fft_lib: " + lib_str);
                }
                lib = FourierTransformLibrary::CUDA;
            } else {
                lib = FourierTransformLibrary::CUDA;
            }
        } else {
            std::string lib_str;
            if (opts.Get("fft_lib", lib_str)) {
                if (lib_str == "KISS") {
                    lib = FourierTransformLibrary::KISS;
                } else if (lib_str == "FFTW") {
                    lib = FourierTransformLibrary::FFTW;
                } else {
                    throw std::invalid_argument("invalid fft_lib: " + lib_str);
                }
            }
        }
        const_cast<FourierTransformLibrary&>(fLib) = lib;
    }
};
