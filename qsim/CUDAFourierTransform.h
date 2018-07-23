#pragma once
#ifdef USE_CUDA

#include "FourierTransform.h"
FourierTransform *CreateCUDAFourierTransform(size_t n, bool inverse);

#endif

