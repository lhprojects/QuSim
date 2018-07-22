#pragma once

#include "FourierTransform.h"

FourierTransform *CreateCUDAFourierTransform(size_t n, bool inverse);
