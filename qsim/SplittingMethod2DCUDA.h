#pragma once

#include "EvolverImpl.h"
EvolverImpl2D *CreateSplittingMethod2DCUDA(std::map<std::string, std::string> const &opts);
