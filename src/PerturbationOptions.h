#pragma once

#include <map>
#include <string>
#include "Perburbation.h"
#include "QuSim.h"



struct PerturbationCommon {
    bool fPreconditional = false;
    BornSerisePreconditioner const fPreconditioner = BornSerisePreconditioner::Identity;
    Real const fSlow = 0;
        
    int const fOrder = 0;
    Real const fEpsilon = 0;

    Complex* const fPsiK = nullptr;
    Complex* const fTmpPsi = nullptr;
    Device* const fDev = nullptr;

    void InitPerturbationCommon(OptionsImpl const& opts, Real epsilon, size_t n, Device* dev);
    ~PerturbationCommon();
};
