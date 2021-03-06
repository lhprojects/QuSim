#pragma once

#include "QuSim.h"
#include "Linear.h"
#include "OptionsImpl.h"

struct IVPSolverImpl {

    IVPSolverImpl();
    virtual ~IVPSolverImpl() {}

    void initSystem(
        Real en,
        Real mass,
        Real hbar,
        SolverMethod met,
        OptionsImpl const &opts);

    virtual void Compute() = 0;

    OptionsImpl fOpts;
    Real const fE;
    Real const fMass;
    Real const fHbar;
    SolverMethod const fMethod;

};

struct IVPSolver1DImpl : IVPSolverImpl {

    IVPSolver1DImpl();

    void initSystem1D(
        std::function<Complex(Real)> const & v,
        Real x0,
        Real x1,
        size_t n,
        Real en,
        Complex initPsi,
        Complex initPsiPrime,
        SolverMethod met,
        Real mass,
        Real hbar,
        OptionsImpl const &opts);

    virtual void Compute() = 0;

    std::function<Complex(Real)> fVFunc;
    size_t fNPoints;
    size_t fNBins;
    Real fDx;
    Real fX0;
    Real fE;
    Real fInitJ;
    Real fFinalJ;
    PsiVector fPsi;
    PsiVector fPsiPrime;

    Real fT;
    Real fR;


};
