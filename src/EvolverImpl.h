#pragma once
//#define _CRT_SECURE_NO_WARNINGS
//#define _SCL_SECURE_NO_WARNINGS

#include "QuSim.h"
#include <functional>
#include "Linear.h"
#include "Device.h"
#include "OptionsImpl.h"
#include "Utils.h"


struct QuEvolverImpl
{
    using RealType = Real;
    using ComplexType = Complex;

    RealType* const fV = nullptr;
    RealType* const fVHost = nullptr;
    ComplexType* const fPsi = nullptr;
    ComplexType* const fPsiHost = nullptr;
    UInt fStep;

    DeviceType const fDeviceType = DeviceType::CPU_SEQ;
    std::unique_ptr<Device> const fDevice;
    Complex const fDt;
    // force normalize
    bool const fFNES = false;

    size_t const fN = 0; // Nx*Ny*Nz
    Real const fDv = 0;  // Dx*Dy*Dz
    Real const fMass = 0;
    Real const fHbar = 0;

    BoundaryCondition const fBoundaryCondition = BoundaryCondition::Period;
    SolverMethod const fSolverMethod = SolverMethod::Unknown;;

    size_t const fBatchSize = 1;
    bool const fFN = 0; // force normalization after initialization
    OptionsImpl const fOpts;

    void InitSystem(bool force_normalization,
        Complex dt, bool force_normalization_each_step,
        BoundaryCondition b, SolverMethod solver,
        Real mass, Real hbar,
        Real dv, size_t n,
        OptionsImpl const &opts);

    Real PotEn() const;
    Real KinEn() const;
    Real Norm2() const;
    Complex Time() const;

    virtual void Step();
    virtual void UpdatePsi() = 0;
    virtual Real CalPotEn() const;
    virtual Real CalKinEn() const = 0;

    virtual ~QuEvolverImpl();

};

inline Complex QuEvolverImpl::Time() const
{
    return 1.0 * fStep * fDt;
}

struct QuEvolver1DImpl : QuEvolverImpl
{

    Real const fX1 = 0;
    Real const fX0 = 0;
    Real const fDx = 0;
    size_t const fNx = 0;

    std::function<Complex(Real)> const fVFunc;
    std::function<Complex(Real)> const fPsi0Func;

    virtual void InitSystem1D(std::function<Complex(Real)> const & psi, bool force_normalization,
        Complex dt, bool force_normalization_each_step,
        std::function<Complex(Real)> const &v, Real x0, Real x1, size_t n,
        BoundaryCondition b, SolverMethod solver,
        Real mass, Real hbar, OptionsImpl const &opts);

    RealType GetX(size_t i) const;
    RealType Xavg() const;

    using QuEvolverImpl::Norm2;
    RealType Norm2(Real x0, Real x1) const;
private:
    void InitPsi();
    void InitPotential();
};

inline Real QuEvolver1DImpl::GetX(size_t i) const
{
    return fDx * i + fX0;
}



struct QuEvolver2DImpl : QuEvolverImpl {

    Real const fX0 = 0;
    Real const fDx = 0;
    Real const fY0 = 0;
    Real const fDy = 0;
    size_t const fNx = 0;
    size_t const fNy = 0;

    std::function<Complex(Real, Real)> fVFunc;
    std::function<Complex(Real, Real)> fPsi0Func;
    // col major
    // y major

    Real GetX(size_t i) const;
    Real GetY(size_t i) const;
    // x is major
    // y is minor
    ptrdiff_t Index(ptrdiff_t i, ptrdiff_t j) const;
    ptrdiff_t IndexFold(ptrdiff_t i, ptrdiff_t j) const;

    virtual void InitSystem2D(std::function<Complex(Real, Real)> const &psi, bool force_normalization,
        Complex dt, bool force_normalization_each_step,
        std::function<Complex(Real, Real)> const &v, Real x0, Real x1, size_t nx,
        Real y0, Real y1, size_t ny,
        BoundaryCondition b, SolverMethod solver,
        Real mass, Real hbar, OptionsImpl const &opts);

private:
    void InitPsi();
    void InitPotential();

};


inline Real QuEvolver2DImpl::GetX(size_t i) const
{
    return fDx * i + fX0;
}

inline Real QuEvolver2DImpl::GetY(size_t i) const
{
    return fDy * i + fY0;
}

inline ptrdiff_t QuEvolver2DImpl::Index(ptrdiff_t i, ptrdiff_t j) const
{
    return QuIdx(i, j, (ptrdiff_t)fNx, (ptrdiff_t)fNy);
}

inline ptrdiff_t QuEvolver2DImpl::IndexFold(ptrdiff_t i, ptrdiff_t j) const
{
    return QuIdxFold(i, j, (ptrdiff_t)fNx, (ptrdiff_t)fNy);
}

