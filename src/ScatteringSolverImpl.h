#pragma once

#include "QuSim.h"
#include "Linear.h"
#include "OptionsImpl.h"
#include "eigen/Eigen/Dense"
#include "Device.h"
#include "Utils.h"

struct ScatteringSolverImpl
{
    void InitScatteringSolver(Real en, Real k0, SolverMethod met, Real mass, Real hbar, size_t n,
        OptionsImpl const& opts);

    virtual void Compute() = 0;

    Real GetMomentum() const;

    using RealType = Real;
    using ComplexType = Complex;

    DeviceType const fDeviceType = DeviceType::CPU_SEQ;
    std::unique_ptr<Device> const fDevice;
    ComplexType const* const fV = nullptr;
    ComplexType const* const fVHost = nullptr;

    // spatial space psi0
    Complex const* const fPsi0X = nullptr;
    // spatial space delta psi
    Complex* const fPsiX = nullptr;

    Complex const* const fPsi0XHost = nullptr;
    Complex* const fPsiXHost = nullptr;

    size_t const fN = 0;
    Real const fE = 0;
    Real const fK0 = 0;
    SolverMethod const fMet = SolverMethod::Unknown;
    Real const fMass = 0;
    Real const fHbar = 0;
    OptionsImpl const fOpts;

    virtual ~ScatteringSolverImpl() {
        if (fDevice) {
            fDevice->SafeFree(mutable_ptr_cast(mutable_cast(fV)));
            fDevice->SafeFree(mutable_ptr_cast(mutable_cast(fPsiX)));
            fDevice->SafeFree(mutable_ptr_cast(mutable_cast(fPsi0X)));
            if (!fDevice->OnMainMem()) {
                SafeFree(mutable_ptr_cast(mutable_cast(fPsi0XHost)));
                SafeFree(mutable_ptr_cast(mutable_cast(fVHost)));
                SafeFree(mutable_ptr_cast(mutable_cast(fPsiXHost)));
            }
        }
    }
};

inline Real ScatteringSolverImpl::GetMomentum() const { 
    return sqrt(2 * fMass * fE);
}


struct ScatteringSolver1DImpl : ScatteringSolverImpl {

    virtual void InitScatteringSolver1D(
        std::function<Complex(Real)> const& v,
        Real x0,
        Real x1,
        size_t n,
        Real en,
        Real k0,
        Real direction,
        SolverMethod met,
        Real mass,
        Real hbar,
        OptionsImpl const &opts);

    Real GetX(ptrdiff_t i) const { return fX0 + fDx * i; }
    void ComputeRT();

    size_t const fNx = 0;
    Real const fX0 = 0;
    Real const fDx = 0;
    std::function<Complex(Real)> const fVFunc;
    Real const fDirection = 0;
    Real const fK0X = 0;

    Real fT;
    Real fR;
private:
    void InitPotential();
    void InitPsi();
};

struct ScatteringSolver2DImpl : ScatteringSolverImpl {

    virtual void InitScatteringSolver2D(
        std::function<Complex(Real, Real)> const & v,
        Real x0,
        Real x1,
        size_t nx,
        Real y0,
        Real y1,
        size_t ny,
        Real en,
        Real k0,
        Real directionx,
        Real directiony,
        SolverMethod met,
        Real mass,
        Real hbar,
        OptionsImpl const &opts);

    Real GetX(ptrdiff_t i) const { return fX0 + fDx * i; }
    Real GetY(ptrdiff_t i) const { return fY0 + fDy * i; }

    virtual Real ComputeXSection(Real cosx, Real cosy);
    virtual Real ComputeTotalXSection(Int n);

    size_t const fNx = 0;
    size_t const fNy = 0;
    Real const fX0 = 0;
    Real const fY0 = 0;
    Real const fDx = 0;
    Real const fDy = 0;
    std::function<Complex(Real, Real)> const fVFunc;
    Real const fK0X = 0;
    Real const fK0Y = 0;
    
    const bool fLastPsiXRecord = false;
    ComplexType* const fLastPsiX = nullptr;
    Real fNormDeltaPsi;
private:
    void InitPotential();
    void InitPsi();
};

struct ScatteringSolver3DImpl : ScatteringSolverImpl {

    virtual void InitScatteringSolver3D(
        std::function<Complex(Real, Real, Real)> const & v,
        Real x0,
        Real x1,
        size_t nx,
        Real y0,
        Real y1,
        size_t ny,
        Real z0,
        Real z1,
        size_t nz,
        Real en,
        Real k0,
        Real directionx,
        Real directiony,
        Real directionz,
        SolverMethod met,
        Real mass,
        Real hbar,
        OptionsImpl const &opts);

    Real GetX(ptrdiff_t i) const { return fX0 + fDx * i; }
    Real GetY(ptrdiff_t i) const { return fY0 + fDy * i; }
    Real GetZ(ptrdiff_t i) const { return fZ0 + fDz * i; }

    // cosx = cos(psi)
    // cosy = sin(psi)
    // cosz = cos(theta)
    virtual Real ComputeXSection(Real cosx, Real cosy, Real cosz);
    // npsi = number of sampling points for psi
    // ntheta = number of sampling points for theta
    virtual Real ComputeTotalXSection(Int npsi, Int ntheta);

    size_t const fNx = 0;
    size_t const fNy = 0;
    size_t const fNz = 0;
    Real const fX0 = 0;
    Real const fY0 = 0;
    Real const fZ0 = 0;
    Real const fDx = 0;
    Real const fDy = 0;
    Real const fDz = 0;
    std::function<Complex(Real, Real, Real)> const fVFunc;
    Real const fK0X = 0;
    Real const fK0Y = 0;
    Real const fK0Z = 0;

    bool fVComputeInTime;
    bool fPsi0XComputeInTime;
    static const bool fLastPsiXRecord = false;
    ComplexType* const fLastPsiX = nullptr;
private:
    void InitPotential();
    void InitPsi();
};
