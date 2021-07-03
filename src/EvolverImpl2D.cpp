#include "EvolverImpl.h"
#include "Device.h"


void QuEvolver2DImpl::InitSystem2D(std::function<Complex(Real, Real)> const& psi, bool force_normalization,
    Complex dt, bool force_normalization_each_step,
    std::function<Complex(Real, Real)> const& v,
    Real x0, Real x1, size_t nx,
    Real y0, Real y1, size_t ny,
    BoundaryCondition b, SolverMethod solver,
    Real mass, Real hbar, OptionsImpl const& opts)
{
    auto dx = (x1 - x0) / nx;
    auto dy = (y1 - y0) / ny;

    InitSystem(force_normalization,
        dt, force_normalization_each_step,
        b, solver,
        mass, hbar,
        dx * dy, nx * ny,
        opts);

    mutable_cast(fPsi0Func) = psi;
    mutable_cast(fVFunc) = v;

    mutable_cast(fX0) = x0;
    mutable_cast(fDx) = dx;
    mutable_cast(fNx) = nx;

    mutable_cast(fY0) = y0;
    mutable_cast(fDy) = dy;
    mutable_cast(fNy) = ny;

    InitPsi();
    InitPotential();

}

void QuEvolver2DImpl::InitPotential()
{
    mutable_cast(fV) = fDevice->Alloc<RealType>(fN);

    if (!fDevice->OnMainMem()) {
        mutable_cast(fVHost) = (RealType*)malloc(sizeof(ComplexType) * fN);
    } else {
        mutable_cast(fVHost) = fV;
    }

    // initialization on cpu
    for (size_t i = 0; i < fNx; ++i) {
        for (size_t j = 0; j < fNy; ++j) {
            Real x = GetX(i);
            Real y = GetY(j);
            size_t idx = Index(i, j);
            fVHost[idx] = QuGetReal(fVFunc(x, y));
            
        }
    }

    if (!fDevice->OnMainMem()) {
        fDevice->ToDevice(fV, fVHost, fN);
    }

}

void QuEvolver2DImpl::InitPsi()
{

    mutable_cast(fPsi) = fDevice->Alloc<ComplexType>(fN);

    if (!fDevice->OnMainMem()) {
        mutable_cast(fPsiHost) = (ComplexType*)malloc(sizeof(ComplexType) * fN);
    } else {
        mutable_cast(fPsiHost) = fPsi;
    }

    // initialization on cpu
    for (size_t i = 0; i < fNx; ++i) {
        for (size_t j = 0; j < fNy; ++j) {
            Real x = GetX(i);
            Real y = GetY(j);
            size_t idx = Index(i, j);
            fPsiHost[idx] = fPsi0Func(x, y);
        }
    }

    if (fFN) {
        auto dseq = Device::Create(DeviceType::CPU_SEQ);
        auto norm2 = dseq->Norm2(fPsiHost, fN) * fDv;
        auto scale = 1. / sqrt(norm2);
        dseq->Scale(fPsiHost, scale, fN);
    }

    if (!fDevice->OnMainMem()) {
        fDevice->ToDevice(fPsi, fPsiHost, fN);
    }

}

