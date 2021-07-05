#include "EvolverImpl.h"
#include "Device.h"

void QuEvolver1DImpl::InitSystem1D(std::function<Complex(Real)> const& psi, bool force_normalization,
    Complex dt, bool force_normalization_each_step,
    std::function<Complex(Real)> const& v, Real x0, Real x1, size_t n,
    BoundaryCondition b, SolverMethod solver,
    Real mass, Real hbar, OptionsImpl const& opts)
{

    auto dx = (x1 - x0) / n;
    InitSystem(force_normalization,
        dt, force_normalization_each_step,
        b, solver,
        mass, hbar,
        dx, n,
        opts);

    mutable_cast(fPsi0Func) = psi;
    mutable_cast(fVFunc) = v;

    mutable_cast(fNx) = n;
    mutable_cast(fX0) = x0;
    mutable_cast(fX1) = x1;
    mutable_cast(fDx) = dx;

    InitPsi();
    InitPotential();

}


void QuEvolver1DImpl::InitPotential()
{
    mutable_cast(fV) = fDevice->Alloc<RealType>(fN);

    if (!fDevice->OnMainMem()) {
        mutable_cast(fVHost) = (RealType*)malloc(sizeof(RealType) * fN);
    } else {
        mutable_cast(fVHost) = fV;
    }

    // initialization on cpu
    for (size_t i = 0; i < fNx; ++i) {
        Real x = GetX(i);
        fVHost[i] = QuGetReal(fVFunc(x));
    }

    if (!fDevice->OnMainMem()) {
        fDevice->ToDevice(fV, fVHost, fNx);
    }
}

void QuEvolver1DImpl::InitPsi()
{
    mutable_cast(fPsi) = fDevice->Alloc<ComplexType>(fN);

    if (!fDevice->OnMainMem()) {
        mutable_cast(fPsiHost) = (ComplexType*)malloc(sizeof(ComplexType) * fN);
    } else {
        mutable_cast(fPsiHost) = fPsi;
    }

    // initialization on cpu
    for (size_t i = 0; i < fN; ++i) {
        Real x = GetX(i);
        fPsiHost[i] = fPsi0Func(x);
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

QuEvolver1DImpl::RealType
QuEvolver1DImpl::Xavg() const
{
    RealType Idx = fDevice->Abs2Idx(fPsi, fN) / fDevice->Norm2(fPsi, fN);
    return fX0 + fDx * Idx;
}

QuEvolver1DImpl::RealType
QuEvolver1DImpl::Norm2(Real x0, Real x1) const
{
    Int Idx1 = Int((x0 - fX0) / fDx);
    Int Idx2 = Int((x1 - fX0) / fDx);
    if (Idx1 <= 0) Idx1 = 0;
    if (Idx2 >= Int(fN)) Idx2 = fN;

    if (Idx1 >= Idx2) {
        return 0;
    }
    return fDevice->Norm2(fPsi + Idx1, Idx2 - Idx1) * fDv;
}
