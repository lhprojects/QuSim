#include "EvolverImpl.h"
#include "Device.h"
#include "Utils.h"
#include "FourierTransformOptions.h"
#include "DeviceType.h"

void QuEvolverImpl::InitSystem(bool force_normalization,
    Complex dt, bool force_normalization_each_step,
    BoundaryCondition b,
    SolverMethod solver,
    Real mass,
    Real hbar,
    Real dv,
    size_t n,
    OptionsImpl const &opts)
{
    fStep = 0;
    
    mutable_cast(fBatchSize) = opts.GetInt("batch", 1);
    mutable_cast(fHbar) = hbar;
    mutable_cast(fFN) = force_normalization;
    mutable_cast(fBoundaryCondition) = b;
    mutable_cast(fDv) = dv;
    mutable_cast(fDt) = dt;
    mutable_cast(fN) = n;
    mutable_cast(fFNES) = force_normalization_each_step;
    mutable_cast(fMass) = mass;
    mutable_cast(fSolverMethod) = solver;
    mutable_cast(fOpts) = opts;
    
    mutable_cast(fDeviceType) = GetDeviceType(fOpts);
    mutable_cast(fDevice) = Device::Create(fDeviceType);
}

Real QuEvolverImpl::PotEn() const
{
    return CalPotEn();
}

Real QuEvolverImpl::KinEn() const
{
    return CalKinEn();
}

Real QuEvolverImpl::CalPotEn() const
{
    Complex norm2 = fDevice->Abs2Mul(fPsi, fV, fN) * fDv;
    return norm2.real() / Norm2();
}

Real QuEvolverImpl::Norm2() const
{
    return fDevice->Norm2(fPsi, fN) * fDv;
}

void QuEvolverImpl::Step()
{
    UpdatePsi();

    if (fFNES) {
        double scale = 1.0 / sqrt(Norm2());
        fDevice->Scale(fPsi, scale, fN);
    }

    if (!fDevice->OnMainMem()) {
        fDevice->ToHost(fPsiHost, fPsi, fN);
    }
}

QuEvolverImpl::~QuEvolverImpl()
{
    if (fDevice) {
        fDevice->SafeFree(mutable_cast(fPsi));
        fDevice->SafeFree(mutable_cast(fV));
        if (!fDevice->OnMainMem()) {
            free(fPsiHost);
            free(fVHost);
        }
    }
}

