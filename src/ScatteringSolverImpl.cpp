#include "ScatteringSolverImpl.h"
#include "Utils.h"
#include "DeviceType.h"
void ScatteringSolverImpl::InitScatteringSolver(Real en, SolverMethod met, Real mass, Real hbar, size_t n,
    OptionsImpl const& opts)
{

    const_cast<Real&>(fE) = en;
    const_cast<SolverMethod&>(fMet) = met;
    const_cast<Real&>(fMass) = mass;
    const_cast<Real&>(fHbar) = hbar;
    const_cast<OptionsImpl&>(fOpts) = opts;
    mutable_cast(fN) = n;
    const_cast<Real&>(fK0) = sqrt(2 * fMass * fE) / fHbar;

    mutable_cast(fDeviceType) = GetDeviceType(opts);
    mutable_cast(fDevice) = Device::Create(fDeviceType);
}

void ScatteringSolver1DImpl::InitScatteringSolver1D(std::function<Complex(Real)> const & v,
    Real x0, Real x1, size_t n, Real en, Real direction,
    SolverMethod met, Real mass, Real hbar, OptionsImpl const & opts)
{
    InitScatteringSolver(en, met, mass, hbar, n, opts);

    mutable_cast(fN) = n;
    const_cast<size_t&>(fNx) = n;
    const_cast<Real&>(fX0) = x0;
    const_cast<Real&>(fDx) = (x1 - x0) / n;
    const_cast<std::function<Complex(Real)>&>(fVFunc) = v;
    if (direction < 0) direction = -1;
    else if (direction > 0) direction = 1;
    const_cast<Real&>(fDirection) = direction > 0 ? 1 : -1;
    const_cast<Real&>(fK0) = sqrt(2 * fMass * fE) / fHbar;

    InitPotential();
    InitPsi();
}

void ScatteringSolver1DImpl::ComputeRT() {
    // calculate in real space

    Complex r = 0;
    Complex t = 0;
    for (size_t i = 0; i < fNx; ++i) {
        r += (fPsi0X[i] + fPsiX[i])*fV[i].real() * exp(I * (-fDirection*fK0)*(-GetX(i)));
        t += (fPsi0X[i] + fPsiX[i])*fV[i].real() * exp(I * (+fDirection*fK0)*(-GetX(i)));
    }
    fR = abs2(r*fDx*fMass / (fHbar*fHbar*std::abs(fK0) * I));
    fT = abs2(t*fDx*fMass / (fHbar*fHbar*std::abs(fK0) * I) + Complex(1, 0));
}

void ScatteringSolver1DImpl::InitPotential()
{
    mutable_cast(fV) = fDevice->Alloc<Complex>(fNx);

    for_each_global_idx(fX0, fDx, fNx,
        [&, v = mutable_ptr_cast(fV)](size_t idx, Real x) {
        v[idx] = fVFunc(x);
    });
}

void ScatteringSolver1DImpl::InitPsi()
{
    mutable_cast(fPsi0X) = fDevice->Alloc<Complex>(fNx);
    mutable_cast(fPsiX) = fDevice->Alloc<Complex>(fNx);

    for_each_global_idx(fX0, fDx, fNx,
        [psi0 = mutable_ptr_cast(fPsi0X), dir=fDirection, kx = fK0](size_t idx, Real x) {
        psi0[idx] = exp(dir *kx * x * I);
    });

}




void ScatteringSolver2DImpl::InitScatteringSolver2D(std::function<Complex(Real, Real)> const & v,
    Real x0, Real x1, size_t nx,
    Real y0, Real y1, size_t ny,
    Real en,
    Real directionx,
    Real directiony,
    SolverMethod met, Real mass, Real hbar, OptionsImpl const & opts)
{
    InitScatteringSolver(en, met, mass, hbar, nx * ny, opts);

    const_cast<size_t&>(fNx) = nx;
    const_cast<size_t&>(fNy) = ny;
    const_cast<Real&>(fX0) = x0;
    const_cast<Real&>(fY0) = y0;
    const_cast<Real&>(fDx) = (x1 - x0) / nx;
    const_cast<Real&>(fDy) = (y1 - y0) / ny;
    const_cast<std::function<Complex(Real, Real)>&>(fVFunc) = v;

    {
        Real const nm = 1 / sqrt(directionx*directionx + directiony * directiony);
        directionx *= nm;
        directiony *= nm;
    }

    mutable_cast(fK0X) = fK0 * directionx;
    mutable_cast(fK0Y) = fK0 * directiony;

    InitPotential();
    InitPsi();

    mutable_cast(fLastPsiX) = fDevice->Alloc<ComplexType>(fN);
    fNormDeltaPsi = 0;
}

void ScatteringSolver2DImpl::InitPotential()
{
    mutable_cast(fV) = fDevice->Alloc<ComplexType>(fN);
    if (fDevice->OnMainMem()) {
        mutable_cast(fVHost) = fV;
    } else {
        mutable_cast(fVHost) = (ComplexType*)malloc(sizeof(ComplexType) * fN);
    }

    for_each_global_idx(fX0, fDx, fNx,
        fY0, fDy, fNy,
        [&, v = mutable_ptr_cast(fVHost)](size_t idx, Real x, Real y) {
        v[idx] = fVFunc(x, y);
    });

    if (!fDevice->OnMainMem()) {
        fDevice->ToDevice(mutable_ptr_cast(fV), fVHost, fN);
    }
}

void ScatteringSolver2DImpl::InitPsi()
{
    mutable_cast(fPsi0X) = fDevice->Alloc<ComplexType>(fN);
    mutable_cast(fPsiX) = fDevice->Alloc<ComplexType>(fN);

    if (!fDevice->OnMainMem()) {
        mutable_cast(fPsi0XHost) = (ComplexType*)malloc(sizeof(ComplexType) * fN);
        mutable_cast(fPsiXHost) = (ComplexType*)malloc(sizeof(ComplexType) * fN);
    } else {
        mutable_cast(fPsi0XHost) = fPsi0X;
        mutable_cast(fPsiXHost) = fPsiX;
    }

    for_each_global_idx(fX0, fDx, fNx,
        fY0, fDy, fNy,
        [psi0 = mutable_ptr_cast(fPsi0XHost), kx = fK0X, ky = fK0Y]
    (size_t idx, Real x, Real y) {
        psi0[idx] = exp((kx * x + ky * y) * I);
    });

    memset(fPsiXHost, 0, sizeof(ComplexType) * fN);

    if (!fDevice->OnMainMem()) {
        fDevice->ToDevice(mutable_ptr_cast(fPsi0X), fPsi0XHost, fN);
        fDevice->ToDevice(fPsiX, fPsiXHost, fN);
    }
}


Real ScatteringSolver2DImpl::ComputeXSection(Real cosx, Real cosy)
{
    Real const nm = 1 / sqrt(cosx*cosx + cosy * cosy);
    cosx *= nm;
    cosy *= nm;

    Complex r = 0;

    // Prolbem: (nabla^2+k0^2) delta psi = 2m/hbar^2 v psi
    // Green's: (nabla^2+k0^2) G = delta(x)
    //     => : delta psi  = \int G (2m / hbar^2) v psi dx dy
    // G = 1/4 (1st kind Hankel, alpha = 0) (k r)
    // (1st kind Hankel, alpha = 0)(z) ~ sqrt(2/(Pi z)) exp(i k z)
    //     => : delta psi  = \int psi v exp(i k r)  [(2m / hbar^2)  1/4 sqrt(2/(Pi k r))] dx dy

    Real const kx = fK0 * cosx;
    Real const ky = fK0 * cosy;

    ForEach2DGrid(fX0, fDx, fNx,
        fY0, fDy, fNy,
        [&, kx, ky](size_t idx, size_t, size_t, Real x, Real y) {

            r += (fPsi0X[idx] + fPsiX[idx]) * fV[idx].real() * exp(-I * (kx * x + ky * y));

    });
    Real dX_dTheta = abs2(r*(2 * fMass) / (fHbar*fHbar)*fDx*fDy*(1. / 4))*(2 / (Pi * fK0));

    return dX_dTheta;

}

Real ScatteringSolver2DImpl::ComputeTotalXSection(Int n)
{
    Real xsec_t = 0;
    for (Int i = 0; i < n; ++i) {
        Real theta = i * 2 * Pi / n;
        Real cosx = cos(theta);
        Real cosy = sin(theta);
        Real xsec = ComputeXSection(cosx, cosy);
        xsec_t += xsec;
    }
    xsec_t *= 2 * Pi / n;

    return xsec_t;
}




void ScatteringSolver3DImpl::InitScatteringSolver3D(std::function<Complex(Real, Real, Real)> const & v, Real x0, Real x1, size_t nx,
    Real y0, Real y1, size_t ny, Real z0, Real z1, size_t nz,
    Real en,
    Real directionx, Real directiony, Real directionz,
    SolverMethod met, Real mass, Real hbar, OptionsImpl const & opts)
{
    InitScatteringSolver(en, met, mass, hbar, nx*ny*nz, opts);

    const_cast<size_t&>(fNx) = nx;
    const_cast<size_t&>(fNy) = ny;
    const_cast<size_t&>(fNz) = nz;
    const_cast<Real&>(fX0) = x0;
    const_cast<Real&>(fY0) = y0;
    const_cast<Real&>(fZ0) = z0;
    const_cast<Real&>(fDx) = (x1 - x0) / nx;
    const_cast<Real&>(fDy) = (y1 - y0) / ny;
    const_cast<Real&>(fDz) = (z1 - z0) / nz;
    const_cast<std::function<Complex(Real, Real, Real)>&>(fVFunc) = v;

    {
        Real const nm = 1 / sqrt(directionx * directionx + directiony * directiony + directionz * directionz);
        directionx *= nm;
        directiony *= nm;
        directionz *= nm;
    }
    const_cast<Real&>(fK0X) = fK0 * directionx;
    const_cast<Real&>(fK0Y) = fK0 * directiony;
    const_cast<Real&>(fK0Z) = fK0 * directionz;

    InitPotential();
    InitPsi();
}

void ScatteringSolver3DImpl::InitPotential()
{
    mutable_cast(fV) = fDevice->Alloc<ComplexType>(fN);
    if (fDevice->OnMainMem()) {
        mutable_cast(fVHost) = fV;
    } else {
        mutable_cast(fVHost) = (ComplexType*)malloc(sizeof(ComplexType) * fN);
    }

    for_each_global_idx(fX0, fDx, fNx,
        fY0, fDy, fNy,
        fZ0, fDz, fNz,
        [&, v = mutable_ptr_cast(fVHost)](size_t idx, Real x, Real y, Real z) {
            v[idx] = fVFunc(x, y, z);
    });

    if (!fDevice->OnMainMem()) {
        fDevice->ToDevice(mutable_ptr_cast(fV), fVHost, fN);
    }
}

void ScatteringSolver3DImpl::InitPsi()
{
#if 1
    mutable_cast(fPsi0X) = fDevice->Alloc<ComplexType>(fN);
#endif
    mutable_cast(fPsiX) = fDevice->Alloc<ComplexType>(fN);

    if (!fDevice->OnMainMem()) {
#if 1
        mutable_cast(fPsi0XHost) = (ComplexType*)malloc(sizeof(ComplexType) * fN);
#endif
        mutable_cast(fPsiXHost) = (ComplexType*)malloc(sizeof(ComplexType) * fN);
    } else {
#if 1
        mutable_cast(fPsi0XHost) = fPsi0X;
#endif
        mutable_cast(fPsiXHost) = fPsiX;
    }

#if 1
    for_each_global_idx(fX0, fDx, fNx,
        fY0, fDy, fNy,
        fZ0, fDz, fNz,
        [&, psi0 = mutable_ptr_cast(fPsi0XHost), kx = fK0X, ky = fK0Y, kz = fK0Z]
        (size_t idx, Real x, Real y, Real z) {
        psi0[idx] = exp((kx*x + ky*y + kz*z)*I);
    });
#endif

    memset(fPsiXHost, 0, sizeof(ComplexType) * fN);

    if (!fDevice->OnMainMem()) {
#if 1
        fDevice->ToDevice(mutable_ptr_cast(fPsi0X), fPsi0XHost, fN);
#endif
        fDevice->ToDevice(fPsiX, fPsiXHost, fN);
    }

    if (fLastPsiXRecord) {
        mutable_cast(fLastPsiX) = fDevice->Alloc<ComplexType>(fN);
        fDevice->SetZero(mutable_ptr_cast(fLastPsiX), fN);
    }
}



Real ScatteringSolver3DImpl::ComputeXSection(Real cosx, Real cosy, Real cosz)
{
    Real const nm = 1 / sqrt(cosx * cosx + cosy * cosy + cosz * cosz);
    cosx *= nm;
    cosy *= nm;
    cosz *= nm;

    Complex r = 0;

    // Prolbem: (nabla^2+k0^2) delta psi = 2m/hbar^2 v psi
    // Green's: (nabla^2+k0^2) G = delta(x)
    //     => : delta psi  = \int G (2m / hbar^2) v psi dx dy dz
    // G = 1/(4 Pi) exp(i k z) / z
    //     => : delta psi  = \int psi v exp(i k r) / r   1/(4 Pi) * [(2m / hbar^2) dx dy dz

    Real kx = fK0 * cosx;
    Real ky = fK0 * cosy;
    Real kz = fK0 * cosz;

#if 1
    r = fDevice->Xsection3D(fPsi0X, fPsiX, fV,
        kx, ky, kz,
        fX0, fY0, fZ0,
        fDx, fDy, fDz,
        fNx, fNy, fNz);
#else
    for_each_global_idx(fX0, fDx, fNx,
        fY0, fDy, fNy,
        fZ0, fDz, fNz,
        [&r, psi0 = fPsi0XHost, psi = fPsiXHost, v = fVHost, kx, ky, kz,
        kx0 = fK0X, ky0 = fK0Y, kz0 = fK0Z]
    (size_t idx, Real x, Real y, Real z) {
#if 0
            auto psi0_ = psi0[idx];
#else
            auto psi0_ = exp(QuI * (x * kx0 + y * ky0 + z * kz0));
#endif
            r += (psi0_ + psi[idx]) * v[idx].real() * exp(-QuI * (kx * x + ky * y + kz * z));
        });
#endif

    Real dX_dOmega = abs2(r*(2 * fMass) / (fHbar*fHbar)*fDx*fDy*fDz / (4 * Pi));
    return dX_dOmega;
}

Real ScatteringSolver3DImpl::ComputeTotalXSection(Int npsi, Int ntheta)
{
    Real t = 0;
    for (Int j = 0; j < ntheta; ++j) {
        for (Int i = 0; i < npsi; ++i) {
            Real phi = i * 2 * Pi / npsi;
            Real theta = i * Pi / ntheta;
            Real cosx = sin(theta)*cos(phi);
            Real cosy = sin(theta)*sin(phi);
            Real cosz = cos(theta);
            Real xsec = ComputeXSection(cosx, cosy, cosz)*sin(theta);
            t += xsec;
        }
    }
    t *= 4 * Pi / (ntheta * npsi);
    return t;
}
