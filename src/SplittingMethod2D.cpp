#include "SplittingMethod2D.h"
#include "SplittingUtils.h"
#include "Utils.h"

void SplittingMethod2D::InitSystem2D(std::function<Complex(Real, Real)> const &psi, bool force_normalization,
	Complex dt, bool force_normalization_each_step,
	std::function<Complex(Real, Real)> const &vs, Real x0, Real x1,
	size_t nx, Real y0, Real y1,
	size_t ny, BoundaryCondition b,
	SolverMethod solver, Real mass, Real hbar,
	OptionsImpl const &opts)
{
	QuEvolver2DImpl::InitSystem2D(psi, force_normalization, dt, force_normalization_each_step,
		vs, x0, x1, nx, y0, y1, ny,
		b, solver, mass, hbar, opts);
	fFourierTransformOptions.Init(opts, fDeviceType);

	InitExpV();
	InitExpT();

	if (b == BoundaryCondition::Period) {
		fFTPsi = fDevice->Alloc<ComplexType>(fN);

		fft.reset(FourierTransform2D::Create(fNx, fNy, false, fFourierTransformOptions.fLib));
		inv_fft.reset(FourierTransform2D::Create(fNx, fNy, true, fFourierTransformOptions.fLib));

	} else {
		throw std::runtime_error("unsupported boundary condition!");
	}

}

void SplittingMethod2D::UpdatePsi()
{
	QuUpdatePsi(this);
}

void SplittingMethod2D::InitExpV()
{
	if (fCacheExp) {
		if (SolverMethod::SplittingMethodO2 == fSolverMethod) {
			mutable_cast(fExpVDt_0D5) = fDevice->Alloc<Complex>(fN);
			ComplexType f = -I / fHbar * fDt * 0.5;
			fDevice->Exp(mutable_ptr_cast(fExpVDt_0D5), fV, f, fN);
		} else if (SolverMethod::SplittingMethodO4 == fSolverMethod) {
			mutable_cast(fExpVDt_C1) = fDevice->Alloc<Complex>(fN);
			mutable_cast(fExpVDt_C2) = fDevice->Alloc<Complex>(fN);
			ComplexType f1 = -I / fHbar * fDt * SplitingConstants<Real>::C1;
			fDevice->Exp(mutable_ptr_cast(fExpVDt_C1), fV, f1, fN);
			ComplexType f2 = -I / fHbar * fDt * SplitingConstants<Real>::C2;
			fDevice->Exp(mutable_ptr_cast(fExpVDt_C2), fV, f2, fN);
		}
	}
}

void SplittingMethod2D::InitExpT()
{
	if (fCacheExp) {
        RealType const Dkx = 2 * Pi / (fDx * fNx);
        RealType const Dky = 2 * Pi / (fDy * fNy);
        RealType const DTx = QuSqr(Dkx * fHbar) / (2 * fMass);
        RealType const DTy = QuSqr(Dky * fHbar) / (2 * fMass);
        ComplexType const alpha = -I * fDt / fHbar;

		if (SolverMethod::SplittingMethodO2 == fSolverMethod) {
			mutable_cast(fExpTDt) = fDevice->Alloc<ComplexType>(fN);
			fDevice->SetOne(mutable_ptr_cast(fExpTDt), fN);
			fDevice->MulExpK2D(mutable_ptr_cast(fExpTDt), alpha, DTx, fNx, DTy, fNy);
		} else if (SolverMethod::SplittingMethodO4 == fSolverMethod) {
			mutable_cast(fExpTDt_D1) = fDevice->Alloc<ComplexType>(fN);
			mutable_cast(fExpTDt_D2) = fDevice->Alloc<ComplexType>(fN);
			fDevice->SetOne(mutable_ptr_cast(fExpTDt_D1), fN);
			fDevice->SetOne(mutable_ptr_cast(fExpTDt_D2), fN);

			fDevice->MulExpK2D(mutable_ptr_cast(fExpTDt_D1), alpha * SplitingConstants<Real>::D1, DTx, fNx, DTy, fNy);
			fDevice->MulExpK2D(mutable_ptr_cast(fExpTDt_D2), alpha * SplitingConstants<Real>::D2, DTx, fNx, DTy, fNy);
		}
	}
}

void SplittingMethod2D::ExpV(SplittingMethod2D::ComplexType* psi,
	SplittingMethod2D::RealType tt) const
{
	if (fCacheExp && tt == 0.5 && fExpVDt_0D5) {
		fDevice->Mul(psi, fExpVDt_0D5, fN);
	} else if (fCacheExp && tt == SplitingConstants<Real>::C1 && fExpVDt_C1) {
		fDevice->Mul(psi, fExpVDt_C1, fN);
	} else if (fCacheExp && tt == SplitingConstants<Real>::C2 && fExpVDt_C2) {
		fDevice->Mul(psi, fExpVDt_C2, fN);
	} else {
		ComplexType f = -I / fHbar * fDt * tt;
		fDevice->MulExp(psi, f, fV, fN);
	}
}

void SplittingMethod2D::ExpT(SplittingMethod2D::ComplexType*psi,
	SplittingMethod2D::RealType tt) const
{
	fft->Transform(psi, fFTPsi);

	if (fCacheExp && tt == 1.0 && fExpTDt) {
		fDevice->Mul(fFTPsi, fExpTDt, fN);
	} else if (fCacheExp && tt == SplitingConstants<Real>::D1 && fExpTDt_D1) {
		fDevice->Mul(fFTPsi, fExpTDt_D1, fN);
	} else if(fCacheExp && tt == SplitingConstants<Real>::D2 && fExpTDt_D2) {
		fDevice->Mul(fFTPsi, fExpTDt_D2, fN);
	} else {

        RealType const Dkx = 2 * Pi / (fDx * fNx);
        RealType const Dky = 2 * Pi / (fDy * fNy);
        RealType const DTx = QuSqr(Dkx * fHbar) / (2 * fMass);
        RealType const DTy = QuSqr(Dky * fHbar) / (2 * fMass);
        ComplexType const alpha = -I * fDt / fHbar;

        fDevice->MulExpK2D(fFTPsi, alpha * tt, DTx, fNx, DTy, fNy);
	}

    inv_fft->Transform(fFTPsi, psi);
    fDevice->Scale(psi, 1. / fN, fN);

}

SplittingMethod2D::RealType
SplittingMethod2D::CalKinEn() const
{
	RealType normnorm = fDevice->Norm2(fPsi, fN);
	fft->Transform(fPsi, fFTPsi);

	RealType Dkx = 2 * Pi / (fDx * fNx);
	RealType Dky = 2 * Pi / (fDy * fNy);
	RealType DTx = QuSqr(Dkx * fHbar) / (2 * fMass);
	RealType DTy = QuSqr(Dky * fHbar) / (2 * fMass);

	RealType kin = fDevice->Abs2K2D(fFTPsi, DTx, fNx, DTy, fNy);
	RealType norm2 = fDevice->Norm2(fFTPsi, fN);
	kin =  kin / norm2;
	return kin;
}
