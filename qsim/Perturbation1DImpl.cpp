#include "Perturbation1DImpl.h"

void QuPerturbation1DImpl::InitQuPerturbation1D(std::function<Complex(Real)> const & v,
	Real x0, Real x1, size_t n, Real en, Real epsilon, Real direction,
	SolverMethod met, Real mass, Real hbar,
	std::map<std::string, std::string> const & opts)
{
	initPerturbation(en, met, mass, hbar, opts);
	const_cast<Real&>(fEpsilon) = epsilon;
	const_cast<size_t&>(fNx) = n;
	const_cast<Real&>(fX0) = x0;
	const_cast<Real&>(fDx) = (x1 - x0) / n;
	const_cast<std::function<Complex(Real)>&>(fVFunc) = v;
	const_cast<Real&>(fK0) = sqrt(2 * fMass * fE) / fHbar;

	fFFT.reset(FourierTransform::Create(fNx, false, FourierTransformLibrary::KISS));
	fInvFFT.reset(FourierTransform::Create(fNx, true, FourierTransformLibrary::KISS));

	initPotential();

	fPsi0X.resize(fNx);
	fPsi1X.resize(fNx);
	fPsi2X.resize(fNx);

	fPsi0K.resize(fNx);
	fPsi1K.resize(fNx);
	fPsi2K.resize(fNx);
	ftmp1.resize(fNx);
	ftmp2.resize(fNx);

}

void QuPerturbation1DImpl::initPotential()
{
	fV.resize(fNx);

	for (size_t i = 0; i < fNx; ++i) {
		Real x = GetX(i);
		Complex com = fVFunc(x);
		fV[i] = com.real();

	}
}


void QuPerturbation1DImpl::Compute()
{

	if (true) {
		for (size_t i = 0; i < fNx; ++i) {
			fPsi0X[i] = exp(fK0 * GetX(i) * I);
		}

		fFFT->Transform(fPsi0X.data(), fPsi0K.data());

		auto VProd = [this](PsiVector const &psix, PsiVector &retx) {
			for (size_t i = 0; i < fNx; ++i) {
				retx[i] = fV[i] * psix[i];
			}
		};
		auto G0Prod = [this](PsiVector const &psik, PsiVector &retk) {
			for (size_t i = 0; i < fNx; ++i) {
				ptrdiff_t ii = i < fNx / 2 ? i : fNx - i;
				Real p = 2 * Pi / (fNx * fDx) * ii * fHbar;
				Real e = p * p / (2 * fMass);
				Complex Green0K = 1. / (fE + I * fEpsilon - e);
				retk[i] = Green0K * psik[i];
			}

		};

		auto X2K = [this](PsiVector const &psix, PsiVector &psik) {
			fFFT->Transform(psix.data(), psik.data());
		};

		auto K2X = [this](PsiVector const &psik, PsiVector &psix) {
			fInvFFT->Transform(psik.data(), psix.data());
			for (size_t i = 0; i < fNx; ++i) {
				psix[i] *= 1. / (fNx);
			}
		};

		VProd(fPsi0X, ftmp1);
		X2K(ftmp1, ftmp2);
		G0Prod(ftmp2, fPsi1K);
		K2X(fPsi1K, fPsi1X);

		VProd(fPsi1X, ftmp1);
		X2K(ftmp1, ftmp2);
		G0Prod(ftmp2, fPsi2K);
		K2X(fPsi2K, fPsi2X);

		Real r = 0;
		Real t = 0;
		for (size_t i = fNx / 2; i < fNx; ++i) {
			r += abs2(fPsi1K[i]);
		}

		for (size_t i = 0; i < fNx / 2; ++i) {
			t += abs2(fPsi1K[i]);
		}

		fR = r * fEpsilon / fE * fDx / fNx;
		fT = t * fEpsilon / fE * fDx / fNx;
	}

}

