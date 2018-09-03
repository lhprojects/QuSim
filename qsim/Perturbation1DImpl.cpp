#define _CRT_SECURE_NO_WARNINGS
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

	fPsiX.resize(fNx);
	fPsiK.resize(fNx);
	fPsi0X.resize(fNx);
	ftmp1.resize(fNx);
	ftmp2.resize(fNx);

	if (fMet == SolverMethod::BornSerise) {

		{
			int the_order = 1;
			auto it = fOpts.find("order");
			if (it != fOpts.end()) {
				auto &order = it->second;
				if (sscanf(order.c_str(), "%d", &the_order) < 1) {
					throw std::runtime_error("not valid order");
				}
			}
			const_cast<int&>(fOrder) = the_order;
		}

		{
			int the_split = 1;
			auto it = fOpts.find("split_n");
			if (it != fOpts.end()) {
				auto &split = it->second;
				if (sscanf(split.c_str(), "%d", &the_split) < 1) {
					throw std::runtime_error("not valid split");
				}
			}
			const_cast<int&>(fSplit) = the_split;
		}
	}
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

	if (fMet == SolverMethod::BornSerise) {
		for (size_t i = 0; i < fNx; ++i) {
			fPsi0X[i] = exp(fK0 * GetX(i) * I);
		}

		auto add = [this](PsiVector const &a, PsiVector &b) {
			for (size_t i = 0; i < fNx; ++i) {
				b[i] += a[i];
			}
		};
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

		if (fSplit == 1) {
			for (int i = 0; i < fOrder; ++i) {
				if (i == 0) {
					VProd(fPsi0X, ftmp1);
				} else {
					add(fPsi0X, fPsiX);
					VProd(fPsiX, ftmp1);
				}
				X2K(ftmp1, ftmp2);
				G0Prod(ftmp2, fPsiK);
				K2X(fPsiK, fPsiX);
			}

		}


		// post calculation
		Real r = 0;
		Real t = 0;
		for (size_t i = fNx / 2; i < fNx; ++i) {
			r += abs2(fPsiK[i]);
		}

		for (size_t i = 0; i < fNx / 2; ++i) {
			t += abs2(fPsiK[i]);
		}

		fR = r * fEpsilon / fE * fDx / fNx;
		fT = t * fEpsilon / fE * fDx / fNx;
	}

}

