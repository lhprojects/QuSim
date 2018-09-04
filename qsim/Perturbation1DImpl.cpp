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

		auto sub = [this](PsiVector const &a, PsiVector &b) {
			for (size_t i = 0; i < fNx; ++i) {
				b[i] -= a[i];
			}
		};

		auto Scale = [this](Real f, PsiVector &a) {
			for (size_t i = 0; i < fNx; ++i) {
				a[i] *= f;
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
					// G psi0
					VProd(fPsi0X, ftmp1);
					X2K(ftmp1, ftmp2);
					G0Prod(ftmp2, fPsiK);
					K2X(fPsiK, fPsiX);
				} else {
					// O2: G ( 1 + G ) psi0 = G ( psi0 + O1)
					// O3: G ( 1 + G ( 1 + G ) ) psi0 = G ( psi0 + O2)
					add(fPsi0X, fPsiX);
					VProd(fPsiX, ftmp1);
					X2K(ftmp1, ftmp2);
					G0Prod(ftmp2, fPsiK);
					K2X(fPsiK, fPsiX);
				}
			}
		} else {

			size_t const NP = (int)pow(fOrder + 1, fSplit);
			std::vector<Real> coef(NP);
			std::vector<Real> ccoef(NP);
			coef[0] = 1;
			std::function<void(int n, int cur, std::vector<Real> const &coef, std::vector<Real> &ccoef)> CalCoef =
				[&](int n, int cur, std::vector<Real> const &coef, std::vector<Real> &ccoef) {

				if (cur == 0) {
					ccoef = coef;
					return;
				}  
				std::vector<Real> g(NP);
				std::vector<Real> cg(NP);

				g = coef;
				for (int o = 0; o < fOrder; ++o) {
					for (size_t i = coef.size()-1; i > 0; --i) {
						g[i] = 1. / fSplit * g[i - 1];
					}
					g[0] = 0;
					CalCoef(n, cur - 1, g, cg);
					for (size_t i = 0; i < g.size(); ++i) {
						g[i] = coef[i] + cg[i];
					}
				}

				CalCoef(n, cur - 1, g, ccoef);


			};
			CalCoef(fSplit, fSplit, coef, ccoef);

			if (0) {
				for (size_t i = 0; i < ccoef.size(); ++i) {
					printf("%lf\n", ccoef[i]);
				}
			}

			for (size_t i = 1; i < ccoef.size()-1; ++i) {
				for (size_t j = i+1; j < ccoef.size(); ++j) {
					ccoef[j] /= ccoef[i];
				}
			}

			if (0) {
				for (size_t i = 0; i < ccoef.size(); ++i) {
					printf("%lf\n", ccoef[i]);
				}
			}

			for (size_t i = ccoef.size()-1; i > 0; --i) {

				if (i == ccoef.size() - 1) {
					VProd(fPsi0X, ftmp1);
					Scale(ccoef[i], ftmp1);
					X2K(ftmp1, ftmp2);
					G0Prod(ftmp2, fPsiK);
					K2X(fPsiK, fPsiX);
				} else {
					add(fPsi0X, fPsiX);
					VProd(fPsiX, ftmp1);
					Scale(ccoef[i], ftmp1);
					X2K(ftmp1, ftmp2);
					G0Prod(ftmp2, fPsiK);
					K2X(fPsiK, fPsiX);
				}
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

