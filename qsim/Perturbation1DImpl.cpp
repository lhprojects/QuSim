#define _CRT_SECURE_NO_WARNINGS
#include "Perturbation1DImpl.h"

void QuPerturbation1DImpl::InitPerturbation1D(std::function<Complex(Real)> const & v,
	Real x0, Real x1, size_t n, Real en, Real epsilon, Real direction,
	SolverMethod met, Real mass, Real hbar,
	std::map<std::string, std::string> const & opts)
{
	InitScatteringSolver1D(v, x0, x1, n, en, direction, met, mass, hbar, opts);

	fFFT.reset(FourierTransform::Create(fNx, false, FourierTransformLibrary::KISS));
	fInvFFT.reset(FourierTransform::Create(fNx, true, FourierTransformLibrary::KISS));
	const_cast<Real&>(fEpsilon) = epsilon;

	fPsiX.resize(fNx);
	fPsiK.resize(fNx);
	//ftmp1.resize(fNx);
	//ftmp2.resize(fNx);

	if (fMet == SolverMethod::BornSerise) {

		{
			int the_order = 0;
			auto it = opts.find("order");
			if (it != opts.end()) {
				auto &order = it->second;
				if (sscanf(order.c_str(), "%d", &the_order) < 1) {
					throw std::runtime_error("not valid order");
				}
			}
			const_cast<int&>(fOrder) = the_order;
		}

		{
			int the_split = 0;
			auto it = fOpts.find("split_n");
			if (it != fOpts.end()) {
				auto &split = it->second;
				if (sscanf(split.c_str(), "%d", &the_split) < 1) {
					throw std::runtime_error("not valid split");
				}
			}
			const_cast<int&>(fSplit) = the_split;
		}

		fPerturbationOptions.Init(opts);
		if (fPerturbationOptions.fPreconditional) {

			fVasb.resize(fNx);
			{
				for (size_t i = 0; i < fNx; ++i) {
					Real lambda = 2 * Pi / (sqrt(2 * fE*fMass) / fHbar);
					Real x = GetX(i);
					if (i < fNx / 2) {
						Real xx = (x - fX0) / (4 * lambda);
						fVasb[i] = -fE * exp(-xx * xx);
					} else {
						Real xx = (x - (fX0 + fDx * fNx)) / (4 * lambda);
						fVasb[i] = -fE * exp(-xx * xx);
					}
				}
			}

			PreconditionalBornSerise pbs;
			Real minEpsilon = pbs.GetMinEpsilon(fNx, fV.data(), fVasb.data());
			const_cast<Real&>(fEpsilon) = (epsilon < minEpsilon ? minEpsilon : epsilon);
		}

	}
}


void QuPerturbation1DImpl::Compute()
{

	auto X2K = [this](Complex const *psix, Complex *psik) {
		fFFT->Transform(psix, psik);
	};

	auto K2X = [this](Complex const *psik, Complex *psix) {
		fInvFFT->Transform(psik, psix);
		for (size_t i = 0; i < fNx; ++i) {
			psix[i] *= 1. / (fNx);
		}
	};


	if (fSplit != 0 && fSplit != 0) { // n-split Born serise

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

		auto VProd = [&](PsiVector const &psix, PsiVector &retx) {
			for (size_t i = 0; i < fNx; ++i) {
				retx[i] = fV[i] * psix[i];
			}
		};

		auto VProdU = [&](PsiVector &psix) {
			for (size_t i = 0; i < fNx; ++i) {
				psix[i] *= fV[i];
			}
		};

		auto G0ProdU = [this](PsiVector &psik) {
			for (size_t i = 0; i < fNx; ++i) {
				ptrdiff_t ii = i < fNx / 2 ? i : fNx - i;
				Real p = 2 * Pi / (fNx * fDx) * ii * fHbar;
				Real e = p * p / (2 * fMass);
				Complex Green0K = 1. / (fE + I * fEpsilon - e);
				psik[i] = Green0K * psik[i];
			}

		};

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
				for (size_t i = coef.size() - 1; i > 0; --i) {
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

		for (size_t i = 1; i < ccoef.size() - 1; ++i) {
			for (size_t j = i + 1; j < ccoef.size(); ++j) {
				ccoef[j] /= ccoef[i];
			}
		}

		if (0) {
			for (size_t i = 0; i < ccoef.size(); ++i) {
				printf("%lf\n", ccoef[i]);
			}
		}

		for (size_t i = ccoef.size() - 1; i > 0; --i) {

			if (i == ccoef.size() - 1) {
				VProd(fPsi0X, fPsiX);
			} else {
				add(fPsi0X, fPsiX);
				VProdU(fPsiX);
			}
			Scale(ccoef[i], fPsiX);
			X2K(fPsiX.data(), fPsiK.data());
			G0ProdU(fPsiK);
			K2X(fPsiK.data(), fPsiX.data());
		}

	} else if (fPerturbationOptions.fPreconditional) { // Preconditional Born serise

		ftmp1.resize(fNx);

		PreconditionalBornSerise pbs;
		for (int i = 0; i < fOrder; ++i) {
			pbs.Update1D(fNx, fPsi0X.data(), fPsiX.data(), fPsiK.data(),
				fV.data(), fVasb.data(), ftmp1.data(), fEpsilon, fE,
				fMass, fHbar, fDx, X2K, K2X,
				fPerturbationOptions.fSlow,
				fPerturbationOptions.fPreconditioner);
		}

	} else { // naive born serise
		BornSerise bs;

		for (int i = 0; i < fOrder; ++i) {
			bs.Update1D(fNx, fPsi0X.data(), fPsiX.data(),
				fPsiK.data(), fV.data(), fEpsilon, fE, fMass, fHbar, fDx,
				X2K, K2X);
		}
	}


	// post calculation
	if (0) {
		// calculate in momentum space
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

	} else {
		// calculate in real space
		Complex r = 0;
		Complex t = 0;
		for (size_t i = 0; i < fNx; ++i) {
			r += (fPsi0X[i] + fPsiX[i])*fV[i] * exp(+I * fK0*GetX(i));
			t += (fPsi0X[i] + fPsiX[i])*fV[i] * exp(-I * fK0*GetX(i));
		}
		fR = abs2(r*fDx*fMass / (fHbar*fHbar*fK0 * I));
		fT = abs2(t*fDx*fMass / (fHbar*fHbar*fK0 * I));

	}

}

