#include "SolverImpl.h"
#include "Matrix2.h"
#include "Linear.h"

SolverImpl1D::SolverImpl1D()
{
}

void SolverImpl1D::initSystem1D(std::function<Complex(Real)> const & v,
	Real x0, Real x1, size_t n,
	Real en, Complex initPsi, Complex initPsiPrime,
	SolverMethod met,
	Real mass, Real hbar,
	OptionsImpl const &opts)
{
	initSystem(en, mass, hbar, met, opts);

	fSmallRoundError = opts.GetBool("small_round_error", true);
	fNPoints = n;
	fNBins = n - 1;
	fX0 = x0;
	fDx = (x1 - x0) / fNBins;
	fVFunc = v;

	if (met == SolverMethod::ImplicitMidpointMethod) {
		fV.resize(fNBins);
		for (size_t i = 0; i < fNBins; ++i) {
			Real x = x0 + (i + 0.5)*fDx;
			fV[i] = fVFunc(x).real();
		}
	} else if (met == SolverMethod::ExplicitRungeKuttaO4Classical) {
		fV.resize(3 * fNBins);
		for (size_t i = 0; i < fNBins; ++i) {
			Real x1 = x0 + (i + 0.0)*fDx;
			Real x2 = x0 + (i + 0.5)*fDx;
			Real x3 = x0 + (i + 1.0)*fDx;
			fV[3 * i] = fVFunc(x1).real();
			fV[3 * i + 1] = fVFunc(x2).real();
			fV[3 * i + 2] = fVFunc(x3).real();
		}
	} else if (met == SolverMethod::GaussLegendreO4) {
		fV.resize(2 * fNBins);
		Real f1 = 1. / 2 - sqrt(3) / 6;
		Real f2 = 1. / 2 + sqrt(3) / 6;
		for (size_t i = 0; i < fNBins; ++i) {
			Real x1 = x0 + (i + f1)*fDx;
			Real x2 = x0 + (i + f2)*fDx;
			fV[2 * i] = fVFunc(x1).real();
			fV[2 * i + 1] = fVFunc(x2).real();
		}
	} else if (met == SolverMethod::ExplicitRungeKuttaO6Luther1967) {
		fV.resize(7 * fNBins);
		Real f1 = 0;
		Real f2 = 1;
		Real f3 = 1. / 2;
		Real f4 = 2. / 3;
		Real f5 = (7 - sqrt(21)) / 14;
		Real f6 = (7 + sqrt(21)) / 14;
		Real f7 = 1;

		for (size_t i = 0; i < fNBins; ++i) {
			Real x1 = x0 + (i + f1)*fDx;
			Real x2 = x0 + (i + f2)*fDx;
			Real x3 = x0 + (i + f3)*fDx;
			Real x4 = x0 + (i + f4)*fDx;
			Real x5 = x0 + (i + f5)*fDx;
			Real x6 = x0 + (i + f6)*fDx;
			Real x7 = x0 + (i + f7)*fDx;
			fV[7 * i] = fVFunc(x1).real();
			fV[7 * i + 1] = fVFunc(x2).real();
			fV[7 * i + 2] = fVFunc(x3).real();
			fV[7 * i + 3] = fVFunc(x4).real();
			fV[7 * i + 4] = fVFunc(x5).real();
			fV[7 * i + 5] = fVFunc(x6).real();
			fV[7 * i + 6] = fVFunc(x7).real();
		}
	} else {
		throw std::runtime_error("Unsupported method");
	}

	fPsi.resize(fNPoints);
	fPsiPrime.resize(fNPoints);

	fPsi[0] = initPsi;
	fPsiPrime[0] = initPsiPrime;

	fTMat.resize(2, 2);
	fTMat(0, 0) = 1;
	fTMat(0, 1) = 0;
	fTMat(1, 0) = 0;
	fTMat(1, 1) = 1;

	fV0 = fVFunc(x0).real();
	fV1 = fVFunc(x1).real();

	fInitJ = (std::conj(fPsi[0])*fPsiPrime[0] - fPsi[0] * std::conj(fPsiPrime[0])).imag();

}

//typedef Eigen::Matrix2d Matrix;
typedef Mat2<Real> Matrix;


#undef LOOP_FUNC_NAME
#define LOOP_FUNC_NAME MainLoop
#include "SolverMainLoop.h"
#undef LOOP_FUNC_NAME
#define LOOP_FUNC_NAME MainLoopSamllRoundError
#define SMALL_ROUND_ERROR
#include "SolverMainLoop.h"

void SolverImpl1D::Compute()
{

	//          |  psi       |
	//  Psi =   |            |
	//          |  psi^prime |
	//
	//   d Psi / dx = A Psi
	//  
	//               | 0    1 |  
	//   A   =       |        |  
	//               | a    0 |  
	//
	//   a = -2m/hbar^2 (E - V)
	//
	//   Midpoint method
	//                       1 + 1/2 A(x + 1/2dx) dx
	//      Psi(x + dx) =  -------------------------- Psi(x)
	//                       1 - 1/2 A(x + 1/2dx) dx


	if (fSmallRoundError)
		MainLoopSamllRoundError();
	else
		MainLoop();

	CalculateFinalJFromPsi();

	do {
		Real T11 = fTMat(0, 0);
		Real T12 = fTMat(0, 1);
		Real T21 = fTMat(1, 0);
		Real T22 = fTMat(1, 1);

		fT = 0;
		fR = 0;
		Real kin1 = fE - fV0;
		if (kin1 <= 0) break;

		Real kin2 = fE - fV1;
		if (kin2 <= 0) break;

		Real k1 = sqrt(2 * fMass*kin1) / fHbar;
		Real k2 = sqrt(2 * fMass*kin2) / fHbar;

		//t (1 , I k2)^T =  TMat (1 + r, i k1 - i k1 r)^T
		// t = T11 + T12 i k1 + r (T11 - T12 i k1)
		// i k2 t = i k2 (T11 + T12 i k1) + r i k2 (T11 - T12 i k1)
		// i k2 t =  T21 + T22 i k1 + r (T21 - T22 i k1)
		Complex r = -(I * k2*(T11 + T12 * I*k1) - (T21 + T22 * I*k1)) / (I*k2*(T11 - T12 * I*k1) - (T21 - T22 * I*k1));

		Complex One = (-T12 * T21 + T11 * T22);
		One = 1;
		Complex t = 2 * k1 * One / (k2*T11 - I * k1*k2*T12 + I * T21 + k1 * T22);

		fT = 4 * k1*k2 / abs2(k2*T11 - I * k1*k2*T12 + I * T21 + k1 * T22);
		fR = abs2(r);
	} while (false);

}

void SolverImpl1D::CalculateFinalJFromPsi()
{
	fFinalJ = (std::conj(fPsi[fNPoints - 1])*fPsiPrime[fNPoints - 1] - fPsi[fNPoints - 1] * std::conj(fPsiPrime[fNPoints - 1])).imag();
}
