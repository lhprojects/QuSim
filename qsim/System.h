#pragma once

#include <vector>
#include <complex>
#include <cmath>
#include <functional>


using Real = double;
using Complex = std::complex<Real>;
using PsiVector = std::vector<Complex>;
using Int = int64_t;
using UInt = uint64_t;
static constexpr Complex I = Complex(0, 1);
static constexpr Real Pi = 3.141592653589793238462643383279502884197169399375;

#include <cmath>

inline Real abs2(Complex const &c)
{
	return real(c)*real(c) + imag(c)*imag(c);
}

enum class BoundaryCondition {
	InfiniteWall,
	Period,
	ExtendInfinity,
};

struct System {


	Complex fDt;
	bool fFNES;

	Real fMass;
	Real fHbar;

	BoundaryCondition fBoundaryCondition;
	std::string fVStr;
	Real fX0;
	Real fDx;
	size_t fN;

	std::string fPsiStr;
	bool fFN;

	std::vector<Real> fV;

	// extend infinity only
	std::vector<Complex> fProp;
	// period only
	std::vector<Complex> fFTPsi;
	// extend 2
	std::vector<Complex> f2Psi;
	std::vector<Complex> fFT2Psi;
	// infinite wall
	std::vector<Complex> fIWPsi;
	std::vector<Complex> fIWKPsi;


	Int fStep;
	std::vector<Complex> fLastLastPsi;
	std::vector<Complex> fLastPsi;
	std::vector<Complex> fPsi;
	std::vector<Complex> fVPsi;
	std::vector<Complex> fTVPsi;
	std::vector<Complex> fVTVPsi;

	System()
	{
		fN = 0;
	}
	Real getX(size_t i)
	{
		return fDx * i + fX0;
	}

	void init(char const *psi, bool force_normalization,
		Complex dt, bool force_normalization_each_step,
		char const *vs, Real x0, Real x1, size_t n, BoundaryCondition b,
		Real mass, Real hbar = 1);


	void initPsi();
	void initPotential();
	void initFreeParticleProp();

	// vpsi = exp(-i/hbar V Dt) psi
	void ExpV(PsiVector &vpsi, PsiVector const &psi, double t)
	{
		Complex f = -1.0 / fHbar * fDt * t;
		for (size_t i = 0; i < fN; ++i) {
			vpsi[i] = psi[i] * exp(f*fV[i] * I);
		}
	}

	void ExpT(PsiVector &tpsi, PsiVector const &psi);

	void Zero(PsiVector &psi)
	{
		for (size_t i = 0; i < fN; ++i) {
			psi[i] = 0.;
		}
	}

	void Scale(PsiVector &psi, Complex const &c)
	{
		for (size_t i = 0; i < fN; ++i) {
			psi[i] = psi[i] * c;
		}
	}

	void Scale(PsiVector &psi, double c)
	{
		for (size_t i = 0; i < fN; ++i) {
			psi[i] = psi[i] * c;
		}
	}

	void Add(PsiVector &psi, PsiVector const &psi1, PsiVector const &psi2)
	{
		for (size_t i = 0; i < fN; ++i) {
			psi[i] = psi1[i] + psi2[i];
		}
	}

	void Copy(PsiVector &psi, PsiVector const &psi1)
	{
		for (size_t i = 0; i < fN; ++i) {
			psi[i] = psi1[i];
		}
	}

	void step();

	Real Time()
	{
		return fStep * abs(fDt);
	}

	Real Xavg()
	{
		Real norm2 = 0;
		for (size_t i = 0; i < fN; ++i) {
			norm2 += abs2(fPsi[i])*getX(i)*fDx;
		}
		return norm2/Norm2();
	}

	Real PotEn()
	{
		Real norm2 = 0;
		for (size_t i = 0; i < fN; ++i) {
			norm2 += abs2(fPsi[i])*fV[i]*fDx;
		}
		return norm2 / Norm2();
	}

	Real KinEn();

	Real EnPartialT();

	Real Norm2(PsiVector const &psi)
	{
		Real norm2 = 0;
		for (size_t i = 0; i < psi.size(); ++i) {
			norm2 += abs2(psi[i])*fDx;
		}
		return norm2;
	}

	Real Norm2()
	{
		return Norm2(fPsi);
	}


	Real NormLeft()
	{
		Real norm2 = 0;
		for (size_t i = 0; i < fN/2; ++i) {
			norm2 += abs2(fPsi[i])*fDx;
		}
		if (fN % 2) {
			norm2 += 0.5*abs2(fPsi[fN/2])*fDx;
		}
		return norm2;
	}

	Real NormRight()
	{
		Real norm2 = 0;
		for (size_t i = (fN+1)/2; i < fN; ++i) {
			norm2 += abs2(fPsi[i])*fDx;
		}
		if (fN % 2) {
			norm2 += 0.5*abs2(fPsi[fN / 2])*fDx;
		}
		return norm2;
	}

};

