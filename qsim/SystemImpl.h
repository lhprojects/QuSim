#pragma once
#define _CRT_SECURE_NO_WARNINGS

#include "System.h"
#include <stdio.h>

template<class PsiVector>
inline void dump_comp(PsiVector const &v, char const *fn)
{
	FILE *f = fopen(fn, " w");
	for (int i = 0; i < (int)v.size(); ++i) {
		auto x = v[i];
		fprintf(f, "(% .20lf, % .20lf)\n", x.real(), x.imag());
	}
	fclose(f);
};

template<class PsiVector>
inline void dump_real(PsiVector const &v, char const *fn)
{
	FILE *f = fopen(fn, " w");
	for (int i = 0; i < (int)v.size(); ++i) {
		fprintf(f, "% .20lf\n", (double)v(i));
	}
	fclose(f);
};

template<class Matrix>
inline void dump_matrix_real(Matrix const &v, char const *fn)
{
	FILE *f = fopen(fn, " w");
	for (int i = 0; i < (int)v.rows(); ++i) {
		for (int j = 0; j < (int)v.cols(); ++j) {
			fprintf(f, "%+.20lf ", (double)v(i, j));
		}
		fprintf(f, "\n");
	}
	fclose(f);
};

template<class Matrix>
inline void dump_matrix_comp(Matrix const &v, char const *fn)
{
	FILE *f = fopen(fn, " w");
	for (int i = 0; i < (int)v.rows(); ++i) {
		for (int j = 0; j < (int)v.cols(); ++j) {
			fprintf(f, "(% .20lf, % .20lf)", v(i, j).real(), v(i, j).imag());
		}
		fprintf(f, "\n");
	}
	fclose(f);
};

struct SystemImpl
{
	Complex fDt;
	bool fFNES;

	Real fMass;
	Real fHbar;

	std::string fVStr;

	BoundaryCondition fBoundaryCondition;
	SolverMethod fSolverMethod;

	std::string fPsiStr;
	bool fFN;
	Int fStep;

	// init fPsi
	// init fV
	// init fN
	virtual void init(char const *psi, bool force_normalization,
		Complex dt, bool force_normalization_each_step,
		char const *vs, BoundaryCondition b, SolverMethod solver,
		Real mass, Real hbar = 1);

	virtual void step() = 0;
	Real Time()
	{
		return fStep * abs(fDt);
	}

	virtual Real PotEn() = 0;
	virtual Real KinEn() = 0;
	virtual Real EnPartialT() = 0;
	virtual Real Norm2() = 0;



};

struct SystemImpl1D : SystemImpl {

	Real fX0;
	Real fDx;
	size_t fN;
	std::vector<Real> fV;



	std::vector<Complex> fLastLastPsi;
	std::vector<Complex> fLastPsi;
	std::vector<Complex> fPsi;
	
	void initPsi();
	void initPotential();

	SystemImpl1D()
	{
		fN = 0;
	}

	Real getX(size_t i)
	{
		return fDx * i + fX0;
	}

	// init fPsi
	// init fV
	// init fN
	virtual void init(char const *psi, bool force_normalization,
		Complex dt, bool force_normalization_each_step,
		char const *vs, Real x0, Real x1, size_t n,
		BoundaryCondition b, SolverMethod solver,
		Real mass, Real hbar = 1);


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

	void step() override
	{
		fLastLastPsi = fLastPsi;
		fLastPsi = fPsi;

		update_psi();

		if (fFNES) {
			Scale(fPsi, 1.0 / sqrt(Norm2()));
		}
		++fStep;

	}

	virtual Real CalPotEn();
	virtual Real CalKinEn();
	// update fPsi
	virtual void update_psi() = 0;

	Real Xavg()
	{
		Real norm2 = 0;
		for (size_t i = 0; i < fN; ++i) {
			norm2 += abs2(fPsi[i])*getX(i)*fDx;
		}
		return norm2 / Norm2();
	}

	Real PotEn();
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
		for (size_t i = 0; i < fN / 2; ++i) {
			norm2 += abs2(fPsi[i])*fDx;
		}
		if (fN % 2) {
			norm2 += 0.5*abs2(fPsi[fN / 2])*fDx;
		}
		return norm2;
	}

	Real NormRight()
	{
		Real norm2 = 0;
		for (size_t i = (fN + 1) / 2; i < fN; ++i) {
			norm2 += abs2(fPsi[i])*fDx;
		}
		if (fN % 2) {
			norm2 += 0.5*abs2(fPsi[fN / 2])*fDx;
		}
		return norm2;
	}

};