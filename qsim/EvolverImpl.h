#pragma once
//#define _CRT_SECURE_NO_WARNINGS
//#define _SCL_SECURE_NO_WARNINGS

#include "QuSim.h"
#include "eigen/Eigen/Dense"
#include <stdio.h>
#include <functional>

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

struct EvolverImpl
{
	Complex fDt;
	bool fFNES;

	Real fMass;
	Real fHbar;

	BoundaryCondition fBoundaryCondition;
	SolverMethod fSolverMethod;

	bool fFN;
	Int fStep;

	std::map<std::string, std::string> const fOpts;

	// init fPsi
	// init fV
	// init fN
	void initSystem(bool force_normalization,
		Complex dt, bool force_normalization_each_step,
		BoundaryCondition b, SolverMethod solver,
		Real mass, Real hbar, std::map<std::string, std::string> const &opts);

	virtual void step() = 0;
	virtual Real PotEn() = 0;
	virtual Real KinEn() = 0;
	virtual Real EnPartialT() = 0;
	// int abs2(psi) dx for 1D
	//( sum_i psi_i Dx)
	// or int abs2(psi) dx dy for 2D
	//( sum_ij psi_ij Dx Dy)
	virtual Real Norm2() = 0;
	virtual ~EvolverImpl() { }

	Real Time()
	{
		return fStep * abs(fDt);
	}

};

struct EvolverImpl1D : EvolverImpl {

	Real fX0;
	Real fDx;
	size_t fN;

	std::function<Complex(Real)> fVFunc;
	std::vector<Real> fV;

	std::function<Complex(Real)> fPsi0Func;
	std::vector<Complex> fLastLastPsi;
	std::vector<Complex> fLastPsi;
	std::vector<Complex> fPsi;
	
	EvolverImpl1D()
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
	virtual void initSystem1D(std::function<Complex(Real)> const & psi, bool force_normalization,
		Complex dt, bool force_normalization_each_step,
		std::function<Complex(Real)> const &v, Real x0, Real x1, size_t n,
		BoundaryCondition b, SolverMethod solver,
		Real mass, Real hbar, std::map<std::string, std::string> const &opts);

	virtual Real CalPotEn();
	virtual Real CalKinEn();
	// update fPsi
	virtual void update_psi() = 0;

	void step() override;
	Real PotEn() override;
	Real KinEn() override;
	Real EnPartialT() override;
	Real Norm2() override;


private:
	void initPsi();
	void initPotential();
public:

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

	Real Xavg()
	{
		Real norm2 = 0;
		for (size_t i = 0; i < fN; ++i) {
			norm2 += abs2(fPsi[i])*getX(i)*fDx;
		}
		return norm2 / Norm2();
	}

	Real Norm2(PsiVector const &psi)
	{
		Real norm2 = 0;
		for (size_t i = 0; i < psi.size(); ++i) {
			norm2 += abs2(psi[i])*fDx;
		}
		return norm2;
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

struct EvolverImpl2D : EvolverImpl {

	Real fX0;
	Real fDx;
	Real fY0;
	Real fDy;
	size_t fNx;
	size_t fNy;
	Eigen::MatrixXd fV;
	std::function<Complex(Real, Real)> fVFunc;

	std::function<Complex(Real, Real)> fPsi0Func;
	// col major
	// y major
	Eigen::MatrixXcd fLastLastPsi;
	Eigen::MatrixXcd fLastPsi;
	Eigen::MatrixXcd fPsi;

	EvolverImpl2D()
	{
		fNx = 0;
		fNy = 0;
	}

	Real getX(size_t i)
	{
		return fDx * i + fX0;
	}

	Real getY(size_t i)
	{
		return fDy * i + fY0;
	}

	// init fPsi
	// init fV
	// init fN
	virtual void initSystem2D(std::function<Complex(Real, Real)> const &psi, bool force_normalization,
		Complex dt, bool force_normalization_each_step,
		std::function<Complex(Real, Real)> const &v, Real x0, Real x1, size_t nx,
		Real y0, Real y1, size_t ny,
		BoundaryCondition b, SolverMethod solver,
		Real mass, Real hbar, std::map<std::string, std::string> const &opts);
	virtual Real CalPotEn() const;
	virtual Real CalKinEn() const;
	// update fPsi
	virtual void update_psi() = 0;

	void step() override;
	Real PotEn() override { return CalPotEn(); }
	Real KinEn() override { return CalKinEn(); }
	Real EnPartialT() override;
	Real Norm2() override;

private:
	void initPsi();
	void initPotential();
public:







};
