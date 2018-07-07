#pragma once

#include <vector>
#include <complex>
#include <cmath>
#include "eigen/Eigen/Dense"

using Real = double;
using Complex = std::complex<Real>;
using PsiVector = std::vector<Complex>;
using Int = int64_t;
using UInt = uint64_t;
static constexpr Complex I = Complex(0, 1);
static constexpr Real Pi = 3.141592653589793238462643383279502884197169399375;

inline Real abs2(Complex const &c)
{
	return real(c)*real(c) + imag(c)*imag(c);
}


enum class BoundaryCondition {
	InfiniteWall,
	Period,
};

enum class SolverMethod {
	SplittingMethodO2,
	SplittingMethodO4,
	Eigen,
	ImplicitMidpointMethod,
	RungeKuttaMethod
};

struct SystemImpl;
struct SystemImpl1D;
struct SystemImpl2D;

struct System {

	void step();
	Real Norm2();
	Real Time();
	Real PotEn();
	Real KinEn();
	Real EnPartialT();

	System();
protected:
	SystemImpl *fImpl;
};


struct System1D : System {

	void init(char const *psi, bool force_normalization,
		Complex dt, bool force_normalization_each_step,
		char const *vs, Real x0, Real x1, size_t n,
		BoundaryCondition b, SolverMethod solver,
		Real mass, Real hbar);

	PsiVector const &GetPsi();
	std::vector<Real> const &GetV();
	Real Xavg();
	UInt GetN();
	Real NormLeft();
	Real NormRight();

	System1D();
protected:
	SystemImpl1D * fImpl;

};


struct System2D : System
{
	System2D();

	void init(char const *psi, bool force_normalization,
		Complex dt, bool force_normalization_each_step,
		char const *vs, Real x0, Real x1, size_t nx,
		Real y0, Real y1, size_t ny,
		BoundaryCondition b, SolverMethod solver,
		Real mass, Real hbar);

	Eigen::MatrixXcd const &GetPsi();
	Eigen::MatrixXd const &GetV();
	UInt GetNx();
	UInt GetNy();

protected:
	SystemImpl2D * fImpl;

};
