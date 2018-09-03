#pragma once

#include <vector>
#include <complex>
#include <cmath>
#include <map>
#include <memory>
#include <functional>
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
	// for evolver only & H is time independent
	SplittingMethodO2,
	SplittingMethodO4,
	Eigen,

	// Runge¨CKutta family
	ImplicitMidpointMethod, // which is same as GaussLegendreO2
	GaussLegendreO4,
	GaussLegendreO6,

	ExplicitRungeKuttaO4Classical,
	ExplicitRungeKuttaO6Luther1967,

};


struct Cal;
struct EvolverImpl;
struct SolverImpl;
struct QuPerturbationImpl;

struct FunctorWrapper
{
	FunctorWrapper(char const *);
	Complex operator()(Real x);
	~FunctorWrapper();
private:
	Complex *fX;
	std::shared_ptr<Cal> fCal;
};

struct Functor2DWrapper {
	Functor2DWrapper(char const *);
	Complex operator()(Real x, Real y);
	~Functor2DWrapper();
private:
	Complex *fX;
	Complex *fY;
	std::shared_ptr<Cal> fCal;
};


struct Evolver {

	void step();
	Real Norm2();
	Real Time();
	Real PotEn();
	Real KinEn();
	Real EnPartialT();

	~Evolver();
protected:
	Evolver();
	std::shared_ptr<EvolverImpl> fImpl;
};

struct Evolver1D : Evolver {

	void init(std::function<Complex(Real)> const &psi, bool force_normalization,
		Complex dt, bool force_normalization_each_step,
		std::function<Complex(Real)> const &v, Real x0, Real x1, size_t n,
		BoundaryCondition b, SolverMethod solver,
		Real mass, Real hbar,
		std::map<std::string, std::string> const &opts);

	PsiVector const &GetPsi();
	std::vector<Real> const &GetV();
	Real Xavg();
	size_t GetN();
	Real NormLeft();
	Real NormRight();

	Evolver1D();

};


struct Evolver2D : Evolver
{
	Evolver2D();

	void init(std::function<Complex(Real, Real)> const &psi, bool force_normalization,
		Complex dt, bool force_normalization_each_step,
		std::function<Complex(Real, Real)> const &v, Real x0, Real x1, size_t nx,
		Real y0, Real y1, size_t ny,
		BoundaryCondition b, SolverMethod solver,
		Real mass, Real hbar,
		std::map<std::string, std::string> const &opts);

	Eigen::MatrixXcd const &GetPsi();
	Eigen::MatrixXd const &GetV();
	size_t GetNx();
	size_t GetNy();


};

struct Solver {

	void Compute();
	~Solver();
protected:
	Solver();
	std::shared_ptr<SolverImpl> fImpl;
};


struct Solver1D : Solver {


	void init(
		std::function<Complex(Real)> const & v,
		Real x0,
		Real x1,
		size_t n,
		Real en,
		Complex initPsi,
		Complex initPsiPrime,
		SolverMethod met,
		Real mass,
		Real hbar,
		std::map<std::string, std::string> const &opts);

	PsiVector const &GetPsi();
	std::vector<Real> const &GetV();

	size_t GetNPoints();
	Eigen::Matrix2cd GetTMat();
	Real GetT();
	Real GetR();
	Real GetEnergy();
	Complex InitPsi();
	Complex InitPsiPrime();
	Complex FinalPsi();
	Complex FinalPsiPrime();

};

struct QuPerturbation {
	~QuPerturbation();
protected:
	QuPerturbation();
	std::shared_ptr<QuPerturbationImpl> fImpl;

};

struct QuPerturbation1D : QuPerturbation {

	QuPerturbation1D();
	void init(
		std::function<Complex(Real)> const & v,
		Real x0,
		Real x1,
		size_t n,
		Real en,
		Real epsilon,
		Real direction,
		SolverMethod met,
		Real mass,
		Real hbar,
		std::map<std::string, std::string> const &opts);

	PsiVector const &GetPsi();
	std::vector<Real> const &GetV();

	size_t GetNPoints();
	void Compute();
	Real GetT();
	Real GetR();
	Real GetEnergy();
	Real GetMomentumGap();
	Real GetEpsilonMomentumWidth();
	Real GetEnergyGap();
	Real GetEpsilon();
	Real GetEpsilonBoundaryError();

};

