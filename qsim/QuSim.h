#pragma once

#include <vector>
#include <complex>
#include <cmath>
#include <map>
#include <memory>
#include <functional>
#include <assert.h>

#ifdef QUEXPORT
#ifdef _MSC_VER
#define EXPORT_FUNC __declspec(dllexport)
#define EXPORT_STRUCT __declspec(dllexport)
#endif
#else
#define EXPORT_FUNC __declspec(dllimport)
#define EXPORT_STRUCT __declspec(dllimport)
#endif

template<class T>
struct VectorView {

	VectorView(T const *d, size_t s) : fData(d), fSize(s) { }
	VectorView(std::vector<T> const &v) : fData(v.data()), fSize(v.size()) { }

	T const *begin() const { return fData; }
	T const *end() const { return fData + fSize; }
	T const *cbegin() const { return fData; }
	T const *cend() const { return fData + fSize; }
	size_t size() const { return fSize; }
	T const & operator[](size_t i) const { return fData[i]; }
	T const & operator()(size_t i) const { return fData[i]; }
private:
	T const * const fData;
	size_t const fSize;
};

template<class T>
struct MatrixView {

	MatrixView(T const *d, size_t rows, size_t cols) : fData(d), fRows(rows), fCols(cols) { }

	T const *begin() const { return fData; }
	T const *end() const { return fData + fRows * fCols; }
	T const *cbegin() const { return fData; }
	T const *cend() const { return fData + fRows * fCols; }
	size_t size() const { return fRows * fCols; }
	size_t cols() const { return fCols; }
	size_t rows() const { return fRows; }
	T const *data() const { return fData; }
	T const & operator[](size_t i) const { assert(i < size()); return fData[i]; }
	T const & operator()(size_t i) const { assert(i < size()); return fData[i]; }
	T const & operator()(size_t i, size_t j) const {
		assert(i < fRows );
		assert(j < fCols);
		return fData[i + j*fCols];
	}
private:
	T const * const fData;
	size_t const fRows;
	size_t const fCols;
};

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

	// Perburbation Method
	BornSerise,
	MatrixInverse
};


struct Cal;
struct EvolverImpl;
struct IVPSolverImpl;
struct ScatteringSolverImpl;

struct EXPORT_STRUCT FunctorWrapper
{
	FunctorWrapper(char const *);
	Complex operator()(Real x);
	~FunctorWrapper();
private:
	Complex *fX;
	std::shared_ptr<Cal> fCal;
};

struct EXPORT_STRUCT Functor2DWrapper {
	Functor2DWrapper(char const *);
	Complex operator()(Real x, Real y);
	~Functor2DWrapper();
private:
	Complex *fX;
	Complex *fY;
	std::shared_ptr<Cal> fCal;
};

struct EXPORT_STRUCT Calculator {
	Calculator(char const *expr);
	Complex &SetVaraible(char const *, Complex v);
	Complex &GetVaraible(char const *);
	Complex Evaluate();
private:
	std::shared_ptr<Cal> fImpl;
};

struct EXPORT_STRUCT Evolver {

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

struct EXPORT_STRUCT Evolver1D : Evolver {

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


struct EXPORT_STRUCT Evolver2D : Evolver
{
	Evolver2D();

	void init(std::function<Complex(Real, Real)> const &psi, bool force_normalization,
		Complex dt, bool force_normalization_each_step,
		std::function<Complex(Real, Real)> const &v, Real x0, Real x1, size_t nx,
		Real y0, Real y1, size_t ny,
		BoundaryCondition b, SolverMethod solver,
		Real mass, Real hbar,
		std::map<std::string, std::string> const &opts);

	MatrixView<Complex> GetPsi();
	MatrixView<Real> GetV();
	size_t GetNx();
	size_t GetNy();


};

struct EXPORT_STRUCT Solver {

	void Compute();
	~Solver();
protected:
	Solver();
	std::shared_ptr<IVPSolverImpl> fImpl;
};


struct EXPORT_STRUCT Solver1D : Solver {


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
	MatrixView<Real> GetTMat();
	Real GetT();
	Real GetR();
	Real GetEnergy();
	Complex InitPsi();
	Complex InitPsiPrime();
	Complex FinalPsi();
	Complex FinalPsiPrime();

};

struct EXPORT_STRUCT QuScatteringProblemSolver {
	~QuScatteringProblemSolver();
	void Compute();
protected:
	QuScatteringProblemSolver();
	std::shared_ptr<ScatteringSolverImpl> fImpl;

};

struct EXPORT_STRUCT QuScatteringProblemSolver1D : QuScatteringProblemSolver {

	PsiVector const &GetPsi();
	std::vector<Real> const &GetV();

	size_t GetNPoints();
	Real GetT();
	Real GetR();
	Real GetEnergy();
	Real GetMomentum();

};

struct EXPORT_STRUCT QuScatteringInverseMatrix1D : QuScatteringProblemSolver1D {

	void init(
		std::function<Complex(Real)> const & v,
		Real x0,
		Real x1,
		size_t n,
		Real en,
		Real direction,
		SolverMethod met,
		Real mass,
		Real hbar,
		std::map<std::string, std::string> const &opts);

};


struct EXPORT_STRUCT QuPerturbation1D : QuScatteringProblemSolver1D {

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

	Real GetMaxEnergy();
	Real GetMaxMomentum();

	Real GetMomentumGap();
	Real GetEpsilonMomentumWidth();
	Real GetEnergyGap();
	Real GetEpsilon();
	Real GetEpsilonBoundaryError();

};

