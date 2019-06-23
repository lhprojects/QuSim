#pragma once

#include <complex>
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
	// data stored col major
	// first is the row (represent y)
	// second is the col (represent x)

	T const & operator()(size_t row, size_t col) const {
		assert(row < fRows );
		assert(col < fCols);
		return fData[row + col*fRows];
	}
private:
	T const * const fData;
	size_t const fRows;
	size_t const fCols;
};

template<class T>
struct Tensor3View {

	// sz3 most inner
	// sz1 most outer
	Tensor3View(T const *d, size_t sz1, size_t sz2, size_t sz3) : fData(d), fSize1(sz1), fSize2(sz2), fSize3(sz3) {}

	T const *begin() const { return fData; }
	T const *end() const { return fData + fSize1 * fSize2 * fSize3; }
	T const *cbegin() const { return fData; }
	T const *cend() const { return fData + fSize1 * fSize2 * fSize3; }
	size_t size() const { return fSize1 * fSize2 * fSize3; }
	size_t size1() const { return fSize1; }
	size_t size2() const { return fSize2; }
	size_t size3() const { return fSize3; }
	T const *data() const { return fData; }
	T const & operator[](size_t i) const { assert(i < size()); return fData[i]; }
	T const & operator()(size_t i) const { assert(i < size()); return fData[i]; }
	T const & operator()(size_t i3, size_t i2, size_t i1) const
	{
		assert(i3 < fSize3);
		assert(i2 < fSize2);
		assert(i1 < fSize1);
		return fData[i3 + fSize3 *(i2  + fSize2 * i1)];
	}
private:
	T const * const fData;
	size_t const fSize1;
	size_t const fSize2;
	size_t const fSize3;
};

using Real = double;
using Complex = std::complex<Real>;
using Int = int64_t;
using UInt = uint64_t;
static constexpr Complex I = Complex(0, 1);
static constexpr Real Pi = 3.141592653589793238462643383279502884197169399375;

struct OptionsImpl;
struct Options {

	inline Options &Order(Int n)
	{
		SetInt("order", n);
		return *this;
	}

	inline Options &VellekoopPreconditioner()
	{
		SetString("preconditioner", "Vellekoop");
		return *this;
	}

	inline Options &Hao1Preconditioner()
	{
		SetString("preconditioner", "Hao1");
		return *this;
	}

	inline Options &BornIdentityPreconditioner()
	{
		SetString("preconditioner", "BornIndentity");
		return *this;
	}

	inline Options &Hao2Preconditioner()
	{
		SetString("preconditioner", "Hao2");
		return *this;
	}

	inline Options &Preconditional(bool v)
	{
		SetBool("preconditional", v);
		return *this;
	}

	inline Options &Slow(Real v)
	{
		SetReal("slow", v);
		return *this;
	}

	inline Options &SplitN(Int n)
	{
		SetInt("split_n", n);
		return *this;
	}

	inline Options &SpaceOrder(Int n)
	{
		SetInt("space_order", n);
		return *this;
	}

	inline Options &Cuda()
	{
		SetString("fft_lib", "cuda");
		return *this;
	}

	inline Options &Batch(Int n)
	{
		SetInt("batch", n);
		return *this;
	}

	inline Options &CudaPrecisionSingle()
	{
		SetString("cuda_precision", "single");
		return *this;
	}

	inline Options &SmallRoundError(bool v)
	{
		SetBool("small_round_error", v);
		return *this;
	}

	inline Options &MatrixSolverBiCGSTAB()
	{
		SetString("matrix_solver", "BiCGSTAB");
		return *this;
	}

	inline Options &MatrixSolverLU()
	{
		SetString("matrix_solver", "LU");
		return *this;
	}

	inline Options &DiagonalPreconditioner()
	{
		SetString("preconditioner", "DiagonalPreconditioner");
		return *this;
	}


	inline Options &IdentityPreconditioner()
	{
		SetString("preconditioner", "IdentityPreconditioner");
		return *this;
	}

	inline Options &IncompleteLUTPreconditioner()
	{
		SetString("preconditioner", "IncompleteLUT");
		return *this;
	}

	EXPORT_FUNC Options();
	EXPORT_FUNC Options(Options const &);
	EXPORT_FUNC Options& operator=(Options const &);
	EXPORT_FUNC ~Options();

	EXPORT_FUNC Int GetInt(char const *);
	EXPORT_FUNC Real GetReal(char const *);
	EXPORT_FUNC Complex GetComplex(char const *);
	EXPORT_FUNC char const *GetString(char const *);

	EXPORT_FUNC void SetBool(char const *, bool v);
	EXPORT_FUNC void SetReal(char const *, Real v);
	EXPORT_FUNC void SetComplex(char const *, Complex v);
	EXPORT_FUNC void SetInt(char const *, Int v);
	EXPORT_FUNC void SetString(char const *, char const *v);

	OptionsImpl *fOpts;
};


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
	FunctorWrapper(FunctorWrapper const &r);
	FunctorWrapper &operator=(FunctorWrapper const &r);
	Complex operator()(Real x);
	~FunctorWrapper();
private:
	Complex *fX;
	Cal *fCal;
	int *fRef;
};

struct EXPORT_STRUCT Functor2DWrapper {
	Functor2DWrapper(char const *);
	Functor2DWrapper(Functor2DWrapper const &r);
	Functor2DWrapper &operator=(Functor2DWrapper const &r);
	Complex operator()(Real x, Real y);
	~Functor2DWrapper();
private:
	Complex *fX;
	Complex *fY;
	Cal *fCal;
	int *fRef;
};

struct EXPORT_STRUCT Calculator {
	Calculator(char const *expr);
	Calculator(Calculator const &) = delete;
	Calculator &operator=(Calculator const &) = delete;
	~Calculator();
	Complex &SetVaraible(char const *, Complex v);
	Complex &GetVaraible(char const *);
	Complex Evaluate();
private:
	Cal *fImpl;
};

struct EXPORT_STRUCT Evolver {
	void step();
	Real Norm2();
	Real Time();
	Real PotEn();
	Real KinEn();
	Real EnPartialT();
	SolverMethod GetMethod();

	Evolver(Evolver const &) = delete;
	Evolver const &operator=(Evolver const &) = delete;
	~Evolver();
protected:
	Evolver();
	EvolverImpl *fImpl;
};

struct EXPORT_STRUCT Evolver1D : Evolver {

	void init(std::function<Complex(Real)> const &psi, bool force_normalization,
		Complex dt, bool force_normalization_each_step,
		std::function<Complex(Real)> const &v, Real x0, Real x1, size_t n,
		BoundaryCondition b, SolverMethod solver,
		Real mass, Real hbar,
		Options const &opts);

	VectorView<Complex> GetPsi();
	VectorView<Real> GetV();
	Real Xavg();
	size_t GetN();
	Real NormLeft();
	Real NormRight();

	Evolver1D();

};


struct EXPORT_STRUCT Evolver2D : Evolver
{
	void init(std::function<Complex(Real, Real)> const &psi, bool force_normalization,
		Complex dt, bool force_normalization_each_step,
		std::function<Complex(Real, Real)> const &v, Real x0, Real x1, size_t nx,
		Real y0, Real y1, size_t ny,
		BoundaryCondition b, SolverMethod solver,
		Real mass, Real hbar,
		Options const &opts);

	MatrixView<Complex> GetPsi();
	MatrixView<Real> GetV();
	size_t GetNx();
	size_t GetNy();


};

struct EXPORT_STRUCT Solver {

	void Compute();
	~Solver();
	Solver(Solver const &) = delete;
	Solver const &operator=(Solver const &) = delete;
	SolverMethod GetMethod();

protected:
	Solver();
	IVPSolverImpl *fImpl;
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
		Options const &opts);

	VectorView<Complex> GetPsi();
	VectorView<Real> GetV();

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
	Real GetEnergy();
	Real GetMomentum();
	SolverMethod GetMethod();

	QuScatteringProblemSolver(QuScatteringProblemSolver const &) = delete;
	QuScatteringProblemSolver const &operator=(QuScatteringProblemSolver const &) = delete;

protected:
	QuScatteringProblemSolver();
	ScatteringSolverImpl *fImpl;

};

struct EXPORT_STRUCT QuScatteringProblemSolver1D : QuScatteringProblemSolver {

	VectorView<Complex> GetPsi();
	VectorView<Complex> GetPsi0();
	VectorView<Real> GetV();

	size_t GetNPoints();
	Real GetT();
	Real GetR();

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
		Options const &opts);

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
		Options const &opts);

	Real GetMaxEnergy();
	Real GetMaxMomentum();

	Real GetMomentumGap();
	Real GetEpsilonMomentumWidth();
	Real GetEnergyGap();
	Real GetEpsilon();
	Real GetEpsilonBoundaryError();

};

struct EXPORT_STRUCT QuScatteringProblemSolver2D : QuScatteringProblemSolver {

	MatrixView<Complex> GetPsi0();
	// scattering part of wave function
	MatrixView<Complex> GetPsi();

	MatrixView<Real> GetV();

	Real ComputeXSection(Real cosx, Real cosy);
	Real ComputeTotalXSection(Int n);

};

struct EXPORT_STRUCT QuScatteringInverseMatrix2D : QuScatteringProblemSolver2D {

	void init(
		std::function<Complex(Real, Real)> const & v,
		Real x0,
		Real x1,
		size_t nx,
		Real y0,
		Real y1,
		size_t ny,
		Real en,
		Real directionx,
		Real directiony,
		SolverMethod met,
		Real mass,
		Real hbar,
		Options const &opts);

};

struct EXPORT_STRUCT QuPerturbation2D : QuScatteringProblemSolver2D {

	void init(
		std::function<Complex(Real, Real)> const & v,
		Real x0,
		Real x1,
		size_t nx,
		Real y0,
		Real y1,
		size_t ny,
		Real en,
		Real epsilon,
		Real directionx,
		Real directiony,
		SolverMethod met,
		Real mass,
		Real hbar,
		Options const &opts);

	MatrixView<Real> GetVabsb();
	Real GetDeltaPsiNorm();
};

struct EXPORT_STRUCT QuScatteringProblemSolver3D : QuScatteringProblemSolver {

	Tensor3View<Complex> GetPsi();
	Tensor3View<Real> GetV();

	Real ComputeXSection(Real cosx, Real cosy, Real cosz);
	Real ComputeTotalXSection(Int npsi, Int ntheta);
};

struct EXPORT_STRUCT QuScatteringInverseMatrix3D : QuScatteringProblemSolver3D {


	void init(
		std::function<Complex(Real, Real, Real)> const & v,
		Real x0,
		Real x1,
		size_t nx,
		Real y0,
		Real y1,
		size_t ny,
		Real z0,
		Real z1,
		size_t nz,
		Real en,
		Real directionx,
		Real directiony,
		Real directionz,
		SolverMethod met,
		Real mass,
		Real hbar,
		Options const &opts);

};

struct EXPORT_STRUCT QuPerturbation3D : QuScatteringProblemSolver3D {


	void init(
		std::function<Complex(Real, Real, Real)> const & v,
		Real x0,
		Real x1,
		size_t nx,
		Real y0,
		Real y1,
		size_t ny,
		Real z0,
		Real z1,
		size_t nz,
		Real en,
		Real epsilon,
		Real directionx,
		Real directiony,
		Real directionz,
		SolverMethod met,
		Real mass,
		Real hbar,
		Options const &opts);

};
