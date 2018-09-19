#pragma once

#define  _CRT_SECURE_NO_WARNINGS
#include "ScatteringSolverImpl.h"
#include "MatrixSolver.h"

struct ScatteringProblemSolverInverseMatrix3D : ScatteringSolver3DImpl {

	ScatteringProblemSolverInverseMatrix3D() : fOrder(), fMatrixSolver() {}

	virtual void InitScatteringSolver3D(
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
		std::map<std::string, std::string> const &opts);


	void Compute() override;

	Real ComputeXSection(Real cosx, Real cosy, Real cosz) override;
	// npsi = number of sampling points for psi
	// ntheta = number of sampling points for theta
	Real ComputeTotalXSection(Int npsi, Int ntheta) override;

	void InitEMinusH();
	SparseMatrixSolver fMatrixSolver;
	int const fOrder;
	Eigen::SparseMatrix<Complex> fEMinusH;



};
