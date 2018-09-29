#pragma once
#include "EvolverImpl.h"
#include "eigen/Eigen/Sparse"

struct GaussLegendreMethod2D : EvolverImpl2D {


	Eigen::SparseMatrix<Complex> fh;
	Eigen::SparseMatrix<Complex> fn;
	Eigen::SparseMatrix<Complex> fd;
	Eigen::SparseLU<Eigen::SparseMatrix<Complex> > fLU;
	Eigen::VectorXcd ffnPsi1;
	Eigen::VectorXcd ffnPsi2;


	ptrdiff_t Index(ptrdiff_t i, ptrdiff_t j) const
	{
		return i* (ptrdiff_t)fNy + j;
	}

	ptrdiff_t IndexP(ptrdiff_t i, ptrdiff_t j) const
	{
		i = i >= (ptrdiff_t)fNx ? i - (ptrdiff_t)fNx : i;
		if (i < 0) i += (ptrdiff_t)fNx;
		j = j >= (ptrdiff_t)fNy ? j - (ptrdiff_t)fNy : j;
		if (j < 0) j += (ptrdiff_t)fNy;
		return i* (ptrdiff_t)fNy + j;
	}

	void initSystem2D(std::function<Complex(Real, Real)> const &psi, bool force_normalization,
		Complex dt, bool force_normalization_each_step,
		std::function<Complex(Real, Real)> const &vs, Real x0, Real x1, size_t nx,
		Real y0, Real y1, size_t ny,
		BoundaryCondition b, SolverMethod solver,
		Real mass, Real hbar, OptionsImpl const &opts) override;

	void update_psi() override;
	Real CalKinEn() const override;



};
