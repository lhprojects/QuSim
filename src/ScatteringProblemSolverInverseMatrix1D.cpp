#define  _CRT_SECURE_NO_WARNINGS

#include "ScatteringProblemSolverInverseMatrix1D.h"
#include "Utils.h"
#include "Perburbation.h"
#include "GaussLegendreMethod.h"
#include "NumericalDiff.h"

void ScatteringProblemSolverInverseMatrix1D::InitScatteringSolver1D(std::function<Complex(Real)> const & v,
	Real x0, Real x1, size_t n, Real en, Real direction, SolverMethod met,
	Real mass, Real hbar, OptionsImpl const & opts)
{
	ScatteringSolver1DImpl::InitScatteringSolver1D(v, x0, x1, n, en, direction, met, mass, hbar, opts);

	InitCommon(fDevice.get(), fN, opts);	

	{		
		Real lambda = PerburbationUtility::CalLambda(fE, fMass, fHbar);
		if (24 * lambda > fNx * fDx) {
			throw std::runtime_error("too small size of to fill absorbtion layer");
		}

		PerburbationUtility::GaussMAsbLayer1D(fNx, fDx, mutable_ptr_cast(fV), fHbar, fMass, fE, 4.);

		std::vector<Eigen::Triplet<Complex> > elems;
        FillElems1D(elems, fNx, fSpaceOrder, fE, fV, Minus1(), fHbar, fMass, fDx);

        fEMinusH.setFromTriplets(elems.begin(), elems.end());
        fEMinusH.makeCompressed();
	}

}

void ScatteringProblemSolverInverseMatrix1D::Compute()
{

	InverseMatrixCompute(this);
	ComputeRT();

}
