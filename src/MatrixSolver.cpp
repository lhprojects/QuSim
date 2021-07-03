#include "MatrixSolver.h"
#include "OptionsImpl.h"
#include "Utils.h"

void SparseMatrixSolver::Init(OptionsImpl const & opts)
{
	{
		MatrixSolverMethod solver = MatrixSolverMethod::LU;
		std::string mat;

		if (opts.Get("matrix_solver", mat)) {
			if (mat == "LU") {
				solver = MatrixSolverMethod::LU;
			} else if (mat == "BiCGSTAB") {
				solver = MatrixSolverMethod::BiCGSTAB;
			} else {
				throw std::runtime_error("can't parse solver");
			}
		}
		const_cast<MatrixSolverMethod&>(fMatrixSolver) = solver;
	}

	{
		Preconditioner prc = Preconditioner::DiagonalPreconditioner;

		std::string prec_str;
		if (opts.Get("matrix_preconditioner", prec_str)) {
			if (prec_str == "DiagonalPreconditioner") {
				prc = Preconditioner::DiagonalPreconditioner;
			} else if (prec_str == "IdentityPreconditioner") {
				prc = Preconditioner::IdentityPreconditioner;
			} else if (prec_str == "IncompleteLUT") {
				prc = Preconditioner::IncompleteLUT;
			} else {
				throw std::runtime_error("can't parse preconditioner");
			}
		}
		const_cast<Preconditioner&>(fPreconditioner) = prc;
	}

	mutable_cast(fMaxIters) = opts.GetDefaut("matrix_solver_max_iters", Real(1));

}

SparseMatrixSolver::SparseMatrixSolver() : fMatrixSolver(), fPreconditioner()
{
}

void SparseMatrixSolver::Solve(QuSparseMatrix const & m,
	Complex const *b_,
	Complex *x, size_t xsz)
{
	if (fMatrixSolver == MatrixSolverMethod::LU) {
		fSparseLU.compute(m);
	} else if (fMatrixSolver == MatrixSolverMethod::BiCGSTAB) {
		if (fPreconditioner == Preconditioner::DiagonalPreconditioner) {
			fBiCGSTAB_diag.compute(m);
		} else if (fPreconditioner == Preconditioner::IncompleteLUT) {
			fBiCGSTAB_ilu.compute(m);
		} else if (fPreconditioner == Preconditioner::IdentityPreconditioner) {
            fBiCGSTAB_ident.compute(m);
		}
	}

	Eigen::Map<Eigen::VectorXcd> res(x, xsz);
	Eigen::Map<const Eigen::VectorXcd> b(b_, xsz);

	if (fMatrixSolver == MatrixSolverMethod::LU) {
		res = fSparseLU.solve(b);
	} else if (fMatrixSolver == MatrixSolverMethod::BiCGSTAB) {
		if (fPreconditioner == Preconditioner::DiagonalPreconditioner) {
			fBiCGSTAB_diag.setMaxIterations((Int)(fMaxIters * m.cols() * 2));
			res = fBiCGSTAB_diag.solve(b);
		} else if (fPreconditioner == Preconditioner::IncompleteLUT) {
			fBiCGSTAB_ilu.setMaxIterations((Int)(fMaxIters * m.cols() * 2));
			res = fBiCGSTAB_ilu.solve(b);
		} else if (fPreconditioner == Preconditioner::IdentityPreconditioner) {
			fBiCGSTAB_ident.setMaxIterations((Int)(fMaxIters * m.cols() * 2));
			res = fBiCGSTAB_ident.solve(b);
		}
	}

}
