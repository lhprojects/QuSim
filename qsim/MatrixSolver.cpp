#include "MatrixSolver.h"
#include "OptionsImpl.h"

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
		if (opts.Get("preconditioner", prec_str)) {
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

}

SparseMatrixSolver::SparseMatrixSolver() : fMatrixSolver(), fPreconditioner()
{
}

void SparseMatrixSolver::Solve(Eigen::SparseMatrix<Complex> const & m,
	Eigen::VectorXcd const &x,
	Eigen::VectorXcd &v1)
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

	if (fMatrixSolver == MatrixSolverMethod::LU) {
		v1 = fSparseLU.solve(x);
	} else if (fMatrixSolver == MatrixSolverMethod::BiCGSTAB) {
		if (fPreconditioner == Preconditioner::DiagonalPreconditioner) {
			fBiCGSTAB_diag.setMaxIterations(500);
			v1 = fBiCGSTAB_diag.solve(x);
		} else if (fPreconditioner == Preconditioner::IncompleteLUT) {
			fBiCGSTAB_ilu.setMaxIterations(500);
			v1 = fBiCGSTAB_ilu.solve(x);
		} else if (fPreconditioner == Preconditioner::IdentityPreconditioner) {
			fBiCGSTAB_ident.setMaxIterations(500);
			v1 = fBiCGSTAB_ident.solve(x);
		}
	}

}
