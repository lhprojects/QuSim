#include "MatrixSolver.h"

void SparseMatrixSolver::Init(std::map<std::string, std::string> const & opts)
{
	{
		MatrixSolverMethod solver = MatrixSolverMethod::LU;
		auto it = opts.find("matrix_solver");
		if (it != opts.end()) {
			if (it->second == "LU") {
				solver = MatrixSolverMethod::LU;
			} else if (it->second == "BiCGSTAB") {
				solver = MatrixSolverMethod::BiCGSTAB;
			} else {
				throw std::runtime_error("can't parse solver");
			}
		}
		const_cast<MatrixSolverMethod&>(fMatrixSolver) = solver;
	}

	{
		Preconditioner prc = Preconditioner::DiagonalPreconditioner;
		auto it = opts.find("preconditioner");
		if (it != opts.end()) {
			if (it->second == "DiagonalPreconditioner") {
				prc = Preconditioner::DiagonalPreconditioner;
			} else if (it->second == "IdentityPreconditioner") {
				prc = Preconditioner::IdentityPreconditioner;
			} else if (it->second == "IncompleteLUT") {
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
			fBiCGSTAB_diag.setMaxIterations(10000);
			v1 = fBiCGSTAB_diag.solve(x);
		} else if (fPreconditioner == Preconditioner::IncompleteLUT) {
			fBiCGSTAB_ilu.setMaxIterations(10000);
			v1 = fBiCGSTAB_ilu.solve(x);
		} else if (fPreconditioner == Preconditioner::IdentityPreconditioner) {
			fBiCGSTAB_ilu.setMaxIterations(10000);
			v1 = fBiCGSTAB_ident.solve(x);
		}
	}

}
