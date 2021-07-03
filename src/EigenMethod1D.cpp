#include "EigenMethod1D.h"
#include "Utils.h"
#include "NumericalDiff.h"

QuEigenMethod1D::QuEigenMethod1D()
{
}

void QuEigenMethod1D::InitSystem1D(std::function<Complex(Real)> const &psi, bool force_normalization,
    Complex dt, bool force_normalization_each_step,
    std::function<Complex(Real)> const &vs, Real x0, Real x1, size_t n,
    BoundaryCondition b, SolverMethod solver,
    Real mass, Real hbar,
    OptionsImpl const &opts)
{
    QuEvolver1DImpl::InitSystem1D(psi, force_normalization,
        dt, force_normalization_each_step,
        vs, x0, x1, n, b, solver,
        mass, hbar, opts);


    if(b == BoundaryCondition::Period) {

        fH.resize(fN, fN);
        fH.setZero();
        Real f = -fHbar * fHbar / (fDx * fDx * 2 * fMass);

        // No method to calculate eigenvalues of sparse matrix
        for (ptrdiff_t i = 0; i < (ptrdiff_t)fN; ++i) {
            fH(i, i - 2 < 0 ? fN + i - 2 : i - 2) = -1. / 12 * f;
            fH(i, i - 1 < 0 ? fN + i - 1 : i - 1) = 4. / 3 * f;
            fH(i, i + 0) = -5. / 2 * f;
            fH(i, i + 1 > (ptrdiff_t)fN - 1 ? i + 1 - fN : i + 1) = 4. / 3 * f;
            fH(i, i + 2 > (ptrdiff_t)fN - 1 ? i + 2 - fN : i + 2) = -1. / 12 * f;
            fH(i, i + 0) += QuGetReal(fV[i]);
        }

        fSolver.compute(fH);

        Eigen::VectorXd const& eigenvalues = fSolver.eigenvalues();
        expDt.resize(fN);


        for (ptrdiff_t i = 0; i < (ptrdiff_t)eigenvalues.size(); ++i) {
            expDt[i] = exp(-I * eigenvalues[i] * fDt / fHbar);
        }
        //dump_comp(expDt, "expDt.txt");

        psi0.resize(fN);

        for (ptrdiff_t i = 0; i < (ptrdiff_t)eigenvalues.size(); ++i) {
            psi0(i) = fPsi[i];
        }

        psi0_eigenspace = (psi0.transpose() * fSolver.eigenvectors()).transpose();
        //dump_matrix_comp(psi0_eigenspace, "psi0_eigenspace.txt");

        this->psi = psi0;
        this->psi_eigenspace = psi0_eigenspace;
    } else  {
        throw std::runtime_error("not supported boundary conditions!");
    }


}

void QuEigenMethod1D::UpdatePsi()
{
    for (ptrdiff_t i = 0; i < (ptrdiff_t)psi.size(); ++i) {
        psi_eigenspace[i] *= expDt(i);
    }
    psi = fSolver.eigenvectors() * psi_eigenspace;
    for (ptrdiff_t i = 0; i < (ptrdiff_t)psi.size(); ++i) {
        fPsi[i] = psi(i);
    }

    fStep += 1;
}

QuEigenMethod1D::RealType QuEigenMethod1D::CalKinEn() const
{
    return CalKinEnergy(fPsi, fNx, 4, fHbar, fMass, fDx, fDevice.get());
}
