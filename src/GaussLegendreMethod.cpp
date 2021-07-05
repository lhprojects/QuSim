
#include "Utils.h"
#include "OptionsImpl.h"
#include "eigen/Eigen/Sparse"
#include "GaussLegendreMethod.h"

void QuGaussLegendreMethod::Init(OptionsImpl const& opts, size_t n)
{
    mutable_cast(fN) = n;
    fh.resize(fN, fN); // = H Dt / hbar
    fTmpPsi.resize(fN);
}

void QuGaussLegendreMethod::Seth(std::vector<Eigen::Triplet<Complex> > elems)
{
    return fh.setFromTriplets(elems.begin(), elems.end());
}

void QuGaussLegendreMethod::ComputeLU(SolverMethod solver_method)
{
    Eigen::SparseMatrix<Complex> id(fN, fN);
    id.setIdentity();


    //            1 - 1/2 i H Dt /hbar - 1/12 (H Dt /hbar)^2 + ... 
    //  psi  ->  --------------------------------------------------   psi 
    //            1 + 1/2 i H Dt /hbar - 1/12 (H Dt /hbar)^2 + ...


    if (SolverMethod::GaussLegendreO6 == solver_method) {
        fn = id - 0.5 * fh * I - 1 / 10.0 * (fh * fh) + 1 / 120.0 * I * (fh * fh * fh);
        fd = id + 0.5 * fh * I - 1 / 10.0 * (fh * fh) - 1 / 120.0 * I * (fh * fh * fh);
    } else if (SolverMethod::GaussLegendreO4 == solver_method) {
        fn = id - 0.5 * fh * I - 1 / 12.0 * (fh * fh);
        fd = id + 0.5 * fh * I - 1 / 12.0 * (fh * fh);
        // fd = fd1 * f2
        // fd1 = id + (3I - sqrt(3))/12 fh
        // fd2 = id + (3I + sqrt(3))/12 fh

    } else if (SolverMethod::ImplicitMidpointMethod == solver_method) {
        fn = id - 0.5 * fh * I;
        fd = id + 0.5 * fh * I;
    } else {
        throw std::invalid_argument("unkown solver method");
    }

    fdLU.compute(fd);

}


void QuGaussLegendreMethod::UpdatePsi(Complex* psi)
{
    Eigen::Map<Eigen::VectorXcd> psiVec(psi, fN);
    fTmpPsi = fn * psiVec;
    psiVec = fdLU.solve(fTmpPsi);
}


Int GetSpaceOrder(OptionsImpl const& opts, SolverMethod sm)
{
    Int space_order = 2;
    if (opts.Get("space_order", space_order)) {
        if (space_order <= 2) {
            space_order = 2;
        } else if (space_order <= 4) {
            space_order = 4;
        } else if (space_order <= 6) {
            space_order = 6;
        } else {
            throw std::invalid_argument("space_order too large");
        }
    } else {
        bool space_O2 = false;
        if (opts.Get("space_O2", space_O2)) {
            space_order = 2;
        } else {
            // not set
            if (SolverMethod::ImplicitMidpointMethod == sm) {
                space_order = 2;
            } else if (SolverMethod::GaussLegendreO4 == sm) {
                space_order = 4;
            } else if (SolverMethod::GaussLegendreO6 == sm) {
                space_order = 6;
            }
        }
    }
    return space_order;
}
