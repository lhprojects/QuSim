
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


    if (SolverMethod::GaussLegendreO6 == solver_method) {
        //fn1 = id - 0.5 * fh * I - 1 / 10.0 * (fh * fh) + 1 / 120.0 * I * (fh * fh * fh);
        //fd1 = id + 0.5 * fh * I - 1 / 10.0 * (fh * fh) - 1 / 120.0 * I * (fh * fh * fh);
        //fdLU1.compute(fd1);
        //fd1(above) = fd1 * fd2 * fd3 (follow)
        // (confusion on notation)
        fn1 = id + (-0.1357999257081538030691227 - 0.1423427884419439108776332 * I) * fh;
        fn2 = id + (-0.2153144231161121782447335 * I) * fh;
        fn3 = id + (0.1357999257081538030691227 - 0.1423427884419439108776332 * I) * fh;

        fd1 = id + (-0.1357999257081538030691227 + 0.1423427884419439108776332 * I) * fh;
        fd2 = id + (+0.2153144231161121782447335 * I) * fh;
        fd3 = id + (0.1357999257081538030691227 + 0.1423427884419439108776332 * I) * fh;

        fdLU1.compute(fd1);   
        fdLU2.compute(fd2);
        fdLU3.compute(fd3);


    } else if (SolverMethod::GaussLegendreO4 == solver_method) {
        //            1 - 1/2 i H Dt /hbar - 1/12 (H Dt /hbar)^2 
        //  psi  ->  --------------------------------------------------   psi 
        //            1 + 1/2 i H Dt /hbar - 1/12 (H Dt /hbar)^2
        //fn1 = id - 0.5 * fh * I - 1 / 12.0 * (fh * fh);
        //fd1 = id + 0.5 * fh * I - 1 / 12.0 * (fh * fh);
        //fdLU1.compute(fd1);

        fn1 = id + (-3. * I - sqrt(3)) / 12. * fh;
        fn2 = id + (-3. * I + sqrt(3)) / 12. * fh;
        fd1 = id + (+3. * I - sqrt(3)) / 12. * fh;
        fd2 = id + (+3. * I + sqrt(3)) / 12. * fh;
        
        fdLU1.compute(fd1);
        if (fdLU1.info() != Eigen::ComputationInfo::Success) {
            throw std::runtime_error("computation failed\n");
        }
        fdLU2.compute(fd2);
        if (fdLU2.info() != Eigen::ComputationInfo::Success) {
            throw std::runtime_error("computation failed\n");
        }

    } else if (SolverMethod::ImplicitMidpointMethod == solver_method) {
        fn1 = id - 0.5 * fh * I;
        fd1 = id + 0.5 * fh * I;
        fdLU1.compute(fd1);
    } else {
        throw std::invalid_argument("unkown solver method");
    }


}


void QuGaussLegendreMethod::UpdatePsi(Complex* psi, SolverMethod solver_method)
{
    Eigen::Map<Eigen::VectorXcd> psiVec(psi, fN);
    if (solver_method == SolverMethod::ImplicitMidpointMethod) {
        fTmpPsi = fn1 * psiVec;
        psiVec = fdLU1.solve(fTmpPsi);
    } else if (solver_method == SolverMethod::GaussLegendreO4) {
        fTmpPsi = fn1 * psiVec;
        psiVec = fn2 * fTmpPsi;
        fTmpPsi = fdLU2.solve(psiVec);
        psiVec = fdLU1.solve(fTmpPsi);
    } else if (solver_method == SolverMethod::GaussLegendreO6) {
        fTmpPsi = fn1 * psiVec;
        psiVec = fn2 * fTmpPsi;
        fTmpPsi = fn2 * psiVec;
        psiVec = fdLU2.solve(fTmpPsi);
        fTmpPsi = fdLU2.solve(psiVec);
        psiVec = fdLU1.solve(fTmpPsi);
    }
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
