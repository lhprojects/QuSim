#pragma once

#include "ScatteringSolverImpl.h"
#include "MatrixSolver.h"

struct InverseMatrixMethodCommon {
    SparseMatrixSolver fMatrixSolver;
    bool const fPreferPreciseSmallWaveFunction = false;
    int const fSpaceOrder = 0;
    QuSparseMatrix fEMinusH;
    Complex* const fTmpPsi = nullptr;
    Device * const fDevice_ = nullptr;

    void InitCommon(Device *dev,size_t n, OptionsImpl const & opts);
    ~InverseMatrixMethodCommon();
};

template<class Solver>
void InverseMatrixCompute(Solver* s)
{
    if (s->fPreferPreciseSmallWaveFunction) {

        s->fDevice->MulR(s->fTmpPsi, s->fPsi0X, s->fV, s->fN);
        s->fMatrixSolver.Solve(s->fEMinusH, s->fTmpPsi, s->fPsiX, s->fN);
        s->fDevice->Sub(s->fPsiX, s->fPsiX, s->fPsi0X, s->fN);

    } else {

        s->fDevice->MulR(s->fTmpPsi, s->fPsi0X, s->fV, s->fN);
        s->fMatrixSolver.Solve(s->fEMinusH, s->fTmpPsi, s->fPsiX, s->fN);

    }

}


