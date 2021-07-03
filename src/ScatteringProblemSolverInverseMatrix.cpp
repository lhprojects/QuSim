#include "ScatteringProblemSolverInverseMatrix.h"
#include "Utils.h"

void InverseMatrixMethodCommon::InitCommon(Device *dev, size_t n, OptionsImpl const & opts)
{
    mutable_cast(fPreferPreciseSmallWaveFunction) = opts.GetBool("PreferPreciseSmallWaveFunction", false);
    mutable_cast(fSpaceOrder) = (int)opts.GetInt("space_order", 2);
    fMatrixSolver.Init(opts);
    fEMinusH.resize(n, n);
    mutable_cast(fTmpPsi) = dev->Alloc<Complex>(n);
    mutable_cast(fDevice_) = dev;
}

InverseMatrixMethodCommon::~InverseMatrixMethodCommon()
{
    fDevice_->SafeFree(mutable_cast(fTmpPsi));
}


