#pragma once

#include "SolverImpl.h"
#include "Linear.h"

struct ComplexPotentialIVPSolver1DImpl : SolverImpl1D {


	PsiVector fV;
	void Compute() override;


};
