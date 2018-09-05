#pragma once

#include "SolverImpl.h"


struct ComplexPotentialIVPSolver1DImpl : SolverImpl1D {


	PsiVector fV;
	void Compute() override;


};
