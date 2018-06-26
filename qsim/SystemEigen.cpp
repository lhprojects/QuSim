#include "SystemEigen.h"

SystemEigen::SystemEigen()
{
}

void SystemEigen::init(char const * psi, bool force_normalization, 
	Complex dt, bool force_normalization_each_step,
	char const * vs, Real x0, Real x1, size_t n,
	BoundaryCondition b, SolverMethod solver,
	Real mass, Real hbar)
{
	SystemImpl::init(psi, force_normalization,
		dt, force_normalization_each_step,
		vs, x0, x1, n, b, solver,
		mass, hbar);

	//fH.resize(fN, fN);
}

void SystemEigen::update_psi()
{

}
