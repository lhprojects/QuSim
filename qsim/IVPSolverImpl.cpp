#include "IVPSolverImpl.h"

IVPSolverImpl::IVPSolverImpl() : fMass(0), fHbar(0), fE(0), fMethod()
{
}

void IVPSolverImpl::initSystem(
	Real en,
	Real mass,
	Real hbar,
	SolverMethod met,
	OptionsImpl const &opts)
{
	const_cast<Real&>(fE) = en;
	const_cast<Real&>(fMass) = mass;
	const_cast<Real&>(fHbar) = hbar;
	const_cast<SolverMethod&>(fMethod) = met;
	const_cast<OptionsImpl&>(fOpts) = opts;
}
