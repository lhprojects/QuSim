#pragma once

#include "QuSim.h"
struct QuPerturbationImpl {

	QuPerturbationImpl() : fHbar(0), fMass(0), fE(0), fMet(), fOpts() { }
	~QuPerturbationImpl() { }

	void initPerturbation(Real en, SolverMethod met, Real mass, Real hbar, std::map<std::string, std::string> const &opts)
	{
		
		const_cast<Real&>(fE) = en;
		const_cast<SolverMethod&>(fMet) = met;
		const_cast<Real&>(fMass) = mass;
		const_cast<Real&>(fHbar) = hbar;
		const_cast<std::map<std::string, std::string>&>(fOpts) = opts;
	}

	Real const fE;
	SolverMethod const fMet;
	Real const fMass;
	Real const fHbar;
	std::map<std::string, std::string> const fOpts;
};
