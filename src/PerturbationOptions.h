#pragma once

#include <map>
#include <string>
#include "Perburbation.h"
#include "QuSim.h"

struct PerturbationOptions {

	bool fPreconditional;
	BornSerisePreconditioner fPreconditioner;
	Real fSlow;
	PerturbationOptions() : fSlow(), fPreconditional(), fPreconditioner() {
	}

	void Init(OptionsImpl const &opts);
};

