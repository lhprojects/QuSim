#define _CRT_SECURE_NO_WARNINGS
#include "PerturbationOptions.h"

void PerturbationOptions::Init(std::map<std::string, std::string> const & opts)
{

	{
		bool prec = false;
		auto it = opts.find("preconditional");
		if (it != opts.end()) {
			prec = it->second != "0";
		}
		const_cast<bool&>(fPreconditional) = prec;

	}

	{
		BornSerisePreconditioner precer = BornSerisePreconditioner::Vellekoop;
		auto it = opts.find("preconditioner");
		if (it != opts.end()) {
			if (it->second == "Vellekoop") {
				precer = BornSerisePreconditioner::Vellekoop;
			} else if (it->second == "Hao1") {
				precer = BornSerisePreconditioner::Hao1;
			} else if (it->second == "Hao2") {
				precer = BornSerisePreconditioner::Hao2;
			} else {
				throw std::runtime_error("Unknown preconditioner");
			}
		}
		const_cast<BornSerisePreconditioner&>(fPreconditioner) = precer;

	}

	{

		auto it = opts.find("slow");
		Real slow = 1;
		if (it != opts.end()) {
			if (sscanf(it->second.c_str(), "%lf", &slow) < 1) {
				throw std::runtime_error("parse slow factor error");
			}

		}
		const_cast<Real&>(fSlow) = slow;

	}

}
