#define _CRT_SECURE_NO_WARNINGS
#include "PerturbationOptions.h"
#include "OptionsImpl.h"

void PerturbationOptions::Init(OptionsImpl const & opts)
{

	const_cast<bool&>(fPreconditional) = opts.GetBool("preconditional");
	{
		BornSerisePreconditioner precer = BornSerisePreconditioner::Vellekoop;

		std::string pre;
		if (opts.Get("preconditioner", pre)) {
			if (pre == "Vellekoop") {
				precer = BornSerisePreconditioner::Vellekoop;
			} else if (pre == "Hao1") {
				precer = BornSerisePreconditioner::Hao1;
			} else if (pre == "Hao2") {
				precer = BornSerisePreconditioner::Hao2;
			} else {
				throw std::runtime_error("Unknown preconditioner");
			}
		}
		const_cast<BornSerisePreconditioner&>(fPreconditioner) = precer;

	}
	const_cast<Real&>(fSlow) = opts.GetReal("slow", 1);

}
