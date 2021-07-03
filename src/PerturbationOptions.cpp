#define _CRT_SECURE_NO_WARNINGS
#include "PerturbationOptions.h"
#include "OptionsImpl.h"
#include "Utils.h"

void PerturbationCommon::InitPerturbationCommon(OptionsImpl const & opts,
    Real epsilon, size_t n, Device *dev)
{

    const_cast<bool&>(fPreconditional) = opts.GetBool("preconditional", false);
    if(fPreconditional) {
        BornSerisePreconditioner precer = BornSerisePreconditioner::Vellekoop;

        std::string pre;
        if (opts.Get("preconditioner", pre)) {
            if (pre == "Vellekoop") {
                precer = BornSerisePreconditioner::Vellekoop;
            } else if (pre == "Hao1") {
                precer = BornSerisePreconditioner::Hao1;
            } else if (pre == "Hao2") {
                precer = BornSerisePreconditioner::Hao2;
            } else if (pre == "BornIndentity") {
                precer = BornSerisePreconditioner::Identity;
            }
            else {
                throw std::runtime_error("Unknown preconditioner");
            }
        }
        const_cast<BornSerisePreconditioner&>(fPreconditioner) = precer;

    }

    mutable_cast(fOrder) = (int)opts.GetInt("order", 1);
    mutable_cast(fSlow) = opts.GetReal("slow", 1);

    mutable_cast(fDev) = dev;
    mutable_cast(fEpsilon) = epsilon;
    mutable_cast(fPsiK) = dev->Alloc<Complex>(n);
    if (fPreconditional) {
        mutable_cast(fTmpPsi) = dev->Alloc<Complex>(n);
    }

}


PerturbationCommon::~PerturbationCommon()
{
    fDev->SafeFree(mutable_cast(fPsiK));
    fDev->SafeFree(mutable_cast(fTmpPsi));
}
