#pragma once

#include "QuSim.h"

template<class R>
struct SplitingConstants {
    static const R C1;
    static const R C2;
    static const R D1;
    static const R D2;
};


template<class SplitingMethod>
void QuUpdatePsi(SplitingMethod* sm)
{
    // exp( (h1+h2£©t) ~ 1 + h1 t + h2 t + 0.5 (h1 + h2 + h1 h2 + h2 h1) t^2
    // exp( h1 t) ~ 1 + h1 t + 0.5 (h1 t)^2
    // exp( h2 t) ~ 1 + h2 t + 0.5 (h2 t)^2
    // exp( 0.5 h1 t) exp( h2 t) exp( 0.5 h1 t) 
    // ~ 1 + h1 t + h2 t + 0.5 (0.5 h1 t)^2 + 0.5 (h2 t)^2 + 0.5 (0.5 h1 t)^2 
    //   + 0.5 h1 h2 t^2 + 0.5 h2 h1 t^2 + ( 0.5  h1 t)^2
    // =  1 + h1 t + h2 t + 0.5 (h1 t)^2 + 0.5 (h2 t)^2 + 0.5 (h1 h2 + h2 h1)t^2 
    // ~ exp( (h1 + h2)t)

    auto psi = sm->fPsi;
    bool const fold_head_tail = false;
    for (size_t i = 0; i < sm->fBatchSize; ++i) {
        if (sm->fSolverMethod == SolverMethod::SplittingMethodO2) {
            if (i == 0 || !fold_head_tail) {
                sm->ExpV(psi, 0.5);
            }

            sm->ExpT(psi, 1.);

            if (i == sm->fBatchSize - 1 || !fold_head_tail) {
                sm->ExpV(psi, 0.5);
            } else {
                sm->ExpV(psi, 1.);
            }

        } else if (sm->fSolverMethod == SolverMethod::SplittingMethodO4) {
            if (i == 0 || !fold_head_tail) {
                sm->ExpV(psi, SplitingConstants<Real>::C1);
            }

            sm->ExpT(psi, SplitingConstants<Real>::D1);
            sm->ExpV(psi, SplitingConstants<Real>::C2);
            sm->ExpT(psi, SplitingConstants<Real>::D2);
            sm->ExpV(psi, SplitingConstants<Real>::C2);
            sm->ExpT(psi, SplitingConstants<Real>::D1);

            if (i == sm->fBatchSize - 1 || !fold_head_tail) {
                sm->ExpV(psi, SplitingConstants<Real>::C1);
            } else {
                sm->ExpV(psi, 2 * SplitingConstants<Real>::C1);
            }
        } else {
            throw std::runtime_error("unknown method");
        }
        sm->fStep += 1;
    }

}
