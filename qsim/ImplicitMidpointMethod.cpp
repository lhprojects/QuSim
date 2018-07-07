#include "ImplicitMidpointMethod.h"
#include "eigen/Eigen/SparseLU"

void ImplicitMidpointMethod::initSystem1D(char const * psi, bool force_normalization, Complex dt,
	bool force_normalization_each_step, char const * vs, Real x0, Real x1, size_t n,
	BoundaryCondition b, SolverMethod solver, Real mass, Real hbar)
{
	SystemImpl1D::initSystem1D(psi, force_normalization,
		dt, force_normalization_each_step,
		vs, x0, x1, n, b, solver,
		mass, hbar);

	if (b != BoundaryCondition::Period) {
		throw std::runtime_error("not supported boundary condition");
	}

	//          (1 - 1/2 i H Dt /hbar )
	//  psi  ->  ----------------------- psi
	//          (1 + 1/2 i H Dt /hbar )

	// H looks like
	//
	// 1 + 1/2 i (V + hbar^2/2m 2/(Dx^2)) Dt / hbar        + 1/2 i (- hbar^2/2m 1/(Dx^2)) Dt / hbar                         0                          ...     
	//   + 1/2 i (- hbar^2/2m 1/(Dx^2)) Dt / hbar       1 + 1/2 i (V + hbar^2/2m 2/(Dx^2)) Dt / hbar      + 1/2 i (- hbar^2/2m 1/(Dx^2)) Dt / hbar     ...
	//                   0                                + 1/2 i (- hbar^2/2m 1/(Dx^2)) Dt / hbar       1 + 1/2 i (V + hbar^2/2m 2/(Dx^2)) Dt / hbar  ...
	//                                                                    ...

	fh.resize(fN, fN);
	for (int i = 0; i < fN; ++i) {
		fh.insert(i, i) = 1.0 + 0.5*I*(fV[i] + hbar*hbar/(2*fMass) * 2/(fDx*fDx) )*fDt/hbar;
		fh.insert(i, i + 1 >= fN ? 0 : i + 1) = 0.5*I*(hbar * hbar / (2 * fMass) * (-1) / (fDx*fDx))*fDt / hbar;
		fh.insert(i, i - 1 < 0 ? fN - 1 : i - 1) = 0.5*I*(hbar * hbar / (2 * fMass) * (-1) / (fDx*fDx))*fDt / hbar;
	}
	
	fLU.compute(fh);

	fhPsi.resize(fN);
}

void ImplicitMidpointMethod::update_psi()
{
	//          
	//  psi  -> (1 + 1/2 i H Dt /hbar )^-1 (1 - 1/2 i H Dt /hbar ) psi
	//          

	Complex c = 0.5*I*(fHbar * fHbar / (2 * fMass) * 1 / (fDx*fDx))*fDt / fHbar;
	for (int i = 0; i < fN; ++i) {
		Complex dg = 1.0 - 0.5*fV[i] * fDt / fHbar*I - 2.0*c;
		Complex off =  c;
		fhPsi(i) = dg * fPsi[i];
		fhPsi(i) += off * fPsi[i + 1 < fN ? i + 1 : 0];
		fhPsi(i) += off * fPsi[i - 1 >= 0 ? i - 1 : fN - 1];
	}

	fhPsi = fLU.solve(fhPsi);

	for (int i = 0; i < fN; ++i) {
		fPsi[i] = fhPsi[i];
	}

	
}
