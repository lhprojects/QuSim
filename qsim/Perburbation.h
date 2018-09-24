#pragma once

typedef double PeReal;
#include <math.h>
PeReal const PI = 3.141592653589793;

#include <complex>
typedef std::complex<PeReal> PeComp;


namespace PerburbationUtility {

	inline PeComp ITime(PeComp const &c)
	{
		// I * (re, im) = (-im, re)
		return PeComp(-imag(c), real(c));
	}
	inline PeComp ITime(PeReal re, PeReal im)
	{
		// I * (re, im) = (-im, re)
		return PeComp(-im, re);
	}

	inline void Add(PeComp const *a, PeComp *b, size_t sz)
	{
		for (size_t i = 0; i < sz; ++i) {
			b[i] += a[i];
		}
	}

	inline void Mul(PeComp const *a, PeComp *b, size_t sz)
	{
		for (size_t i = 0; i < sz; ++i) {
			b[i] *= a[i];
		}
	}

	inline void Mul(PeReal const *a, PeComp *b, size_t sz)
	{
		for (size_t i = 0; i < sz; ++i) {
			b[i] *= a[i];
		}
	}

	inline void Mul(PeReal const *re, PeReal const *im, PeComp *b, size_t sz)
	{
		for (size_t i = 0; i < sz; ++i) {
			b[i] *= PeComp(re[i], im[i]);
		}
	}

};



struct BornSerise {

	// S = V psi0
	// deltaPsi = psi - psi0
	// ReV = Re(V)
	// ImV = Im(V)
	template<class F1, class F2, class F3>
	void Update(size_t n,
		PeComp const *psi0x,
		PeComp *deltaPsix,
		PeComp *deltaPsik,
		PeReal const *V,
		F1 const &G0,
		F2 const &X2K,
		F3 const &K2X)
	{
		using namespace PerburbationUtility;

		// psi  =  G (V psi + S)
		Add(psi0x, deltaPsix, n);
		Mul(V, deltaPsix, n);
		X2K(deltaPsix, deltaPsik);
		G0(deltaPsik);
		K2X(deltaPsik, deltaPsix);
	}

	template<class F1, class F2>
	void Update1D(size_t n,
		PeComp const *psi0x,
		PeComp *deltaPsix,
		PeComp *deltaPsik,
		PeReal const *V,
		PeReal epsilon,
		PeReal e,
		PeReal mass,
		PeReal hbar,
		PeReal dx,
		F1 const &X2K,
		F2 const &K2X)
	{
		auto G0 = [&](PeComp *psik) {
			PeReal dp = hbar * 2 * PI / (n * dx);
			PeReal divByTwoMass = 1 / (2 * mass);
			for (size_t i = 0; i < n; ++i) {
				PeReal p = (i < n / 2) ? (i * dp) : ((n - i) *dp);
				PeReal t = p * p * divByTwoMass;
				PeComp g0 = PeReal(1) / (e - t + PeComp(0, epsilon));
				psik[i] = g0 * psik[i];
			}

		};
		Update(n, psi0x, deltaPsix, deltaPsik, V, G0, X2K, K2X);
	}

};

enum class BornSerisePreconditioner {
	Vellekoop,
	Hao1,
	Hao2,
};

struct PreconditionalBornSerise {

	PeReal GetMinEpsilon(size_t n, PeReal const *reV, PeReal const *imV)
	{
		PeReal minEpsilon = 0;
		for (int i = 0; i < n; ++i) {
			if (abs(PeComp(reV[i], imV[i])) > minEpsilon) {
				minEpsilon = abs(PeComp(reV[i], imV[i]));
			}
		}
		return minEpsilon;
	}

	// S = V psi0
	// deltaPsi = psi - psi0
	// ReV = Re(V)
	// ImV = Im(V)
	template<class F1, class F2, class F3>
	void Update(size_t n,
		PeComp const *psi0x,
		PeComp *deltaPsix,
		PeComp *deltaPsik,
		PeReal const *reV,
		PeReal const *imV,
		PeComp *tmp,
		PeReal epsilon,
		F1 const &G0,
		F2 const &X2K,
		F3 const &K2X,
		PeReal const slow,
		BornSerisePreconditioner preconditioner)
	{
		// new psi = (gamma G V - gamma + 1) psi  + gamma G S
		// new psi = gamma G (V psi + S) + (1 - gamma) psi

		/////////////////////////////////////////
		// In Ref: A convergent Born series for solving the inhomogeneous Helmholtz equation in arbitrarily large media
		// V = V0 - I epsilon
		// epsilon >= abs(V0)
		// gamma = I V / espilon
		// new psi = gamma G (V psi + S) + (- I V0 / epsilon) psi

		// Transform:
		// V0 = -2 V0_
		// epsilon = 2 epsilon_

		// In This implementation:
		// V_ =  V0_ + I epsilon_
		// epsilon_ >= abs(V0_)
		// gamma_ = 1 - I V0_ /epsilon_

		// new psi = gamma_ G_ ((V0_ + I epsilon_) psi + S) + (I V0_ / epsilon_) psi
		/////////////////////////////////////////


		/////////////////////////////////////////
		// My preconditioner
		// first:
		//     gamma =  1 - I conj(V0)/ epsilon
		// second:
		//     0 = 1 -  1/2 gamma (1 + I V0 / epsilon)
		/////////////////////////////////////////


		using namespace PerburbationUtility;
		for (size_t i = 0; i < n; ++i) {
			tmp[i] = PeComp(reV[i], imV[i] + epsilon) * deltaPsix[i] + reV[i] * psi0x[i];
		}
		X2K(tmp, deltaPsik);
		G0(deltaPsik);
		K2X(deltaPsik, tmp);
		if (preconditioner == BornSerisePreconditioner::Vellekoop) {
			for (size_t i = 0; i < n; ++i) {
				PeComp delta = ITime(reV[i], imV[i]) / epsilon;
				PeComp gamma = slow * (PeReal(1) - delta);
				PeComp oneMinusGamma = (PeReal(1) - slow) + slow * delta;
				deltaPsix[i] = gamma * tmp[i] + oneMinusGamma * deltaPsix[i];
			}
		} else if (preconditioner == BornSerisePreconditioner::Hao1) {
			for (size_t i = 0; i < n; ++i) {
				PeComp delta = ITime(reV[i], -imV[i]) / epsilon;
				PeComp gamma = slow * (PeReal(1) - delta);
				PeComp oneMinusGamma = (PeReal(1) - slow) + slow * delta;
				deltaPsix[i] = gamma * tmp[i] + oneMinusGamma * deltaPsix[i];
			}
		} else if (preconditioner == BornSerisePreconditioner::Hao2) {
			for (size_t i = 0; i < n; ++i) {
				PeComp f = PeReal(1) + ITime(reV[i], imV[i]) / epsilon;
				PeComp gamma = PeReal(2) * slow / f;
				PeComp oneMinusGamma = PeReal(1) - gamma;
				deltaPsix[i] = gamma * tmp[i] + oneMinusGamma * deltaPsix[i];
			}
		}
	}

	template<class F1, class F2>
	void Update1D(size_t n,
		PeComp const *psi0x,
		PeComp *deltaPsix,
		PeComp *deltaPsik,
		PeReal const *reV,
		PeReal const *imV,
		PeComp *tmp,
		PeReal epsilon,
		PeReal e,
		PeReal mass,
		PeReal hbar,
		PeReal dx,
		F1 const &X2K,
		F2 const &K2X,
		PeReal const slow,
		BornSerisePreconditioner preconditioner)
	{
		auto G0 = [&](PeComp *psik) {
			PeReal dp = hbar * 2 * PI / (n * dx);
			PeReal divByTwoMass = 1 / (2 * mass);
			for (size_t i = 0; i < n; ++i) {
				PeReal p = (i < n / 2) ? (i * dp) : ((n - i) *dp);
				PeReal t = p * p * divByTwoMass;

				PeComp g0 = PeReal(1) / (e - t + PeComp(0, epsilon));
				psik[i] = g0 * psik[i];
			}

		};
		Update(n, psi0x, deltaPsix, deltaPsik, reV, imV, tmp, epsilon, G0, X2K, K2X, slow, preconditioner);
	}

};


