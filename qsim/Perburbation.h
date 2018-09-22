#pragma once

#include <complex>

typedef double PeReal;
typedef std::complex<PeReal> PeComp;
PeReal const PI = 3.141592653589793;

namespace PerburbationUtility {

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
