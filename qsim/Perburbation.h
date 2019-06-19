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

	// psik = G0 * psik
	inline void Green01D(PeComp *psik, PeReal hbar, PeReal mass, PeReal e, PeReal epsilon,
		size_t n, PeReal dx)
	{
		PeReal dp = hbar * 2 * PI / (n * dx);
		PeReal divByTwoMass = PeReal(1) / (PeReal(2) * mass);
		for (size_t i = 0; i < n; ++i) {
			PeReal p = (i < n / 2) ? (i * dp) : ((n - i) *dp);
			PeReal t = p * p * divByTwoMass;
			PeComp g0 = PeReal(1) / (e - t + PeComp(0, epsilon));
			psik[i] = g0 * psik[i];
		}
	}

	// psik = G0 * psik
	inline void Green02D(PeComp *psik, PeReal hbar, PeReal mass, PeReal e, PeReal epsilon,
		size_t nx, size_t ny, PeReal dx, PeReal dy)
	{
		PeReal dpx = hbar * 2 * PI / (nx * dx);
		PeReal dpy = hbar * 2 * PI / (ny * dy);
		PeReal divByTwoMass = PeReal(1) / (PeReal(2) * mass);
		for (size_t i = 0; i < nx; ++i) {
			for (size_t j = 0; j < ny; ++j) {
				PeReal px = (i < nx / 2) ? (i * dpx) : ((nx - i) *dpx);
				PeReal py = (j < ny / 2) ? (j * dpy) : ((ny - j) *dpy);
				PeReal t = (px * px + py * py) * divByTwoMass;

				PeComp g0 = PeReal(1) / (e - t + PeComp(0, epsilon));
				psik[i*ny + j] = g0 * psik[i*ny + j];
			}
		}
	}

	// psik = G0 * psik
	inline void Green03D(PeComp *psik, PeReal hbar, PeReal mass, PeReal e, PeReal epsilon,
		size_t nx, size_t ny, size_t nz,
		PeReal dx, PeReal dy, PeReal dz)
	{
		PeReal dpx = hbar * 2 * PI / (nx * dx);
		PeReal dpy = hbar * 2 * PI / (ny * dy);
		PeReal dpz = hbar * 2 * PI / (nz * dz);
		PeReal divByTwoMass = PeReal(1) / (PeReal(2) * mass);
		for (size_t i = 0; i < nx; ++i) {
			PeReal px = (i < nx / 2) ? (i * dpx) : ((nx - i) *dpx);

			for (size_t j = 0; j < ny; ++j) {
				PeReal py = (j < ny / 2) ? (j * dpy) : ((ny - j) *dpy);

				for (size_t k = 0; k < ny; ++k) {
					PeReal pz = (k < nz / 2) ? (k * dpz) : ((nz - k) *dpz);

					PeReal t = (px * px + py * py + pz * pz) * divByTwoMass;

					PeComp g0 = PeReal(1) / (e - t + PeComp(0, epsilon));
					psik[(i*ny + j)*nz + k] = g0 * psik[(i*ny + j)*nz + k];
				}
			}
		}
	}

	inline void GaussAsbLayer1D(size_t nx, PeReal dx,
		PeReal *v, PeReal hbar, PeReal mass, PeReal e,
		PeReal alpha)
	{
		PeReal lambda = 2 * PI / (sqrt(2 * e * mass) / hbar);
		PeReal betax = dx / (alpha * lambda);
		for (size_t i = 0; i < nx; ++i) {
			PeReal xx;

			if (i < nx / 2) {
				xx = i * betax;
			} else {
				xx = (nx - i) * betax;
			}
			v[i] = -e * exp(-(xx * xx));
		}
	}

	inline void GaussAsbLayer2D(size_t nx, size_t ny, PeReal dx, PeReal dy,
		PeReal *v, PeReal hbar, PeReal mass, PeReal e,
		PeReal alpha)
	{
		PeReal lambda = 2 * PI / (sqrt(2 * e * mass) / hbar);
		PeReal betax = dx / (alpha * lambda);
		PeReal betay = dy / (alpha * lambda);
		for (size_t i = 0; i < nx; ++i) {
			for (size_t j = 0; j < ny; ++j) {
				PeReal xx;
				PeReal yy;

				if (i < nx / 2) {
					xx = i * betax;
				} else {
					xx = (nx -  i)*betax;
				}

				if (j < ny / 2) {
					yy = j * betay;
				} else {
					yy = (ny - j) * betay;
				}

				v[i*ny + j] = -e * (exp(-xx * xx) + exp(-yy * yy));
			}
		}

	}

	inline void GaussAsbLayer3D(size_t nx, size_t ny, size_t nz,
		PeReal dx, PeReal dy, PeReal dz,
		PeReal *v, PeReal hbar, PeReal mass, PeReal e,
		PeReal alpha)
	{
		PeReal lambda = 2 * PI / (sqrt(2 * e * mass) / hbar);
		PeReal betax = dx / (alpha * lambda);
		PeReal betay = dy / (alpha * lambda);
		PeReal betaz = dz / (alpha * lambda);
		for (size_t i = 0; i < nx; ++i) {
			for (size_t j = 0; j < ny; ++j) {
				for (size_t k = 0; k < nz; ++k) {
					PeReal xx;
					PeReal yy;
					PeReal zz;

					if (i < nx / 2) {
						xx = i * betax;
					} else {
						xx = (nx - i)*betax;
					}

					if (j < ny / 2) {
						yy = j * betay;
					} else {
						yy = (ny - j) * betay;
					}

					if (k < nz / 2) {
						zz = k * betaz;
					} else {
						zz = (nz - k) * betaz;
					}

					v[(i*ny + j)*nz + k] = -e * (exp(-xx * xx) + exp(-yy * yy) + exp(-zz * zz));
				}
			}
		}

	}

};

// native born serise
struct BornSerise {


	// (T + V - E) (deltaPsi + psi0) = 0 with scattering condition
	// (E - T) deltaPsi = V (psi0 + deltaPsi)
	// (E - T + i eps) deltaPsi = V (psi0 + deltaPsi) replace `with scattering condition` with `eps -> +0`
	// deltaPsi = (E - T + i eps)^-1 V (psi0 + deltaPsi)
	// deltaPsi = G0 V (psi0 + deltaPsi)
	template<class F1, class F2, class F3>
	void Update(size_t n,    //	   in: x dimension
		PeComp const *psi0x, //	   in: psi0
		PeComp *deltaPsix,   // inout: deltaPsi = psi - psi0
		PeComp *deltaPsik,   //   out: deltaPsi
		PeReal const *V,     //    in: real potential
		F1 const &G0,        //    in: (E - T + i eps)^-1
		F2 const &X2K,       //    in:
		F3 const &K2X)       //    in:
	{
		using namespace PerburbationUtility;

		Add(psi0x, deltaPsix, n);
		Mul(V, deltaPsix, n);
		X2K(deltaPsix, deltaPsik);
		G0(deltaPsik);
		K2X(deltaPsik, deltaPsix);
	}

	template<class F1, class F2>
	void Update1D(size_t n,    //    in:
		PeComp const *psi0x,   //    in:
		PeComp *deltaPsix,     // inout: psi - psi0, could be zeros
		PeComp *deltaPsik,     //    out: 
		PeReal const *V,       //    in: real potential
		PeReal epsilon,        //    in: should be tiny
		PeReal e,              //    in: energy
		PeReal mass,           //    in: mass
		PeReal hbar,           //    in: hbar
		PeReal dx,             //    in: delta x
		F1 const &X2K,         //    in:
		F2 const &K2X)         //    in:
	{
		auto G0 = [&](PeComp *psik) {
			PerburbationUtility::Green01D(psik, hbar, mass, e, epsilon, n, dx);
		};
		Update(n, psi0x, deltaPsix, deltaPsik, V, G0, X2K, K2X);
	}

	template<class F1, class F2>
	void Update2D(size_t nx, size_t ny,
		PeComp const *psi0x,
		PeComp *deltaPsix,
		PeComp *deltaPsik,
		PeReal const *V,
		PeReal epsilon,
		PeReal e,
		PeReal mass,
		PeReal hbar,
		PeReal dx,
		PeReal dy,
		F1 const &X2K,
		F2 const &K2X)
	{
		auto G0 = [&](PeComp *psik) {
			PerburbationUtility::Green02D(psik, hbar, mass, e, epsilon, nx, ny, dx, dy);
		};
		Update(nx*ny, psi0x, deltaPsix, deltaPsik, V, G0, X2K, K2X);
	}

	template<class F1, class F2>
	void Update3D(size_t nx, size_t ny, size_t nz,
		PeComp const *psi0x,
		PeComp *deltaPsix,
		PeComp *deltaPsik,
		PeReal const *V,
		PeReal epsilon,
		PeReal e,
		PeReal mass,
		PeReal hbar,
		PeReal dx,
		PeReal dy,
		PeReal dz,
		F1 const &X2K,
		F2 const &K2X)
	{
		auto G0 = [&](PeComp *psik) {
			PerburbationUtility::Green03D(psik, hbar, mass, e, epsilon, nx, ny, nz, dx, dy, dz);
		};
		Update(nx*ny*nz, psi0x, deltaPsix, deltaPsik, V, G0, X2K, K2X);
	}

};

enum class BornSerisePreconditioner {
	Identity,
	Vellekoop,
	Hao1,
	Hao2,
};

struct PreconditionalBornSerise {

	void GetEpsilon(PeReal &epsilon, BornSerisePreconditioner pre, size_t n, PeReal const *reV, PeReal const *imV)
	{

		PeReal minEpsilon = 0;
		if (pre == BornSerisePreconditioner::Identity) {
			PeReal A = 0;
			PeReal A2 = 0;
			PeReal V2 = 0;
			for (int i = 0; i < n; ++i) {
				A += imV[i];
				A2 += imV[i] * imV[i];
				V2 += reV[i] * reV[i];
			}
			A2 /= n;
			V2 /= n;
			A /= n;

			minEpsilon = 0.5*(A2 + V2) / (-A);
			//input > 0.5 (A2 + V2) / (-A) 
			//epsilon = (A2 + V2) / (-A);
			//printf("%f %f %f\n", A2, V2, epsilon);
			//printf("%f\n", 1 - A*A/(A2+V2));
		}
		else {
			for (int i = 0; i < n; ++i) {
				if (abs(PeComp(reV[i], imV[i])) > minEpsilon) {
					minEpsilon = abs(PeComp(reV[i], imV[i]));
				}
			}
		}
		if (epsilon < minEpsilon) {
			epsilon = 2 * minEpsilon;
		}
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
		if (preconditioner == BornSerisePreconditioner::Identity) {
			for (size_t i = 0; i < n; ++i) {
				PeReal gamma = slow;
				PeReal oneMinusGamma = PeReal(1) - gamma;
				deltaPsix[i] = gamma * tmp[i] + oneMinusGamma * deltaPsix[i];
			}
		} else if (preconditioner == BornSerisePreconditioner::Vellekoop) {
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
			PerburbationUtility::Green01D(psik, hbar, mass, e, epsilon, n, dx);
		};
		Update(n, psi0x, deltaPsix, deltaPsik, reV, imV, tmp, epsilon, G0, X2K, K2X, slow, preconditioner);
	}

	template<class F1, class F2>
	void Update2D(size_t nx, size_t ny,
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
		PeReal dy,
		F1 const &X2K,
		F2 const &K2X,
		PeReal const slow,
		BornSerisePreconditioner preconditioner)
	{
		auto G0 = [&](PeComp *psik) {
			PerburbationUtility::Green02D(psik, hbar, mass, e, epsilon, nx, ny, dx, dy);
		};
		Update(nx * ny, psi0x, deltaPsix, deltaPsik, reV, imV, tmp, epsilon, G0, X2K, K2X, slow, preconditioner);
	}

	template<class F1, class F2>
	void Update3D(size_t nx, size_t ny, size_t nz,
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
		PeReal dy,
		PeReal dz,
		F1 const &X2K,
		F2 const &K2X,
		PeReal const slow,
		BornSerisePreconditioner preconditioner)
	{
		auto G0 = [&](PeComp *psik) {
			PerburbationUtility::Green03D(psik, hbar, mass, e, epsilon,
				nx, ny, nz,
				dx, dy, dz);
		};
		Update(nx * ny, psi0x, deltaPsix, deltaPsik, reV, imV, tmp, epsilon, G0, X2K, K2X, slow, preconditioner);
	}

};


