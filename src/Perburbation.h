#pragma once

typedef double PeReal;
#include <math.h>
PeReal const PI = 3.141592653589793;

#include <complex>
typedef std::complex<PeReal> PeComp;

#include "Device.h"
#include "QuSim.h"
#include "FourierTransform.h"

namespace PerburbationUtility {

	template<class R>
	static R CalDp2Div2M(R dx, size_t n, R hbar, R mass)
	{
		R Dpx = 2 * QuPi / (dx * n) * hbar;
		R Dpx2Div2M = Dpx * Dpx / (2 * mass);
		return Dpx2Div2M;
	}

	template<class R>
	static R CalLambda(R e, R mass, R hbar)
	{
		return 2 * PI / (sqrt(2 * e * mass) / hbar);
	}

	inline size_t fold_half(size_t i, size_t n) {
		if (i < n / 2) return i;
		else return n - i;
	};

    // total width = 20 / k0 = 3.3 lambda
    inline PeReal cal_v10(PeReal x, PeReal k0, PeReal e0)
    {
        // for e=0.5 e0:  relfaction=2.0E-6, through=1.9E-4
        // for e=0.2 e0:  relfaction=5.3E-4, through=9.7E-6
        x *= k0;
        return (-e0) * ((1 - 0.01 * x * x) > 0 ? exp(-0.03 * x * x * 1. / sqrt(1 - 0.01 * x * x)) : 0);
    };

    // total width = 30 / k0 = 5.0 lambda
    inline PeReal cal_v15(PeReal x, PeReal k0, PeReal e0)
    {
        // for e=0.5 e0:  relfaction=1.1E-7, through=4.6E-5
        // for e=0.2 e0:  relfaction=1.5E-4, through=8.7E-7
        x *= k0;
        return (-0.8 * e0) * ((1 - x * x / (15. * 15)) > 0 ? exp(-0.015 * x * x * 1. / sqrt(1 - x * x / (15. * 15))) : 0);
    };

    // total width = 40 / k0 = 6.4 lambda
    inline PeReal cal_v20(PeReal x, PeReal k0, PeReal e0)
    {
        // for e=0.5 e0:  relfaction=6.4E-8, through=3.1E-6
        // for e=0.2 e0:  relfaction=1.8E-4, through=1.5E-8
        x *= k0;
        return (-0.7 * e0) * ((1 - x * x / (20. * 20)) > 0 ? exp(-0.007 * x * x * 1. / sqrt(1 - x * x / (20. * 20))) : 0);
    };

	inline void GaussAsbLayer1D(size_t nx, PeReal dx,
		PeComp *v, PeReal hbar, PeReal mass, PeReal e,
		PeReal alpha)
	{
		PeReal lambda = CalLambda(e, mass, hbar);
		PeReal betax = dx / (alpha * lambda);
		for (size_t i = 0; i < nx; ++i) {
			PeReal xx;

			if (i < nx / 2) {
				xx = i * betax;
			} else {
				xx = (nx - i) * betax;
			}

            v[i].imag(v[i].imag() - e * exp(-(xx * xx)));
		}
	}

	inline void GaussMAsbLayer1D(size_t nx, PeReal dx,
        PeComp* v, PeReal hbar, PeReal mass, PeReal e,
        PeReal alpha)
    {
        PeReal k0 = sqrt(2 * mass * e) / hbar;
        for (size_t i = 0; i < nx; ++i) {
            PeReal xx;
            PeReal x = fold_half(i, nx) * dx;
            auto idx = CalGlobalIdx(i, nx);

            v[idx].imag(v[idx].imag() +
                (cal_v15(x, k0, e)));

		}
	}

	inline void GaussAsbLayer2D(size_t nx, size_t ny, PeReal dx, PeReal dy,
		PeComp *v, PeReal hbar, PeReal mass, PeReal e,
		PeReal alpha)
	{
		PeReal lambda = CalLambda(e, mass, hbar);
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

                auto idx = i * ny + j;
                v[idx].imag(v[idx].imag() - e * (exp(-xx * xx) + exp(-yy * yy)));
			}
		}

	}

	inline void GaussMAsbLayer2D(size_t nx, size_t ny, PeReal dx, PeReal dy,
		PeComp* v, PeReal hbar, PeReal mass, PeReal e,
		PeReal alpha)
	{
		PeReal k0 = sqrt(2 * mass * e) / hbar;

		for (size_t i = 0; i < nx; ++i) {
			for (size_t j = 0; j < ny; ++j) {

				PeReal x = fold_half(i, nx) * dx;
				PeReal y = fold_half(j, ny) * dy;
				auto idx = CalGlobalIdx(i, j, nx, ny);

                v[idx].imag(v[idx].imag() +
                    (cal_v15(x, k0, e) + cal_v15(y, k0, e)));
			}
		}
	}

	inline void GaussMAsbLayer3D(size_t nx, size_t ny, size_t nz,
		PeReal dx, PeReal dy, PeReal dz,
		PeComp* v, PeReal hbar, PeReal mass, PeReal e,
		PeReal alpha)
	{
		PeReal k0 = sqrt(2 * mass * e) / hbar;

        for (size_t i = 0; i < nx; ++i) {
            for (size_t j = 0; j < ny; ++j) {
                for (size_t k = 0; k < nz; ++k) {

                    PeReal x = fold_half(i, nx) * dx;
                    PeReal y = fold_half(j, ny) * dy;
                    PeReal z = fold_half(k, nz) * dz;

                    auto idx = CalGlobalIdx(i, j, k, nx, ny, nz);

                    if (0) {
                        v[idx].imag(v[idx].imag() +
                            (cal_v10(x, k0, e) + cal_v10(y, k0, e) + cal_v10(z, k0, e)));
                    } else {
                        v[idx].imag(v[idx].imag() +
                            (cal_v15(x, k0, e) + cal_v15(y, k0, e) + cal_v15(z, k0, e)));
                    }
                }
            }
        }
	}

	inline void GaussAsbLayer3D(size_t nx, size_t ny, size_t nz,
		PeReal dx, PeReal dy, PeReal dz,
		PeComp*v, PeReal hbar, PeReal mass, PeReal e,
		PeReal alpha)
	{
		PeReal lambda = CalLambda(e, mass, hbar);
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
					auto idx = (i * ny + j) * nz + k;
                    v[idx].imag(v[idx].imag() - e * (exp(-xx * xx) + exp(-yy * yy) + exp(-zz * zz)));
				}
			}
		}

	}

};

// naive born series
struct BornSeries {


	// (T + V - E) (deltaPsi + psi0) = 0 with scattering condition
	// (E - T) deltaPsi = V (psi0 + deltaPsi)
	// (E - T + i eps) deltaPsi = V (psi0 + deltaPsi) replace `with scattering condition` with `eps -> +0`
	// deltaPsi = (E - T + i eps)^-1 V (psi0 + deltaPsi)
	// deltaPsi = G0 V (psi0 + deltaPsi)

	void Update1D(size_t n,    //    in:
		PeComp const* psi0x,   //    in:
		PeComp* deltaPsix,     //inout: psi - psi0, could be zeros
		PeComp* deltaPsik,     //    out: 
		PeComp const* V,       //    in: real potential
		PeReal epsilon,        //    in: should be tiny
		PeReal e,              //    in: energy
		PeReal mass,           //    in: mass
		PeReal hbar,           //    in: hbar
		PeReal dx,             //    in: delta x
		FourierTransform1D* fft,
		FourierTransform1D* invfft,
		Device *dev)           //    in:
	{
		using namespace PerburbationUtility;
		PeReal Dpx2Div2M = CalDp2Div2M(dx, n, hbar, mass);

        dev->AddMul(deltaPsix, psi0x, V, n);
        fft->Transform(deltaPsix, deltaPsik);
        dev->G01D(deltaPsik, e, epsilon, Dpx2Div2M, n);
        dev->Scale(deltaPsik, 1. / n, n);
        invfft->Transform(deltaPsik, deltaPsix);

	}

	void Update2D(size_t nx, size_t ny,
		PeComp const *psi0x,
		PeComp *deltaPsix,
		PeComp *deltaPsik,
		PeComp const *V,
		PeReal epsilon,
		PeReal e,
		PeReal mass,
		PeReal hbar,
		PeReal dx,
		PeReal dy,
		FourierTransform2D* fft,
		FourierTransform2D* invfft,
		Device* dev)
	{
		using namespace PerburbationUtility;
		PeReal Dpx2Div2M = CalDp2Div2M(dx, nx, hbar, mass);
		PeReal Dpy2Div2M = CalDp2Div2M(dy, ny, hbar, mass);

        dev->AddMul(deltaPsix, psi0x, V, nx * ny);
        fft->Transform(deltaPsix, deltaPsik);
        dev->G02D(deltaPsik, e, epsilon, Dpx2Div2M, Dpy2Div2M, nx, ny);
        dev->Scale(deltaPsik, 1. / (nx * ny), nx * ny);
		invfft->Transform(deltaPsik, deltaPsix);
	}

	void Update3D(size_t nx, size_t ny, size_t nz,
		PeComp const *psi0x,
		PeComp *deltaPsix,
		PeComp *deltaPsik,
		PeComp const *V,
		PeReal epsilon,
		PeReal e,
		PeReal mass,
		PeReal hbar,
		PeReal dx,
		PeReal dy,
		PeReal dz,
		FourierTransform3D* fft,
		FourierTransform3D* invfft,
		Device* dev)
	{
		using namespace PerburbationUtility;
		PeReal Dpx2Div2M = CalDp2Div2M(dx, nx, hbar, mass);
        PeReal Dpy2Div2M = CalDp2Div2M(dy, ny, hbar, mass);
        PeReal Dpz2Div2M = CalDp2Div2M(dz, nz, hbar, mass);

        dev->AddMul(deltaPsix, psi0x, V, nx * ny * nz);
        fft->Transform(deltaPsix, deltaPsik);
		dev->G03D(deltaPsik, e, epsilon, Dpx2Div2M, Dpy2Div2M, Dpz2Div2M, nx, ny, nz);
        dev->Scale(deltaPsik, 1. / (nx * ny * nz), nx * ny * nz);
		invfft->Transform(deltaPsik, deltaPsix);
	}

};

enum class BornSerisePreconditioner {
	Identity,
	Vellekoop,
	Hao1,
	Hao2,
};


struct PreconditionalBornSeries {

	template<class U>
	PeReal Abs2(U u)
	{
		return u.real() * u.real() + u.imag() * u.imag();
	}

	void GetEpsilon(PeReal &epsilon, BornSerisePreconditioner pre,
		size_t n, PeComp const *v, Device *dev)
	{
        PeReal minEpsilon = 0;
        if (pre == BornSerisePreconditioner::Identity) {
            PeReal A = dev->SumImag(v, n);
            PeReal A2 = dev->Norm2(v, n);
            minEpsilon = 0.5 * A2 / (-A);
        } else {
			minEpsilon = dev->Max(v, n);
        }

        if (epsilon < minEpsilon) {
            epsilon = minEpsilon;
        }
	}

	// S = V psi0
	// deltaPsi = psi - psi0
	// ReV = Re(V)
	// ImV = Im(V)
    template<class F1, class F2, class F3>
    void Update(size_t n,
        PeComp const* psi0x,
        PeComp* deltaPsix,
        PeComp* deltaPsik,
		PeComp const* v,
        PeComp* tmp,
        PeReal epsilon,
        F1 const& G0,
        F2 const& X2K,
        F3 const& K2X,
        PeReal const slow,
        BornSerisePreconditioner preconditioner,
        Device* dev
		)
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

        // tmp = (deltaPsix + psi0x)*(v + I epsilon) + v.real()*psi0x
        dev->Add3(tmp, deltaPsix, psi0x, v, epsilon, n);
		
        X2K(tmp, deltaPsik);
		G0(deltaPsik);
		K2X(deltaPsik, tmp);

        if (preconditioner == BornSerisePreconditioner::Identity) {
			if (slow == 1.) {
				dev->Scale(deltaPsix, tmp, 1. / n, n);
			} else {
				dev->Scale(tmp, 1. / n, n);
				dev->LinearUpdate(deltaPsix, tmp, slow, n);
			}
        } else if (preconditioner == BornSerisePreconditioner::Vellekoop) {
			dev->Scale(tmp, 1. / n, n);
			dev->Vellekoop(deltaPsix, tmp, v, slow, epsilon, n);
        } else if (preconditioner == BornSerisePreconditioner::Hao1) {
			dev->Scale(tmp, 1. / n, n);
			dev->Hao1(deltaPsix, tmp, v, slow, epsilon, n);
        } else if (preconditioner == BornSerisePreconditioner::Hao2) {
			dev->Scale(tmp, 1. / n, n);
			dev->Hao2(deltaPsix, tmp, v, slow, epsilon, n);
        } else {
            throw std::invalid_argument("Uknown preconditioner");
        }
    }

	void Update1D(size_t n,
		PeComp const *psi0x,
		PeComp *deltaPsix,
		PeComp *deltaPsik,
		PeComp const *v,
		PeComp *tmp,
		PeReal epsilon,
		PeReal e,
		PeReal mass,
		PeReal hbar,
		PeReal dx,
		PeReal const slow,
		BornSerisePreconditioner preconditioner,
		FourierTransform1D* fft,
		FourierTransform1D* invfft,
		Device* dev)
	{
		using namespace PerburbationUtility;
		auto G0 = [&](PeComp *psik) {
			PeReal Dp2Div2M = CalDp2Div2M(dx, n, hbar, mass);
			dev->G01D(psik, e, epsilon, Dp2Div2M, n);
		};
		
		auto X2K = [&](PeComp const* psix, PeComp* psik) {
			fft->Transform(psix, psik);
		};

		auto K2X = [&](PeComp const* psik, PeComp* psix) {
			invfft->Transform(psik, psix);
		};

		Update(n, psi0x, deltaPsix, deltaPsik, v, tmp, epsilon, G0, X2K, K2X, slow, preconditioner, dev);
	}

	void Update2D(size_t nx, size_t ny,
		PeComp const *psi0x,
		PeComp *deltaPsix,
		PeComp *deltaPsik,
		PeComp const* v,
		PeComp *tmp,
		PeReal epsilon,
		PeReal e,
		PeReal mass,
		PeReal hbar,
		PeReal dx,
		PeReal dy,
		PeReal const slow,
		BornSerisePreconditioner preconditioner,
		FourierTransform2D* fft,
		FourierTransform2D* invfft,
		Device *dev)
	{
		using namespace PerburbationUtility;
		auto G0 = [&](PeComp* psik) {
            PeReal Dpx2Div2M = CalDp2Div2M(dx, nx, hbar, mass);
            PeReal Dpy2Div2M = CalDp2Div2M(dy, ny, hbar, mass);
            dev->G02D(psik, e, epsilon, Dpx2Div2M, Dpy2Div2M, nx, ny);
        };

        auto X2K = [&](PeComp const* psix, PeComp* psik) {
            fft->Transform(psix, psik);
        };

        auto K2X = [&](PeComp const* psik, PeComp* psix) {
            invfft->Transform(psik, psix);
        };

        Update(nx * ny, psi0x, deltaPsix, deltaPsik, v, tmp, epsilon, G0, X2K, K2X, slow, preconditioner, dev);
	}

	void Update3D(size_t nx, size_t ny, size_t nz,
		PeComp const *psi0x,
		PeComp *deltaPsix,
		PeComp *deltaPsik,
		PeComp const* v,
		PeComp *tmp,
		PeReal epsilon,
		PeReal e,
		PeReal mass,
		PeReal hbar,
		PeReal dx,
		PeReal dy,
		PeReal dz,
		PeReal const slow,
		BornSerisePreconditioner preconditioner,
		FourierTransform3D* fft,
		FourierTransform3D* invfft,
		Device *dev)
	{
		using namespace PerburbationUtility;
		auto G0 = [&](PeComp* psik) {
			PeReal Dpx2Div2M = CalDp2Div2M(dx, nx, hbar, mass);
			PeReal Dpy2Div2M = CalDp2Div2M(dy, ny, hbar, mass);
			PeReal Dpz2Div2M = CalDp2Div2M(dz, nz, hbar, mass);
			dev->G03D(psik, e, epsilon, Dpx2Div2M, Dpy2Div2M, Dpz2Div2M, nx, ny, nz);
		};

		auto X2K = [&](PeComp const* psix, PeComp* psik) {
			fft->Transform(psix, psik);
		};

		auto K2X = [&](PeComp const* psik, PeComp* psix) {
			invfft->Transform(psik, psix);
		};

        Update(nx * ny * nz, psi0x, deltaPsix, deltaPsik, v, tmp, epsilon, G0, X2K, K2X, slow, preconditioner, dev);
	}

};


