#pragma once

#include "eigen/Eigen/Sparse"
#include "QuSim.h"
#include "Utils.h"
#include "Device.h"

template<class R>
struct NumericalDiffUtils {

	// this numers checkout from
	// https://web.media.mit.edu/~crtaylor/calculator.html

    static R C1_6_L2Half() { return -9. / 1920; }
    static R C1_6_L1Half() { return +125. / 1920; }
    static R C1_6_L0Half() { return -2250. / 1920; }
	static R C1_6_R0Half() { return -C1_6_L0Half();  }
	static R C1_6_R1Half() { return -C1_6_L1Half(); }
	static R C1_6_R2Half() { return -C1_6_L2Half(); }

    static R C1_4_L1Half() { return +1. / 24.; }
    static R C1_4_L0Half() { return -27. / 24; }
    static R C1_4_R0Half() { return -C1_4_L0Half(); }
    static R C1_4_R1Half() { return -C1_4_L1Half(); }

	static Minus1 C1_2_L0Half() { return Minus1(); }
	static One C1_2_R0Half() { return -C1_2_L0Half(); }

	// d^2 f(x) / dx^2 =  (C2_1 f(x-dx) + C2_0 f(x) + C2_1 f(x+dx)) / dx^2
	static One C2_2_1() { return  One(); }
	static R C2_2_0() { return  -2.; }

	static R C2_4_2() { return  -1. / 12; }
	static R C2_4_1() { return 16. / 12; }
	static R C2_4_0() { return -30. / 12; }

	static R C2_6_3() { return 2. / 180; }
	static R C2_6_2() { return -27. / 180; }
	static R C2_6_1() { return 270. / 180; }
	static R C2_6_0() { return -490. / 180; }
};

template<class R>
R CalKinEnergy(Complex const *psi, size_t nx, int order, R hbar, R mass, R dx, Device *dev)
{
	// T = -hbar^2/2m nabla
	R const tx = hbar * hbar / (2 * mass) / (dx * dx);
	auto foldx = [=](size_t i) -> size_t {
		if ((ptrdiff_t)i < 0) i += nx;
		else if (i >= nx) i -= nx;
		return i;
	};

	R kin = 0;
	for (size_t i = 0; i < nx; ++i) {
		if (order <= 2) {
            auto psi_prime_x = NumericalDiffUtils<Real>::C1_2_L0Half() * psi[i] 
				+ NumericalDiffUtils<Real>::C1_2_R0Half() * psi[foldx(i + 1)];
            kin += Abs2(psi_prime_x) * tx;
		} else if (order <= 4) {
            auto psi_prime_x = NumericalDiffUtils<Real>::C1_4_L1Half() * psi[i]
                + NumericalDiffUtils<Real>::C1_4_L0Half() * psi[foldx(i + 1)]
            +NumericalDiffUtils<Real>::C1_4_R0Half() * psi[foldx(i + 2)]
            +NumericalDiffUtils<Real>::C1_4_R1Half() * psi[foldx(i + 3)];

			kin += Abs2(psi_prime_x) * tx;
		} else if (order <= 6) {
            auto psi_prime_x = NumericalDiffUtils<Real>::C1_6_L2Half() * psi[i]
                + NumericalDiffUtils<Real>::C1_6_L1Half() * psi[foldx(i + 1)]
            +NumericalDiffUtils<Real>::C1_6_L0Half() * psi[foldx(i + 2)]
            +NumericalDiffUtils<Real>::C1_6_R0Half() * psi[foldx(i + 3)]
            +NumericalDiffUtils<Real>::C1_6_R1Half() * psi[foldx(i + 4)]
            +NumericalDiffUtils<Real>::C1_6_R2Half() * psi[foldx(i + 5)];

			kin += Abs2(psi_prime_x) * tx;
		} else {
			throw std::runtime_error("too higher order");
		}
	}

	kin /= dev->Norm2(psi, nx);
	return kin;
}

// M = dig + f * (V + T)
template<class R, class F, class Potential, class Dig>
void FillElems1D(std::vector<Eigen::Triplet<Complex> >& elems,
	size_t nx,
	int order,
	Dig dig,
	Potential const* v,
	F f,
	R hbar,
	R mass,
	R dx)
{

	// T = -hbar^2/2m nabla
	R const tx = -hbar * hbar / (2 * mass) / (dx * dx);
	auto const ftx = f * tx;

	auto foldx = [=](size_t i) -> size_t {
		if ((ptrdiff_t)i < 0) i += nx;
		else if (i >= nx) i -= nx;
		return i;
	};

	for (size_t i = 0; i < nx; ++i) {
        auto dig_ = dig + f * v[i];
		if (order <= 2) {
			elems.push_back(Eigen::Triplet<Complex>((int)i, (int)foldx(i + 0), NumericalDiffUtils<Real>::C2_2_0() * ftx + dig_));
			elems.push_back(Eigen::Triplet<Complex>((int)i, (int)foldx(i - 1), NumericalDiffUtils<Real>::C2_2_1() * ftx));
			elems.push_back(Eigen::Triplet<Complex>((int)i, (int)foldx(i + 1), NumericalDiffUtils<Real>::C2_2_1() * ftx));
		} else if (order <= 4) {
			elems.push_back(Eigen::Triplet<Complex>((int)i, (int)foldx(i - 2), NumericalDiffUtils<Real>::C2_4_2() * ftx));
			elems.push_back(Eigen::Triplet<Complex>((int)i, (int)foldx(i - 1), NumericalDiffUtils<Real>::C2_4_1() * ftx));
			elems.push_back(Eigen::Triplet<Complex>((int)i, (int)foldx(i + 0), NumericalDiffUtils<Real>::C2_4_0() * ftx + dig_));
			elems.push_back(Eigen::Triplet<Complex>((int)i, (int)foldx(i + 1), NumericalDiffUtils<Real>::C2_4_1() * ftx));
			elems.push_back(Eigen::Triplet<Complex>((int)i, (int)foldx(i + 2), NumericalDiffUtils<Real>::C2_4_2() * ftx));
		} else if (order <= 6) {
			elems.push_back(Eigen::Triplet<Complex>((int)i, (int)foldx(i - 3), NumericalDiffUtils<Real>::C2_6_3() * ftx));
			elems.push_back(Eigen::Triplet<Complex>((int)i, (int)foldx(i - 2), NumericalDiffUtils<Real>::C2_6_2() * ftx));
			elems.push_back(Eigen::Triplet<Complex>((int)i, (int)foldx(i - 1), NumericalDiffUtils<Real>::C2_6_1() * ftx));
			elems.push_back(Eigen::Triplet<Complex>((int)i, (int)foldx(i + 0), NumericalDiffUtils<Real>::C2_6_0() * ftx + dig_));
			elems.push_back(Eigen::Triplet<Complex>((int)i, (int)foldx(i + 1), NumericalDiffUtils<Real>::C2_6_1() * ftx));
			elems.push_back(Eigen::Triplet<Complex>((int)i, (int)foldx(i + 2), NumericalDiffUtils<Real>::C2_6_2() * ftx));
			elems.push_back(Eigen::Triplet<Complex>((int)i, (int)foldx(i + 3), NumericalDiffUtils<Real>::C2_6_3() * ftx));
		} else {
			throw std::runtime_error("too higher order");
		}
	}

};


// M = dig + f * (V + T)
template<class R, class F, class Potential, class Dig>
void FillElems2D(std::vector<Eigen::Triplet<Complex> >& elems,
	size_t nx,
	size_t ny,
	int order,
	Dig dig,
	Potential const* v,
	F f,
	R hbar,
	R mass,
	R dx,
	R dy)
{

	// T = -hbar^2/2m nabla
    R const tx = -hbar * hbar / (2 * mass) / (dx * dx);
    R const ty = -hbar * hbar / (2 * mass) / (dy * dy);
	auto const ftx = f * tx;
	auto const fty = f * tx;

	auto foldx = [=](size_t i) -> size_t {
		if ((ptrdiff_t)i < 0) i += nx;
		else if (i >= nx) i -= nx;
		return i;
	};

	auto foldy = [=](size_t i) -> size_t {
		if ((ptrdiff_t)i < 0) i += ny;
		else if (i >= ny) i -= ny;
		return i;
	};

	auto global = [=](size_t x, size_t y)->size_t {
		return foldy(y) + foldx(x) * ny;
	};


	ForEach2DIdx(nx, ny, [=, &elems](size_t idx, size_t i, size_t j){
        auto dig_ = dig + f * v[idx];

		if (order <= 2) {

            elems.push_back(Eigen::Triplet<Complex>(idx, global(i, j),
                NumericalDiffUtils<double>::C2_2_0() * ftx
                + NumericalDiffUtils<double>::C2_2_0() * fty + dig_));
            elems.push_back(Eigen::Triplet<Complex>(idx, global(i - 1, j), NumericalDiffUtils<double>::C2_2_1() * ftx));
            elems.push_back(Eigen::Triplet<Complex>(idx, global(i + 1, j), NumericalDiffUtils<double>::C2_2_1() * ftx));
            elems.push_back(Eigen::Triplet<Complex>(idx, global(i, j - 1), NumericalDiffUtils<double>::C2_2_1() * fty));
            elems.push_back(Eigen::Triplet<Complex>(idx, global(i, j + 1), NumericalDiffUtils<double>::C2_2_1() * fty));

		} else if (order <= 4) {

			elems.push_back(Eigen::Triplet<Complex>(idx, global(i, j),
				NumericalDiffUtils<double>::C2_4_0() * ftx
				+ NumericalDiffUtils<double>::C2_4_0() * fty + dig_));
			
            elems.push_back(Eigen::Triplet<Complex>(idx, global(i - 2, j), NumericalDiffUtils<double>::C2_4_2() * ftx));
            elems.push_back(Eigen::Triplet<Complex>(idx, global(i - 1, j), NumericalDiffUtils<double>::C2_4_1() * ftx));
			elems.push_back(Eigen::Triplet<Complex>(idx, global(i + 1, j), NumericalDiffUtils<double>::C2_4_1() * ftx));
            elems.push_back(Eigen::Triplet<Complex>(idx, global(i + 2, j), NumericalDiffUtils<double>::C2_4_2() * ftx));

			elems.push_back(Eigen::Triplet<Complex>(idx, global(i, j - 2), NumericalDiffUtils<double>::C2_4_2() * fty));
			elems.push_back(Eigen::Triplet<Complex>(idx, global(i, j - 1), NumericalDiffUtils<double>::C2_4_1() * fty));
			elems.push_back(Eigen::Triplet<Complex>(idx, global(i, j + 1), NumericalDiffUtils<double>::C2_4_1() * fty));
			elems.push_back(Eigen::Triplet<Complex>(idx, global(i, j + 2), NumericalDiffUtils<double>::C2_4_2() * fty));
		} else if (order <= 6) {

			elems.push_back(Eigen::Triplet<Complex>(idx, global(i, j),
				NumericalDiffUtils<double>::C2_6_0() * ftx
				+ NumericalDiffUtils<double>::C2_6_0() * fty + dig_));
			elems.push_back(Eigen::Triplet<Complex>(idx, global(i - 3, j), NumericalDiffUtils<double>::C2_6_3() * ftx));
			elems.push_back(Eigen::Triplet<Complex>(idx, global(i - 2, j), NumericalDiffUtils<double>::C2_6_2() * ftx));
			elems.push_back(Eigen::Triplet<Complex>(idx, global(i - 1, j), NumericalDiffUtils<double>::C2_6_1() * ftx));
			elems.push_back(Eigen::Triplet<Complex>(idx, global(i + 1, j), NumericalDiffUtils<double>::C2_6_1() * ftx));
			elems.push_back(Eigen::Triplet<Complex>(idx, global(i + 2, j), NumericalDiffUtils<double>::C2_6_2() * ftx));
			elems.push_back(Eigen::Triplet<Complex>(idx, global(i + 3, j), NumericalDiffUtils<double>::C2_6_3() * ftx));

            elems.push_back(Eigen::Triplet<Complex>(idx, global(i, j - 3), NumericalDiffUtils<double>::C2_6_3() * fty));
            elems.push_back(Eigen::Triplet<Complex>(idx, global(i, j - 2), NumericalDiffUtils<double>::C2_6_2() * fty));
            elems.push_back(Eigen::Triplet<Complex>(idx, global(i, j - 1), NumericalDiffUtils<double>::C2_6_1() * fty));
            elems.push_back(Eigen::Triplet<Complex>(idx, global(i, j + 1), NumericalDiffUtils<double>::C2_6_1() * fty));
            elems.push_back(Eigen::Triplet<Complex>(idx, global(i, j + 2), NumericalDiffUtils<double>::C2_6_2() * fty));
            elems.push_back(Eigen::Triplet<Complex>(idx, global(i, j + 3), NumericalDiffUtils<double>::C2_6_3() * fty));

		} else {
			throw std::runtime_error("too higher order");
		}
	});
};


// M = dig + f * (V + T)
template<class R, class F, class Potential, class Dig>
void FillElems3D(std::vector<Eigen::Triplet<Complex> >& elems,
	size_t nx,
	size_t ny,
	size_t nz,
	int order,
	Dig dig,
	Potential const* v,
    F f,
    R hbar,
    R mass,
    R dx,
    R dy,
    R dz)
{

	// T = -hbar^2/2m nabla
	R const tx = -hbar * hbar / (2 * mass) / (dx * dx);
	R const ty = -hbar * hbar / (2 * mass) / (dy * dy);
	R const tz = -hbar * hbar / (2 * mass) / (dz * dz);
	auto const ftx = f * tx;
	auto const fty = f * tx;
	auto const ftz = f * tz;

	auto foldx = [=](size_t i) -> size_t {
		if ((ptrdiff_t)i < 0) i += nx;
		else if (i >= nx) i -= nx;
		return i;
	};

	auto foldy = [=](size_t i) -> size_t {
		if ((ptrdiff_t)i < 0) i += ny;
		else if (i >= ny) i -= ny;
		return i;
	};
	
	auto foldz = [=](size_t i) -> size_t {
		if ((ptrdiff_t)i < 0) i += nz;
		else if (i >= nz) i -= nz;
		return i;
	};

	auto global = [=](size_t x, size_t y, size_t z)->size_t {
        return foldz(z) + (foldy(y) + foldx(x) * ny) * nz;
	};

    // reserve memory to reduce memory waste
    if (order <= 2)
        elems.reserve(7 * nx * ny * nz);
    else if (order <= 4)
        elems.reserve(13 * nx * ny * nz);
    else if (order <= 6)
        elems.reserve(19 * nx * ny * nz);

    ForEach3DIdx(nx, ny, nz, [=, &elems](size_t idx, size_t i, size_t j, size_t k) {
        auto dig_ = dig + f * v[idx];

		if (order <= 2) {

            elems.push_back(Eigen::Triplet<Complex>(idx, global(i, j, k),
                NumericalDiffUtils<double>::C2_2_0() * ftx
                + NumericalDiffUtils<double>::C2_2_0() * fty
                + NumericalDiffUtils<double>::C2_2_0() * ftz
                + dig_));
            elems.push_back(Eigen::Triplet<Complex>(idx, global(i - 1, j, k), NumericalDiffUtils<double>::C2_2_1() * ftx));
            elems.push_back(Eigen::Triplet<Complex>(idx, global(i + 1, j, k), NumericalDiffUtils<double>::C2_2_1() * ftx));

            elems.push_back(Eigen::Triplet<Complex>(idx, global(i, j - 1, k), NumericalDiffUtils<double>::C2_2_1() * fty));
            elems.push_back(Eigen::Triplet<Complex>(idx, global(i, j + 1, k), NumericalDiffUtils<double>::C2_2_1() * fty));

            elems.push_back(Eigen::Triplet<Complex>(idx, global(i, j, k - 1), NumericalDiffUtils<double>::C2_2_1() * ftz));
            elems.push_back(Eigen::Triplet<Complex>(idx, global(i, j, k + 1), NumericalDiffUtils<double>::C2_2_1() * ftz));

		} else if (order <= 4) {

            elems.push_back(Eigen::Triplet<Complex>(idx, global(i, j, k),
                NumericalDiffUtils<double>::C2_4_0() * ftx
                + NumericalDiffUtils<double>::C2_4_0() * fty
                + NumericalDiffUtils<double>::C2_4_0() * ftz
                + dig_));

			elems.push_back(Eigen::Triplet<Complex>(idx, global(i - 2, j, k), NumericalDiffUtils<double>::C2_4_2() * ftx));
			elems.push_back(Eigen::Triplet<Complex>(idx, global(i - 1, j, k), NumericalDiffUtils<double>::C2_4_1() * ftx));
			elems.push_back(Eigen::Triplet<Complex>(idx, global(i + 1, j, k), NumericalDiffUtils<double>::C2_4_1() * ftx));
			elems.push_back(Eigen::Triplet<Complex>(idx, global(i + 2, j, k), NumericalDiffUtils<double>::C2_4_2() * ftx));

			elems.push_back(Eigen::Triplet<Complex>(idx, global(i, j - 2, k), NumericalDiffUtils<double>::C2_4_2() * fty));
			elems.push_back(Eigen::Triplet<Complex>(idx, global(i, j - 1, k), NumericalDiffUtils<double>::C2_4_1() * fty));
			elems.push_back(Eigen::Triplet<Complex>(idx, global(i, j + 1, k), NumericalDiffUtils<double>::C2_4_1() * fty));
			elems.push_back(Eigen::Triplet<Complex>(idx, global(i, j + 2, k), NumericalDiffUtils<double>::C2_4_2() * fty));

			elems.push_back(Eigen::Triplet<Complex>(idx, global(i, j, k - 2), NumericalDiffUtils<double>::C2_4_2() * ftz));
			elems.push_back(Eigen::Triplet<Complex>(idx, global(i, j, k - 1), NumericalDiffUtils<double>::C2_4_1() * ftz));
			elems.push_back(Eigen::Triplet<Complex>(idx, global(i, j, k + 1), NumericalDiffUtils<double>::C2_4_1() * ftz));
			elems.push_back(Eigen::Triplet<Complex>(idx, global(i, j, k + 2), NumericalDiffUtils<double>::C2_4_2() * ftz));
		} else if (order <= 6) {

            elems.push_back(Eigen::Triplet<Complex>(idx, global(i, j, k),
                NumericalDiffUtils<double>::C2_6_0() * ftx
                + NumericalDiffUtils<double>::C2_6_0() * fty
                + NumericalDiffUtils<double>::C2_6_0() * ftz
                + dig_));
			elems.push_back(Eigen::Triplet<Complex>(idx, global(i - 3, j, k), NumericalDiffUtils<double>::C2_6_3() * ftx));
			elems.push_back(Eigen::Triplet<Complex>(idx, global(i - 2, j, k), NumericalDiffUtils<double>::C2_6_2() * ftx));
			elems.push_back(Eigen::Triplet<Complex>(idx, global(i - 1, j, k), NumericalDiffUtils<double>::C2_6_1() * ftx));
			elems.push_back(Eigen::Triplet<Complex>(idx, global(i + 1, j, k), NumericalDiffUtils<double>::C2_6_1() * ftx));
			elems.push_back(Eigen::Triplet<Complex>(idx, global(i + 2, j, k), NumericalDiffUtils<double>::C2_6_2() * ftx));
			elems.push_back(Eigen::Triplet<Complex>(idx, global(i + 3, j, k), NumericalDiffUtils<double>::C2_6_3() * ftx));

			elems.push_back(Eigen::Triplet<Complex>(idx, global(i, j - 3, k), NumericalDiffUtils<double>::C2_6_3() * fty));
			elems.push_back(Eigen::Triplet<Complex>(idx, global(i, j - 2, k), NumericalDiffUtils<double>::C2_6_2() * fty));
			elems.push_back(Eigen::Triplet<Complex>(idx, global(i, j - 1, k), NumericalDiffUtils<double>::C2_6_1() * fty));
			elems.push_back(Eigen::Triplet<Complex>(idx, global(i, j + 1, k), NumericalDiffUtils<double>::C2_6_1() * fty));
			elems.push_back(Eigen::Triplet<Complex>(idx, global(i, j + 2, k), NumericalDiffUtils<double>::C2_6_2() * fty));
			elems.push_back(Eigen::Triplet<Complex>(idx, global(i, j + 3, k), NumericalDiffUtils<double>::C2_6_3() * fty));

            elems.push_back(Eigen::Triplet<Complex>(idx, global(i, j, k - 3), NumericalDiffUtils<double>::C2_6_3() * ftz));
            elems.push_back(Eigen::Triplet<Complex>(idx, global(i, j, k - 2), NumericalDiffUtils<double>::C2_6_2() * ftz));
            elems.push_back(Eigen::Triplet<Complex>(idx, global(i, j, k - 1), NumericalDiffUtils<double>::C2_6_1() * ftz));
            elems.push_back(Eigen::Triplet<Complex>(idx, global(i, j, k + 1), NumericalDiffUtils<double>::C2_6_1() * ftz));
            elems.push_back(Eigen::Triplet<Complex>(idx, global(i, j, k + 2), NumericalDiffUtils<double>::C2_6_2() * ftz));
            elems.push_back(Eigen::Triplet<Complex>(idx, global(i, j, k + 3), NumericalDiffUtils<double>::C2_6_3() * ftz));
		} else {
			throw std::runtime_error("too higher order");
		}
		});
};
