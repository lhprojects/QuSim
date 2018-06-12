#pragma once

#include <complex>
#include <cmath>
#include <math.h>

inline std::complex<double> erfz(std::complex<double> const &c)
{
	double const Pi = 3.141592653589793238462643383279502884197169399375;
	double const Re = c.real();
	double const R2 = Re * Re;
	double Im = c.imag();

	if (Im == 0) return std::complex<double>(erf(Re), 0);

	bool conj = Im < 0;
	if (conj) Im = -Im;

	double erfR = erf(Re);

	double E_R = exp(-R2) / (2 * Pi * Re) * (1 - cos(-2 * Re * Im));
	double E_I = exp(-R2) / (2 * Pi * Re) * (0 - sin(-2 * Re * Im));

	int const N = 15;
	double F_R = 0;
	for (int n = N; n >= 1; --n) {
		F_R += exp(-0.25*n*n) / (0.25*n*n + R2);
	}
	F_R *= exp(-R2) * Re / Pi;

	int M = (int)(2 * Im);
	double H_R = 0;
	double H_I = 0;
	double G_R = 0;
	double G_I = 0;
	int M_N = M - N > 1 ? M - N : 1;
	for (int n = M + N; n >= M_N; --n) {
		H_R += exp(-0.25*n*n - n * Im - R2) / (0.25*n*n + R2) * Re;
		H_I += exp(-0.25*n*n - n * Im - R2) / (0.25*n*n + R2) * (0.5 * n);

		G_R += exp(-0.25*(n - 2 * Im)*(n - 2 * Im) + (Im*Im - R2)) / (0.25*n*n + R2) * Re;
		G_I += exp(-0.25*(n - 2 * Im)*(n - 2 * Im) + (Im*Im - R2)) / (0.25*n*n + R2) * (-0.5 * n);
	}
	H_R *= 1 / (2 * Pi);
	H_I *= 1 / (2 * Pi);
	G_R *= 1 / (2 * Pi);
	G_I *= 1 / (2 * Pi);

	std::complex<double> a(cos(-2 * Re * Im), sin(-2 * Re * Im));
	std::complex<double> b = a * std::complex<double>(G_R + H_R, G_I + H_I);
	double real = erfR + E_R + F_R - b.real();
	double imag = 0 + E_I + 0 - b.imag();

	if (conj) return std::complex<double>(real, -imag);
	return std::complex<double>(real, imag);
}
