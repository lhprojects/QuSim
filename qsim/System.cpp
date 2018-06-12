#include "System.h"

using std::abs;

int g = 0;
Complex simpsons_rule_rec(std::function<Complex(Real)> const &f, Real a, Real b,
	Complex fa, Complex fb, Complex fc,
	Real eps, Complex whole, Int n)
{
	g++;
	Real c = 0.5*(a + b);
	Real c1 = 0.5*(a + c);
	Real c2 = 0.5*(c + b);
	Real h = 1 / 12.0 * (b - a);
	Complex fc1 = f(c1);
	Complex fc2 = f(c2);
	Complex l = h * (fa + 4.0*fc1 + fc);
	Complex r = h * (fc + 4.0*fc2 + fb);
	//printf("n: %d\n", (int)n);
	//printf("left: %.20f\n", l.real());
	//printf("right: %.20f\n", r.real());
	//printf("whole: %.20f\n", whole.real());
	//printf("err: %.20f\n", abs(l + r - whole)/whole);
	if (n < 0 && abs(l + r - whole) < 15 * eps* abs(whole)) {
		return l + r + 1 / 15.0*(l + r - whole);
	}
	Complex nl = simpsons_rule_rec(f, a, c, fa, fc, fc1, eps, l, n - 1);
	Complex nr = simpsons_rule_rec(f, c, b, fc, fb, fc2, eps, r, n - 1);
	return nl + nr;
}


Complex Integrate(std::function<Complex(Real)> const &f, Real a, Real b)
{
	Complex fa = f(a);
	Complex fb = f(b);
	Complex fc = f(0.5*(a + b));
	Real h = 1 / 6.0 * (b - a);
	return simpsons_rule_rec(f, a, b, fa, fb, fc, 1E-14, h*(fa + 4.0*fc + fb), 8);
}


