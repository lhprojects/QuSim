#define _CRT_SECURE_NO_WARNINGS

#include "../qsim/System.h"
int main()
{

	System p;
	p.init("abs(x) < 0.0000001 ? 1 : 0",
		false, 0.01, "0", -50, 50, 1001,
		BoundaryCondition::Period, 1, 1);

	System l;
	l.init("abs(x) < 0.0000001 ? 1 : 0",
		false, 0.01, "0", -50, 50, 1001,
		BoundaryCondition::ExtendInfinity, 1, 1);

	for (int i = 0; i < 1; ++i) {
		p.step();
		l.step();
	}
	FILE *f = fopen("a.txt", "w");
	for (size_t i = 0; i < l.fN; ++i) {
		fprintf(f, "(% f, % f, % f) (% f, % f, % f)\n", l.fPsi[i].real(), l.fPsi[i].imag(), abs(l.fPsi[i]),
			p.fPsi[i].real(), p.fPsi[i].imag(), abs(p.fPsi[i]));
	}
	fclose(f);

	return 0;
}

