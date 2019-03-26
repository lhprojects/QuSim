
#include "../qsim/Cal.h"
#include <chrono>
#include <functional>


void test(char const *expr, std::function<CCom(CCom x, CCom y)> const &func)
{
	printf("benchmark for %s\n", expr);
	int niter = 10000000;

	{
		Cal cal(expr);
		double c = 0;
		auto t0 = std::chrono::high_resolution_clock::now();
		cal.SetVarVal("x", 0);
		cal.SetVarVal("y", 0);
		auto &x = cal.GetVarVal("x");
		auto &y = cal.GetVarVal("y");
		cal.GenPseudoCode();
		for (int i = 0; i < niter; ++i) {
			x = i;
			y = i;
			c += cal.RunPseudoCode().real();
		}
		//printf("%f\n", c);
		auto t1 = std::chrono::high_resolution_clock::now();
		auto d = std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0);
		printf("Cal: %f MCLPS\n", 1.0 * niter / d.count());

	}

	{
		double c = 0;
		auto t0 = std::chrono::high_resolution_clock::now();
		for (int i = 0; i < niter; ++i) {
			c += func(i, i).real();
		}
		//printf("%f\n", c);
		auto t1 = std::chrono::high_resolution_clock::now();
		auto d = std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0);
		printf("Native Code: %f MCLPS\n", 1.0*niter / d.count());

	}

}
int main()
{
	printf("MCLPS = Million call per second\n");
	test("x*x+y*y", [](CCom x, CCom y) {return x * x + y * y; });
	test("sqrt(x*x+y*y)", [](CCom x, CCom y) {return sqrt(x * x + y * y); });
	test("exp(I*y)*exp(I*x)", [&](CCom x, CCom y) {return exp(CCom(0,1)*y)*exp(CCom(0, 1)*x); });

}

