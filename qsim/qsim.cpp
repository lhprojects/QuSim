// qsim.cpp: 定义控制台应用程序的入口点。
//

#include "System.h"

#include <vector>
#include <complex>
#include <cmath>
#include <functional>
#include <gl/GL.h>
#include <gl/GLU.h>


int main()
{
	System sys;
	sys.init(0.01, 0.01, 1000, 1, 1);
	for (int i = 0; i < 100; ++i) {
		sys.step();
		printf("%f\n", sys.Norm());
	}

	return 0;
}

