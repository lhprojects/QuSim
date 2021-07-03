#include "Cal.h"
#include <assert.h>

#define ASS(x) if(!(x)) { printf("%s failed\n", #x);}

double const Pi = 3.1415926;
void testCal()
{
	{
		Cal cal(" 1");
		ASS(abs(cal.Val().real() - 1) < 0.001);
	}
	{
		Cal cal(" -1");
		ASS(abs(cal.Val().real() - -1) < 0.001);
	}

	{
		Cal cal(" (3.1415926)");
		ASS(abs(cal.Val().real() - 3.1415926) < 0.001);
	}

	{
		Cal cal(" 1 ? 2 : 3");
		ASS(abs(cal.Val().real() - 2) < 0.001);
	}
	{
		Cal cal(" 0 ? 2 : 3");
		ASS(abs(cal.Val().real() - 3) < 0.001);
	}
	{
		Cal cal(" 1 && 1");
		ASS(abs(cal.Val().real() - 1) < 0.001);
	}
	{
		Cal cal(" 1 && 0");
		assert(abs(cal.Val().real() - 0) < 0.001);
	}
	{
		Cal cal(" 1 || 1");
		ASS(abs(cal.Val().real() - 1) < 0.001);
	}
	{
		Cal cal(" 1 || 0");
		ASS(abs(cal.Val().real() - 1) < 0.001);
	}
	{
		Cal cal(" 0 || 0");
		ASS(abs(cal.Val().real() - 0) < 0.001);
	}

	{
		Cal cal(" 0 || 1 ? 2 : 3");
		ASS(abs(cal.Val().real() - 2) < 0.001);
	}
	{
		Cal cal(" 0 && 1 ? 2 : 3");
		ASS(abs(cal.Val().real() - 3) < 0.001);
	}
	{
		Cal cal(" 1+2*3^4");
		ASS(abs(cal.Val().real() - 163) < 0.001);
	}

	{
		Cal cal(" sin(3.1415926/2)");
		ASS(abs(cal.Val().real() - 1) < 0.001);
	}

	{
		Cal cal(" cos(3.1415926)");
		ASS(abs(cal.Val().real() - -1) < 0.001);
	}
	{
		Cal cal(" tan(3.1415926/4)");
		ASS(abs(cal.Val().real() - 1) < 0.001);
	}

	{
		Cal cal(" asin(1)");
		ASS(abs(cal.Val().real() - Pi / 2) < 0.001);
	}

	{
		Cal cal(" acos(1)");
		ASS(abs(cal.Val().real() - 0) < 0.001);
	}

	{
		Cal cal(" atan(1)");
		ASS(abs(cal.Val().real() - 3.1415926 / 4) < 0.001);
	}

	{
		Cal cal(" exp(1)");
		ASS(abs(cal.Val().real() - 2.718281) < 0.001);
	}

	{
		Cal cal(" log(1)");
		ASS(abs(cal.Val().real() - 0.0) < 0.001);
	}

	{
		Cal cal(" sqrt(4)");
		ASS(abs(cal.Val().real() - 2) < 0.001);
	}
	{
		Cal cal(" 1 + 1 ");
		ASS(abs(cal.Val().real() - 2) < 0.001);
	}

}
