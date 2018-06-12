#include "Cal.h"
#include <assert.h>

void testCal()
{
	{
		Cal cal(" -1");
		assert(abs(cal.Val().real() - -1) < 0.001);
	}
	{
		Cal cal(" 1+2*3^4");
		assert(abs(cal.Val().real() - 163) < 0.001);
	}

	{
		Cal cal(" (3.1415926)");
		assert(abs(cal.Val().real() - 3.1415926) < 0.001);
	}

	{
		Cal cal(" sin(3.1415926)");
		assert(abs(cal.Val().real() - 0) < 0.001);
	}
	{
		Cal cal(" cos(3.1415926)");
		assert(abs(cal.Val().real() - -1) < 0.001);
	}

	{
		Cal cal(" sqrt(4)");
		assert(abs(cal.Val().real() - 2) < 0.001);

	}
	{
		Cal cal(" 1 + 1 ");
		assert(abs(cal.Val().real() - 2) < 0.001);
	}

}