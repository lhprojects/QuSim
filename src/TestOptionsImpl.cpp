#include "OptionsImpl.h"
#define ASS(x) if(!(x)) { printf("%s failed\n", #x);}

void testOptionsImpl()
{
    OptionsImpl opt;

    opt.SetBool("a", true);
    opt.SetInt("a", 1);

    bool b = false;
    ASS(!opt.Get("a", b));
    Int i = 0;
    ASS(opt.Get("a", i));

}
