#include <cstdio>
#include "qd/dd_real.h"
#include "str_dd.h"
#include "top.h"

int main(int argc, char *argv[])
{
    if (argc != 3) throw gen_err("usage: ./cmp_dd_double fname1 fname2");
    ifstream in1(argv[1]);
    ifstream in2(argv[2]);
    dd_real x1, x2;
    
    while (in1.good() && in2.good()) {
	in1 >> x1;
	in2 >> x2;
	printf("%s\n", str(x1-x2));
    }
}
