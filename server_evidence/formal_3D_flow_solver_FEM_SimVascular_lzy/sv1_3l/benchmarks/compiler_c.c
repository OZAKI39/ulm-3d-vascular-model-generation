#include <stdio.h>
#if __GNUC__ != 12
#error expected GCC12
#endif
int main(void){int a=21;printf("C result=%d\n",2*a);return 2*a!=42;}
