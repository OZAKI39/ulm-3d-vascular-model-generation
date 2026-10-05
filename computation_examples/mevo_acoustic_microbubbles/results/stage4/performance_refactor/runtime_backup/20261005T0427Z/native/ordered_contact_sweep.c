
#include <math.h>
#include <stdint.h>
typedef double (*ddot64)(int64_t, const double *, int64_t, const double *, int64_t);
int ordered_contact_sweep(double *x, const double *r,
                          const int64_t *pairs, int64_t count, ddot64 dot) {
    for (int64_t k=0; k<count; ++k) {
        int64_t i=pairs[2*k], j=pairs[2*k+1];
        double delta[3];
        for (int a=0; a<3; ++a) delta[a]=x[3*j+a]-x[3*i+a];
        double squared = dot ? dot(3,delta,1,delta,1) :
            (delta[0]*delta[0]+delta[1]*delta[1])+delta[2]*delta[2];
        double d=sqrt(squared);
        if (d==0.0) return 1;
        double amount=fmax(0.0,(r[i]+r[j])-d);
        double fraction=r[j]/(r[i]+r[j]);
        for (int a=0; a<3; ++a) x[3*i+a]-=((fraction*amount)*delta[a])/d;
        for (int a=0; a<3; ++a) x[3*j+a]+=(((1.0-fraction)*amount)*delta[a])/d;
    }
    return 0;
}
