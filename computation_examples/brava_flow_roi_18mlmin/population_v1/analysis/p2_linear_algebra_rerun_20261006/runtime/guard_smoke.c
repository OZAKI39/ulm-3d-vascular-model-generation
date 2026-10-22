#include <stddef.h>
#include <stdlib.h>
extern int cudaMalloc(void**,size_t);extern int cudaFree(void*);int main(int c,char**v){void*p=0;size_t n=c>1?((size_t)1<<42):4096;int e=cudaMalloc(&p,n);if(e)return e;return cudaFree(p);}
