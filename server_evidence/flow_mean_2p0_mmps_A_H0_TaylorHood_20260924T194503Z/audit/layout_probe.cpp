#include "ComMod.h"
#include <cstddef>
#include <cstdio>
int main(){printf("{\"msh_eNoN\":%zu,\"msh_nFs\":%zu,\"msh_eType\":%zu,\"msh_fs\":%zu,\"fs_size\":%zu,\"fs_eNoN\":%zu,\"fs_eType\":%zu}\n",offsetof(mshType,eNoN),offsetof(mshType,nFs),offsetof(mshType,eType),offsetof(mshType,fs),sizeof(fsType),offsetof(fsType,eNoN),offsetof(fsType,eType));}
