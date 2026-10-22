/* Validation-only allocation admission; no library/system modifications. */
#define _GNU_SOURCE
#include <dlfcn.h>
#include <sys/file.h>
#include <fcntl.h>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
static int before(size_t n,const char* api){
 const char* path=getenv("BRAVA_CUDA_ADMISSION_LOCK");
 if(!path){fprintf(stderr,"CUDA_ADMISSION_UNCONFIGURED\n");_exit(86);}
 int fd=open(path,O_CREAT|O_RDWR,0600);if(fd<0||flock(fd,LOCK_EX)){perror("CUDA_ADMISSION_LOCK");_exit(86);}
 int(*info)(size_t*,size_t*)=dlsym(RTLD_NEXT,"cudaMemGetInfo");size_t freebytes=0,total=0;
 if(!info||info(&freebytes,&total)){fprintf(stderr,"CUDA_ADMISSION_QUERY_FAILED\n");_exit(86);}
 size_t reserve=(size_t)2048*1024*1024;
 int allowed=freebytes>=reserve && n<=freebytes-reserve;
 const char* log=getenv("BRAVA_CUDA_ADMISSION_LOG");
 if(log){FILE* f=fopen(log,"a");if(!f)_exit(86);fprintf(f,"{\"pid\":%ld,\"api\":\"%s\",\"request_bytes\":%zu,\"free_before_bytes\":%zu,\"reserve_bytes\":%zu,\"allowed\":%s}\n",(long)getpid(),api,n,freebytes,reserve,allowed?"true":"false");fclose(f);}
 if(!allowed){fprintf(stderr,"CUDA_ADMISSION_DENIED_BEFORE_ALLOCATION api=%s bytes=%zu free=%zu reserve=%zu\n",api,n,freebytes,reserve);flock(fd,LOCK_UN);close(fd);_exit(86);}
 return fd;
}
static void after(int fd){flock(fd,LOCK_UN);close(fd);}
int cudaMalloc(void** p,size_t n){int fd=before(n,"cudaMalloc");int(*fn)(void**,size_t)=dlsym(RTLD_NEXT,"cudaMalloc");if(!fn)_exit(86);int r=fn(p,n);after(fd);return r;}
int cudaMallocAsync(void** p,size_t n,void* stream){int fd=before(n,"cudaMallocAsync");int(*fn)(void**,size_t,void*)=dlsym(RTLD_NEXT,"cudaMallocAsync");if(!fn)_exit(86);int r=fn(p,n,stream);after(fd);return r;}
int cudaMallocManaged(void** p,size_t n,unsigned flags){int fd=before(n,"cudaMallocManaged");int(*fn)(void**,size_t,unsigned)=dlsym(RTLD_NEXT,"cudaMallocManaged");if(!fn)_exit(86);int r=fn(p,n,flags);after(fd);return r;}
