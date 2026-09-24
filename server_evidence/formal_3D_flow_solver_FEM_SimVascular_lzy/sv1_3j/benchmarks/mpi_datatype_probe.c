#include <mpi.h>
#include <stdio.h>
#include <string.h>
int main(int argc,char **argv) {
  MPI_Init(&argc,&argv);
  MPI_Comm_set_errhandler(MPI_COMM_WORLD,MPI_ERRORS_RETURN);
  MPI_Comm_set_errhandler(MPI_COMM_SELF,MPI_ERRORS_RETURN);
  int rank,size,failed=0;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);MPI_Comm_size(MPI_COMM_WORLD,&size);
  MPI_Datatype types[]={MPI_INT,MPI_DOUBLE,MPI_CHAR,MPI_CXX_BOOL,MPI_INTEGER,MPI_DOUBLE_PRECISION,MPI_CHARACTER,MPI_LOGICAL};
  const char *names[]={"MPI_INT","MPI_DOUBLE","MPI_CHAR","MPI_CXX_BOOL","MPI_INTEGER","MPI_DOUBLE_PRECISION","MPI_CHARACTER","MPI_LOGICAL"};
  for(int i=0;i<8;i++) {
    int bytes=-1,errclass=-1;union { double alignment;unsigned char bytes[64]; } buffer;
    memset(&buffer,0,sizeof(buffer));
    int ts=MPI_Type_size(types[i],&bytes);
    int bc=MPI_Bcast(&buffer,1,types[i],0,MPI_COMM_WORLD);
    MPI_Error_class(bc,&errclass);
    printf("DATATYPE rank=%d ranks=%d name=%s size=%d type_size_rc=%d bcast_rc=%d error_class=%d\n",rank,size,names[i],bytes,ts,bc,errclass);
    failed+=(ts!=MPI_SUCCESS || bc!=MPI_SUCCESS);
  }
  MPI_Finalize();return failed?1:0;
}
