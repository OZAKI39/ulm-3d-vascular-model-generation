// Same MPI handles used by Stage J and CmMod.h, with real content checks added.
#include <mpi.h>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <array>
int main(int argc, char **argv) {
  MPI_Init(&argc,&argv);
  MPI_Comm_set_errhandler(MPI_COMM_WORLD,MPI_ERRORS_RETURN);
  MPI_Comm_set_errhandler(MPI_COMM_SELF,MPI_ERRORS_RETURN);
  int rank,ranks,failed=0,total=0;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank); MPI_Comm_size(MPI_COMM_WORLD,&ranks);
  if(argc!=5){MPI_Finalize();return 2;}
  MPI_Datatype types[]={MPI_INT,MPI_DOUBLE,MPI_CHAR,MPI_CXX_BOOL,
                       MPI_INTEGER,MPI_DOUBLE_PRECISION,MPI_CHARACTER,MPI_LOGICAL};
  const char *names[]={"MPI_INT","MPI_DOUBLE","MPI_CHAR","MPI_CXX_BOOL",
                      "MPI_INTEGER","MPI_DOUBLE_PRECISION","MPI_CHARACTER","MPI_LOGICAL"};
  int expected_sizes[]={sizeof(int),sizeof(double),sizeof(char),sizeof(bool),
    std::atoi(argv[1]),std::atoi(argv[2]),std::atoi(argv[3]),std::atoi(argv[4])};
  for(int i=0;i<8;++i){
    int bytes=-1,errclass=-1;
    int ts=MPI_Type_size(types[i],&bytes);
    std::array<unsigned char,64> reference{},buffer{};
    int integer=42,logical=1; double real=1.25; char character='Z'; bool boolean=true;
    const void *values[]={&integer,&real,&character,&boolean,&integer,&real,&character,&logical};
    int storage[]={sizeof(integer),sizeof(real),sizeof(character),sizeof(boolean),
                   sizeof(integer),sizeof(real),sizeof(character),sizeof(logical)};
    bool size_ok=ts==MPI_SUCCESS && bytes>0 && bytes==expected_sizes[i] && bytes==storage[i];
    int bc=MPI_ERR_TYPE;bool content_ok=false;
    if(size_ok){
      std::memcpy(reference.data(),values[i],bytes);
      if(rank==0)buffer=reference;
      bc=MPI_Bcast(buffer.data(),1,types[i],0,MPI_COMM_WORLD);
      content_ok=bc==MPI_SUCCESS && std::memcmp(reference.data(),buffer.data(),bytes)==0;
    }
    MPI_Error_class(bc,&errclass);
    std::printf("DATATYPE rank=%d ranks=%d name=%s size=%d expected_size=%d type_size_rc=%d bcast_rc=%d error_class=%d data_ok=%d\n",
                rank,ranks,names[i],bytes,expected_sizes[i],ts,bc,errclass,content_ok?1:0);
    failed+=!size_ok || !content_ok;
  }
  MPI_Allreduce(&failed,&total,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
  MPI_Finalize();return total?1:0;
}
