program mpi_fortran_smoke
  use mpi
  implicit none
  integer :: ierr, rank, ranks, value, failures, total, sizes(4)
  double precision :: real_value
  character :: character_value
  logical :: logical_value
  call MPI_Init(ierr)
  call MPI_Comm_set_errhandler(MPI_COMM_WORLD, MPI_ERRORS_RETURN, ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, ranks, ierr)
  sizes = [storage_size(value)/8, storage_size(real_value)/8, &
           storage_size(character_value)/8, storage_size(logical_value)/8]
  value = -1
  real_value = -1.0d0
  character_value = 'x'
  logical_value = .false.
  if (rank == 0) then
    value = 42
    real_value = 1.25d0
    character_value = 'Z'
    logical_value = .true.
  end if
  failures = 0
  call MPI_Bcast(value, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  if (ierr /= MPI_SUCCESS .or. value /= 42) failures = failures+1
  call MPI_Bcast(real_value, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  if (ierr /= MPI_SUCCESS .or. real_value /= 1.25d0) failures = failures+1
  call MPI_Bcast(character_value, 1, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
  if (ierr /= MPI_SUCCESS .or. character_value /= 'Z') failures = failures+1
  call MPI_Bcast(logical_value, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
  if (ierr /= MPI_SUCCESS .or. .not. logical_value) failures = failures+1
  write(*,'(A,I0,A,I0,A,4(I0,1X),A,I0)') 'FORTRAN rank=',rank, &
    ' ranks=',ranks,' sizes=',sizes,'failures=',failures
  call MPI_Allreduce(failures,total,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
  call MPI_Finalize(ierr)
  if (total /= 0) stop 1
end program
