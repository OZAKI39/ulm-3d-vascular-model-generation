program smoke
implicit none
integer :: a
a=21
if (2*a /= 42) stop 1
print *, "FORTRAN result=42"
end program
