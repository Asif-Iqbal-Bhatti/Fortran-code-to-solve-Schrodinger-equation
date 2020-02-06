MODULE SE_mod
    implicit none
 contains
 subroutine hamil(H,N,l) 

 double precision                                :: t
 double precision, parameter                     :: rmax = 100 
 integer                                         :: i, j, k, M, g
 double precision, dimension(:), allocatable     :: r
 double precision, dimension(:,:), allocatable, intent(inout)   :: H
 integer, intent(in)             :: N, l

 allocate(r(N-1),H(N-1,N-1))
!****************Intialization 

 t = rmax/dble(N)
 r = 0d0 
 H = 0d0 
!****************Constructing the corresponding Hamiltonian**************
!*********************for 1D schrodinger equation************************  
!************Theoretical Energy level -1/2n**2***************************
 do i = 1, N-2
    r(i)     = i * t
    H(i,i)   = 2*(1 - (1/r(i))*t*t + (l*(l+1))*(1/r(i)**2)*((t*t)/2) )
    H(i+1,i) = -1d0
    H(i,i+1) = -1d0
 enddo
 H(N-1,N-1)   = 2*(1 - (1/(r(N-2)+t)*t*t) + (l*(l+1))*((1/r(N-2)+t)**2)*((t*t)/2) )
 H = H/(2*t**2) 
!************************************************************************

 end subroutine hamil
 
END MODULE SE_mod 
