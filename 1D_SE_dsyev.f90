program SE
 implicit none
 double precision                              :: t, start_time, stop_time
 double precision, parameter                   :: rmax = 50 
 integer                                       :: i, j, k, N, M, INFO, LWORK, l
 double precision, dimension(:), allocatable   :: r, W, WORK
 double precision, dimension(:,:), allocatable :: H
 intrinsic  INT, MIN, MAX

 print*, "Enter the size of a grid N > 2 " ; read*, N
 print*, "Enter the value for angular momemtum l = " ; read*, l
 
 LWORK = 3*N-1
 allocate(WORK(LWORK),r(N-1),H(N-1,N-1),W(N))

!****************************Intialization******************************* 
 t = rmax/dble(N)
 r = 0d0 
 H = 0d0 
!****************Constructing the corresponding Hamiltonian**************
!*********************for 1D schrodinger equation************************  
 do i = 1, N-2
    r(i)     = i*t
    H(i,i)   = 2*(1 - (1/r(i))*t*t + (dble(l)*(dble(l)+1))*(1/(r(i)**2))*((t*t)/2) )
    H(i+1,i) = -1d0
    H(i,i+1) = -1d0
 enddo
 H(N-1,N-1) = 2*(1-(1/(r(N-2)+t)*t*t) + (dble(l)*(dble(l)+1))*(1/(r(N-2)+t)**2) &
               *((t*t)/2) )
 H = H/(2*t**2) 
!************************************************************************

 open (unit = 10, file="matrix", action="write", status="replace")
 do i = 1, N-1
   write(10,'(1X, 10000F12.6)') (H(i,j), j= 1, N-1)
 enddo
 close(10)
 M = N - 1
!****************Eigenvalues & Eigen vectors***************************** 
 call cpu_time(start_time)
 call dsyev('V', 'L', M, H, M, W, WORK, LWORK, INFO)
 call cpu_time(stop_time)
 if( INFO.gt.0 ) then
   write(*,*)'The algorithm failed to compute eigenvalues.'
   stop
 endif 
!*************************************************************************  
 print*, ":: The Eigenvalues for l = ", l
 write(*,'(9X,100F8.4)') (W(i), i = 1, 5)
 write (*,*) "The corresponding Eigenvectors are"
 open(unit = 20, file = 'l= .dat', action = 'write', status = 'replace')
 write (20,'(9X,100F8.4)') (W(i), i = 1, 5)
 write(20,*) "For l = ", l
 write(20,'(2X,I4)') N-2  
 do i = 1, N-2
    write(20,'(1X, 1000F12.6)') r(i), H(i,1)*(1/r(i)), H(i,2)*(1/r(i)), H(i,3)*(1/r(i)), H(i,4)*(1/r(i))
 enddo
 close(20)
 print *, "Original loop time:", stop_time - start_time, "seconds"
end program SE 
