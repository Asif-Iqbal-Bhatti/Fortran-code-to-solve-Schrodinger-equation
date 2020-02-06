program SE
    implicit none
 double precision                                :: t
 double precision, parameter                     :: rmax = 40 
 integer                                         :: i, j, k, N, M, INFO, kd, LDAB, l
 double precision, dimension(:), allocatable     :: r, W, WORK
 double precision, dimension(:,:), allocatable   :: H, Z, ab
 character(1), parameter                         :: uplo = 'U'
 intrinsic  INT, MIN, MAX
 double precision                                :: st_time, sp_time 

 print*, "Enter the size of a grid N > 5: " ;  read*, N
 print*, "Enter the value for angular momentum l = " ; read*, l 
 
 kd = 1 !******* For H we have only one super/sub diagonal vector*******
 LDAB = kd + 1
 allocate(r(N-1),H(N-1,N-1),W(N),Z(N,N),ab(LDAB,N), WORK(3*N-2))
!***********************Intialization*********************************** 

 t = rmax/dble(N)
 r = 0d0 
 H = 0d0 
!****************Constructing the corresponding Hamiltonian**************
!*********************for 1D schrodinger equation************************  
 do i = 1, N-2
    r(i)     = i*t
    H(i,i)   = 2*(1 - (1/r(i))*t*t + (l*(l+1))*(1/r(i)**2)*((t*t)/2) )
    H(i+1,i) = -1d0
    H(i,i+1) = -1d0
 enddo
 H(N-1,N-1)  = 2*(1 - (1/(r(N-2)+t)*t*t) + (l*(l+1))*(1/(r(N-2)+t)**2)*((t*t)/2) )
 H = H/(2*t**2) 
!************************************************************************
 open(unit= 10, file = 'hamilmatrix', action = 'write', status = 'replace')
 do i = 1 , N-1
    write(10,'(1X,10000F12.8)') (H(i,j), j = 1, N-1)
 enddo
 close(10)
!********Upper or lower SYM BAND matrix is stored in ******************
 if (uplo == 'U') then
    do i = 1, N-1
       do j = i, MIN(N-1, i+kd)   
         ab(kd+1+i-j,j) = H(i,j)
       enddo
    enddo 
 else if (uplo == 'L') then
    do i = 1, N-1
       do j = i, MAX(1, i-kd)
         ab(1+i-j,j) = H(i,j)
       enddo
    enddo
 endif
!*******************************************************************************
N = N -1 !******************************** (N-1)x(N-1) matrix
!*******calling a LAPACK dsbev subroutine*******************************
 call cpu_time(st_time)
 call dsbev('V',uplo, N, kd, ab, LDAB, W, Z, N, WORK, INFO) 
 call cpu_time(sp_time)
!***********************printing the eigenvalues***********************

 print*, ":: The Eigen values for l =  ", l 
 WRITE(*,'(1x,1000f12.5)'), (W(i), i = 1, 5)
 print*, ":: The correspoding Eigenvector are"

!*********Saving the radial part in file to plot************************ 

 open(unit = 20, file = 'radialWF.dat', action = 'write', status = 'replace')
 write(20,*) N-1
 do i = 1, N-1
    write(20,'(1X, 1000F12.5)') r(i), -Z(i,1)*(1/r(i)) , -Z(i,2)*(1/r(i)) 
 enddo
 close(20)
 write(*,*) "Wall Time", sp_time - st_time, "seconds"

end program SE 
