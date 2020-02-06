program SE
 implicit none
 double precision                                :: t
 double precision, parameter                     :: rmax = 10 
 integer                                         :: i, j, k, N, M, INFO, LWORK
 integer, parameter                              :: l = 0, LWMAX = 10000
 double precision, dimension(:), allocatable     :: r, W
 double precision, dimension(:,:), allocatable   :: H
 double precision, dimension(LWMAX)              :: WORK
 intrinsic  INT, MIN
 print*, "Enter the size of a grid N > 5: "
 read*, N 
 LWORK = -1
 allocate(r(N))
 allocate(H(N,N))  
 allocate(W(N))

!****************Intialization 

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
 H(N-1,N-1)   = 2*(1 - (1/(r(N-2)+t)*t*t) + (l*(l+1))*(1/r(i)**2)*((t*t)/2) )
 H = H/(2*t**2) 
!************************************************************************
! do i = 1, N
!    write(*,'(8F8.3)') (H(i,j), j= 1, N)
! enddo
 M = N
! Query for the optimal workspace
 call dsyev('V', 'L', M, H, M, W, WORK, LWORK, INFO)
 LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )
 call dsyev('V', 'L', M, H, M, W, WORK, LWORK, INFO)
 if( INFO.gt.0 ) then
    write(*,*)'The algorithm failed to compute eigenvalues.'
    stop
 endif 
  
 print*, ":: The Eigen values for l = 0 are ..."
 do i = 1, 10
 print*, W(i)
 enddo
 do i = 1, N
    write(*,'(8F8.3, 8F8.3)') i*t, H(i,1)
 enddo  
end program SE 
