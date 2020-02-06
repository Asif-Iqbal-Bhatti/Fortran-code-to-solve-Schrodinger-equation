program eigendstebz
    USE SE_mod
    implicit none
 
 character(1), parameter                   :: RANGE = 'A', ORDER = 'E'
 integer                                   :: N, M, NSPLIT, INFO, i, j, l
 integer, parameter                        :: IL = 0, IU = 0
 integer, dimension(:), allocatable        :: IBLOCK, ISPLIT,  IWORK
 double precision, parameter               :: VL = 0, VU = 0 
 double precision, parameter               :: ABSTOL = -1
 double precision, dimension (:), allocatable    :: W, WORK, D, E
 double precision, dimension (:,:), allocatable  :: A
 double precision                          :: start_time, stop_time     

!********Requesting matrix dimension*************************************** 

 print*, "Enter the matrix dimension > 2:  "; read*, N
 print*, "Enter the corresponding l = "; read*, l
 
 call hamil(A,N,l)
 allocate (IBLOCK(N), ISPLIT(N), IWORK(3*N), W(N), WORK(4*N), D(N), E(N-1))
 N = N - 1
 open (unit = 10, file = 'hamilmatrix', action = 'write', status = 'replace')
  do i = 1, N
   write(10,'(1x,1000000F12.6)') (A(i,j), j= 1, N)
  enddo
 close(10)
!*****assigning diagonal and sub/super diagonal terms******************
 do i = 1, N
   D(i) = A(i,i)
 enddo
 do j = 1, N-1
   E(j) = A(j,j+1)
 enddo  
 
!*******************All Eigenvalues is requested***********************

 call cpu_time(start_time)
 call DSTEBZ(RANGE, ORDER, N, VL, VU, IL, IU, ABSTOL, D, E, M, &
             NSPLIT, W, IBLOCK, ISPLIT, WORK, IWORK, INFO)
 call cpu_time(stop_time) 
 if (INFO.NE.0) then
   write(*,*) 'DSTEBZ', INFO
 else
   write(*,*)
   write(*,*) "Number of Eigen values found", M, NSPLIT
   Write (*, *) 'Eigenvalues are::'
   Write (*, *) (W(j), j= 1, 5)
 endif    
 print*, "Original loop time", start_time - stop_time, "seconds"

deallocate(IBLOCK,ISPLIT,IWORK,W,WORK,D,E,A)
end program eigendstebz
