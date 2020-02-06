program readMATRIX 
    implicit none
    integer   :: i, j, k, l, n, kd
    character(2) :: g
    double precision, dimension(:,:), allocatable  :: A, B, AB
    
    print*, LEN(g)
    kd = 2
    n = 0
    open (unit = 10, file = 'matrix', action = 'read', status = 'old')
      do
        read (10,*,end=10)
        n = n + 1
      enddo
10  close(10)
    print*, n
    allocate(B(n,n))
    allocate(A(n,n))
    allocate(AB(kd+1,n))
    AB=0d0 
    open(unit=20, file = 'matrix', status = 'old', action = 'read')
       read(20,*) ((A(i,j),j = 1, n), i = 1, n) 
    close(20)
    open(unit=30, file = 'saved.dat', action = 'write')
      do i = 1, n
         write(30,'(F5.3, 2F8.3)') (A(i,j), j = 1, n) 
      enddo
    close(30)

   
   open(unit=40, file = 'tri.dat', action = 'write',status = 'replace')
!      do i = 1, n
!         do j = i, min(n, i+kd)   
    open(unit = 30, file = 'saved.dat', action="read", recl = n)    
      do i = 1,n
         read (30, FMT=*, END=30) (AB(kd+1+i-j,j),j = i, min(n,i+kd))
      enddo 
30    close(30)
!            AB(kd+1+i-j,j) = A(i,j)
!!         enddo
!      enddo
        write(40,'(3F8.3)') AB
     close(40)
      do i = 1, 3
      write(*,'(3F8.3)') (AB(i,j),j=1,n)
      enddo
      
!!!Triming a matrix for zeroes!!!!!
       open(unit=40, file = 'tri.dat', action = 'read')
k = 0       
       do i = 1, kd +1
          do j= 1,n
             if (AB(i,j) == 0.000) then
             k=k+1
                   AB(i,j) = 0
             endif
          enddo
       enddo   
             
       close(40)
       print*,k
      do i = 1, 3
         write(*,'(3F8.3)') (AB(i,j),j=1,n)
      enddo
      
deallocate(A,B, AB)
end program
