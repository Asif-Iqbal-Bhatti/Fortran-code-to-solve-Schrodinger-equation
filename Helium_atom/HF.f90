program HF_Helium
     	
    use matrixmultiplication ! using module for matrix multi. & calling DSYEV
    implicit none

  integer    				        :: p, n_alpha, iter
  integer   				        :: i, j, k, l, maxiter = 200
  double precision 			        :: start_time, stop_time
  double precision, parameter 		:: pi = 4.0*ATAN(1.0d0), Z = 2, x = 1
  double precision 					:: aa, eold, enew
  double precision, dimension(:), 	    allocatable	    ::  alpha, e, c
  double precision, dimension(:,:), 	allocatable 	::  h, f, s, v
  double precision, dimension(:,:,:,:), allocatable  	::  q
  character (len=100)                                   :: filename

!*****************Reading parameters from file********************************

  write (*,*) "Enter the Filename for Gaussian parameters ::   "
  read  (*,*) filename
  write (*,*) "Reading from a file ... "
  open  (unit=10, file=trim(filename),status='old',action='read')
  write (*,*) "Number of Gaussians detected ~~>> "
  read  (10,*) n_alpha
  allocate (alpha (n_alpha))
  do p = 1, n_alpha
     read (10,*) alpha(p)
     write (*,*) "For Gaussian # ",p,"Coefficient detected ~~>> ", alpha(p)
     if ( alpha(p) <= 0.0d0 ) stop '-ve coefficient non acceptable !!! '
  end do
!**********************************************************************************
  
!!!!Inserting the Coulomb integral with the matrix elements of the e-e interaction:
!!!q(i,j,k,l) = Integral[ Integral{ chi_i(r) chi_j(r') 1/|r'-r| chi_k(r) chi_l(r') dr dr'} ]
!!where chi_i(r) = Exponential(-alpha(i)*(r**2)) are the Gaussians  

  allocate ( q(n_alpha, n_alpha, n_alpha, n_alpha) )

  do i=1,n_alpha
     do j=1,n_alpha
        do k=1,n_alpha
           do l=1,n_alpha
              q(i,j,k,l) = x * (2.0d0*pi**2.5d0)/((alpha(i)+alpha(j))*(alpha(k)+alpha(l)) &                            * sqrt(alpha(i)+alpha(j)+alpha(k)+alpha(l)))
           end do
        end do
     end do
  end do

!***********Overlap integrals S and one-electron hamiltonian H on the GB*****
  
  allocate ( s(n_alpha, n_alpha), h(n_alpha, n_alpha), f(n_alpha, n_alpha), &
             v(n_alpha, n_alpha), c(n_alpha), e (n_alpha)  )
  
  do i=1,n_alpha
     do j=1,n_alpha
        aa = alpha(i)+alpha(j)
        s(i,j) = (pi/aa)**1.5d0
        h(i,j) = (s(i,j)*3.0d0*alpha(i)*alpha(j))/aa - (Z*2*pi)/aa 
     end do
  end do
!**********************writing matrix to a file*****************************

     open(unit = 555, file = 'fock.dat',status = 'replace', action = 'write')
     do i = 1, n_alpha
        write (555,'(F15.8)') (h(i,j), j = 1, n_alpha)
     enddo
     close(555)
  
!*************Starting Initial Guess for solving HF equations***************
  
  do i=1,n_alpha
     c(i) = 0.0d0
  end do
  c(1) = 1.1d0; c(2) = 1.50d0
!**************************Self-consistency iteration***********************

  enew = 0.0d0
  print*, "Coefficients initialised with a guess ..." 
  print*, "Entering Main SCF Loop ... "
  write (*,*) "SCF #     HF-Eigenvalue       Energy         Enew-Eold"
  call cpu_time(start_time) 
  
  do iter = 1, maxiter
!******************Fill the Fock matrix*************************************
     
  do i=1,n_alpha
     do j=1,n_alpha
        f(i,j) = h(i,j)
        do k=1,n_alpha
           do l=1,n_alpha
              f(i,j) = f(i,j) + q(i,j,k,l) * c(k)*c(l)
           end do
        end do
     end do
  end do
!*****************Solution [expansion coefficients are stored into v(j,i)
!*******************j=basis function index, i= eigenvalue index]
!********Reference:: Jos Thijssen "Comp. physic"*****************************     
     call diag( n_alpha, n_alpha, f, s, e, v )
     
     c(:) =  v(:,1)
     
     eold = enew
     enew = 0.0d0
     do i = 1, n_alpha
        do j = 1, n_alpha
           enew = enew + 2.0*h(i,j)*c(i)*c(j)
           do k = 1, n_alpha
              do l = 1, n_alpha
                 enew = enew + q(i, j, k, l)*c(i)*c(j)*c(k)*c(l)
              end do
           end do
        end do
     end do
     write(*, 100) iter, e(1), enew, enew-eold
100  format(2x,I4,2(6x,F10.6),6x,F15.12)     
     if ( abs (enew-eold) < 1.0d-8 ) then
     print*, ''//achar(27)//'[94m Convergence criterion satisfied exiting ... '//achar(27)//'[0m'
        call cpu_time(stop_time)
        print*, "Loop time", stop_time - start_time, "seconds" 
        deallocate ( e, c, v, f, h, s, alpha )
        stop 
     end if
     
  end do
  print*, ''//achar(27)//'[31m Convergence failed ... '//achar(27)//'[0m'
  deallocate ( e, c, v, f, h, s, alpha )
  stop 
  
end program HF_Helium
