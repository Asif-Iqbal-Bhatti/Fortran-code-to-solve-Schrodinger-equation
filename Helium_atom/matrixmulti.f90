module matrixmultiplication
	implicit none
	contains
  subroutine diag ( n, ldh, h, s, e, v )
  
  integer, intent(in)           :: n, ldh
  double precision, intent(in)  :: h(ldh,n), s(ldh,n)
  double precision, intent(out) :: e(n), v(ldh,n)
  integer lwork, info, i, j, nn
  double precision, parameter :: small=1.d-10
  double precision, allocatable :: work(:), b(:,:), h1(:,:)
  info = 0
  lwork = 3*n
  allocate (work(lwork), b(ldh,n))
!*************Copy S into an auxiliary matrix because dsyev destroys the matrix
  b = s
  
!****************Diagonalize S matrix*******************************************

  call dsyev ( 'V', 'U', n, b, ldh, e, work, lwork, info )
  if (info /= 0) stop 'S-matrix diagonalization failed '
!*******************************************************************************

!***Keep only linearly independent combinations (within a given threshold)
!***store into matrix "aux" the eigenvectors of S divided by the squares
!***of the eigenvalues
  nn = 0
  do i=1,n
     if (e(i) > small) then
        nn = nn + 1
        b(:,nn) = b(:,i) / sqrt(e(i))
     end if
!********print*, aux
  end do
  if ( nn < n) write(*,*) " # of linearly independent vectors =", nn, n
!*******************************************************************************
!***************Transform H using the "B" matrix******************************

!**********V(i,j) = \sum_{k=1}^{n}{H(i,k)x B(k,j)},  i=1,n, j=1,nn

  call dgemm ( 'N', 'N', n, nn, n, 1.0d0, h, ldh, b, ldh, 0.0d0, v, ldh )

!***********h1(i,j) = \sum_{k=1}^{n} {B(k,i) x v(k,j)},  i=1,nn, j=1,nn
!*********************************H' = b^T*H*b******************************

  allocate (h1(nn,nn) )
  call dgemm ( 'T', 'N', nn, nn, n, 1.0d0, b, ldh, v, ldh, 0.0d0, h1, nn )
!****************************Diagonalize H again********************************

  info = 0
  call dsyev ('V', 'U', nn, h1, nn, e, work, lwork, info)
  if (info /= 0) stop 'H-matrix diagonalization failed '
!*********************Back-transform eigenvectors*******************************

  call dgemm ( 'N', 'N', n, nn, nn, 1.0d0, b, ldh, h1, nn, 0.0d0, v, ldh )

  deallocate (h1, b, work)
  end subroutine diag

end module matrixmultiplication
