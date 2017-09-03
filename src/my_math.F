MODULE my_math
!
USE typcnst, only:  dp
!  
CONTAINS
FUNCTION dot_prod(a,b)
   REAL(KIND=DP) :: dot_prod
   REAL(KIND=DP), DIMENSION(:), INTENT(IN) :: a,b
   dot_prod = SUM(a*b)
   !dot_prod=a(1)*b(1)+a(2)*b(2)+a(3)*b(3)
ENDFUNCTION dot_prod
!
FUNCTION x_prod(a,b)
   REAL(KIND=DP), DIMENSION(3) :: x_prod
   REAL(KIND=DP), DIMENSION(3), INTENT(IN) :: a,b
   x_prod(1)=a(2)*b(3)-a(3)*b(2) 
   x_prod(2)=a(3)*b(1)-a(1)*b(3) 
   x_prod(3)=a(1)*b(2)-a(2)*b(1) 
ENDFUNCTION x_prod
!
FUNCTION mix_prod(a,b,c)
   REAL(KIND=DP) :: mix_prod
   REAL(KIND=DP), DIMENSION(3), INTENT(IN) :: a,b,c
   mix_prod=dot_prod(a,x_prod(b,c))
ENDFUNCTION mix_prod
!
FUNCTION vect_len(a)
   REAL(KIND=DP) :: vect_len
   REAL(KIND=DP), DIMENSION(3), INTENT(IN) :: a
   vect_len = SQRT(SUM(a*a))
ENDFUNCTION
!
function cos_ab(a,b)
   real(kind=DP) :: cos_ab
   real(kind=DP), dimension(3), intent(in) :: a,b
   cos_ab = dot_prod(a,b)/vect_len(a)/vect_len(b)
endfunction
!
function square_twovec(a,b)
   real(kind=DP) :: square_twovec
   real(kind=DP), dimension(3), intent(in) :: a, b
   square_twovec = sqrt( ( vect_len(a)*vect_len(b) )**2 - dot_prod(a,b)**2 )
endfunction
!
function angle(x,y)
    real(kind=dp) :: angle, x, y, r
    r = sqrt(x*x + y*y)
    angle = acos(x/(r+1.0e-20_dp))
    if (y<0)  angle = -angle
endfunction
!
! function othmat(): whether D is an othogonal matrix
function othmat(D,n)
    integer :: n
    real(dp), dimension(n,n) :: D
    logical :: othmat
    !
    real(dp) :: A
    integer :: ii,jj,kk
    !
    othmat = .true.
    do ii=1,n
    do jj=ii,n
        A = 0.0_dp
        do kk=1,n
            A=A + D(kk,ii)*D(kk,jj)
        enddo
        !
        if ( ii==jj .and. abs(A-1.0_dp)>1.0e-3_dp ) othmat = .false.
        if ( ii/=jj .and. abs(A)>1.0e-3_dp )        othmat = .false.
        !
        if (.not. othmat) exit
    enddo
        if (.not. othmat) exit
    enddo
    !
endfunction
!
! function unitmat(): whether D is an unitary matrix
function unimat(D,n)
    integer :: n
    complex(dp),dimension(n,n) :: D
    logical :: unimat
    !
    complex(dp) :: A
    integer :: ii, jj, kk
    !
    unimat = .true.
    do ii=1,n
    do jj=ii,n
        A = 0.0_dp
        do kk=1,n
            A = A + conjg(D(kk,ii))*D(kk,jj)
        enddo
        !
        if( ii==jj .and. abs(A-1.0_dp)>1.0e-3_dp ) unimat = .false.
        if( ii/=jj .and. abs(A)>1.0e-3_dp )        unimat = .false.
        !
        if (.not. unimat) exit
    enddo
        if (.not. unimat) exit
    enddo
    !
endfunction

subroutine init_random_seed()
    integer :: i, n, clock
    integer, dimension(:), allocatable :: seed
    !
    call random_seed(size = n)
    allocate(seed(n))
    call system_clock(count=clock)
    seed = clock + 37 * (/ (i - 1, i = 1, n) /)
    call random_seed(put = seed)
    deallocate(seed)
end subroutine init_random_seed

ENDMODULE my_math
