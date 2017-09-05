MODULE typcnst
!
IMPLICIT NONE
PUBLIC
!
! Constants =========================================================
!
INTEGER, PARAMETER :: DP = selected_real_kind(14,200)
!INTEGER, PARAMETER :: DP = 8
!INTEGER, PARAMETER :: sgl = selected_real_kind(6,30)
!INTEGER, PARAMETER :: i4b = selected_int_kind(9)
INTEGER, PARAMETER :: i12 = selected_int_kind(12)
!
real(kind=dp), parameter    :: c_pi = 3.141592653589793_dp
complex(kind=dp), parameter :: c_im = cmplx(0.0_dp,1.0_dp)
!
integer, parameter       :: lsize = 30   ! max list size
integer, parameter       :: Lmax  = 3    ! max L quantum number for orbital and D matrix
integer, parameter       :: LTmax = 10   ! max lattice distance of TB hopping 
!
!
integer, parameter :: mxRint=2      ! max lattice distabce of interaction
!
! Type ==============================================================
!
! atom type
type :: atom_type
    integer :: elmt
    real(kind=dp), dimension(3) :: pos 
    real(kind=dp) :: alph, beta, gamm
    integer, dimension(0:Lmax) :: nobtL
    integer, dimension(1:2*Lmax+1,0:Lmax) :: M   !  m(i,L) gives the M number of i-th L orbitals
    integer, dimension(1:2*Lmax+1,0:Lmax) :: idx !ind(i,L) gives the orbital number (up) in Hr
    !                                             ind(i,L)+nobt gives the spin down one
    logical, dimension(1:2*Lmax+1,0:Lmax) :: cr  ! whether this orbital correlated
endtype
!
! atomic orbital type
type :: obt_type
    integer :: elmt                     ! Element
    integer :: at                       ! Atom
    real(kind=dp), dimension(3) :: pos  ! Position
    real(kind=dp) :: alph, beta, gamm   ! Euler's angles
    integer :: L                        ! Total angular momentum
    integer :: nthL                     ! This is the nth L orbital
    integer :: M                        ! Lz
    integer :: S                        ! Spin
    logical :: cr                       ! whether this orbital correlated
endtype
!
!
type :: sym_op
    character(len=30) :: str                ! The name of this operation
    real(kind=dp), dimension(3,3) :: R      ! Rotation matrix (in cart)
    integer :: det                          ! Determinant
    real(kind=dp), dimension(3) :: tau      ! Translational part (in cart)
    real(kind=dp) :: alph, beta, gamm       ! Euler's angles
    real(kind=dp), dimension(3) :: axis     ! Rotation axis
    real(kind=dp)               :: thet     ! Rotation angle
    integer                     :: rotm     ! Rotational order = 2pi/thet
endtype
!
! Irreps ------------------------------------------------------------
!
integer, parameter :: maxdim = 4    ! max dimension of irrep
integer, parameter :: maxbcp = 3    ! max basis components 
integer, parameter :: maxrep = 18   ! max number of irrep
!
type :: irepbs_type
    integer     :: J         ! Total angular momentum
    integer     :: M(maxbcp) ! Jz angular momentums  
    complex(dp) :: C(maxbcp) ! Coefficients
endtype
!
type :: irep_type
    character(len=30)   :: name         ! name of this irrep
    integer             :: D            ! dimension
    integer             :: P            ! parity
    type(irepbs_type)   :: bs(maxdim)   ! bases
    ! irrep matrices
    complex(kind=dp), dimension(maxdim,maxdim,48) :: mat  ! representation matrices 
    complex(kind=dp), dimension(48)       :: cha          ! character table
endtype

!
! pseudo-atom and wannier functions ---------------------------------
!
integer, parameter :: maxrep_pat=30
type :: pat_type
    integer              :: kd
    real(dp)             :: pos(3)
    character(len=30)    :: SSG
    real(dp)             :: alph, beta, gamm
    integer              :: nirep
    integer              :: idx1(maxrep_pat), idx2(maxrep_pat)
endtype
!
type :: wan_type
    integer     :: kd                   ! kind
    integer     :: at                   ! Atom
    integer     :: iirep                ! belongs to the iirep-th irep at atom at
    real(dp)    :: pos(3)               ! Position
    real(dp)    :: alph, beta, gamm     ! Euler's angles
    integer     :: J                    ! Total angular momentum
    integer     :: M(maxbcp)            ! Jz angular momentums  
    complex(dp) :: C(maxbcp)            ! Coefficients
endtype
!
ENDMODULE typcnst
