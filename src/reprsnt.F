module reprsnt
!
use typcnst,    only: dp, c_pi, c_im, maxdim, maxbcp, maxrep, irep_type, sym_op, &
                      irepbs_type
use my_math,    only: vect_len, unimat, x_prod, othmat
use pntgrp,     only: euler_angle
!
implicit none
public
private :: fac, bnum, bst
!
real(dp)            :: fac(40)
integer,parameter   :: bnum=50
type(irepbs_type)   :: bst(bnum)
!
contains
!
!**** init_pntgrp_bst() *********************************************
!
! Define all the bases needed to construct irrep of point groups
!   For now, only single group irreps are considered
!
!********************************************************************
subroutine init_pntgrp_bst()
    real(kind=dp),parameter :: sq12=sqrt(0.5_dp), &
                               sq16=sqrt(1.0_dp/6.0_dp), &
                               sq56=sqrt(5.0_dp/6.0_dp)
    integer :: ii
    !
    do ii=1,bnum
        bst(ii)%J=0
        bst(ii)%M(:)=0
        bst(ii)%C(:)=0.0_dp
    enddo
    !
    ! NSOC bases ----------------------------------------------------
    !
    ! |00>
    bst(1)%J=0
    bst(1)%M(1)=0
    bst(1)%C(1)=1.0_dp
    !
    ! sqrt(1/2)|32> - sqrt(1/2)|3-2>
    bst(2)%J=3
    bst(2)%M(1:2)=(/ 2,-2/) 
    bst(2)%C(1:2)=(/ sq12,-sq12/)
    !
    ! sqrt(1/2)|20> - 0.5i*|22> - 0.5i*|2-2>
    bst(3)%J=2
    bst(3)%M(1:3)=(/ 0, 2, -2/)
    bst(3)%C(1)= sq12; bst(3)%C(2)=-0.5*c_im; bst(3)%C(3)=-0.5*c_im 
    !
    ! sqrt(1/2)|20> + 0.5i*|22> + 0.5i*|2-2>
    bst(4)%J=2
    bst(4)%M(1:3)=(/ 0, 2, -2/)
    bst(4)%C(1)= sq12; bst(4)%C(2)= 0.5*c_im; bst(4)%C(3)= 0.5*c_im 
    !
    ! sqrt(1/2)|11> + sqrt(1/2)|1-1>
    bst(5)%J=1
    bst(5)%M(1:2)=(/ 1,-1/)
    bst(5)%C(1:2)=(/ sq12, sq12/)
    !
    ! |10>
    bst(6)%J=1
    bst(6)%M(1)=0
    bst(6)%C(1)=1.0_dp
    !
    ! sqrt(1/2)|11> - sqrt(1/2)|1-1>
    bst(7)%J=1
    bst(7)%M(1:2)=(/ 1,-1/)
    bst(7)%C(1:2)=(/ sq12,-sq12/)
    !
    ! sqrt(1/2)|21> - sqrt(1/2)|2-1>
    bst(8)%J=2
    bst(8)%M(1:2)=(/ 1,-1/)
    bst(8)%C(1:2)=(/ sq12,-sq12/)
    !
    ! - sqrt(1/2)|22> + sqrt(1/2)|2-2>
    bst(9)%J=2
    bst(9)%M(1:2)=(/ 2,-2/)
    bst(9)%C(1:2)=(/-sq12, sq12/)
    !
    ! - sqrt(1/2)|21> - sqrt(1/2)|2-1>
    bst(10)%J=2
    bst(10)%M(1:2)=(/ 1,-1/)
    bst(10)%C(1:2)=(/-sq12,-sq12/)
    !
    ! sqrt(1/2)|33> - sqrt(1/2)|3-3>
    bst(11)%J=3
    bst(11)%M(1:2)=(/3,-3/)
    bst(11)%C(1:2)=(/sq12,-sq12/)
    !
    ! sqrt(1/2)|33> + sqrt(1/2)|3-3>
    bst(12)%J=3
    bst(12)%M(1:2)=(/3,-3/)
    bst(12)%C(1:2)=(/sq12,sq12/)
    !
    ! |11>
    bst(13)%J=1
    bst(13)%M(1)=1
    bst(13)%C(1)=1.0_dp
    !
    ! |1-1>
    bst(14)%J=1
    bst(14)%M(1)=-1
    bst(14)%C(1)=1.0_dp
    !
    ! |2-2>
    bst(15)%J=2
    bst(15)%M(1)=-2
    bst(15)%C(1)=1.0_dp
    !
    ! -|22>
    bst(16)%J=2
    bst(16)%M(1)=2
    bst(16)%C(1)=-1.0_dp
    !
    ! sqrt(1/2)|22> + sqrt(1/2)|2-2>
    bst(17)%J=2
    bst(17)%M(1:2)=(/2,-2/)
    bst(17)%C(1:2)=(/sq12,sq12/)
    !
    ! sqrt(1/2)|22> - sqrt(1/2)|2-2>
    bst(18)%J=2
    bst(18)%M(1:2)=(/2,-2/)
    bst(18)%C(1:2)=(/sq12,-sq12/)
    !
    ! SOC bases -----------------------------------------------------
    !
    ! O: 1E12 = |1/2 1/2>
    bst(30)%J=-1
    bst(30)%M(1)=1
    bst(30)%C(1)=1.0_dp
    !
    ! O: 2E12 = |1/2 -1/2>
    bst(31)%J=-1
    bst(31)%M(1)=-1
    bst(31)%C(1)=1.0_dp
    !
    ! O: 1E52 = sqrt(1/6)|5/2 5/2> - sqrt(5/6)|5/2 -3/2>
    bst(32)%J=-5
    bst(32)%M(1:2)=(/5, -3/)
    bst(32)%C(1:2)=(/sq16, -sq56/)
    !
    ! O: 2E52 = sqrt(1/6)|5/2 3/2> - sqrt(5/6)|5/2 -5/2>
    bst(33)%J=-5
    bst(33)%M(1:2)=(/3, -5/)
    bst(33)%C(1:2)=(/-sq56, sq16/)
    !
    ! O: 1F32 = sqrt(1/2)|3/2 1/2> - i*sqrt(1/2)|3/2 -3/2>
    bst(34)%J=-3
    bst(34)%M(1:2)=(/1, -3/)
    bst(34)%C(1:2)=(/cmplx(sq12,0.0_dp), -cmplx(0.0_dp,sq12)/)
    !
    ! O: 2F32 = i*sqrt(1/2)|3/2 3/2> - sqrt(1/2)|3/2 -1/2>
    bst(35)%J=-3
    bst(35)%M(1:2)=(/3, -1/)
    bst(35)%C(1:2)=(/cmplx(0.0_dp,sq12), -cmplx(sq12,0.0_dp)/)
    !
    ! O: 3F32 = sqrt(1/2)|3/2 1/2> + i*sqrt(1/2)|3/2 -3/2>
    bst(36)%J=-3
    bst(36)%M(1:2)=(/1, -3/)
    bst(36)%C(1:2)=(/cmplx(sq12,0.0_dp), cmplx(0.0_dp,sq12)/)
    !
    ! O: 4F32 =-i*sqrt(1/2)|3/2 3/2> - sqrt(1/2)|3/2 -1/2>
    bst(37)%J=-3
    bst(37)%M(1:2)=(/3, -1/)
    bst(37)%C(1:2)=(/-cmplx(0.0_dp,sq12), -cmplx(sq12,0.0_dp)/)
    !
    ! D6: 1E32
    bst(38)%J=-3
    bst(38)%M(1)=3
    bst(38)%C(1)=1.0_dp
    !
    ! D6: 2E32
    bst(39)%J=-3
    bst(39)%M(1)=-3
    bst(39)%C(1)=1.0_dp
    !
    ! D6: 1E52
    bst(40)%J=-5
    bst(40)%M(1)=-5
    bst(40)%C(1)=1.0_dp
    !
    ! D6: 2E52
    bst(41)%J=-5
    bst(41)%M(1)=5
    bst(41)%C(1)=1.0_dp
    !
    ! D3: 1E3/2 = sqrt(1/2)|3/2 3/2> + sqrt(1/2)|3/2 -3/2>
    bst(42)%J=-3
    bst(42)%M(1:2)=(/3, -3/)
    bst(42)%C(1:2)=(/sq12, sq12/)
    !
    ! D3: 2E3/2 = sqrt(1/2)|3/2 3/2> - sqrt(1/2)|3/2 -3/2>
    bst(43)%J=-3
    bst(43)%M(1:2)=(/3, -3/)
    bst(43)%C(1:2)=(/sq12,-sq12/)
    !
endsubroutine init_pntgrp_bst
!
!**** set_pntgrp_irep() *********************************************
!
! Set irep matrix and character for a given point group
!   For now, only single group irreps are supported
!
!********************************************************************
subroutine set_pntgrp_irep(nop,op,gname,nspin,cax,aax,nirep,irep)
    integer,intent(in) :: nop                   ! Number of group operations
    type(sym_op),intent(in) :: op(nop)          ! The group operations
    character(len=*),intent(in) :: gname        ! The point group name
    integer,intent(in) :: nspin                 ! NSOC irreps or SOC irreps
    real(kind=dp),intent(in) :: cax(3),aax(3)   ! The coordinates of the basis
    integer,intent(out) :: nirep                ! Number of irrep
    type(irep_type),intent(out) :: irep(:)      ! Contains the matrix representations
                                                !         and characters
    !
    integer :: ii, jj, nn, mm, gg, xi, Jt
    real(dp) :: bax(3), tmat(3,3), rmat(3,3)
    real(dp) :: alph(nop), beta(nop), gamm(nop)
    complex(dp)  :: dmatg(maxdim,maxdim), chtmp(48), matmp(maxdim,maxdim,48), &
                    dmat(-5:5,-5:5), cha2sum
    !
    ! Set basis =====================================================
    !
matmp=0.0_dp
chtmp=0.0_dp
if (nspin==1) then
    select case(trim(gname))
    case ('O')
        nirep=5
        irep( 1)=irep_type(   'A1',1, 0, (/bst( 1),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 2)=irep_type(   'A2',1, 0, (/bst( 2),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 3)=irep_type(    'E',2, 0, (/bst( 3),bst( 4),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 4)=irep_type(   'T1',3, 0, (/bst( 5),bst( 6),bst( 7),bst( 1)/), matmp,chtmp)
        irep( 5)=irep_type(   'T2',3, 0, (/bst( 8),bst( 9),bst(10),bst( 1)/), matmp,chtmp)
    case ('Oh')
        nirep=10
        irep( 1)=irep_type(  'A1g',1, 1, (/bst( 1),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 2)=irep_type(  'A2g',1, 1, (/bst( 2),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 3)=irep_type(   'Eg',2, 1, (/bst( 3),bst( 4),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 4)=irep_type(  'T1g',3, 1, (/bst( 5),bst( 6),bst( 7),bst( 1)/), matmp,chtmp)
        irep( 5)=irep_type(  'T2g',3, 1, (/bst( 8),bst( 9),bst(10),bst( 1)/), matmp,chtmp)
        irep( 6)=irep_type(  'A1u',1,-1, (/bst( 1),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 7)=irep_type(  'A2u',1,-1, (/bst( 2),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 8)=irep_type(   'Eu',2,-1, (/bst( 3),bst( 4),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 9)=irep_type(  'T1u',3,-1, (/bst( 5),bst( 6),bst( 7),bst( 1)/), matmp,chtmp)
        irep(10)=irep_type(  'T2u',3,-1, (/bst( 8),bst( 9),bst(10),bst( 1)/), matmp,chtmp)
    case ('T')
        nirep=4
        irep( 1)=irep_type(    'A',1, 0, (/bst( 1),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 2)=irep_type(   '1E',1, 0, (/bst( 3),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 3)=irep_type(   '2E',1, 0, (/bst( 4),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 4)=irep_type(    'T',3, 0, (/bst( 5),bst( 6),bst( 7),bst( 1)/), matmp,chtmp)
    case ('Th')
        nirep=8
        irep( 1)=irep_type(   'Ag',1, 1, (/bst( 1),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 2)=irep_type(  '1Eg',1, 1, (/bst( 3),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 3)=irep_type(  '2Eg',1, 1, (/bst( 4),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 4)=irep_type(   'Tg',3, 1, (/bst( 5),bst( 6),bst( 7),bst( 1)/), matmp,chtmp)
        irep( 5)=irep_type(   'Au',1,-1, (/bst( 1),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 6)=irep_type(  '1Eu',1,-1, (/bst( 3),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 7)=irep_type(  '2Eu',1,-1, (/bst( 4),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 8)=irep_type(   'Tu',3,-1, (/bst( 5),bst( 6),bst( 7),bst( 1)/), matmp,chtmp)
    case ('Td')
        nirep=5
        irep( 1)=irep_type(   'A1',1, 1, (/bst( 1),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 2)=irep_type(   'A2',1, 1, (/bst( 2),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 3)=irep_type(    'E',2, 1, (/bst( 3),bst( 4),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 4)=irep_type(   'T1',3,-1, (/bst( 8),bst( 9),bst(10),bst( 1)/), matmp,chtmp)
        irep( 5)=irep_type(   'T2',3,-1, (/bst( 5),bst( 6),bst( 7),bst( 1)/), matmp,chtmp)
    case ('D6')
        nirep=6
        irep( 1)=irep_type(   'A1',1, 0, (/bst( 1),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 2)=irep_type(   'A2',1, 0, (/bst( 6),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 3)=irep_type(   'B1',1, 0, (/bst(11),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 4)=irep_type(   'B2',1, 0, (/bst(12),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 5)=irep_type(   'E1',2, 0, (/bst(13),bst(14),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 6)=irep_type(   'E2',2, 0, (/bst(15),bst(16),bst( 1),bst( 1)/), matmp,chtmp)
    case ('D6h')
        nirep=12
        irep( 1)=irep_type(  'A1g',1, 1, (/bst( 1),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 2)=irep_type(  'A2g',1, 1, (/bst( 6),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 3)=irep_type(  'B1g',1, 1, (/bst(11),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 4)=irep_type(  'B2g',1, 1, (/bst(12),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 5)=irep_type(  'E1g',2, 1, (/bst(13),bst(14),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 6)=irep_type(  'E2g',2, 1, (/bst(15),bst(16),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 7)=irep_type(  'A1u',1,-1, (/bst( 1),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 8)=irep_type(  'A2u',1,-1, (/bst( 6),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 9)=irep_type(  'B1u',1,-1, (/bst(11),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
        irep(10)=irep_type(  'B2u',1,-1, (/bst(12),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
        irep(11)=irep_type(  'E1u',2,-1, (/bst(13),bst(14),bst( 1),bst( 1)/), matmp,chtmp)
        irep(12)=irep_type(  'E2u',2,-1, (/bst(15),bst(16),bst( 1),bst( 1)/), matmp,chtmp)
    case ('D4')
        nirep=5
        irep( 1)=irep_type(   'A1',1, 0, (/bst( 1),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 2)=irep_type(   'A2',1, 0, (/bst( 6),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 3)=irep_type(   'B1',1, 0, (/bst(17),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 4)=irep_type(   'B2',1, 0, (/bst(18),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 5)=irep_type(    'E',2, 0, (/bst(13),bst(14),bst( 1),bst( 1)/), matmp,chtmp)
    case ('D4h')
        nirep=10
        irep( 1)=irep_type(  'A1g',1, 1, (/bst( 1),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 2)=irep_type(  'A2g',1, 1, (/bst( 6),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 3)=irep_type(  'B1g',1, 1, (/bst(17),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 4)=irep_type(  'B2g',1, 1, (/bst(18),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 5)=irep_type(   'Eg',2, 1, (/bst(13),bst(14),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 6)=irep_type(  'A1u',1,-1, (/bst( 1),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 7)=irep_type(  'A2u',1,-1, (/bst( 6),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 8)=irep_type(  'B1u',1,-1, (/bst(17),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 9)=irep_type(  'B2u',1,-1, (/bst(18),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
        irep(10)=irep_type(   'Eu',2,-1, (/bst(13),bst(14),bst( 1),bst( 1)/), matmp,chtmp)
    case ('D3')
        nirep=3
        irep( 1)=irep_type(   'A1',1, 0, (/bst( 1),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 2)=irep_type(   'A2',1, 0, (/bst( 6),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 3)=irep_type(    'E',2, 0, (/bst(13),bst(14),bst( 1),bst( 1)/), matmp,chtmp)
    case ('D3d')
        nirep=6
        irep( 1)=irep_type(  'A1g',1, 1, (/bst( 1),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 2)=irep_type(  'A2g',1, 1, (/bst( 6),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 3)=irep_type(   'Eg',2, 1, (/bst(13),bst(14),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 4)=irep_type(  'A1u',1,-1, (/bst( 1),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 5)=irep_type(  'A2u',1,-1, (/bst( 6),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 6)=irep_type(   'Eu',2,-1, (/bst(13),bst(14),bst( 1),bst( 1)/), matmp,chtmp)
    case ('D3h')
        nirep=6
        irep( 1)=irep_type(  'A1`',1, 0, (/bst( 1),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 2)=irep_type(  'A2`',1, 0, (/bst(12),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 3)=irep_type(   'E`',2, 0, (/bst(13),bst(14),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 4)=irep_type(  'A1"',1, 1, (/bst(11),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 5)=irep_type(  'A2"',1, 0, (/bst( 6),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 6)=irep_type(   'E"',2, 1, (/bst(13),bst(14),bst( 1),bst( 1)/), matmp,chtmp)
    case ('D2')
        nirep=4
        irep( 1)=irep_type(    'A',1, 0, (/bst( 1),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 2)=irep_type(   'B1',1, 0, (/bst( 6),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 3)=irep_type(   'B2',1, 0, (/bst( 5),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 4)=irep_type(   'B3',1, 0, (/bst( 7),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
    case ('D2h')
        nirep=8
        irep( 1)=irep_type(   'Ag',1, 1, (/bst( 1),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 2)=irep_type(  'B1g',1, 1, (/bst( 6),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 3)=irep_type(  'B2g',1, 1, (/bst( 5),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 4)=irep_type(  'B3g',1, 1, (/bst( 7),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 5)=irep_type(   'Au',1,-1, (/bst( 1),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 6)=irep_type(  'B1u',1,-1, (/bst( 6),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 7)=irep_type(  'B2u',1,-1, (/bst( 5),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 8)=irep_type(  'B3u',1,-1, (/bst( 7),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
    case ('D2d')
        nirep=5
        irep( 1)=irep_type(   'A1',1, 0, (/bst( 1),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 2)=irep_type(   'A2',1,-1, (/bst(18),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 3)=irep_type(   'B1',1, 0, (/bst(17),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 4)=irep_type(   'B2',1, 0, (/bst( 6),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 5)=irep_type(    'E',2, 0, (/bst(14),bst(13),bst( 1),bst( 1)/), matmp,chtmp)
    case ('C6')
        nirep=6
        irep( 1)=irep_type(    'A',1, 0, (/bst( 1),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 2)=irep_type(    'B',1, 0, (/bst(11),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 3)=irep_type(  '1E1',1, 0, (/bst(13),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 4)=irep_type(  '2E1',1, 0, (/bst(14),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 5)=irep_type(  '1E2',1, 0, (/bst(15),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 6)=irep_type(  '2E2',1, 0, (/bst(16),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
    case ('C6h')
        nirep=12
        irep( 1)=irep_type(   'Ag',1, 1, (/bst( 1),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 2)=irep_type(   'Bg',1, 1, (/bst(11),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 3)=irep_type( '1E1g',1, 1, (/bst(13),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 4)=irep_type( '2E1g',1, 1, (/bst(14),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 5)=irep_type( '1E2g',1, 1, (/bst(15),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 6)=irep_type( '2E2g',1, 1, (/bst(16),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 7)=irep_type(   'Au',1,-1, (/bst( 1),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 8)=irep_type(   'Bu',1,-1, (/bst(11),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 9)=irep_type( '1E1u',1,-1, (/bst(13),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
        irep(10)=irep_type( '2E1u',1,-1, (/bst(14),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
        irep(11)=irep_type( '1E2u',1,-1, (/bst(15),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
        irep(12)=irep_type( '2E2u',1,-1, (/bst(16),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
    case ('C6v')
        nirep=6
        irep( 1)=irep_type(   'A1',1, 0, (/bst( 1),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 2)=irep_type(   'A2',1, 1, (/bst( 6),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 3)=irep_type(   'B1',1,-1, (/bst(11),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 4)=irep_type(   'B2',1,-1, (/bst(12),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 5)=irep_type(   'E1',2, 0, (/bst(13),bst(14),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 6)=irep_type(   'E2',2, 0, (/bst(15),bst(16),bst( 1),bst( 1)/), matmp,chtmp)
    case ('C4')
        nirep=4
        irep( 1)=irep_type(    'A',1, 0, (/bst( 1),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 2)=irep_type(    'B',1, 0, (/bst(16),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 3)=irep_type(   '1E',1, 0, (/bst(13),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 4)=irep_type(   '2E',1, 0, (/bst(14),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
    case ('C4h')
        nirep=8
        irep( 1)=irep_type(   'Ag',1, 1, (/bst( 1),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 2)=irep_type(   'Bg',1, 1, (/bst(16),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 3)=irep_type(  '1Eg',1, 1, (/bst(13),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 4)=irep_type(  '2Eg',1, 1, (/bst(14),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 5)=irep_type(   'Au',1,-1, (/bst( 1),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 6)=irep_type(   'Bu',1,-1, (/bst(16),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 7)=irep_type(  '1Eu',1,-1, (/bst(13),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 8)=irep_type(  '2Eu',1,-1, (/bst(14),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
    case ('C4v')
        nirep=5
        irep( 1)=irep_type(   'A1',1, 0, (/bst( 1),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 2)=irep_type(   'A2',1, 1, (/bst( 6),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 3)=irep_type(   'B1',1, 0, (/bst(17),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 4)=irep_type(   'B2',1, 0, (/bst(18),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 5)=irep_type(    'E',2, 0, (/bst(13),bst(14),bst( 1),bst( 1)/), matmp,chtmp)
    case ('C3')
        nirep=3
        irep( 1)=irep_type(    'A',1, 0, (/bst( 1),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 2)=irep_type(   '1E',1, 0, (/bst(13),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 3)=irep_type(   '2E',1, 0, (/bst(14),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
    case ('C3h')
        nirep=6
        irep( 1)=irep_type(   'A`',1, 0, (/bst( 1),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 2)=irep_type(  '1E`',1, 0, (/bst(13),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 3)=irep_type(  '2E`',1, 0, (/bst(14),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 4)=irep_type(   'A"',1, 0, (/bst( 6),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 5)=irep_type( '1E1"',1, 1, (/bst(13),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 6)=irep_type( '2E1"',1, 1, (/bst(14),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
    case ('C3v')    ! C3v-B set
        nirep=3
        irep( 1)=irep_type(   'A1',1, 0, (/bst( 1),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 2)=irep_type(   'A2',1, 0, (/bst(11),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 3)=irep_type(    'E',2, 0, (/bst(13),bst(14),bst( 1),bst( 1)/), matmp,chtmp)
    case ('C2')
        nirep=2
        irep( 1)=irep_type(    'A',1, 0, (/bst( 1),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 2)=irep_type(    'B',1, 0, (/bst(13),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
    case ('C2h')
        nirep=4
        irep( 1)=irep_type(   'Ag',1, 1, (/bst( 1),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 2)=irep_type(   'Bg',1, 1, (/bst(13),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 3)=irep_type(   'Au',1,-1, (/bst( 1),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 4)=irep_type(   'Bu',1,-1, (/bst(13),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
    case ('C2v')
        nirep=4
        irep( 1)=irep_type(   'A1',1, 0, (/bst( 1),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 2)=irep_type(   'A2',1, 0, (/bst(18),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 3)=irep_type(   'B1',1, 0, (/bst( 7),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 4)=irep_type(   'B2',1, 0, (/bst( 5),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
    case ('S6')
        nirep=6
        irep( 1)=irep_type(   'Ag',1, 0, (/bst( 1),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 2)=irep_type(  '1Eg',1, 0, (/bst(15),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 3)=irep_type(  '2Eg',1, 0, (/bst(16),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 4)=irep_type(   'Au',1, 0, (/bst( 6),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 5)=irep_type(  '1Eu',1, 0, (/bst(13),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 6)=irep_type(  '2Eu',1, 0, (/bst(14),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
    case ('S4')
        nirep=4
        irep( 1)=irep_type(    'A',1, 0, (/bst( 1),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 2)=irep_type(    'B',1, 0, (/bst( 6),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 3)=irep_type(   '1E',1, 0, (/bst(14),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 4)=irep_type(   '2E',1, 0, (/bst(13),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
    case ('Ci')
        nirep=2
        irep( 1)=irep_type(   'Ag',1, 0, (/bst( 1),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 2)=irep_type(   'Au',1, 0, (/bst( 6),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
    case ('Cs')
        nirep=2
        irep( 1)=irep_type(   'A`',1, 0, (/bst( 1),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 2)=irep_type(   'A"',1, 0, (/bst( 6),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
    case ('C1')
        nirep=1
        irep( 1)=irep_type(    'A',1, 0, (/bst( 1),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
    case default
    endselect
elseif (nspin==2) then
    select case(trim(gname))
    case ('O')
        nirep=3
        irep( 1)=irep_type(  'E12',2, 0, (/bst(30),bst(31),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 2)=irep_type(  'E52',2, 0, (/bst(32),bst(33),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 3)=irep_type(  'F32',4, 0, (/bst(34),bst(35),bst(36),bst(37)/), matmp,chtmp)
    case ('Oh')
        nirep=6
        irep( 1)=irep_type( 'E12g',2, 1, (/bst(30),bst(31),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 2)=irep_type( 'E52g',2, 1, (/bst(32),bst(33),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 3)=irep_type( 'F32g',4, 1, (/bst(34),bst(35),bst(36),bst(37)/), matmp,chtmp)
        irep( 4)=irep_type( 'E12u',2,-1, (/bst(30),bst(31),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 5)=irep_type( 'E52u',2,-1, (/bst(32),bst(33),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 6)=irep_type( 'F32u',4,-1, (/bst(34),bst(35),bst(36),bst(37)/), matmp,chtmp)
    case ('T')
        nirep=3
        irep( 1)=irep_type(  'E12',2, 0, (/bst(30),bst(31),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 2)=irep_type( '1F32',2, 0, (/bst(34),bst(35),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 3)=irep_type( '2F32',2, 0, (/bst(36),bst(37),bst( 1),bst( 1)/), matmp,chtmp)
    case ('Th')
        nirep=6
        irep( 1)=irep_type( 'E12g',2, 1, (/bst(30),bst(31),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 2)=irep_type('1F32g',2, 1, (/bst(34),bst(35),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 3)=irep_type('2F32g',2, 1, (/bst(36),bst(37),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 4)=irep_type( 'E12u',2,-1, (/bst(30),bst(31),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 5)=irep_type('1F32u',2,-1, (/bst(34),bst(35),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 6)=irep_type('2F32u',2,-1, (/bst(36),bst(37),bst( 1),bst( 1)/), matmp,chtmp)
    case ('Td')
        nirep=3
        irep( 1)=irep_type(  'E12',2, 0, (/bst(30),bst(31),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 2)=irep_type(  'E52',2, 0, (/bst(32),bst(33),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 3)=irep_type(  'F32',4, 0, (/bst(34),bst(35),bst(36),bst(37)/), matmp,chtmp)
    case ('D6')
        nirep=3
        irep( 1)=irep_type(  'E12',2, 0, (/bst(30),bst(31),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 2)=irep_type(  'E32',2, 0, (/bst(38),bst(39),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 3)=irep_type(  'E52',2, 0, (/bst(40),bst(41),bst( 1),bst( 1)/), matmp,chtmp)
    case ('D6h')
        nirep=6
        irep( 1)=irep_type( 'E12g',2, 1, (/bst(30),bst(31),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 2)=irep_type( 'E32g',2, 1, (/bst(38),bst(39),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 3)=irep_type( 'E52g',2, 1, (/bst(40),bst(41),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 4)=irep_type( 'E12u',2,-1, (/bst(30),bst(31),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 5)=irep_type( 'E32u',2,-1, (/bst(38),bst(39),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 6)=irep_type( 'E52u',2,-1, (/bst(40),bst(41),bst( 1),bst( 1)/), matmp,chtmp)
    case ('D4')
        nirep=2
        irep( 1)=irep_type(  'E12',2, 0, (/bst(30),bst(31),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 2)=irep_type(  'E32',2, 0, (/bst(38),bst(39),bst( 1),bst( 1)/), matmp,chtmp)
    case ('D4h')
        nirep=4
        irep( 1)=irep_type( 'E12g',2, 1, (/bst(30),bst(31),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 2)=irep_type( 'E32g',2, 1, (/bst(38),bst(39),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 3)=irep_type( 'E12u',2,-1, (/bst(30),bst(31),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 4)=irep_type( 'E32u',2,-1, (/bst(38),bst(39),bst( 1),bst( 1)/), matmp,chtmp)
    case ('D3')
        nirep=3
        irep( 1)=irep_type(  'E12',2, 0, (/bst(30),bst(31),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 2)=irep_type( '1E32',1, 0, (/bst(42),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 3)=irep_type( '2E32',1, 0, (/bst(43),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
    case ('D3d')
        nirep=6
        irep( 1)=irep_type( 'E12g',2, 1, (/bst(30),bst(31),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 2)=irep_type('1E32g',1, 1, (/bst(42),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 3)=irep_type('2E32g',1, 1, (/bst(43),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 4)=irep_type( 'E12u',2,-1, (/bst(30),bst(31),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 5)=irep_type('1E32u',1,-1, (/bst(42),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 6)=irep_type('2E32u',1,-1, (/bst(43),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
    case ('D3h')
        nirep=3
        irep( 1)=irep_type(  'E12',2, 0, (/bst(30),bst(31),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 2)=irep_type(  'E32',2, 0, (/bst(38),bst(39),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 3)=irep_type(  'E52',2, 0, (/bst(40),bst(41),bst( 1),bst( 1)/), matmp,chtmp)
    case ('D2')
        nirep=1
        irep( 1)=irep_type(  'E12',2, 0, (/bst(30),bst(31),bst( 1),bst( 1)/), matmp,chtmp)
    case ('D2h')
        nirep=2
        irep( 1)=irep_type( 'E12g',2, 1, (/bst(30),bst(31),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 2)=irep_type( 'E12u',2,-1, (/bst(30),bst(31),bst( 1),bst( 1)/), matmp,chtmp)
    case ('D2d')
        nirep=2
        irep( 1)=irep_type(  'E12',2, 0, (/bst(30),bst(31),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 2)=irep_type(  'E32',2, 0, (/bst(38),bst(39),bst( 1),bst( 1)/), matmp,chtmp)
    case ('C6')
        nirep=6
        irep( 1)=irep_type( '1E12',1, 0, (/bst(31),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 2)=irep_type( '2E12',1, 0, (/bst(30),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 3)=irep_type( '1E32',1, 0, (/bst(38),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 4)=irep_type( '2E32',1, 0, (/bst(39),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 5)=irep_type( '1E52',1, 0, (/bst(40),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 6)=irep_type( '2E52',1, 0, (/bst(41),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
    case ('C6h')
        nirep=12
        irep( 1)=irep_type('1E12g',1, 1, (/bst(31),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 2)=irep_type('2E12g',1, 1, (/bst(30),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 3)=irep_type('1E32g',1, 1, (/bst(38),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 4)=irep_type('2E32g',1, 1, (/bst(39),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 5)=irep_type('1E52g',1, 1, (/bst(40),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 6)=irep_type('2E52g',1, 1, (/bst(41),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 7)=irep_type('1E12u',1,-1, (/bst(31),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 8)=irep_type('2E12u',1,-1, (/bst(30),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 9)=irep_type('1E32u',1,-1, (/bst(38),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
        irep(10)=irep_type('2E32u',1,-1, (/bst(39),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
        irep(11)=irep_type('1E52u',1,-1, (/bst(40),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
        irep(12)=irep_type('2E52u',1,-1, (/bst(41),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
    case ('C6v')
        nirep=3
        irep( 1)=irep_type(  'E12',2, 0, (/bst(30),bst(31),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 2)=irep_type(  'E32',2, 0, (/bst(38),bst(39),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 3)=irep_type(  'E52',2, 0, (/bst(40),bst(41),bst( 1),bst( 1)/), matmp,chtmp)
    case ('C4')
        nirep=4
        irep( 1)=irep_type( '1E12',1, 0, (/bst(30),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 2)=irep_type( '2E12',1, 0, (/bst(31),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 3)=irep_type( '1E32',1, 0, (/bst(39),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 4)=irep_type( '2E32',1, 0, (/bst(38),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
    case ('C4h')
        nirep=8
        irep( 1)=irep_type('1E12g',1, 1, (/bst(30),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 2)=irep_type('2E12g',1, 1, (/bst(31),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 3)=irep_type('1E32g',1, 1, (/bst(39),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 4)=irep_type('2E32g',1, 1, (/bst(38),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 5)=irep_type('1E12u',1,-1, (/bst(30),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 6)=irep_type('2E12u',1,-1, (/bst(31),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 7)=irep_type('1E32u',1,-1, (/bst(39),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 8)=irep_type('2E32u',1,-1, (/bst(38),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
    case ('C4v')
        nirep=2
        irep( 1)=irep_type(  'E12',2, 0, (/bst(30),bst(31),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 2)=irep_type(  'E32',2, 0, (/bst(38),bst(39),bst( 1),bst( 1)/), matmp,chtmp)
    case ('C3')
        nirep=3
        irep( 1)=irep_type( '1E12',1, 0, (/bst(31),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 2)=irep_type( '2E12',1, 0, (/bst(30),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 3)=irep_type(  'A32',1, 0, (/bst(38),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
    case ('C3h')
        nirep=6
        irep( 1)=irep_type( '1E12',1, 0, (/bst(31),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 2)=irep_type( '2E12',1, 0, (/bst(30),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 3)=irep_type( '1E32',1, 0, (/bst(38),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 4)=irep_type( '2E32',1, 0, (/bst(39),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 5)=irep_type( '1E52',1, 0, (/bst(40),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 6)=irep_type( '2E52',1, 0, (/bst(41),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
    case ('C3v')    ! C3v-B set
        nirep=3
        irep( 1)=irep_type(  'E12',2, 0, (/bst(30),bst(31),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 2)=irep_type( '1E32',1, 0, (/bst(42),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 3)=irep_type( '2E32',1, 0, (/bst(43),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
    case ('C2')
        nirep=2
        irep( 1)=irep_type( '1E12',1, 0, (/bst(31),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 2)=irep_type( '2E12',1, 0, (/bst(30),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
    case ('C2h')
        nirep=4
        irep( 1)=irep_type('1E12g',1, 1, (/bst(31),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 2)=irep_type('2E12g',1, 1, (/bst(30),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 3)=irep_type('1E12u',1,-1, (/bst(31),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 4)=irep_type('2E12u',1,-1, (/bst(30),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
    case ('C2v')
        nirep=1
        irep( 1)=irep_type(  'E12',2, 0, (/bst(30),bst(31),bst( 1),bst( 1)/), matmp,chtmp)
    case ('S6')
        nirep=6
        irep( 1)=irep_type('1E12g',1, 1, (/bst(31),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 2)=irep_type('2E12g',1, 1, (/bst(30),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 3)=irep_type( 'A32g',1, 1, (/bst(38),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 4)=irep_type('1E12u',1,-1, (/bst(31),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 5)=irep_type('2E12u',1,-1, (/bst(30),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 6)=irep_type( 'A32u',1,-1, (/bst(38),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
    case ('S4')
        nirep=4
        irep( 1)=irep_type( '1E12',1, 0, (/bst(30),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 2)=irep_type( '2E12',1, 0, (/bst(31),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 3)=irep_type( '1E32',1, 0, (/bst(39),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 4)=irep_type( '2E32',1, 0, (/bst(38),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
    case ('Ci')
        nirep=2
        irep( 1)=irep_type( 'A12g',1, 1, (/bst(30),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 2)=irep_type( 'A12u',1,-1, (/bst(30),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
    case ('Cs')
        nirep=2
        irep( 1)=irep_type( '1E12',1, 0, (/bst(31),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
        irep( 2)=irep_type( '2E12',1, 0, (/bst(30),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
    case ('C1')
        nirep=1
        irep( 1)=irep_type(  'A12',1, 0, (/bst(30),bst( 1),bst( 1),bst( 1)/), matmp,chtmp)
    case default
    endselect
endif
    !
    do xi=1,nirep
        !
        ! Beware that: 
        !   1. For a given point group irrep, all bases have same Jt !!!
        !   2. Jt>=0: integer momentum; Jt<0: half-integer momentum
        !
        if( irep(xi)%P==0 ) then
            if ( irep(xi)%bs(1)%J >=0 ) then
                irep(xi)%P = fun_gu(irep(xi)%bs(1)%J )
            else
                irep(xi)%P = 1
            endif
        endif
        !
        do ii=1,irep(xi)%D
            irep(xi)%bs(ii)%P = irep(xi)%P
        enddo
    enddo
    !
    ! Calculate Euler angles in cax aax ====================
    !
    bax=x_prod(cax,aax)
    tmat(:,1)=aax; tmat(:,2)=bax; tmat(:,3)=cax
    if( .not. othmat(tmat,3) ) then
        write(*,*) 'Error in set_irrep(): wrong cax and aax !!!'
        stop
    endif
    !
    do gg=1,nop
        rmat=matmul( op(gg)%R, tmat )
        rmat=matmul(transpose(tmat),rmat)
        rmat=rmat*op(gg)%det    ! make rmat positive definite 
        call euler_angle( rmat(:,3),rmat(:,1),alph(gg),beta(gg),gamm(gg) )
    enddo
    !
    ! Set matrix ====================================================
    !
    do xi=1,nirep   ! xi loop 
        !
        !************************************************************
        ! Begin of xi Loop
        !************************************************************
        !
        Jt=irep(xi)%bs(1)%J
        !
        cha2sum=0.0_dp
        irep(xi)%cha(:)=0.0_dp
        do gg=1,nop
            !
            ! Matrix Representation ..................................
            !
            call dmat_cmplx(alph(gg),beta(gg),gamm(gg), &
                            dmat(-abs(Jt):abs(Jt),-abs(Jt):abs(Jt)), Jt)
            !
            dmatg(:,:) = 0.0_dp
            do ii=1,irep(xi)%D
            do jj=1,irep(xi)%D
                do nn=1,maxbcp
                do mm=1,maxbcp
                    dmatg(ii,jj) = dmatg(ii,jj) &
                      + dmat( irep(xi)%bs(ii)%M(nn), irep(xi)%bs(jj)%M(mm) ) &
                      * conjg(irep(xi)%bs(ii)%C(nn))*irep(xi)%bs(jj)%C(mm)
                enddo
                enddo
            enddo
            enddo
            !
            if (op(gg)%det==-1 .and. irep(xi)%P==-1) dmatg = -dmatg
            !
            if (.not. unimat( dmatg(1:irep(xi)%D,1:irep(xi)%D), irep(xi)%D ) ) then
                write(*,*) 'Error in set_irep(): irep%mat of', trim(irep(xi)%name),&
                           'is not uniary'
                stop
            endif
            !
            irep(xi)%mat(1:maxdim,1:maxdim,gg)=dmatg(1:maxdim,1:maxdim)
            !
            ! Character =================================================
            !
            do ii=1,irep(xi)%D
                irep(xi)%cha(gg)= irep(xi)%cha(gg) + dmatg(ii,ii) 
            enddo
            !
            cha2sum = cha2sum + conjg(irep(xi)%cha(gg))*irep(xi)%cha(gg)
            !
        enddo
        !print*, cha2sum, nop
        !
        ! Check if irrep is irreducible==================================
        !
        if ( nint(abs(cha2sum)) /= nop ) then
            write(*,*) 'Error in set_irep(): ', trim(irep(xi)%name), '  is reducible !!!'
            stop
        endif
        !
        !************************************************************
        ! End of xi Loop
        !************************************************************
    enddo       
    !
    ! Check Orthogonality ===========================================
    !
    do ii=1,nirep
    do jj=1,ii-1
        cha2sum=0.0_dp
        do gg=1,nop
            cha2sum=cha2sum+conjg(irep(ii)%cha(gg))*irep(jj)%cha(gg)
        enddo
        if( abs(cha2sum)> nop*1.0e-3 ) then
            write(*,*) 'Error in set_irep(): ', ii,jj, '-th irreps are not orthogonal !!!' 
            write(*,*) cha2sum
500         format(12(2F8.4,3X))
            write(*,500) irep(ii)%cha(1:nop)
            write(*,500) irep(jj)%cha(1:nop)
            write(*,'(12F8.4)') alph(1:nop)
            write(*,'(12F8.4)') beta(1:nop)
            write(*,'(12F8.4)') gamm(1:nop)
            write(*,'(12I8)')   op(1:nop)%det
            stop
        endif
    enddo
    enddo
    !
    ! Complicity ====================================================
    !
    nn=0
    do ii=1,nirep
        nn = nn + irep(ii)%D*irep(ii)%D
    enddo
    if (nn/=nop) then
        write(*,*) 'Error in set_irep(): irreps set is not complete !!!'
        stop
    endif
    !
endsubroutine set_pntgrp_irep
!
!
subroutine Dmat_rl_op(op,Dmat,L)
    type(sym_op) :: op
    real(kind=dp), dimension(-L:L,-L:L) :: Dmat
    integer      :: L
    !   
    call Dmat_rl(op%alph, op%beta, op%gamm, Dmat, L ) 
    if ( op%det==-1 .and. mod(L,2)==1) Dmat=-Dmat
    !   
endsubroutine
!
! real D matrix, on real spherical harmonics
subroutine Dmat_rl(alph,beta,gamm,Dmat,L)
    !   
    real(kind=dp) :: alph, beta, gamm
    real(kind=dp), dimension(-L:L,-L:L) :: Dmat, dmat_
    integer :: L
    !   
    integer :: ii, jj
    complex(kind=dp), parameter :: c_q2 = cmplx(sqrt(0.5_dp),0.0_dp), &
                                   c_r1 = cmplx(1.0_dp,0.0_dp)
    complex(kind=dp), dimension(-3:3), save :: &
       aa=(/ c_im*c_q2, c_im*c_q2, c_im*c_q2, 0.5_dp*c_r1, -c_q2, c_q2, -c_q2/), &
       bb=(/ c_im*c_q2,-c_im*c_q2, c_im*c_q2, 0.5_dp*c_r1,  c_q2, c_q2,  c_q2/)
    !   
    call dmat_lt(beta, dmat_, L)
    do ii=-L,L
    do jj=-L,L
        Dmat(ii,jj) = real( conjg(aa(ii))*aa(jj)*exp(-c_im*(alph*ii+gamm*jj)) ) & 
                      *2*dmat_(ii,jj) + & 
                      real( conjg(aa(ii))*bb(jj)*exp(-c_im*(alph*ii-gamm*jj)) ) & 
                      *2*dmat_(ii,-jj)
    enddo
    enddo
    !   
endsubroutine Dmat_rl
!
!**** dmat_cmplx() **************************************************
!
! If L>=0, returns D-matrix of L
! If L<0,  returns D-matrix of abs(L)/2
!
!********************************************************************
subroutine dmat_cmplx(alph,beta,gamm,dmat,L)
    real(kind=dp) :: alph, beta, gamm
    complex(kind=dp), dimension(-abs(L):abs(L),-abs(L):abs(L)) :: dmat
    integer :: L
    !
    integer :: nn, mm
    real(kind=dp),  dimension(-abs(L):abs(L),-abs(L):abs(L)) :: dmat_
    !
    call dmat_lt(beta,dmat_,L)
    if (L>=0) then
        do nn=-L,L
        do mm=-L,L
            dmat(nn,mm)=dmat_(nn,mm)*exp(-c_im*alph*nn)*exp(-c_im*gamm*mm)
        enddo
        enddo
    else
        do nn=-abs(L),abs(L),2
        do mm=-abs(L),abs(L),2
            dmat(nn,mm)=dmat_(nn,mm)*exp(-0.5_dp*c_im*alph*nn)*exp(-0.5_dp*c_im*gamm*mm)
        enddo
        enddo
    endif
    !
endsubroutine
!
!**** dmat_lt() *****************************************************
!
! If L>=0, returns d-matrix of L
! If L<0,  returns d-matrix of abs(L)/2
!
!********************************************************************
subroutine dmat_lt(beta, dmat, L)
    real(kind=dp) :: beta, s,c 
    real(kind=dp), dimension(-abs(L):abs(L),-abs(L):abs(L)) :: dmat
    integer     :: L, ii, jj
    !   
    real(kind=dp) :: c6s0, c5s1, c4s2, c3s3, c2s4, c1s5, c0s6, &
                     c5s0, c4s1, c3s2, c2s3, c1s4, c0s5, &
                     c3s0, c2s1, c1s2, c0s3
    !   
    s=sin(beta)
    c=cos(beta)
    !   
    if (L==0) then
        dmat(0,0)=1.0_dp
    elseif (L==1) then
        dmat(1, 1)= 0.5_dp*(1+c)
        dmat(1, 0)=-s/sqrt(2.0_dp)
        dmat(1,-1)= 0.5_dp*(1-c)
        dmat(0, 0)= c
        !
        dmat(-1,-1)= dmat(1,1)
        dmat(-1, 0)=-dmat(1,0)
        dmat( 0, 1)=-dmat(1,0)
        dmat( 0,-1)= dmat(1,0)
        dmat(-1, 1)= dmat(1,-1)
    elseif (L==2) then
        dmat(2, 2)= 0.25_dp*(1+c)*(1+c)
        dmat(2, 1)=-0.5_dp*s*(1+c)
        dmat(2, 0)= sqrt(3.0/8.0_dp)*s*s
        dmat(2,-1)=-0.5_dp*s*(1-c)
        dmat(2,-2)= 0.25_dp*(1-c)*(1-c)
        dmat(1, 1)= 0.5_dp*(2*c*c+c-1)
        dmat(1, 0)=-sqrt(3.0_dp/8.0_dp)*2*s*c
        dmat(1,-1)= 0.5_dp*(-2*c*c+c+1)
        dmat(0, 0)= 0.5_dp*(3*c*c-1)
        !
        dmat(-2,-2)= dmat(2,2)
        dmat(-2,-1)=-dmat(2,1)
        dmat( 1, 2)=-dmat(2,1)
        dmat(-1,-2)= dmat(2,1)
        dmat(-2, 0)= dmat(2,0)
        dmat( 0, 2)= dmat(2,0)
        dmat( 0,-2)= dmat(2,0)
        dmat(-2, 1)=-dmat(2,-1)
        dmat(-1, 2)=-dmat(2,-1)
        dmat( 1,-2)= dmat(2,-1)
        dmat(-2, 2)= dmat(2,-2)
        dmat(-1,-1)= dmat(1,1)
        dmat(-1, 0)=-dmat(1,0)
        dmat( 0, 1)=-dmat(1,0)
        dmat( 0,-1)= dmat(1,0)
        dmat(-1, 1)= dmat(1,-1)
    elseif( L==3 ) then
        !
        c=cos(beta/2)
        s=sin(beta/2)
        c6s0=c*c*c*c*c*c
        c5s1=c*c*c*c*c*s
        c4s2=c*c*c*c*s*s
        c3s3=c*c*c*s*s*s
        c2s4=c*c*s*s*s*s
        c1s5=c*s*s*s*s*s
        c0s6=s*s*s*s*s*s
        !
        dmat( 3, 3)=c6s0
        dmat( 2, 3)=sqrt(6.0_dp)*c5s1
        dmat( 2, 2)=c6s0-5.0_dp*c4s2
        dmat( 1, 3)=sqrt(15.0_dp)*c4s2
        dmat( 1, 2)=sqrt(10.0_dp)*(c5s1-2.0*c3s3)
        dmat( 1, 1)=c6s0-8.0*c4s2+6.0*c2s4
        dmat( 0, 3)=2.0*sqrt(5.0_dp)*c3s3
        dmat( 0, 2)=sqrt(30.0_dp)*(c4s2-c2s4)
        dmat( 0, 1)=2.0*sqrt(3.0_dp)*(c5s1-3.0*c3s3+c1s5)
        dmat( 0, 0)=c6s0-9.0*c4s2+9.0*c2s4-c0s6
        !
        dmat( 3,-3)=c0s6
        dmat( 2,-3)=-sqrt(6.0_dp)*c1s5
        dmat( 2,-2)=-c0s6+5.0_dp*c2s4
        dmat( 1,-3)=sqrt(15.0_dp)*c2s4
        dmat( 1,-2)=sqrt(10.0_dp)*(c1s5-2.0*c3s3)
        dmat( 1,-1)=c0s6-8.0*c2s4+6.0*c4s2
        dmat( 0,-3)=-2.0*sqrt(5.0_dp)*c3s3
        dmat( 0,-2)=-sqrt(30.0_dp)*(c2s4-c4s2)
        dmat( 0,-1)=-2.0*sqrt(3.0_dp)*(c1s5-3.0*c3s3+c5s1)
        !
        do ii=1,3
            do jj=-ii+1,0
                dmat(ii,jj)=dmat(-jj,-ii)
            enddo
            do jj=1,ii-1
                dmat(ii,jj)=fun_gu(ii-jj)*dmat(jj,ii)
            enddo
        enddo
        !
        do ii=-3,-1
        do jj=-3,3
            dmat(ii,jj)=fun_gu(ii-jj)*dmat(-ii,-jj)
        enddo
        enddo
        !
    elseif (L==-1) then     ! spin 1/2
        c=cos(beta/2); s=sin(beta/2)
        dmat( 1, 1)=c
        dmat( 1,-1)=-s
        dmat(-1, 1)=s
        dmat(-1,-1)=c
    elseif (L==-3) then     ! spin 3/2
        c=cos(beta/2); s=sin(beta/2)
        c3s0=c*c*c
        c2s1=c*c*s
        c1s2=c*s*s
        c0s3=s*s*s
        !
        dmat( 3, 3)=c3s0
        dmat( 1, 3)=sqrt(3.0_dp)*c2s1
        dmat( 1, 1)=c3s0-2.0*c1s2
        !
        dmat( 3,-3)=-c0s3
        dmat( 1,-3)=sqrt(3.0_dp)*c1s2
        dmat( 1,-1)=c0s3-2.0*c2s1
        !
        dmat( 3,-1)=dmat(1,-3)
        dmat( 3, 1)=-dmat(1,3)
        !
        do ii=-3,-1,2
        do jj=-3,3,2
            dmat(ii,jj) = fun_gu((ii-jj)/2)*dmat(-ii,-jj)
        enddo
        enddo
        !
    elseif (L==-5) then     ! spin 5/2
        c=cos(beta/2); s=sin(beta/2)
        c5s0=c*c*c*c*c
        c4s1=c*c*c*c*s
        c3s2=c*c*c*s*s
        c2s3=c*c*s*s*s
        c1s4=c*s*s*s*s
        c0s5=s*s*s*s*s
        !
        dmat( 5, 5)=c5s0
        dmat( 3, 5)=sqrt(5.0_dp)*c4s1
        dmat( 3, 3)=c5s0-4.0*c3s2
        dmat( 1, 5)=sqrt(10.0_dp)*c3s2
        dmat( 1, 3)=sqrt(2.0_dp)*(2.0*c4s1-3.0*c2s3)
        dmat( 1, 1)=c5s0-6.0*c3s2+3.0*c1s4
        !
        dmat( 5,-5)=-c0s5
        dmat( 3,-5)=sqrt(5.0_dp)*c1s4
        dmat( 3,-3)=c0s5-4.0*c2s3
        dmat( 1,-5)=-sqrt(10.0_dp)*c2s3
        dmat( 1,-3)=-sqrt(2.0_dp)*(2.0*c1s4-3.0*c3s2)
        dmat( 1,-1)=-c0s5+6.0*c2s3-3.0*c4s1
        !
        do ii=3,5,2
            do jj=-ii+2,-1,2
                dmat(ii,jj)=dmat(-jj,-ii)
            enddo
            do jj=1,ii-2,2
                dmat(ii,jj)=fun_gu((ii-jj)/2)*dmat(jj,ii)
            enddo
        enddo
        !
        do ii=-5,-1,2
        do jj=-5,5,2
            dmat(ii,jj)=fun_gu((ii-jj)/2)*dmat(-ii,-jj)
        enddo
        enddo
        !
    else
        write(*,*) 'Eorror in dmat_lt(): L is out of range !!!'
    endif
    ! 
endsubroutine dmat_lt
!
! dmat_spin
subroutine dmat_spin(alph,beta,gamm,dmat)
    real(kind=dp)   :: alph, gamm, beta
    complex(kind=dp), dimension(2,2) :: dmat
    !   
    dmat(1,:) = (/ exp(0.5*c_im*(-alph-gamm))*cos(beta/2), -exp(0.5*c_im*(-alph+gamm))*sin(beta/2)  /)  
    dmat(2,:) = (/ exp(0.5*c_im*( alph-gamm))*sin(beta/2),  exp(0.5*c_im*( alph+gamm))*cos(beta/2)  /)  
    !   
endsubroutine
!
! function gu: (-1)^a
function fun_gu(a)
    integer :: a, fun_gu
    if( mod(abs(a),2)==0 ) then
        fun_gu=1
    else
        fun_gu=-1
    endif
endfunction fun_gu
!
! dmat_check() is used for debug
subroutine dmat_check(alph,beta,gamm)
    real(kind=dp) :: alph, beta, gamm
    !
    integer :: LL, J2, nn, mm
    complex(kind=dp), dimension(-5:5,-5:5) :: dmat
    !
    write(*,'("check D matrix:", 3F8.1)') alph/c_pi*180, beta/c_pi*180, gamm/c_pi*180
    !
    do LL=0,3
        !
        call dmat_cmplx(alph,beta,gamm,dmat(-LL:LL,-LL:LL),LL)
        !
        write(*,'("L=",I3)') LL
        do nn=-LL,LL
        do mm=-LL,LL
            if (mm==LL) then
                write(*,'(2F9.4)') dmat(nn,mm)
            else
                write(*,'(2F9.4)', advance='no') dmat(nn,mm)
            endif
        enddo
        enddo
    enddo
    !
    do J2=1,5,2
        !
        call dmat_cmplx(alph,beta,gamm,dmat(-J2:J2,-J2:J2),-J2)
        !
        write(*,'("L=",I3,"/2")') J2
        do nn=-J2,J2,2
        do mm=-J2,J2,2
            if (mm==J2) then
                write(*,'(2F9.4)') dmat(nn,mm)
            else
                write(*,'(2F9.4)',advance='no') dmat(nn,mm)
            endif
        enddo
        enddo
    enddo
    !
endsubroutine
!
!**** extrct_irep() *************************************************
!
! The subroutine to reduce a reducible representation to irreps
!
!********************************************************************
subroutine extrct_irep(nop,drep,rep,irep,nirep,ncgc,cgc)
    integer,intent(in) :: nop                    ! Number of symmetry operations
    integer,intent(in) :: drep                   ! Dimension of rep
    complex(dp),intent(in) :: rep(drep,drep,nop) ! Representation matrices
    type(irep_type),intent(in) :: irep           ! The standard irrep
    integer,intent(out) :: nirep                 ! The multiplicity of irep in rep
    integer,intent(in)  :: ncgc                  !
    complex(dp),intent(out) :: cgc(drep,ncgc)    ! The CG coefficient
    !
    integer :: ii, jj, kk, gg, nn, mm, mmp, cnt, aa
    complex(dp) :: chasum, mat(6,6), phi(drep,6), alph(drep)
    real(dp)    :: rr
    logical :: suc
    !
    cgc(:,:)=0.0_dp
    !
    if (irep%bs(1)%J>=0) then
        !
        ! calculate nirep ===============================================
        !
        chasum=0.0_dp
        do gg=1,nop
        do ii=1,drep
            chasum = chasum + conjg(rep(ii,ii,gg))*irep%cha(gg)
        enddo
        enddo
        chasum=chasum/nop
        nirep = nint(real(chasum))
        !
        if ( abs(nirep-chasum)>1.0e-2 ) then
            write(*,*) 'Error in extct_irep(): rep or irep is not well defined !!!'
            write(*,'("nirep=",I5," chasum=",2F12.6)') nirep, chasum
            stop
        endif
        !
        if (nirep==0) return
        !
        ! extract ireps =================================================
        !
        cnt=0
        do ii=1,drep        ! The outer loop about ii and nn: for each pair of
        do nn=1,irep%D      !   nn and ii, an irrep can be projected out
            suc=.true.
            !
            ! Project out an irreducible representation
            phi(:,:)=0.0_dp
            do mm=1,irep%D
                chasum=0.0_dp
                do jj=1,drep
                    do gg=1,nop
                        phi(jj,mm)=phi(jj,mm)+conjg(rep(ii,jj,gg))*irep%mat(nn,mm,gg)
                    enddo
                    chasum=chasum + conjg(phi(jj,mm))*phi(jj,mm)
                enddo
                !
                if( sqrt(abs(chasum))<1.0e-2 ) then
                    suc=.false.
                    exit
                else
                    phi(:,mm) = phi(:,mm)/sqrt(abs(chasum))
                endif
            enddo
            !
            if (.not. suc) cycle
            !
            ! Inner product 
            alph(:)=0.0_dp
            do aa=1,cnt
                do jj=1,drep
                    alph(aa)=alph(aa) + conjg(cgc(jj,(aa-1)*irep%D+1))*phi(jj,1)
                enddo
            enddo
            !
            ! Orthogonalize
            do aa=1,cnt
            do mm=1,irep%D
                phi(:,mm) = phi(:,mm) - cgc(:,(aa-1)*irep%D+mm)*alph(aa)
            enddo
            enddo
            !
            ! Normalization
            do mm=1,irep%D
                !
                chasum=0.0_dp
                do jj=1,drep
                    chasum = chasum + conjg(phi(jj,mm))*phi(jj,mm)
                enddo
                !
                if (abs(chasum)<1.0e-1_dp) then
                    suc=.false.
                    exit
                endif
                !
                phi(:,mm)=phi(:,mm)/sqrt(abs(chasum))
                !
            enddo
            !
            ! Append it to cgc
            if (suc) then
                !
                ! check
                do gg=1,nop
                    rr=0.0_dp
                    mat(:,:)=0.0_dp
                    do mm =1,irep%D
                    do mmp=1,irep%D
                        do jj=1,drep
                        do kk=1,drep
                            mat(mm,mmp) = mat(mm,mmp) + conjg(phi(jj,mm))*rep(jj,kk,gg)*phi(kk,mmp)
                        enddo
                        enddo
                        !
                        if (abs(mat(mm,mmp)-irep%mat(mm,mmp,gg))>rr) &
                            rr=abs(mat(mm,mmp)-irep%mat(mm,mmp,gg))
                        !
                    enddo
                    enddo
                    !
                    if ( rr>1.0e-2_dp ) then
                        write(*,*) 'Error in extct_irep(): projection operator dose not work !!!'
                        print*,trim(irep%name),cnt,gg
                        write(*,'(4(2F8.4,"  "))') mat(1,1:irep%D)
                        write(*,'(4(2F8.4,"  "))') mat(2,1:irep%D)
                        print*, '****'
                        write(*,'(4(2F8.4,"  "))') irep%mat(1,1:irep%D,gg)
                        write(*,'(4(2F8.4,"  "))') irep%mat(2,1:irep%D,gg)
                        print*, '****'
                        write(*,'(4(2F8.4,"  "))') phi(:,1)
                        write(*,'(4(2F8.4,"  "))') phi(:,2)
                        stop
                    endif
                    !
                enddo
                !
                cgc(:,cnt*irep%D+1:cnt*irep%D+irep%D) = phi(:,1:irep%D)
                cnt=cnt+1
                !print*, cnt, nirep
                !
            endif
            !
            if (cnt==nirep) exit    ! If we have collected enough irreps, quit
        enddo
            if (cnt==nirep) exit    ! If we have collected enough irreps, quit
        enddo
        !
        if (cnt<nirep) then
            write(*,*) 'Error in extct_irep(): can not find all irreps !!!'
            stop
        endif
    else
        write(*,*) 'Error: For now extrct_irep() does not support double group representation !!!'
        stop
    endif
    !
endsubroutine extrct_irep

!**** int3Ylm *********************************************************
!
! \int d(theta)d(phi)sin(theta)Y(j1m1)Y(j2m2)Y(j3m3) =
!  /          \ /          \
!  | J1 J2 J3 | | J1 J2 J3 |    sqrt( (2J1+1)(2J2+1)(2J3+1) )
!  | 0  0  0  | | M1 M2 M3 |  ---------------------------------
!  \          / \          /             sqrt(4*pi)
!********************************************************************
function int3Ylm(j1,j2,j3,m1,m2,m3)
    integer :: j1, j2, j3, m1, m2, m3
    real(dp) :: int3Ylm
    !
    int3Ylm = sqrt((2*j1+1)*(2*j2+1)*(2*j3+1)/(4*c_pi)) &
            * wgn3j(j1,j2,j3,0,0,0)*wgn3j(j1,j2,j3,m1,m2,m3)
    !
endfunction

!**** wgn3j *********************************************************
!  /          \
! | J1  J2  J3 |
! | M1  M2  M3 |
!  \          /
!********************************************************************
function wgn3j(J1,J2,J3,M1,M2,M3)
    integer :: J1, J2, J3, M1, M2, M3
    REAL(dp) :: wgn3j
    !
    wgn3j = fun_gu(J1-J2-M3)*clebg(J1,J2,J3,M1,M2,-M3)/sqrt(2.0_dp*J3+1)
    !   
endfunction wgn3j

!**************** FUNCTION CLEBGO ***********************************
!
! caculate Clebsch-Gordan-coeff. <J1 J2 M1 M2 | J3 M3>
! using racah-formel
!
!********************************************************************
function clebg(j1,j2,j3,m1,m2,m3)
    integer :: J1, J2, J3, M1, M2, M3
    REAL(dp) :: CLEBG
    !   
    integer :: K1,K2,K3,K4,K5,K6,N1,N2, M, N
    real(dp) :: T,T1
    !   
    IF(M3/=M1+M2) GO TO 2
    K1=J1+J2-J3+1
    K2=J3+J1-J2+1
    K3=J3+J2-J1+1
    K4=J1+J2+J3+2
    T= (2*J3+1)*FAC(K1)*FAC(K2)*FAC(K3)/FAC(K4)
    K1=J1+M1+1
    K2=J1-M1+1
    K3=J2+M2+1
    K4=J2-M2+1
    K5=J3+M3+1
    K6=J3-M3+1
    T=SQRT(T*FAC(K1)*FAC(K2)*FAC(K3)*FAC(K4)*FAC(K5)*FAC(K6))
    N1=MAX0(J2-J3-M1,J1-J3+M2,0)+1
    N2=MIN0(J1+J2-J3,J1-M1,J2+M2)+1
    IF(N1>N2) GOTO 2
    T1=0.0_dp
    DO M=N1,N2
       N=M-1
       K1=J1+J2-J3-N+1
       K2=J1-M1-N+1
       K3=J2+M2-N+1
       K4=J3-J2+M1+N+1
       K5=J3-J1-M2+N+1
       T1=T1+ (1+4*(N/2)-2*N)/(FAC(M)*FAC(K1)*FAC(K2)*FAC(K3) &
            &  *FAC(K4)*FAC(K5))
    ENDDO
    CLEBG=T*T1
    RETURN
    ! coefficient is zero, drop back
2   CONTINUE
    CLEBG=0.0_dp
    RETURN
endfunction

!**** fac_init ******************************************************
!
!   Initialize the factorials
!
!********************************************************************
subroutine fac_init()
    integer :: ii
    fac(1) = 1.0_dp
    do ii=1,39
        fac(ii+1) = fac(ii)*ii      ! fac(n)=(n-1)!
    enddo
endsubroutine

endmodule reprsnt
