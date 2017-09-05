module slicewann
!
use typcnst,    only: DP, pat_type, wan_type, lsize, sym_op
use pntgrp,     only: TonA, find_ltgrp_r, pgrp_name, euler_angle
!
use tghtbd,     only: cell_a, cell_b
use symm,       only: nspgop, spgop
!
implicit none
!
public
!
! kz and bands
real(dp)                        :: kz_sls(3)        ! kz of the 2D slice
integer                         :: bstart, bend     ! bands to construct WFs
!
! pseudo-atoms and swannier orbitals
integer                         :: nswan            ! number of swannier functions
integer                         :: nsat             ! number of pseudo-atoms
type(pat_type)                  :: sat(lsize)       ! pseudo-atoms
type(wan_type),allocatable      :: swan(:)          ! swannier functions
!
contains
!
subroutine get_swan()
    !
    ! read-in variables
    character(len=8)    :: k_conv, sat_pos_conv
    integer             :: sat_kind(lsize)
    real(dp)            :: sat_pos(3,lsize)
    character(len=8)    :: swan_list(3,lsize)
    !
    namelist /SliceWannier/ k_conv, kz_sls, bstart, bend, sat_kind, &
                            sat_pos_conv, sat_pos, swan_list
    !
    integer      :: ii, jj
    real(dp)     :: vtmp(3), cax(3), aax(3)
    integer      :: ltmp(3)
    type(sym_op) :: ltgop(48)   ! Little group operations
    integer      :: nltgop      ! Number of Little group operations
    integer      :: ltg2spg(48) ! Index of little group operations in spgop()
    !
    ! read structure ================================================
    !
    sat_kind  = 0
    swan_list = ''
    open(unit=100, file='job.in', status='old')
    read(100, SliceWannier)
    close(100)
    !
    ! count pseudo-atoms
    nsat = 0
    do ii=1, lsize
        if (sat_kind(ii)==0) exit
        nsat = nsat + 1
        sat(nsat)%kd = sat_kind(ii)
    enddo
    !
    ! sat positions and SSGs
    !
    if (sat_pos_conv(1:1)=='d') then
        do ii=1,nsat
            sat_pos(:,ii) = sat_pos(1,ii)*cell_a(:,1) + sat_pos(2,ii)*cell_a(:,2) &
                          + sat_pos(3,ii)*cell_a(:,3)
        enddo
    endif
    !
    do ii=1,nsat
        call TonA( sat_pos(:,ii), ltmp, vtmp, cell_a, cell_b)
        if (ltmp(1)/=0 .or. ltmp(2)/=0 .or. ltmp(3)/=0) then
            write(*,*) 'Error in get_obt(): atom is not in 0 0 0 lattice !!!'
            stop
        endif
        sat(ii)%pos = sat_pos(:,ii)
        !
        ! find the SSG of sat(ii)%pos
        call find_ltgrp_r(sat(ii)%pos, nspgop, spgop, cell_a, cell_b, nltgop, ltgop, ltg2spg)
        call pgrp_name(nltgop, ltgop, sat(ii)%SSG, cax, aax)
        call euler_angle(cax, aax, sat(ii)%alph, sat(ii)%beta, sat(ii)%gamm)
        !
    enddo
    !
endsubroutine get_swan
endmodule
