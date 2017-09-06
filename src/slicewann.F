module slicewann
!
use typcnst,    only: DP, pat_type, irepbs_type, lsize, sym_op, irep_type
use pntgrp,     only: TonA, find_ltgrp_r, pgrp_name, euler_angle
use reprsnt,    only: init_pntgrp_bst, set_pntgrp_irep
!
use tghtbd,     only: cell_a, cell_b, nspin
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
integer,parameter               :: maxnswan=100     ! max number of nswan
integer                         :: nswan            ! number of swannier functions
integer                         :: nsat             ! number of pseudo-atoms
type(pat_type)                  :: sat(lsize)       ! pseudo-atoms
type(irepbs_type)               :: swan(maxnswan)   ! swannier functions
!
contains
!
subroutine get_swan()
    !
    ! read-in variables
    character(len=30)   :: k_conv, sat_pos_conv
    integer             :: sat_kind(lsize)
    real(dp)            :: sat_pos(3,lsize)
    character(len=30)   :: swan_list(3,lsize)
    !
    namelist /SliceWannier/ k_conv, kz_sls, bstart, bend, sat_kind, &
                            sat_pos_conv, sat_pos, swan_list
    !
    integer             :: ii, jj, kk, nn, xi
    real(dp)            :: vtmp(3), cax(3), aax(3)
    integer             :: ltmp(3)
    type(sym_op)        :: ltgop(48)   ! Little group operations
    integer             :: nltgop      ! Number of Little group operations
    integer             :: ltg2spg(48) ! Index of little group operations in spgop()
    character(len=30)   :: kdtmp
    type(irep_type)     :: irep(20)
    integer             :: nirep
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
    ! sat positions
    !
    if (sat_pos_conv(1:1)=='d') then
        do ii=1,nsat
            sat_pos(:,ii) = sat_pos(1,ii)*cell_a(:,1) + sat_pos(2,ii)*cell_a(:,2) &
                          + sat_pos(3,ii)*cell_a(:,3)
        enddo
    endif
    !
    ! sat SSGs and slice-wannier projectors
    !
    nswan = 0;
    call init_pntgrp_bst();
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
        ! set SSG-ireps
        call set_pntgrp_irep(nltgop,ltgop,sat(ii)%SSG,nspin,cax,aax,nirep,irep)
        !
        ! slice-wannier projector
        !
        sat(ii)%nirep=0
        do jj=1,lsize
            if ( trim(adjustl(swan_list(1,jj)))=='' ) exit
            write(kdtmp,*) sat(ii)%kd
            if ( trim(adjustl(kdtmp)) == trim(adjustl(swan_list(1,jj))) ) then
                !
                if ( trim(adjustl(swan_list(2,jj))) /= trim(adjustl(sat(ii)%SSG)) ) then
                    write(*,*) 'Error in get_swan(): SSG does mismatch !!!'
                    stop
                endif
                !
                ! set the projector
                do xi=1,nirep
                if ( trim(adjustl(irep(xi)%name)) == trim(adjustl(swan_list(3,jj))) ) then
                    !
                    sat(ii)%nirep = sat(ii)%nirep + 1    
                    nn = sat(ii)%nirep
                    sat(ii)%idx1(nn) = nswan - irep(xi)%D + 1
                    sat(ii)%idx2(nn) = nswan
                    sat(ii)%irep(nn) = irep(xi)%name
                    !
                    nswan = nswan + irep(xi)%D
                    do kk=1, irep(xi)%D
                        swan(nswan-irep(xi)%D+kk) = irep(xi)%bs(kk)
                    enddo
                    !
                    exit
                endif
                enddo
                !
            endif
        enddo
    enddo
    !
    ! output projector information
    !
    write(*,'(1X, "Slice wannier projector:")')
    write(*,'(6X, "Position (in cart)", 11X, "SSG", 4X, "Euler angle", 17X, "Projector")')
    do ii=1,nsat
    do xi=1,sat(ii)%nirep
        write(*,'(6X, "(", 3F8.4, ")   ",  A5, "  (", 3F8.4, ")  ", A5 )')  &
        sat(ii)%pos, sat(ii)%SSG, sat(ii)%alph, sat(ii)%beta, sat(ii)%gamm, sat(ii)%irep(xi)
    enddo
    enddo
endsubroutine get_swan
endmodule
