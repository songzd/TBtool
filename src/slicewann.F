module slicewann
!
use typcnst,    only: DP, pat_type, wan_type, lsize
use pntgrp,     only: TonA
use tghtbd,     only: cell_a, cell_b
!
implicit none
!
public
!
! kz and bands
real(dp)                        :: kz_sls(3)        ! kz of the 2D slice
integer                         :: bstart, bend     ! bands to construct WFs
!
! pseudo-atoms and wannier orbitals
integer                         :: nwan             ! number of wannier functions
integer                         :: npat             ! number of pseudo-atoms
type(pat_type)                  :: pat(lsize)       ! pseudo-atoms
type(wan_type),allocatable      :: wan(:)           ! wannier functions
!
contains
!
subroutine get_obt()
    !
    ! read-in variables
    character(len=30)   :: k_conv, pat_pos_conv
    integer             :: pat_kind(lsize)
    real(dp)            :: pat_pos(3,lsize)
    character(len=30)   :: wan_list(3,lsize)
    !
    namelist /SliceWannier/ k_conv, kz_sls, bstart, bend, pat_kind, &
                            pat_pos
    !
    integer     :: ii, jj
    real(dp)    :: vtmp
    integer     :: ltmp
    !
    ! read structure ================================================
    !
    pat_kind  = 0
    wan_list = ''
    open(unit=100, file='job.in', status='old')
    read(100, Slicewannier)
    close(100)
    !
    ! count pseudo-atoms
    npat = 0
    do ii=1, lsize
        if (pat_kind(ii)==0) exit
        npat = npat + 1
        pat(npat)%kd = pat_kind(ii)
    enddo
    !
    ! pat positions
    if (pat_pos_conv(1:1)=='d') then
        do ii=1,npat
            pat_pos(:,ii) = pat_pos(1,ii)*cell_a(:,1) + pat_pos(2,ii)*cell_a(:,2) &
                          + pat_pos(3,ii)*cell_a(:,3)
        enddo
    endif
    !
    do ii=1,npat
        call TonA( pat_pos(:,ii), ltmp, vtmp, cell_a, cell_b)
        if (ltmp(1)/=0 .or. ltmp(2)/=0 .or. ltmp(3)/=0) then
            write(*,*) 'Error in get_obt(): atom is not in 0 0 0 lattice !!!'
            stop
        endif
        pat(ii)%pos = pat_pos(:,ii)
    enddo
endmodule
