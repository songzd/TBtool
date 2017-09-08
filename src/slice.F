!**** slice *********************************************************
!
! This module offers slice-wannier projection information and slice-
!  wannierise controling parameters
!
!********************************************************************
module slice
!
use typcnst,    only: DP, sat_type, irepbs_type, lsize, sym_op, irep_type, &
                      maxdim, Lmax, maxbcp, c_pi
use my_math,    only: vect_len, unimat
use pntgrp,     only: TonA, find_ltgrp_r, pgrp_name, euler_angle
use reprsnt,    only: init_pntgrp_bst, set_pntgrp_irep, dmat_cmplx
!
use tghtbd,     only: cell_a, cell_b, nspin
use symm,       only: nspgop, spgop, symtol, time_rev
!
implicit none
!
public
!
! kz and bands
real(dp)                        :: kz_sls(3)        ! kz of the 2D slice
integer                         :: bstart, bend     ! bands to construct WFs
!
! pseudo-atoms and slice-wannier projector
integer,parameter               :: maxnsproj=100     ! max number of nsproj
integer                         :: nsproj            ! number of sprojnier functions
integer                         :: maxrep_sat        ! number of sprojnier functions
integer                         :: nsat              ! number of pseudo-atoms
type(sat_type)                  :: sat(lsize)        ! pseudo-atoms
type(irepbs_type)               :: sproj(maxnsproj)  ! sprojnier functions
real(dp)                        :: dist_cut          ! distance cut off for overlap
!
! rotation matrix
integer, allocatable        :: satmap(:,:)          ! satoms mapping, gg*aa is the atmap(aa,gg)
integer, allocatable        :: satmapR(:,:,:)       ! satom in the atmapR(:,aa,gg) lattice
complex(dp), allocatable    :: rep_sat(:,:,:,:,:)   ! rotation matrix
complex(dp), allocatable    :: Ut_sat(:,:,:,:)      ! Time reversal matrix
!
!
contains
!
!**** read_swan() ****************************************************
! read the slice-wannier projectors
!********************************************************************
subroutine read_swan()
    !
    ! read-in variables
    character(len=30)   :: k_conv, sat_pos_conv
    integer             :: sat_kind(lsize)
    real(dp)            :: sat_pos(3,lsize)
    character(len=30)   :: sproj_list(3,lsize)
    !
    namelist /SliceWannier/ k_conv, kz_sls, bstart, bend, sat_kind, &
                            sat_pos_conv, sat_pos, sproj_list, dist_cut
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
    logical             :: findproj
    !
    ! read structure ================================================
    !
    sat_kind  = 0
    sproj_list = ''
    dist_cut = 3
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
    nsproj = 0;
    call init_pntgrp_bst();
    do ii=1,nsat
        call TonA( sat_pos(:,ii), ltmp, vtmp, cell_a, cell_b)
        if (ltmp(1)/=0 .or. ltmp(2)/=0 .or. ltmp(3)/=0) then
            write(*,*) 'Error in read_swan(): sat is not in 0 0 0 lattice !!!'
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
            if ( trim(adjustl(sproj_list(1,jj)))=='' ) exit
            write(kdtmp,*) sat(ii)%kd
            if ( trim(adjustl(kdtmp)) == trim(adjustl(sproj_list(1,jj))) ) then
                !
                if ( trim(adjustl(sproj_list(2,jj))) /= trim(adjustl(sat(ii)%SSG)) ) then
                    write(*,*) 'Error in read_swann(): SSG does mismatch !!!'
                    stop
                endif
                !
                ! set the projector
                !
                findproj = .false.
                !
                nn = sat(ii)%nirep + 1
                sat(ii)%deg(nn)  = 0
                sat(ii)%idx(nn)  = nsproj + 1
                sat(ii)%irep(nn) = ''
                !
                do xi=1,nirep
                    if ( index( sproj_list(3,jj), trim(adjustl(irep(xi)%name)) )>0 ) then
                        !
                        findproj = .true.
                        !
                        sat(ii)%deg(nn) = sat(ii)%deg(nn) + irep(xi)%D
                        nsproj = nsproj + irep(xi)%D
                        if ( trim(adjustl(sat(ii)%irep(nn)))=='' ) then
                            sat(ii)%irep(nn) = trim(adjustl(irep(xi)%name))
                        else
                            sat(ii)%irep(nn) = trim(adjustl(sat(ii)%irep(nn))) // '+' // trim(adjustl(irep(xi)%name))
                        endif
                        !
                        do kk=1, irep(xi)%D
                            sproj(nsproj - irep(xi)%D + kk) = irep(xi)%bs(kk)
                        enddo
                        !
                    endif
                enddo
                !
                if ( findproj ) then
                    sat(ii)%nirep = nn
                else
                    write(*,*) 'Error in get_sproj(): sproj does mismatch !!!' 
                    stop
                endif
                !
            endif
        enddo
    enddo
    !
    maxrep_sat=0
    do ii=1,nsat
        if( sat(ii)%nirep>maxrep_sat ) maxrep_sat = sat(ii)%nirep
    enddo
    !
    ! output projector information
    !
    write(*,'(1X, "Slice wannier projector:")')
    write(*,'(6X, "Position (in cart)", 11X, "SSG", 4X, "Euler angle", 17X, "Projector")')
    do ii=1,nsat
    do xi=1,sat(ii)%nirep
        write(*,'(6X, "(", 3F8.4, ")   ",  A5, "  (", 3F8.4, ")  ", A12, I4, " - ", I4)')  &
        sat(ii)%pos, sat(ii)%SSG, sat(ii)%alph, sat(ii)%beta, sat(ii)%gamm, sat(ii)%irep(xi), &
        sat(ii)%idx(xi), sat(ii)%idx(xi) + sat(ii)%deg(xi) - 1
    enddo
    enddo
    !
endsubroutine read_swan

!**** cstrct_rep_slice() ********************************************
! Construct the rotation matrix of slice projectors
!********************************************************************
subroutine cstrct_rep_slice()
    integer         :: aa, bb, aap, ii, jj, LL, gg, xi, idx1, idx1p, Jt, &
                       deg, nn, mm, ii_, jj_
    real(dp)        :: aapos(3)    ! used for atom mapping
    !
    ! D matrix of slice-atoms
    complex(dp)     :: dmat_at(-Lmax:Lmax,-Lmax:Lmax,-Lmax:Lmax,nsat)
    !
    ! D matrix of operation
    complex(dp)     :: dmat_op(-Lmax:Lmax,-Lmax:Lmax,-Lmax:Lmax)
    complex(dp)     :: dmat(-Lmax:Lmax,-Lmax:Lmax)
    !
    ! Local D matrices ==============================================
    !
    dmat_at(:,:,:,:) = 0.0_dp
    do aa=1,nsat
        !
        ! the integer part
        do LL=0, Lmax
            call dmat_cmplx( sat(aa)%alph, sat(aa)%beta, sat(aa)%gamm, &
                             dmat_at(-LL:LL,-LL:LL, LL, aa), LL)
        enddo
        !
        ! the half-integer part
        do LL=-1,-Lmax,-2
            call dmat_cmplx( sat(aa)%alph, sat(aa)%beta, sat(aa)%gamm, &
                             dmat_at(-abs(LL):abs(LL),-abs(LL):abs(LL), LL, aa), LL)
        enddo
    enddo
    !
    ! sats mapping and rep_sat ======================================
    !
    open(unit=100,file='out.rep_sat',status='replace')
    !
    allocate(satmap(nsat,nspgop) )
    allocate(satmapR(3,nsat,nspgop))
    allocate(rep_sat(maxdim, maxdim, maxrep_sat, nsat, nspgop))
    !
    dmat_op(:,:,:) = 0.0_dp
    satmap(:,:)=0
    satmapR(:,:,:)=0
    rep_sat(:,:,:,:,:) = 0.0_dp
    !
    ! Space group operation =========================================
    !
    do gg=1,nspgop
        !
        ! sat mapping -----------------------------------------------
        !
        do aa=1,nsat
            call TonA( matmul(spgop(gg)%R,sat(aa)%pos) + spgop(gg)%tau, &
                       satmapR(:,aa,gg), aapos, cell_a, cell_b )
            do bb=1,nsat
                if ( vect_len( aapos-sat(bb)%pos )<symtol .and. &
                     sat(bb)%kd==sat(aa)%kd ) exit
            enddo
            satmap(aa,gg)=bb
            !
            !write(100,'(5X,7I5)') gg, aa, satmap(aa,gg), satmapR(:,aa,gg)
            !
        enddo
        !
        ! calculate D-matrix of gg : only the rotational part -------
        !
        ! the integer part
        do LL=0,Lmax
            call dmat_cmplx( spgop(gg)%alph, spgop(gg)%beta, spgop(gg)%gamm, &
                             dmat_op(-LL:LL,-LL:LL,LL), LL )
        enddo
        !
        ! the half-integer part
        do LL=-1,-Lmax,-2
            call dmat_cmplx( spgop(gg)%alph, spgop(gg)%beta, spgop(gg)%gamm, &
                             dmat_op(-abs(LL):abs(LL),-abs(LL):abs(LL),LL), LL )
        enddo
        !
        ! rep_at
        do aa=1, nsat
            aap = satmap(aa,gg)
            do xi=1, sat(aa)%nirep
                !
                ! get projector index and total angular momentum
                idx1  = sat( aa)%idx(xi)
                idx1p = sat(aap)%idx(xi)
                deg   = sat( aa)%deg(xi)
                Jt = sproj(idx1)%J
                LL = abs(Jt)
                ! 
                dmat(-LL:LL,-LL:LL) = &
                    matmul(dmat_op(-LL:LL, -LL:LL, Jt), dmat_at(-LL:LL, -LL:LL, Jt, aa))
                dmat(-LL:LL,-LL:LL) = & 
                    matmul(conjg(transpose(dmat_at(-LL:LL, -LL:LL, Jt, aap))), dmat(-LL:LL, -LL:LL))
                if ( spgop(gg)%det==-1 .and. sproj(idx1)%P==-1) dmat = - dmat
                !
                do ii=1,deg
                do jj=1,deg
                    do nn=1,maxbcp
                    do mm=1,maxbcp
                        rep_sat(ii,jj,xi,aa,gg) = rep_sat(ii,jj,xi,aa,gg) + &
                            dmat(sproj(idx1p+ii-1)%M(nn), sproj(idx1+jj-1)%M(mm)) &
                           *conjg( sproj(idx1p+ii-1)%C(nn) ) * sproj(idx1+jj-1)%C(mm)
                    enddo
                    enddo
                enddo
                enddo
                !
                if ( .not. unimat( rep_sat(1:deg,1:deg,xi,aa,gg), deg ) ) then
                    write(*,*) 'Error in cstrct_rep(): orbitals are not compatible with symmetry !!!'
                    stop
                endif
                !
                ! write to out.rep_sat 
                write(100, '("iop=", I4, "  aa=", I4, "  aap=", I4, "  Rp=",3I4, "  ********************")') &
                    gg, aa, aap, satmapR(:,aa,gg)
                write(100, '(4X, "Projector on aa :  ", A12, 2X, I4, " - ",I4 )' ) sat( aa)%irep(xi),  idx1, idx1+deg-1
                write(100, '(4X, "Projector on aap:  ", A12, 2X, I4, " - ",I4 )' ) sat(aap)%irep(xi), idx1p, idx1p+deg-1 
                !
                do ii=1,deg
                    write(100,'(4X,  4(2F10.5, 5X) )') rep_sat(ii,1:deg,xi,aa,gg)
                enddo
            enddo
        enddo
        !
    enddo
    !
    ! Time reversal =================================================
    !
    allocate(Ut_sat(maxdim,maxdim,maxrep_sat,nsat))
    Ut_sat = 0.0_dp
    !
    do aa=1,nsat
        do xi=1,sat(aa)%nirep
            !
            idx1 = sat(aa)%idx(xi)
            deg  = sat(aa)%deg(xi)
            !
            Jt = sproj(idx1)%J
            LL = abs(Jt)
            dmat = 0.0_dp
            call dmat_cmplx( 0.0_dp, c_pi, 0.0_dp, dmat(-LL:LL,-LL:LL), Jt )
            ! 
            do ii=1,deg
            do jj=1,deg
                ii_ = ii - 1 + idx1
                jj_ = jj - 1 + idx1
                do nn=1,maxbcp
                do mm=1,maxbcp
                    Ut_sat(ii,jj,xi,aa) = Ut_sat(ii,jj,xi,aa)  &
                    + conjg(sproj(ii_)%C(nn))*dmat(sproj(ii_)%M(nn), sproj(jj_)%M(mm) )*conjg(sproj(jj_)%C(mm))
                enddo
                enddo
                !
            enddo
            enddo
            !
            ! Check Ut
            !
            if ( .not. unimat( Ut_sat(1:deg,1:deg,xi,aa), deg ) ) then
                write(*,*) 'Error in cstrct_rep_slice(): Ut is not unitary !!!'
                stop
            endif
            !
            ! write to out.rep_sat 
            write(100, '("Time reversal on at aa=", I4, "  *************************************")') aa
            write(100, '(4X, "Projector on aa :  ", A12, 2X, I4, " - ",I4 )' ) sat( aa)%irep(xi),  idx1, idx1+deg-1
            !
            do ii=1,deg
                write(100,'(4X,  4(2F10.5, 5X) )') Ut_sat(ii,1:deg,xi,aa)
            enddo
            !
        enddo
    enddo
    !
    close(100)
endsubroutine cstrct_rep_slice

endmodule slice
