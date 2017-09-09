!**** slice *********************************************************
!
! This module offers slice-wannier projection information and slice-
!  wannierise controling parameters
!
!********************************************************************
module slice
!
use typcnst,    only: dp, sat_type=>wat_type, irepbs_type, lsize, sym_op, irep_type, &
                      maxdim, Lmax, maxbcp, c_pi, LTmax
use my_math,    only: vect_len, unimat, dot_prod, inner_prod, angle
use pntgrp,     only: TonA, find_ltgrp_r, pgrp_name, euler_angle
use reprsnt,    only: init_pntgrp_bst, set_pntgrp_irep, dmat_cmplx
use band,       only: bndspc_ini, bndspc_rls, bnd_slv
use wannier,    only: random_ovlp, symm_ovlp, proj_WF, proj_WF_ini, proj_WF_rls
!
use tghtbd,     only: cell_a, cell_b, nat=>natom, at=>atom, nobt, nspin, nobt2, &
                      nlatt=>h0_lat_num, latt=>h0_lat, LTidx=>h0_LTidx, h0_dat
use symm,       only: nspgop, spgop, symtol, time_rev, atmap, atmapR, rep_at, reps_at
!
implicit none
!
public
!
! kz and bands
real(dp)                        :: kz_sls           ! kz of the 2D slice, in direct 
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
! slice group: the group leaves kz_sls invariant
integer                         :: nslgop            ! number of slice group operations
type(sym_op)                    :: slgop(48)         ! slice group operations
character(len=30)               :: slgnm             ! slice group name
integer                         :: slg2spg(48)       ! index in space group spgop
logical                         :: TR_slice          ! whether this slice is TR-invariant
!
! rotation matrix
integer, allocatable        :: satmap(:,:)          ! satoms mapping, gg*aa is the atmap(aa,gg)
integer, allocatable        :: satmapR(:,:,:)       ! satom in the atmapR(:,aa,gg) lattice
complex(dp), allocatable    :: rep_sat(:,:,:,:,:)   ! rotation matrix
complex(dp), allocatable    :: Ut_sat(:,:,:,:)      ! Time reversal matrix
!
! kx ky mesh
integer                     :: nk1, nk2        ! number of kx/ky grid
real(dp),allocatable        :: kmesh(:,:,:)         ! k points mesh
!
contains
!
!**** read_swan() ****************************************************
! read the slice-wannier projectors
!********************************************************************
subroutine read_swan()
    !
    ! read-in variables
    character(len=30)   :: sat_pos_conv
    integer             :: sat_kind(lsize)
    real(dp)            :: sat_pos(3,lsize)
    character(len=30)   :: sproj_list(3,lsize)
    !
    namelist /SliceWannier/ kz_sls, bstart, bend, sat_kind, &
                            sat_pos_conv, sat_pos, sproj_list, dist_cut, &
                            nk1, nk2
    !
    integer             :: ii, jj, kk, nn, xi, gg
    real(dp)            :: vtmp(3), cax(3), aax(3)
    integer             :: ltmp(3)
    type(sym_op)        :: ltgop(48)   ! Little group operations
    integer             :: nltgop      ! Number of Little group operations
    integer             :: ltg2spg(48) ! Index of little group operations in spgop()
    character(len=30)   :: kdtmp
    type(irep_type)     :: irep(20)
    integer             :: nirep
    logical             :: findproj, kz0p
    real(dp)            :: caxis(3)
    !
    ! There should be a3 \perp a1,a2 ================================
    !
    if ( abs(dot_prod(cell_a(:,3),cell_a(:,1)))>1.0d-4 .or. &
         abs(dot_prod(cell_a(:,3),cell_a(:,2)))>1.0d-4 ) then
        write(*,*) 'Error in read_swan(): the lattice is not supported !!!'
    endif
    caxis = cell_a(:,3)/vect_len(cell_a(:,3))
    !
    ! read  =========================================================
    !
    sat_kind  = 0
    sproj_list = ''
    dist_cut = 3
    open(unit=100, file='job.in', status='old')
    read(100, SliceWannier)
    close(100)
    !
    ! slice group ===================================================
    !
    kz0p = mod(kz_sls,0.5)<1.0d-8
    TR_slice = kz0p .and. time_rev
    !
    nslgop=0
    slg2spg=0
    do gg=1,nspgop
        ! Identity & inversion
        if ( abs(spgop(gg)%thet) < 1.0d-4 ) then
            if ( spgop(gg)%det == 1 .or. ( spgop(gg)%det ==-1 .and. kz0p ) ) then
                nslgop = nslgop + 1
                slgop(nslgop) = spgop(gg)
                slg2spg(nslgop) = gg
            endif
        ! Axis parallel with 001
        elseif ( abs( abs( dot_prod(spgop(gg)%axis, caxis) ) - 1.0_dp )<1.0d-4 ) then
            if ( spgop(gg)%det == 1 .or. ( spgop(gg)%det ==-1 .and. kz0p ) ) then
                nslgop = nslgop + 1
                slgop(nslgop) = spgop(gg)
                slg2spg(nslgop) = gg
            endif
        endif
        ! Axis perpendicular
        if ( abs(dot_prod( spgop(gg)%axis, caxis )) < 1.0d-4 .and. &
             abs(abs(spgop(gg)%thet) - c_pi) < 1.0d-4 ) then
            if( spgop(gg)%det == -1 .or. (spgop(gg)%det == 1 .and. kz0p )) then
                nslgop = nslgop + 1
                slgop(nslgop) = spgop(gg)
                slg2spg(nslgop) = gg
            endif
        endif
        !
    enddo
    !
    call pgrp_name(nslgop, slgop, slgnm, cax, aax)
    !
    write(*,'(1X, "Slice wannier:")')
    write(*,'(5X, "kz = ", F10.5, " (in cell_b(3))")') kz_sls
    write(*,'(5X, "Slice point group:", A9)') trim(slgnm)
    write(*,'(5X, " iop", " det", "    alph", 2X, "axis (in cart)", &
              15X, "tau (in cart)", 15X, "in SPG")')
    do ii=1, nslgop
58      format(5X,I4, I4, F8.1, '  (',3F8.4, ')','   (',3F8.4,')', I8)
        write(*, 58) ii, slgop(ii)%det, slgop(ii)%thet/c_pi*180, &
                     slgop(ii)%axis, slgop(ii)%tau,slg2spg(ii) 
    enddo
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
        call find_ltgrp_r(sat(ii)%pos, nslgop, slgop, cell_a, cell_b, nltgop, ltgop, ltg2spg)
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
                    write(*,*) 'Error in read_swann(): SSG mismatch !!!'
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
                    !
                    ! Physical irrep (TR-invariant irrep), for example: E12, or 1E12+2E12
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
    ! k-mesh --------------------------------------------------------
    !
    allocate( kmesh(3,nk1,nk2) )
    do ii=1,nk1
    do jj=1,nk2
        kmesh(:,ii,jj) = (ii-1)*cell_b(:,1)/nk1 + (jj-1)*cell_b(:,2)/nk2 + kz_sls*cell_b(:,3)
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
    logical         :: found
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
    allocate(satmap(nsat,nslgop) )
    allocate(satmapR(3,nsat,nslgop))
    allocate(rep_sat(maxdim, maxdim, maxrep_sat, nsat, nslgop))
    !
    dmat_op(:,:,:) = 0.0_dp
    satmap(:,:)=0
    satmapR(:,:,:)=0
    rep_sat(:,:,:,:,:) = 0.0_dp
    !
    ! Slice group operation =========================================
    !
    do gg=1,nslgop
        !
        ! sat mapping -----------------------------------------------
        !
        do aa=1,nsat
            found = .false.
            call TonA( matmul(slgop(gg)%R,sat(aa)%pos) + slgop(gg)%tau, &
                       satmapR(:,aa,gg), aapos, cell_a, cell_b )
            do bb=1,nsat
                if ( vect_len( aapos-sat(bb)%pos )<symtol .and. sat(bb)%kd==sat(aa)%kd ) then
                    found = .true.
                    exit
                endif
            enddo
            !
            if (found) then
                satmap(aa,gg)=bb
            else
                write(*,*) 'Error in cstrct_rep_slice(): The slice-atoms does not satisfy &
                            the slice group'
                stop
            endif
            !
            !write(100,'(5X,7I5)') gg, aa, satmap(aa,gg), satmapR(:,aa,gg)
            !
        enddo
        !
        ! calculate D-matrix of gg : only the rotational part -------
        !
        ! the integer part
        do LL=0,Lmax
            call dmat_cmplx( slgop(gg)%alph, slgop(gg)%beta, slgop(gg)%gamm, &
                             dmat_op(-LL:LL,-LL:LL,LL), LL )
        enddo
        !
        ! the half-integer part
        do LL=-1,-Lmax,-2
            call dmat_cmplx( slgop(gg)%alph, slgop(gg)%beta, slgop(gg)%gamm, &
                             dmat_op(-abs(LL):abs(LL),-abs(LL):abs(LL),LL), LL )
        enddo
        !
        ! rep_sat
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
                if ( slgop(gg)%det==-1 .and. sproj(idx1)%P==-1) dmat = - dmat
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
    if (TR_slice) then
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
    endif
    !
    close(100)
endsubroutine cstrct_rep_slice

subroutine build_swan()
    !
    ! overlap
    complex(dp)         :: ovlp(nsproj,nobt2,nlatt)
    !
    ! transformation
    integer             :: atmap_(nat,nslgop), atmapR_(3,nat,nslgop) 
    real(dp)            :: rep_at_(1:2*Lmax+1, 1:2*Lmax+1, 0:Lmax, nat, nslgop)
    complex(dp)         :: reps_at_(nspin,nspin,nat,nslgop)
    !
    ! wave functions
    complex(dp)         :: wave(nobt2,nobt2), Uwan(nobt2,nsproj,nk1,nk2)
    real(dp)            :: eng(nobt2), Sdet(nk1,nk2)
    !
    ! wannier center
    real(dp)            :: WC1(nsproj,nk2), WC2(nsproj,nk1), WC1t(nsproj), WC2t(nsproj)
    complex(dp)         :: phase
    !
    integer             :: RR, aa, ii, jj, gg, nn
    !
    !if ( bend-bstart+1 /= nsproj ) then
    !    write(*,*) 'Error in build_swan(): number of bands and projectos mismatch !!!'
    !    stop
    !endif
    !
    ! set random ovlp ===============================================
    !
    call random_ovlp( nsat, sat, nsproj, nat, at, nobt, nobt2, &
                      nspin, nlatt, latt, cell_a, dist_cut, ovlp)
    !
    ! Get the transformation ========================================
    !
    atmap_  = 0
    atmapR_ = 0
    rep_at_ = 0.0_dp
    reps_at_= 0.0_dp
    do gg=1,nslgop
        do aa=1,nat
            atmap_(aa,gg)    =  atmap(aa,slg2spg(gg))
            atmapR_(:,aa,gg) = atmapR(:,aa,slg2spg(gg))
            rep_at_(:,:,:,aa,gg) = rep_at(:,:,:,aa,slg2spg(gg))
            reps_at_(:,:,aa,gg)  = reps_at(:,:,aa,slg2spg(gg))
        enddo
    enddo
    !
    ! Symmetrize the overlap ========================================
    !
    call symm_ovlp( nsat, sat, maxrep_sat, nsproj, satmap, satmapR, rep_sat, Ut_sat,  & 
                    nat,  at,  nobt, nobt2, nspin, atmap_, atmapR_, rep_at_, reps_at_, &
                    nlatt, latt, LTmax, LTidx, cell_a, cell_b, &
                    TR_slice, nslgop, slgop, ovlp )
    !
    ! Output the overlap ============================================
    !
    open(100,file='out.ovlp_rand',status='replace')
    do RR=1,nlatt
        do ii=1,nsproj
        do jj=1,nobt2
            if ( abs(ovlp(ii,jj,RR))>1.0d-7 ) then
                write(100,'(5I5, 2F20.10)')  latt(:,RR), ii, jj, ovlp(ii,jj,RR)
            endif
        enddo
        enddo
    enddo
    close(100)
    !
    write(*,'(6X, "Symmetrized random overlap matrix written to out.ovlp_rand ")')
    !
    ! Get the wannier function ======================================
    !
    call bndspc_ini(nobt2)
    call proj_WF_ini(nsproj)
    do ii=1,nk1
    do jj=1,nk2
        call bnd_slv( kmesh(:,ii,jj), wave, eng, nobt, nobt2, nspin, nlatt, &
                      h0_dat, latt, nat, at, cell_a )
        call proj_WF( kmesh(:,ii,jj), nobt2, wave, bstart, bend, &
                      nsproj, nlatt, latt, cell_a, ovlp, Uwan(:,:,ii,jj), Sdet(ii,jj) )
    enddo
    enddo
    call bndspc_rls()
    call proj_WF_rls()
    print*, minval( Sdet ), maxval(Sdet)
    !
    ! Wannier center ================================================
    !
    do nn=1,nsproj
        !
        ! WC1
        do ii=1,nk1
            phase=inner_prod( Uwan(:,nn,nk1,ii), Uwan(:,nn,1,ii), nobt2 )
            do jj=1, nk2-1
                phase = phase*inner_prod( Uwan(:,nn,jj,1), Uwan(:,nn,jj+1,1), nobt2 )
            enddo
            WC1(nn,ii) = acos(real(phase)/abs(phase))/(2*c_pi)
            write(40,*) WC1(nn,ii)
        enddo
        write(40,*)' '
    enddo
endsubroutine build_swan
!
endmodule slice
