!**** slicewann *****************************************************
! 
! This module offers the subroutines to construct slice-wannier functions
!
!********************************************************************
module slicewann
!
use typcnst,    only: dp, sat_type, irepbs_type, maxdim, maxrep_sat, Lmax, sym_op, &
                      atom_type, c_im
use my_math,    only: vect_len, init_random_seed
use pntgrp,     only: TonA
!
implicit none
public
!
contains
!
!**** random_ovlp() *************************************************
!
! Generate a random overlap matrix: 
!   O(ss,ii;aa,jj;RR) = < 0,ss,ii | RR,aa,jj >
!
!********************************************************************
subroutine random_ovlp( nsat, sat, nsproj, nat, at, nobt, nobt2, nspin, nlatt, latt, & 
                        cell_a, cell_b, dist_cut, ovlp)
    integer,intent(in)              :: nsat, nsproj, nat, nobt, nobt2, nspin, nlatt
    type(sat_type),intent(in)       :: sat(nsat)                    ! slice-atoms
    type(atom_type),intent(in)      :: at(nat)                      ! atoms
    integer,intent(in)              :: latt(3,nlatt)                ! lattice
    real(dp),intent(in)             :: cell_a(3,3), cell_b(3,3)     ! primitive cell in real/momentum space
    real(dp),intent(in)             :: dist_cut                     ! distance cut off for overlap matrix
    complex(dp),intent(out)         :: ovlp(nsproj, nobt2, nlatt)   ! the overlap matrix
    !
    !
    integer     :: RR, aa, ss, xi, LL, ii, ii_, jj, jj_, nn
    real(dp)    :: bond, rand(2)
    !
    if( nobt2/=nobt*nspin ) then
        write(*,*) 'Error in random_ovlp(): wrong nobt2 !!!'
        stop
    endif
    !
    call init_random_seed()
    !
    ovlp = 0.0_dp
    !
    do RR=1,nlatt
        do aa=1,nat
        do ss=1,nsat
            !
            bond = vect_len(  cell_a(:,1)*latt(1,RR) + cell_a(:,2)*latt(2,RR) &
                             +cell_a(:,3)*latt(3,RR) + at(aa)%pos - sat(ss)%pos )
            if ( bond > dist_cut ) exit
            !
            do LL=0,Lmax
            do ii=1,at(aa)%nobtL(LL)
            do nn=1,nspin
                ii_ = at(aa)%idx(ii,LL) + (nn-1)*nobt
                !
                do xi=1,sat(ss)%nirep
                do jj=1,sat(ss)%deg(xi)
                    jj_ = jj - 1 + sat(ss)%idx(xi)
                    call random_number(rand)
                    rand = rand*10.0_dp*exp(-bond/dist_cut)
                    ovlp(jj_, ii_, RR) = rand(1) + rand(2)*c_im
                enddo
                enddo
                !
            enddo
            enddo
            enddo
            !
        enddo
        enddo
    enddo
    !
endsubroutine random_ovlp

subroutine symm_ovlp( nsat, sat, nsproj, satmap, satmapR, rep_sat,  & 
                      nat,  at,  nobt, nobt2, nspin, atmap,  atmapR, rep_at, reps_at, &
                      nlatt, latt, Rmax, LTidx, cell_a, cell_b, &
                      time_rev, nop, op, &
                      ovlp )
    integer,intent(in)              :: nsat, nsproj, nat, nobt, nobt2, nspin, nlatt, Rmax, nop
    type(sat_type),intent(in)       :: sat(nsat)                ! slice-atoms
    integer,intent(in)              :: satmap(nsat,nop)         ! slice-atom mapping
    integer,intent(in)              :: satmapR(3,nsat,nop)      
    complex(dp),intent(in)          :: rep_sat(maxdim,maxdim,maxrep_sat,nsat,nop)   ! slice rotation matrix
    type(atom_type),intent(in)      :: at(nat)                  ! atoms
    integer,intent(in)              :: atmap(nsat,nop)          ! atom mapping
    integer,intent(in)              :: atmapR(3,nsat,nop)
    real(dp),intent(in)             :: rep_at(1:2*Lmax+1,1:2*Lmax+1,0:Lmax,nat,nop)
    complex(dp),intent(in)          :: reps_at(nspin,nspin,nat,nop)
    integer,intent(in)              :: latt(3,nlatt)
    integer,intent(in)              :: LTidx(-Rmax:Rmax, -Rmax:Rmax, -Rmax:Rmax)
    real(dp),intent(in)             :: cell_a(3,3), cell_b(3,3)
    logical,intent(in)              :: time_rev
    type(sym_op)                    :: op(nop)
    complex(dp),intent(inout)       :: ovlp(nobt2,nsproj,nlatt)
    !
    integer :: aa, ss, aap, ssp, ii, jj, iip, jjp, nn, mm, nnp, mmp, &
               gg, RR, RRp, LL, ii_, jj_, iip_, jjp_, xi
    real(kind=dp), dimension(3) :: aapos    ! used for atom mapping
    !
    ! the symmetrized hr
    complex(dp)         :: ovlp_sym(nsproj, nobt2, nlatt)
    real(dp)            :: Ut(nspin,nspin)
    logical             :: badlat(nsproj, nobt2, nlatt)
    !
    ! g0R  
    real(kind=dp), dimension(3) :: g0R_
    integer, dimension(3) :: g0R
    !
    if( nobt2/=nobt*nspin ) then
        write(*,*) 'Error in impose_sym(): wrong nobt2 !!!'
        stop
    endif
    !
    ! symmetrize the lattice ========================================
    !
    badlat = .false.
    do RR=1,nlatt       ! RR loop
    do aa=1,nat         ! atom loop
    do ss=1,nsat        ! slice-atom loop
        do gg=1,nop     ! operation loop
            g0R_ = latt(1,RR)*cell_a(:,1) + latt(2,RR)*cell_a(:,2) + &
                   latt(3,RR)*cell_a(:,3)
            g0R_ = matmul( op(gg)%R, g0R_ )
            call TonA( g0R_, g0R, aapos, cell_a, cell_b )
            !
            ! find RRp
            RRp = LTidx( g0R(1)+atmapR(1,aa,gg)-satmapR(1,ss,gg), &
                         g0R(2)+atmapR(2,aa,gg)-satmapR(2,ss,gg), &
                         g0R(3)+atmapR(3,aa,gg)-satmapR(3,ss,gg) )
            !
            if (RRp==0) badlat(ss,aa,RR) = .true.
        enddo
        !
        ! TR
        if (time_rev) then
            g0R = - latt(:,RR)
            RRp = LTidx(g0R(1), g0R(2), g0R(3))
            if (RRp==0) badlat(:,:,RR) = .true.
        endif
    enddo
    enddo
    enddo
    !
    ! Point group ===================================================
    !
    ovlp_sym = 0.0_dp
    do RR=1,nlatt
    do aa=1,nat
    do ss=1,nsat
        do gg=1,nop
            !
            g0R_ = latt(1,RR)*cell_a(:,1) + latt(2,RR)*cell_a(:,2) + &
                   latt(3,RR)*cell_a(:,3)
            g0R_ = matmul( op(gg)%R, g0R_ )
            call TonA( g0R_, g0R, aapos, cell_a, cell_b )
            !
            ! find RRp
            RRp = LTidx( g0R(1)+atmapR(1,aa,gg)-satmapR(1,ss,gg), &
                         g0R(2)+atmapR(2,aa,gg)-satmapR(2,ss,gg), &
                         g0R(3)+atmapR(3,aa,gg)-satmapR(3,ss,gg) )
            !
            ! matrix element
            aap = atmap(aa,gg)
            ssp = satmap(ss,gg)
            !
            ! remove bad lattice
            if (badlat(ss,aa,RR )) cycle   ! <=> RRp=0
            if (badlat(ssp,aap,RRp)) cycle
            !
            do LL=0,Lmax
            do jj=1,at(aa)%nobtL(LL)
            do nn=1,nspin
                jj_ = at(aa)%idx(jj,LL) + (nn-1)*nobt
            do xi=1,sat(ss)%nirep
            do ii=1,sat(ss)%deg(xi)
                ii_ = ii + sat(ss)%idx(xi) - 1
                !
                do jjp=1,at(aap)%nobtL(LL)
                do nnp=1,nspin
                    jjp_ = at(aap)%idx(jjp,LL) + (nnp-1)*nobt
                do iip=1,sat(ssp)%deg(xi)
                    iip_ = iip + sat(ssp)%idx(xi) - 1
                    !
                    ovlp_sym(ii_,jj_,RR) = ovlp_sym(ii_,jj_,RR)  &
                        + ovlp(iip_,jjp_,RRp)*rep_at(jjp,jj,LL,aa,gg)*reps_at(nnp,nn,aa,gg) &
                        * conjg( rep_sat(iip,ii,xi,ss,gg) )
                    !
                enddo
                enddo
                enddo
                !
            enddo
            enddo
            !
            enddo
            enddo
            enddo
            !
        enddo
    enddo
    enddo
    enddo
    !
    ovlp = ovlp_sym/nop
    !
    ! Time Reversal =================================================
    !
    if (time_rev) then
        !
        ovlp_sym = 0.0_dp
        !
        if (nspin==2) then
            Ut(1,:) = (/ 0.0_dp,  -1.0_dp /)
            Ut(2,:) = (/ 1.0_dp,   0.0_dp /)
        else
            Ut(1,1) = 1.0_dp
        endif
        !
        do RR=1,nlatt
            !
        enddo
        !
    endif

endsubroutine symm_ovlp

endmodule slicewann
