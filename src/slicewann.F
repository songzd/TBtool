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
    integer     :: RR, aa, ss, xi, LL, ii, ii_, jj, nn
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
                do jj=sat(ss)%idx1(xi), sat(ss)%idx2(xi)
                    call random_number(rand)
                    rand = rand*10.0_dp*exp(-bond/dist_cut)
                    ovlp(jj, ii_, RR) = rand(1) + rand(2)*c_im
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
               gg, RR, RRp, La, Lb, ii_, jj_, iip_, jjp_
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
    enddo
    enddo
    enddo
endsubroutine symm_ovlp

endmodule slicewann
