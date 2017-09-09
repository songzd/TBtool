!**** projwann *****************************************************
! 
! This module offers the subroutines to construct wannier functions
!   by projection method
!
!********************************************************************
module wannier
!
use typcnst,    only: dp, wat_type, irepbs_type, maxdim, Lmax, sym_op, &
                      atom_type, c_im
use my_math,    only: vect_len, dot_prod, init_random_seed
use pntgrp,     only: TonA
!
implicit none
public
!
! for mkl
integer :: lwork
complex(dp), allocatable,  dimension(:) :: work
real(dp), allocatable, dimension(:) :: rwork
!
contains
!
!**** random_ovlp() *************************************************
!
! Generate a random overlap matrix: 
!   O(ss,ii;aa,jj;RR) = < 0,ss,ii | RR,aa,jj >
!
!********************************************************************
subroutine random_ovlp( nwat, wat, nsproj, nat, at, nobt, nobt2, nspin, nlatt, latt, & 
                        cell_a, dist_cut, ovlp)
    integer,intent(in)              :: nwat, nsproj, nat, nobt, nobt2, nspin, nlatt
    type(wat_type),intent(in)       :: wat(nwat)                    ! slice-atoms
    type(atom_type),intent(in)      :: at(nat)                      ! atoms
    integer,intent(in)              :: latt(3,nlatt)                ! lattice
    real(dp),intent(in)             :: cell_a(3,3)                  ! primitive cell
    real(dp),intent(in)             :: dist_cut                     ! distance cut off for overlap matrix
    complex(dp),intent(out)         :: ovlp(nsproj, nobt2, nlatt)   ! the overlap matrix
    !
    !
    integer     :: RR, aa, ss, xi, LL, ii, ii_, jj, jj_, nn
    real(dp)    :: bond, norm
    ! 
    ! random number
    real                 :: rand(2)
    !
    if( nobt2/=nobt*nspin ) then
        write(*,*) 'Error in random_ovlp(): wrong nobt2 !!!'
        stop
    endif
    !
    call init_random_seed()
    !
    ovlp= 0.0_dp
    !
    do ss=1,nwat
    do xi=1,wat(ss)%nirep
    do jj=1,wat(ss)%deg(xi)
        !
        jj_ = jj - 1 + wat(ss)%idx(xi)
        norm = 0.0_dp
        !
        do RR=1,nlatt
        do aa=1,nat
            !
            bond = vect_len(  cell_a(:,1)*latt(1,RR) + cell_a(:,2)*latt(2,RR) &
                             +cell_a(:,3)*latt(3,RR) + at(aa)%pos - wat(ss)%pos )
            if ( bond > dist_cut ) cycle
            !
            do LL=0,Lmax
            do ii=1,at(aa)%nobtL(LL)
            do nn=1,nspin
                ii_ = at(aa)%idx(ii,LL) + (nn-1)*nobt
                !
                call random_number(rand)
                rand = rand*exp(-bond/dist_cut)
                ovlp(jj_, ii_, RR) = rand(1) + rand(2)*c_im
                !
                norm = norm + conjg(ovlp(jj_, ii_, RR))*ovlp(jj_, ii_, RR)
            enddo
            enddo
            enddo
            !
        enddo
        enddo
        !
        ! normalize
        ovlp(jj_,:,:) = ovlp(jj_,:,:)/sqrt(norm)
        !
    enddo
    enddo
    enddo
    !
endsubroutine random_ovlp

!**** symm_ovlp() ***************************************************
!
! Symmetrize the overlap matrix 
!   O(ss,ii;aa,jj;RR) = < 0,ss,ii | RR,aa,jj > = < g * 0,ss,ii | g * RR,aa,jj > 
!
!********************************************************************
subroutine symm_ovlp( nwat, wat, maxrep_wat, nsproj, watmap, watmapR, rep_wat, Ut_wat, & 
                      nat,  at,  nobt, nobt2, nspin, atmap,  atmapR, rep_at, reps_at, &
                      nlatt, latt, Rmax, LTidx, cell_a, cell_b, &
                      time_rev, nop, op, ovlp )
    integer,intent(in)         :: nwat, nsproj, nat, nobt, nobt2, nspin, nlatt, Rmax, nop, maxrep_wat
    type(wat_type),intent(in)  :: wat(nwat)                 ! slice-atom
    integer,intent(in)         :: watmap(nwat,nop)          ! slice-atom mapping
    integer,intent(in)         :: watmapR(3,nwat,nop)       ! slice-atom mapping lattice
    complex(dp),intent(in)     :: rep_wat(maxdim,maxdim,maxrep_wat,nwat,nop)    ! slice rotation matrix
    complex(dp),intent(in)     :: Ut_wat(maxdim,maxdim,maxrep_wat,nwat)         ! TR matrix
    type(atom_type),intent(in) :: at(nat)                   ! atoms
    integer,intent(in)         :: atmap(nat,nop)            ! atom mapping
    integer,intent(in)         :: atmapR(3,nat,nop)         ! atom mapping lattice
    real(dp),intent(in)        :: rep_at(1:2*Lmax+1,1:2*Lmax+1,0:Lmax,nat,nop)  ! rotation matrix
    complex(dp),intent(in)     :: reps_at(nspin,nspin,nat,nop)                  ! spin-rotation matrix
    integer,intent(in)         :: latt(3,nlatt)             ! lattice
    integer,intent(in)         :: LTidx(-Rmax:Rmax, -Rmax:Rmax, -Rmax:Rmax)     ! lattice index
    real(dp),intent(in)        :: cell_a(3,3), cell_b(3,3)  ! primitive cell
    logical,intent(in)         :: time_rev                  ! time reversal symmetry
    type(sym_op),intent(in)    :: op(nop)                   ! symmetry operations
    complex(dp),intent(inout)  :: ovlp(nsproj,nobt2,nlatt)  ! overlap matrix
    !
    integer :: aa, ss, aap, ssp, ii, jj, iip, jjp, nn, nnp, gg, RR, RRp, &
               LL, ii_, jj_, iip_, jjp_, xi
    real(kind=dp), dimension(3) :: aapos    ! used for atom mapping
    !
    ! symmetrize
    complex(dp)         :: ovlp_sym(nsproj, nobt2, nlatt)
    real(dp)            :: Ut(nspin,nspin), norm
    logical             :: badlat(nwat, nat, nlatt)
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
    do ss=1,nwat        ! slice-atom loop
        do gg=1,nop     ! operation loop
            g0R_ = latt(1,RR)*cell_a(:,1) + latt(2,RR)*cell_a(:,2) + &
                   latt(3,RR)*cell_a(:,3)
            g0R_ = matmul( op(gg)%R, g0R_ )
            call TonA( g0R_, g0R, aapos, cell_a, cell_b )
            !
            ! find RRp
            RRp = LTidx( g0R(1)+atmapR(1,aa,gg)-watmapR(1,ss,gg), &
                         g0R(2)+atmapR(2,aa,gg)-watmapR(2,ss,gg), &
                         g0R(3)+atmapR(3,aa,gg)-watmapR(3,ss,gg) )
            !
            if (RRp==0)  then
                !print*, latt(:,RR), ss, aa
                badlat(ss,aa,RR) = .true.
            endif
        enddo
    enddo
    enddo
    enddo
    !
    ! Time Reversal =================================================
    !
    if (time_rev) then
    !if (.false. ) then
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
            !RRp = LTidx(-latt(1,RR),-latt(2,RR),-latt(3,RR)) 
            !
            do aa=1,nat
            do LL=0,Lmax
            do jj=1,at(aa)%nobtL(LL)
            do nn=1,nspin
                jj_ = at(aa)%idx(jj,LL) + (nn-1)*nobt
            do ss=1,nwat
            do xi=1,wat(ss)%nirep
            do ii=1,wat(ss)%deg(xi)
                ii_ = wat(ss)%idx(xi) + ii - 1
                !print*, aa, LL, jj, nn, ii_, jj_
                !
                if(badlat(ss,aa,RR)) cycle
                !
                do nnp=1,nspin
                    jjp_ = at(aa)%idx(jj,LL) + (nnp-1)*nobt
                do iip=1,wat(ss)%deg(xi)
                    iip_ = wat(ss)%idx(xi) + iip - 1 
                    !
                    ! TR: <0si|Raj> -->  Ut(j'j)^* Ut(i'i) <0si'|Raj'>^* 
                    ovlp_sym(ii_,jj_,RR) = ovlp_sym(ii_,jj_,RR) &
                    + conjg(ovlp(iip_,jjp_,RR))*Ut(nnp,nn)*Ut_wat(iip,ii,xi,ss)
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
        enddo
        !
    endif
    !
    ovlp = 0.5*(ovlp + ovlp_sym)
    !
    ! Point group ===================================================
    !
    ovlp_sym = 0.0_dp
    do RR=1,nlatt
    do aa=1,nat
    do ss=1,nwat
        do gg=1,nop
            !
            g0R_ = latt(1,RR)*cell_a(:,1) + latt(2,RR)*cell_a(:,2) + &
                   latt(3,RR)*cell_a(:,3)
            g0R_ = matmul( op(gg)%R, g0R_ )
            call TonA( g0R_, g0R, aapos, cell_a, cell_b )
            !
            ! find RRp
            RRp = LTidx( g0R(1)+atmapR(1,aa,gg)-watmapR(1,ss,gg), &
                         g0R(2)+atmapR(2,aa,gg)-watmapR(2,ss,gg), &
                         g0R(3)+atmapR(3,aa,gg)-watmapR(3,ss,gg) )
            !
            ! matrix element
            aap = atmap(aa,gg)
            ssp = watmap(ss,gg)
            !
            ! remove bad lattice
            if (badlat(ss,aa,RR )) cycle   ! <=> RRp=0
            if (badlat(ssp,aap,RRp)) cycle
            !
            do LL=0,Lmax
            do jj=1,at(aa)%nobtL(LL)
            do nn=1,nspin
                jj_ = at(aa)%idx(jj,LL) + (nn-1)*nobt
            do xi=1,wat(ss)%nirep
            do ii=1,wat(ss)%deg(xi)
                ii_ = ii + wat(ss)%idx(xi) - 1
                !
                do jjp=1,at(aap)%nobtL(LL)
                do nnp=1,nspin
                    jjp_ = at(aap)%idx(jjp,LL) + (nnp-1)*nobt
                do iip=1,wat(ssp)%deg(xi)
                    iip_ = iip + wat(ssp)%idx(xi) - 1
                    !
                    ovlp_sym(ii_,jj_,RR) = ovlp_sym(ii_,jj_,RR)  &
                        + ovlp(iip_,jjp_,RRp)*rep_at(jjp,jj,LL,aa,gg)*reps_at(nnp,nn,aa,gg) &
                        * conjg( rep_wat(iip,ii,xi,ss,gg) )
                    !
                    !if ( ii_==1 .and. jj_==2 .and. abs(ovlp(iip_,jjp_,RRp))>1.0d-6 &
                    !     .and. abs(rep_wat(iip,ii,xi,ss,gg)*rep_at(jjp,jj,LL,aa,gg)*reps_at(nnp,nn,aa,gg))>1.0d-6  ) then
                    !    print*, '**** gg=',gg
                    !    print*, ii_,jj_, '<--', iip_, jjp_
                    !    print*, ovlp(iip_,jjp_,RRp)
                    !    print*, conjg(rep_wat(iip,ii,xi,ss,gg))*rep_at(jjp,jj,LL,aa,gg)*reps_at(nnp,nn,aa,gg)
                    !    print*, ovlp_sym(ii_,jj_,RR)
                    !endif
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
    !print*, ovlp_sym(1,2,LTidx(0,0,0))
    ovlp = ovlp_sym/nop
    !
    ! normalization ===========================================================
    !
    do jj=1,nsproj
        !
        norm = 0.0_dp
        do RR=1,nlatt
        do ii=1,nobt2
            !
            norm = norm + conjg(ovlp(jj, ii, RR))*ovlp(jj, ii, RR)
            !
        enddo
        enddo
        !
        ! normalize
        ovlp(jj,:,:) = ovlp(jj,:,:)/sqrt(norm)
        !
    enddo
    
endsubroutine symm_ovlp

!**** proj_wann() ***************************************************
!
! Get the wannier function by projection method 
!   input : Bloch wave-functions on a k-point mesh
!   output: The wannier gauge Uwan, and the S matrix determination
!
! Ref: Phys. Rev. B 83, 035108 ( A. A. Soluyanov & D. Vanderbilt )
!
!********************************************************************
subroutine proj_WF_ini(nwan)
    integer,intent(in)  :: nwan
    lwork = 6*nwan
    allocate( work(lwork), rwork(lwork) )
endsubroutine

subroutine proj_WF_rls()
    if( allocated(work) ) deallocate(work)
    if( allocated(rwork) ) deallocate(rwork)
endsubroutine

subroutine proj_WF(kp, nobt2, wave, bstart, bend, &
                   nwan, nlatt, latt, cell_a, ovlp, Uwan, Sdet)
    integer,intent(in)      :: nobt2, bstart, bend, nwan, nlatt
    real(dp),intent(in)     :: kp(3)                    ! k-point mesh
    complex(dp),intent(in)  :: wave(nobt2,nobt2)        ! Periodic part of Bloch wave-function
    integer,intent(in)      :: latt(3,nlatt)            ! lattice
    real(dp),intent(in)     :: cell_a(3,3)              ! lattice vector
    complex(dp),intent(in)  :: ovlp(nwan,nobt2,nlatt)   ! overlap matrix
    complex(dp),intent(out) :: Uwan(nobt2,nwan)         ! The wannier gauge
    real(dp),intent(out)    :: Sdet                     ! The determination of S matrix
    !
    integer         :: ii, jj, nn, aa, RR, infor
    real(dp)        :: rvec(3), Seig(nwan)
    complex(dp)     :: Smat(nwan,nwan), Wmat(nwan,nwan), Ugam(nobt2,nwan)
    !
    if ( bend-bstart+1/=nwan ) then
        write(*,*) 'Error in proj_WF(): number of wannier functions and bands mismatch !!!'
        stop
    endif
    !
    !---- W ---------------------------------------------------------
    ! W is the overlap between projector function and the Bloch 
    !   wave-function  
    !   W(i,n,k) = <0,i|k,n> = \sum_{R,a} <0,i|R,a>*U(a,n,k)*exp(iR*k)
    ! The projected Bloch states are
    !   |gam,k,i> = \sum_n |k,n>  W^*(i,n,k)
    !----------------------------------------------------------------
    !
    Wmat = 0.0_dp
    do ii=1,nwan
    do jj=1,nwan
        nn = jj - 1 + bstart
        !
        do RR=1,nlatt
            rvec = cell_a(:,1)*latt(1,RR) + cell_a(:,2)*latt(2,RR) + cell_a(:,3)*latt(3,RR)
            do aa=1,nobt2
                Wmat(ii,jj) = Wmat(ii,jj) + ovlp(ii,aa,RR)*wave(aa,nn) &
                                          * exp(c_im*dot_prod(rvec,kp))
            enddo
        enddo
        !
    enddo
    enddo
    !
    ! The periodic part of |gam,k,i>
    Ugam = matmul(wave(:,bstart:bend),transpose(conjg(Wmat)))
    !
    !---- S ----------------------------------------------------------
    ! S(i,j,k) = <G,k,i|G,k,j> = sum_n W(i,n,k) * W^*(j,n,k)
    !
    ! |wan,k,i> = sum_j |gam,k,j> S^{-1/2}(j,i)
    !-----------------------------------------------------------------
    !
    Smat = matmul(Wmat, transpose(conjg(Wmat)) )
    !
    ! Beware: zheev() saves the eigen-vectors in Smat
    call zheev( 'V','U', nwan, Smat, nwan,  &   ! jobz, uplo, n, a, lda
                Seig, work, lwork, rwork,   &   ! w, work, lwork, rwork
                infor )
    if (infor/=0) then
        write(*,*) 'Error in proj_WF(): zheev fails !!!'
        stop
    endif
    !
    !if (minval(Seig)<1.0e-4) then
    !    write(*,*) 'Error in proj_WF(): Smat singular !!!'
    !    stop
    !endif
    !
    ! Det(S)
    Sdet = 1.0_dp
    do ii=1,nwan
        Sdet = Sdet*Seig(ii)
    enddo
    !
    ! Smat <- Smat^{-1/2}
    do ii=1,nwan
        Wmat(:,ii) = Smat(:,ii)/sqrt(Seig(ii))
    enddo
    Smat = matmul(Wmat,transpose(conjg(Smat)))
    !
    ! Uwan is the periodic part of |wan,k,i>
    Uwan(:,:) = matmul( Ugam, Smat )
    !
endsubroutine proj_WF

endmodule wannier
