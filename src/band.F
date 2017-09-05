module band
!
use typcnst,    only : DP, c_pi, c_im, lsize, atom_type, Lmax, sym_op
use pntgrp,     only: TonA
use my_math,    only : vect_len, dot_prod
!
implicit none
public
!
! workspace for lapack
integer :: lwork
complex(kind=DP), allocatable,  dimension(:) :: work
real(kind=DP), allocatable, dimension(:) :: rwork
!
!
contains

subroutine bndspc_ini(nobt2)
    integer,intent(in)  :: nobt2
    !
    lwork = 6*nobt2
    allocate( work(lwork),rwork(9*nobt2))
endsubroutine

subroutine bndspc_rls()
    if(allocated(work)) deallocate(work)
    if(allocated(rwork)) deallocate(rwork)
endsubroutine

subroutine bnd_slv( kp, wave, eng, nobt, nobt2, nspin, hr_lat_num, hr_dat, hr_lat, natom, atom, cell_a)
    integer,intent(in)              :: nobt, nobt2, nspin
    integer,intent(in)              :: hr_lat_num
    integer,intent(in)              :: natom
    real(dp),intent(in)             :: kp(3)
    complex(dp),intent(out)         :: wave(nobt2,nobt2)
    real(dp),intent(out)            :: eng(nobt2)
    complex(dp),intent(in)          :: hr_dat(nobt2,nobt2,hr_lat_num)
    integer,intent(in)              :: hr_lat(3,hr_lat_num)
    type(atom_type)                 :: atom(natom)
    real(dp),intent(in)             :: cell_a(3,3)
    !
    integer :: ll, aa, bb, La, Lb, ii, jj, nn, mm, ii_, jj_
    !
    real(kind=dp), dimension(3)              :: latt
    complex(kind=dp), dimension(nobt2,nobt2) :: hamk
    real(kind=dp)       :: kdr
    integer :: infor
    !
    if(nobt2/=nobt*nspin) then
        write(*,*) 'Error in bnd_slv(): wrong nobt2 !!!'
        stop
    endif
    !
    ! Fourier transformation
    !
    hamk(:,:) = 0.0_dp
    do ll=1,hr_lat_num
        latt = cell_a(:,1)*hr_lat(1,ll) + cell_a(:,2)*hr_lat(2,ll) + &
                                          cell_a(:,3)*hr_lat(3,ll)
        !
        do aa=1,natom
        do bb=1,natom
            ! Here we adopt the periodic gauge
            !kdr = dot_prod( latt + atom(bb)%pos - atom(aa)%pos, kp )
            kdr = dot_prod( latt, kp ) 
            do La=0,Lmax
            do Lb=0,Lmax
                do ii=1, atom(aa)%nobtL(La)
                do jj=1, atom(bb)%nobtL(Lb)
                    do nn=0,nspin-1
                    do mm=0,nspin-1
                        ii_ = atom(aa)%idx(ii,La) + nn*nobt
                        jj_ = atom(bb)%idx(jj,Lb) + mm*nobt
                        hamk(ii_,jj_) = hamk(ii_,jj_) + hr_dat(ii_,jj_,ll)*exp(c_im*kdr)
                    enddo
                    enddo
                enddo
                enddo
            enddo
            enddo
        enddo
        enddo
        !
    enddo
    ! hamk = 0.5_dp*( hamk+transpose(conjg(hamk)) )
    ! 
    call zheev( 'V','U', &                        ! jobz, uplo
                nobt2, hamk, nobt2, &             ! n, a, lda
                eng, work, lwork, rwork, &        ! w, work, lwork, rwork
                infor )
    !
    wave = hamk
    !
    if (infor/=0) then
        write(*,'(1X, "Error in bnd_slv() when calling to zheevr, info =", I6)') infor
    endif
    !
endsubroutine

!
!**** impose_sym() **************************************************
! Impose symmetry on the real space Hamiltonian
!********************************************************************
subroutine impose_sym(hr_lat_num, nobt, nobt2, nspin, hr_dat, hr_lat, Rmax, LTidx, nop, op, &
                      time_rev, natom, atom, atmap, atmapR, rep_at, reps_at, &
                      cell_a, cell_b)
    !
    integer,intent(in)          :: hr_lat_num
    integer,intent(in)          :: nobt, nobt2, nspin
    complex(dp),intent(inout)   :: hr_dat(nobt2,nobt2,hr_lat_num)
    integer,intent(in)          :: hr_lat(3,hr_lat_num)
    integer,intent(in)          :: Rmax
    integer,intent(in)          :: LTidx(-Rmax:Rmax,-Rmax:Rmax,-Rmax:Rmax)
    integer,intent(in)          :: nop
    type(sym_op),intent(in)     :: op(nop)
    logical,intent(in)          :: time_rev
    integer,intent(in)          :: natom
    type(atom_type),intent(in)  :: atom(natom)
    integer,intent(in)          :: atmap(natom,nop)
    integer,intent(in)          :: atmapR(3,natom,nop)
    real(dp),intent(in)         :: rep_at(1:2*Lmax+1, 1:2*Lmax+1, 0:Lmax, natom, nop)
    complex(dp),intent(in)      :: reps_at(1:nspin,1:nspin,natom,nop)
    real(dp),intent(in)         :: cell_a(3,3), cell_b(3,3)
    !
    integer :: aa, bb, aap, bbp, ii, jj, iip, jjp, nn, mm, nnp, mmp, &
               gg, RR, RRp, La, Lb, ii_, jj_, iip_, jjp_
    real(kind=dp), dimension(3) :: aapos    ! used for atom mapping
    !
    ! the symmetrized hr
    complex(kind=dp), dimension(nobt2,nobt2,hr_lat_num) :: hr_sym
    real(kind=dp), dimension(nspin,nspin) :: Ut 
    logical :: badlat(lsize,lsize,5000)
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
    ! symmetrize the lattice
    !
    badlat = .false.
    do RR=1,hr_lat_num   ! RR loop
    do aa=1,natom        ! atom loop
    do bb=1,natom
        do gg=1,nop      ! gg loop
            g0R_ = hr_lat(1,RR)*cell_a(:,1) + hr_lat(2,RR)*cell_a(:,2) + &
                   hr_lat(3,RR)*cell_a(:,3)
            g0R_ = matmul( op(gg)%R, g0R_ )
            call TonA( g0R_, g0R, aapos, cell_a, cell_b )
            !
            ! find RRp
            RRp = LTidx( g0R(1)+atmapR(1,bb,gg)-atmapR(1,aa,gg), &
                         g0R(2)+atmapR(2,bb,gg)-atmapR(2,aa,gg), &
                         g0R(3)+atmapR(3,bb,gg)-atmapR(3,aa,gg) )
            !
            if (RRp==0) badlat(aa,bb,RR) = .true.
        enddo
    enddo
    enddo
    enddo
    !
    ! Hermitian =====================================================
    !
    !hr_sym(:,:,:) = 0.0_dp
    !do RR=1,hr_lat_num
        !
    !    RRp = LTidx(-hr_lat(1,RR),-hr_lat(2,RR),-hr_lat(3,RR) )
    !    if (RRp==0) cycle 
    !    !
    !    do ii=1,nobt2
    !    do jj=1,nobt2
    !        hr_sym(ii,jj,RR) = conjg(hr_dat(jj,ii,RRp))
    !    enddo
    !    enddo
    !    !
    !enddo
    !hr_dat(:,:,:) = 0.5_dp*( hr_dat(:,:,:)+hr_sym(:,:,:) )
    !
    ! Point Group====================================================
    !
    hr_sym(:,:,:) = 0.0_dp
    !
    do RR=1,hr_lat_num      ! RR loop
    do aa=1,natom           ! atom loop
    do bb=1,natom
        do gg=1,nop      ! gg loop
            g0R_ = hr_lat(1,RR)*cell_a(:,1) + hr_lat(2,RR)*cell_a(:,2) + &
                   hr_lat(3,RR)*cell_a(:,3)
            g0R_ = matmul( op(gg)%R, g0R_ )
            call TonA( g0R_, g0R, aapos, cell_a, cell_b )
            !
            ! find RRp
            RRp = LTidx( g0R(1)+atmapR(1,bb,gg)-atmapR(1,aa,gg), &
                         g0R(2)+atmapR(2,bb,gg)-atmapR(2,aa,gg), &
                         g0R(3)+atmapR(3,bb,gg)-atmapR(3,aa,gg) )
            !
            ! matrix element
            aap = atmap(aa,gg)
            bbp = atmap(bb,gg)
            !
            ! remove bad lattice
            if (badlat(aa,bb,RR )) cycle   ! <=> RRp=0
            if (badlat(aap,bbp,RRp)) cycle
            !
            do La=0,Lmax
            do ii=1,atom(aa)%nobtL(La)
            do nn=1,nspin
                ii_ = atom(aa)%idx(ii,La) + (nn-1)*nobt
            do Lb=0,Lmax
            do jj=1,atom(bb)%nobtL(Lb)
            do mm=1,nspin
                jj_ = atom(bb)%idx(jj,Lb) + (mm-1)*nobt
                !
                do iip=1,atom(aap)%nobtL(La)
                do nnp=1,nspin
                    iip_=atom(aap)%idx(iip,La)+(nnp-1)*nobt
                do jjp=1,atom(bbp)%nobtL(Lb)
                do mmp=1,nspin
                    jjp_=atom(bbp)%idx(jjp,Lb)+(mmp-1)*nobt
                    !
                    hr_sym(ii_,jj_,RR) = hr_sym(ii_,jj_,RR) + &
                        hr_dat(iip_,jjp_,RRp)*rep_at(iip,ii,La,aa,gg)*rep_at(jjp,jj,Lb,bb,gg)&
                                             *conjg(reps_at(nnp,nn,aa,gg))*reps_at(mmp,mm,bb,gg)
                    !
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
            !
        enddo
    enddo       ! atom loop
    enddo       ! gg loop
    enddo       ! RR loop
    !
    hr_dat(:,:,:) = hr_sym(:,:,:)/nop
    !
    ! Time Reversal =================================================
    !
    if (time_rev) then
        !
        hr_sym(:,:,:) = 0.0_dp
        !
        if (nspin==2) then
            Ut(1,:) = (/ 0.0_dp, -1.0_dp /)
            Ut(2,:) = (/ 1.0_dp,  0.0_dp /)
        else
            Ut(1,1) = 1.0_dp
        endif
        !
        do RR=1,hr_lat_num
            !
            RRp = LTidx(-hr_lat(1,RR),-hr_lat(2,RR),-hr_lat(3,RR))
            if(RRp==0) cycle
            !
            do aa=1,natom
            do La=0,Lmax
            do ii=1,atom(aa)%nobtL(La)
            do nn=1,nspin
                ii_ = atom(aa)%idx(ii,La) + (nn-1)*nobt
            do bb=1,natom
            do Lb=0,Lmax
            do jj=1,atom(bb)%nobtL(Lb)
            do mm=1,nspin
                jj_ = atom(bb)%idx(jj,Lb) + (mm-1)*nobt
                !
                do nnp=1,nspin
                    iip_ = atom(aa)%idx(ii,La) + (nnp-1)*nobt
                do mmp=1,nspin
                    jjp_ = atom(bb)%idx(jj,Lb) + (mmp-1)*nobt
                    !
                    hr_sym(ii_,jj_,RR) = hr_sym(ii_,jj_,RR) &
                     + hr_dat(jjp_,iip_,RRp)*Ut(nnp,nn)*Ut(mmp,mm)
                    !
                enddo
                enddo
                !
            enddo
            enddo
            enddo
            enddo
            !
            enddo
            enddo
            enddo
            enddo
        enddo
        !
        hr_dat(:,:,:) = 0.5_dp*(hr_dat(:,:,:)+hr_sym(:,:,:))
        !
    endif
    !
    !
endsubroutine impose_sym

!**** get_sewing() **************************************************
! get the sewing matrix of Bloch bases:
!   g|k,n> = \sum_m S_mn(k) |gk,m>
!********************************************************************
subroutine get_sewing(nobt, nobt2, nspin, nop, op, opidx, natom, atom, atmap, atmapR, &
                      rep_at, reps_at, cell_a, kp, sewing)
    !
    integer,intent(in)          :: nobt, nobt2, nspin
    integer,intent(in)          :: nop
    type(sym_op),intent(in)     :: op(nop)
    integer,intent(in)          :: opidx(48)
    integer,intent(in)          :: natom
    type(atom_type),intent(in)  :: atom(natom)
    integer,intent(in)          :: atmap(natom,48)
    integer,intent(in)          :: atmapR(3,natom,48)
    real(dp),intent(in)         :: rep_at(1:2*Lmax+1, 1:2*Lmax+1, 0:Lmax, natom, 48)
    complex(dp),intent(in)      :: reps_at(1:nspin,1:nspin,natom,48)
    real(dp),intent(in)         :: cell_a(3,3), kp(3)
    complex(dp),intent(out)     :: sewing(nobt2,nobt2,nop)
    !
    integer   :: aa, bb, gg, gidx, La, ii, jj, ii_, jj_, nn, mm
    real(dp)  :: dR(3), gkp(3)
    !
    if( nobt2/=nobt*nspin ) then
        write(*,*) 'Error in get_sewing(): wrong nobt2 !!!'
        stop
    endif
    !
    sewing = 0.0_dp
    do gg=1,nop
        gkp = matmul( op(gg)%R, kp )
        gidx = opidx(gg)
        !
        do aa=1,natom
            dR = cell_a(:,1)*atmapR(1,aa,gidx) + cell_a(:,2)*atmapR(2,aa,gidx) &
               + cell_a(:,3)*atmapR(3,aa,gidx)
            bb = atmap(aa,gidx)
            !
            do La=0,Lmax
            do ii=1,atom(aa)%nobtL(La)
            do nn=1,nspin
                ii_ = atom(aa)%idx(ii,La) + (nn-1)*nobt
                !
                do jj=1,atom(bb)%nobtL(La)
                do mm=1,nspin
                    jj_ = atom(bb)%idx(jj,La) + (mm-1)*nobt
                    sewing(jj_,ii_,gg) = exp( -c_im*dot_prod(gkp,dR) ) &
                          *rep_at(jj,ii,La,aa,gidx)*reps_at(mm,nn,aa,gidx)
                enddo
                enddo
                !
            enddo
            enddo
            enddo
            !
        enddo
    enddo
    !
endsubroutine get_sewing


function fun_Efm( nkp, nobt2, eng, nele, smear)
    integer,intent(in)      :: nkp, nobt2
    real(dp),intent(in)     :: eng(nobt2,nkp),nele,smear
    real(dp)                :: fun_Efm
    !
    real(dp) :: mu1, mu2, mu0, mu0_, N1, N2, N0
    !
    if (nele>=nobt2 .or. nele<=0.0_dp) then
        write(*,*) 'Error in fun_Efm(): nele must satisfy 0<nele<nobt2 !!!'
        stop
    endif
    !
    mu1=minval(eng(1,:))
    mu2=maxval(eng(nobt2,:))
    !
    N1=fun_nele(nkp, nobt2,  eng, mu1, smear)-nele
    N2=fun_nele(nkp, nobt2,  eng, mu2, smear)-nele
    !
    if (  N1>0.0_dp .or. N2<0.0_dp ) then
        write(*,*) 'Error in fun_Efm() !!!'
        stop
    endif
    !
    mu0_=99999._dp
    do 
        mu0=(mu1*N2-mu2*N1)/(N2-N1)
        N0=fun_nele(nkp, nobt2,  eng, mu0, smear)-nele
        if (N0>0) then
            mu2=mu0
            N2 =N0
        elseif (N0<0) then
            mu1=mu0
            N1 =N0
        elseif (N0==0) then
            exit
        endif
        !
        if( abs(mu0-mu0_)<1.0e-4 ) exit
        mu0_ = mu0
        !
    enddo
    !
    fun_Efm=mu0
    !
endfunction

!**** fun_nele ******************************************************
! Electrons number as function of chemical potential
!********************************************************************
function fun_nele(nkp, nobt2, eng, mu, smear)
    integer,intent(in)      :: nkp, nobt2
    real(dp),intent(in)     :: eng(nobt2,nkp),mu,smear
    real(dp)                :: fun_nele
    !
    integer :: ii,jj
    !
    fun_nele = 0.0_dp
    do ii=1,nkp
    do jj=1,nobt2
        fun_nele = fun_nele + fun_fmd( eng(jj,ii)-mu, smear )
    enddo
    enddo
    !
    fun_nele = fun_nele/nkp
endfunction


!**** fun_fmd *******************************************************
! Fermi Distribution
!********************************************************************
function fun_fmd( xi, kT )
    real(dp) :: xi, kT, fun_fmd
    fun_fmd=1.0_dp/(1.0_dp + exp(xi/kT) )
endfunction

!**** fun_pfmd *******************************************************
! Fermi Distribution' derivatives
!********************************************************************
function fun_pfmd( xi, kT )
    real(dp) :: xi, kT, fun_pfmd
    fun_pfmd=-1.0_dp/( kT * (2.0_dp+exp(xi/kT)+exp(-xi/kT)) )
endfunction

endmodule band
