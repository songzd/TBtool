module tghtbd
!
use typcnst,    only : DP, c_pi, c_im, lsize, Lmax, LTmax, atom_type, obt_type
use my_math,    only : x_prod, dot_prod, mix_prod, angle, vect_len
use pntgrp,     only : TonA, euler_angle
use reprsnt,    only : dmat_spin
!
implicit none
!
private 
public :: get_tb, cell_a, cell_b, h0_dat, h0_lat, h0_lat_num, & 
          nobt, nobt2, nspin, natom, atom,  &
          h0_LTidx, obt
!
! Input parameters
real(kind=dp), dimension(3,3)     :: cell_a      ! Lattice and Reciprocal Lattice
real(kind=dp), dimension(3,3)     :: cell_b      ! Lattice and Reciprocal Lattice
integer, dimension(lsize)         :: element     ! Elements of atoms
real(kind=dp), dimension(3,lsize) :: atom_pos    ! Atom postions
real(kind=dp), dimension(3,lsize) :: loc_coord_z ! Local coordinates on atoms 
real(kind=dp), dimension(3,lsize) :: loc_coord_x ! Local coordinates on atoms 
character(len=30)              :: atom_pos_conv  ! Coordinates convention 
character(len=30)              :: loc_coord_conv ! Coorindates convention
integer, dimension(4,lsize)       :: obt_list    ! Orbitals input
real(kind=dp), dimension(3,lsize) :: soc_list    ! Soc input
character(len=30)                 :: soc_type    ! Soc type: nsoc, global, or local
!
namelist /TightBinding/ cell_a,  atom_pos, loc_coord_z, loc_coord_x, obt_list, &
                        element, atom_pos_conv, loc_coord_conv, &
                        soc_type, soc_list
!
! atoms and orbitals
integer                             :: nobt         ! number of nsoc orbitals
integer                             :: nobt2        ! number of  soc orbitals
integer                             :: nspin        ! 1: nsoc,  2: soc
integer                             :: natom        ! number of atoms
type(atom_type), dimension(lsize)   :: atom         ! atom arrays
character(len=30), dimension(-Lmax:Lmax,0:Lmax) :: obtsym ! orbital symbols
type(obt_type),allocatable,dimension(:)         :: obt  ! Orbitals
!
! Hr.dat
complex(kind=dp), allocatable, dimension(:,:,:) :: h0_dat
integer, allocatable, dimension(:,:)            :: h0_lat
integer                                         :: h0_lat_num
!
! soc matrix
complex(kind=dp), dimension(-Lmax:Lmax,1:2,-Lmax:Lmax,1:2, 0:Lmax)   ::  socmat
!
! Lattice index
integer, dimension(-LTmax:LTmax,-LTmax:LTmax,-LTmax:LTmax) :: h0_LTidx
contains
!
!
subroutine get_tb()
    call get_obt()
    call get_hr()
endsubroutine
!
subroutine get_obt()
    integer    :: ii, jj, ll, mm, nn, ss
    real(kind=dp)       :: vol
        real(dp),dimension(3) :: vtmp
        integer, dimension(3) :: ltmp
        !
        ! read structure =================================================
        !
        ! initialize
        element(:)=0
        obt_list(:,:)=0
        soc_list(:,:)=0.0_dp
        !
        open(unit=100, file='system.in', status='old')
        read(100, TightBinding) 
        close(100)
        !
        ! Reciprocal Lattice
        vol = mix_prod(cell_a(:,1),cell_a(:,2),cell_a(:,3))
        cell_b(:,1) = x_prod(cell_a(:,2),cell_a(:,3))/vol *2*c_pi
        cell_b(:,2) = x_prod(cell_a(:,3),cell_a(:,1))/vol *2*c_pi
        cell_b(:,3) = x_prod(cell_a(:,1),cell_a(:,2))/vol *2*c_pi
        !
        ! count atoms
        natom=0
        do ii=1,lsize
            if (element(ii)==0) exit
            natom = natom + 1
            atom(natom)%elmt=element(ii)
        enddo
        !
        ! atom positions
        if (atom_pos_conv(1:1)=='d') then
            do ii=1,natom
                atom_pos(:,ii) = atom_pos(1,ii)*cell_a(:,1) + atom_pos(2,ii)*cell_a(:,2) &
                                                            + atom_pos(3,ii)*cell_a(:,3)
            enddo
        endif
        !
        do ii=1,natom
            call TonA( atom_pos(:,ii), ltmp, vtmp, cell_a, cell_b )
            if (ltmp(1)/=0 .or. ltmp(2)/=0 .or. ltmp(3)/=0) then
                write(*,*) 'Error in get_obt(): atom is not in 0 0 0 lattice !!!'
                stop
            endif
            atom(ii)%pos = atom_pos(:,ii)
        enddo
        !
        ! local coordinates
        if (loc_coord_conv(1:1)=='d') then
            do ii=1,natom
                loc_coord_z(:,ii) =  loc_coord_z(1,ii)*cell_a(:,1) + loc_coord_z(2,ii)*cell_a(:,2) &
                                                                   + loc_coord_z(3,ii)*cell_a(:,3)
                loc_coord_x(:,ii) =  loc_coord_x(1,ii)*cell_a(:,1) + loc_coord_x(2,ii)*cell_a(:,2) &
                                                                   + loc_coord_x(3,ii)*cell_a(:,3)
            enddo
        endif
        do ii=1,natom
            call euler_angle(loc_coord_z(:,ii), loc_coord_x(:,ii), &
                         atom(ii)%alph, atom(ii)%beta, atom(ii)%gamm)
        enddo
        !
        ! Atoms
        nobt=0
        do ii=1,natom       ! atoms
            !
            atom(ii)%nobtL=0
            do jj=1,lsize
                if ( obt_list(1,jj) == 0 ) exit
                if ( obt_list(1,jj) == atom(ii)%elmt ) then
                    nobt = nobt + 1
                    !
                    ll = obt_list(2,jj)
                mm = obt_list(3,jj)
                !
                nn=atom(ii)%nobtL(ll) + 1
                atom(ii)%nobtL(ll)  = nn
                atom(ii)%m(nn,ll)   = mm
                atom(ii)%idx(nn,ll) = nobt
                !
                if ( obt_list(4,jj)==1 ) then
                    atom(ii)%cr(nn,ll) = .true.
                else
                    atom(ii)%cr(nn,ll) = .false.
                endif
                !
            endif
        enddo
        !
    enddo
    !
    if ( soc_type(1:1)=='n' ) then
        nspin=1
    else
        nspin=2
    endif
    nobt2=nobt*nspin
    !
    if (nobt2==0) then
        write(*,*) "Error in get_obt(): no orbitals found !!!"
        stop
    endif
    !
    ! Orbitals
    !
    allocate( obt(nobt2) )
    do nn=1,natom
    do LL=0,Lmax
    do ii=1,atom(nn)%nobtL(LL)
    do ss=1,nspin
        jj=atom(nn)%idx(ii,LL)+(ss-1)*nobt
        !
        obt(jj)%Elmt=atom(nn)%Elmt
        obt(jj)%at  =nn
        obt(jj)%pos =atom(nn)%pos
        obt(jj)%alph=atom(nn)%alph
        obt(jj)%beta=atom(nn)%beta
        obt(jj)%gamm=atom(nn)%gamm
        obt(jj)%L   =LL
        obt(jj)%nthL=ii
        obt(jj)%M   =atom(nn)%M(ii,LL)
        obt(jj)%S   =ss
        obt(jj)%cr  =atom(nn)%cr(ii,LL)
        !
    enddo
    enddo
    enddo
    enddo
    !
    ! write Lattice, Atom and Orbital information
    !
    write(*, '( 1X, "Lattice:")')
    write(*, '( 5X, "  a1=", 3F10.5  )') cell_a(:,1)
    write(*, '( 5X, "  a2=", 3F10.5  )') cell_a(:,2)
    write(*, '( 5X, "  a3=", 3F10.5  )') cell_a(:,3)
    write(*, '( 5X, "  b1=", 3F10.5, " pi")') cell_b(:,1)/(2*c_pi)
    write(*, '( 5X, "  b2=", 3F10.5, " pi")') cell_b(:,2)/(2*c_pi)
    write(*, '( 5X, "  b3=", 3F10.5, " pi")') cell_b(:,3)/(2*c_pi)
    !
    write(*, '( 1X, "Atoms:" )')
    write(*, '( 5X, " Atom", " Elmt", "  Position (cart)", 19X, "Euler Angle"  )')
168 format(5X, 2I5, ' (', 3F10.5, ')  (', 3F8.1, ')' )
    write(*, 168) ( ii, atom(ii)%elmt, atom(ii)%pos, atom(ii)%alph/c_pi*180, &
                  atom(ii)%beta/c_pi*180, atom(ii)%gamm/c_pi*180, ii=1,natom )
    !
    write(*, '( 1X, "Orbitals:")')
    write(*, '( 5X, "  obt", " Atom", "    L", "    M", "    Symbol", " Spin", " Correlated")' )
    !
    call init_obtsym()
    !
#if 1
    do ss=1,nspin
    do ii=1, natom
    do ll=0, Lmax
    if (atom(ii)%nobtL(ll)>0) then
        write(*,175) ( atom(ii)%idx(jj,ll)+(ss-1)*nobt, ii, ll, atom(ii)%m(jj,ll), &
                       trim(obtsym(atom(ii)%m(jj,ll), ll)), ss, &
                       jj=1,atom(ii)%nobtL(ll) )
    endif
    enddo
    enddo
    enddo
175     format(5X,4I5, A10, I5) 
    !
#else
    do ii=1,nobt2
        write(*,176) ii, obt(ii)%at, obt(ii)%L, obt(ii)%M, &
                     trim(obtsym(obt(ii)%M,obt(ii)%L)), obt(ii)%S, obt(ii)%cr
    enddo
#endif
176     format(5X,4I5, A10, I5, L5) 
    !
    !
endsubroutine get_obt

subroutine get_hr
    integer :: ii, jj, ll, aa, bb, La, Lb, nn, mm, nnp, mmp, &
               ii_, jj_, iip_, jjp_, kk
    !
    complex(kind=dp), allocatable, dimension(:,:,:) :: h0_dat_tmp
    integer, allocatable, dimension(:,:)            :: h0_lat_tmp
    integer, dimension(3)       :: lat_tmp
    real(kind=dp)               :: rr, ss
    real(dp),dimension(3)       :: latt
    integer                     :: iost
    !
    complex(kind=dp), dimension(2,2,natom) :: dspin
    !   
    h0_lat_num=5000
    allocate(h0_dat_tmp(nobt2,nobt2,h0_lat_num), h0_lat_tmp(3,h0_lat_num))
    h0_dat_tmp(:,:,:) = 0.0_dp
    h0_lat_tmp(:,:) = 0.0_dp
    !   
    open(unit=100, file='hr.dat', status='old')
    !   
    ll=0
    h0_LTidx(:,:,:)=0
    do  
        read(100,*,iostat=iost) lat_tmp(:), ii, jj, rr, ss
        if (iost/=0) exit
        !
        latt = cell_a(:,1)*lat_tmp(1) + cell_a(:,2)*lat_tmp(2) + cell_a(:,3)*lat_tmp(3)
        !   
        !if ( lat_tmp(1)/=lat_tmp_last(1) .or. lat_tmp(2)/=lat_tmp_last(2) .or. &
        !     lat_tmp(3)/=lat_tmp_last(3) ) then
        if ( h0_LTidx(lat_tmp(1),lat_tmp(2),lat_tmp(3))==0 ) then
            !
            ll=ll+1
            h0_lat_tmp(:,ll) = lat_tmp
            h0_LTidx(lat_tmp(1),lat_tmp(2),lat_tmp(3))=ll
            h0_dat_tmp(ii,jj,ll) = cmplx(rr,ss,kind=dp)
            !
        else
            !
            h0_dat_tmp(ii,jj, h0_LTidx(lat_tmp(1),lat_tmp(2),lat_tmp(3)) ) &
                = cmplx(rr,ss,kind=dp)
            !
        endif
        !   
        !lat_tmp_last = lat_tmp
        !   
        !   
    enddo
    if ( h0_LTidx(0,0,0)==0 ) then
        !
        ! This is used for SOC elements
        ll=ll+1
        h0_lat_tmp(:,ll) = (/0,0,0/)
        h0_LTidx(0,0,0)=ll
        h0_dat_tmp(:,:,ll) = 0.0_dp
    endif
    !   
    close(100)
    !
    h0_lat_num=ll
    write(*,'(1X,"hr.dat:")')
    write(*,'(5X,"h0_lat_num =",I10)') h0_lat_num
    write(*,'(5X,"soc_type   =",A10)') trim(soc_type)
    write(*,'(5X,"soc_list   = Elmt    L    Lamb")')
    do ii=1,lsize
        if ( nint(soc_list(1,ii))==0 ) exit
        write(*,'(17X,2I5,F10.4)') nint(soc_list(1,ii)), nint(soc_list(2,ii)), &
                                   soc_list(3,ii)
    enddo
    !
    ! SOC ===========================================================
    !
    ! if the input is nsoc, then copy up to down
    !if ( soc_type(1:1)=='n' ) then
    !    h0_dat_tmp(nobt+1:nobt2, nobt+1:nobt2, 1:h0_lat_num) = &
    !                h0_dat_tmp(1:nobt, 1:nobt, 1:h0_lat_num)
    !endif
    !
    allocate( h0_dat(nobt2,nobt2,h0_lat_num), h0_lat(3,h0_lat_num) )
    h0_dat(:,:,:) = 0.0_dp
    h0_lat(1:3,1:h0_lat_num) = h0_lat_tmp(1:3,1:h0_lat_num)
    !
    ! if the spin axes are defined in global coordinates, rotate them!
    if ( soc_type(1:1)=='m' .or. soc_type(1:1)=='g' ) then
        do aa=1,natom
            call dmat_spin( atom(aa)%alph, atom(aa)%beta, atom(aa)%gamm, dspin(:,:,aa))
        enddo
        !
        do ll=1,h0_lat_num
        do aa=1,natom
        do La=0,Lmax
        if (atom(aa)%nobtL(La)>0) then
            do bb=1,natom
            do Lb=0,Lmax
            if (atom(bb)%nobtL(Lb)>0) then
                do ii=1,atom(aa)%nobtL(La)
                do jj=1,atom(bb)%nobtL(Lb)
                    do nn =1,nspin
                    do mm =1,nspin
                    do nnp=1,nspin
                    do mmp=1,nspin
                        ii_ = atom(aa)%idx(ii,La) + (nn-1)*nobt
                        jj_ = atom(bb)%idx(jj,Lb) + (mm-1)*nobt
                        iip_= atom(aa)%idx(ii,La) + (nnp-1)*nobt
                        jjp_= atom(bb)%idx(jj,Lb) + (mmp-1)*nobt
                        !
                        h0_dat(ii_,jj_,ll) = h0_dat(ii_,jj_,ll) + h0_dat_tmp(iip_,jjp_,ll) &
                                               *conjg( dspin(nnp,nn,aa) )*dspin(mmp,mm,bb)
                    enddo
                    enddo
                    enddo
                    enddo
                enddo
                enddo
            endif
            enddo
            enddo
        endif
        enddo
        enddo
        enddo
        !
    else
        h0_dat(1:nobt2,1:nobt2,1:h0_lat_num) = h0_dat_tmp(1:nobt2,1:nobt2,1:h0_lat_num)
    endif
    !
    deallocate(h0_dat_tmp,h0_lat_tmp)
    !
#if 0 
    open(100,file='debug_hr.dat',status='replace')
    do ll=1,h0_lat_num
    do jj=1,nobt2
    do ii=1,nobt2
260 format(5I5,2F12.6)
        write(100,260) h0_lat(:,ll), ii, jj, h0_dat(ii,jj,ll)
    enddo
    enddo
    enddo
    close(100)
#endif
    !
    ! if the input hr is nsoc, then add it mannuly
    if ( soc_type(1:1)=='m' ) then
        !
        call init_soc()
        !
        ! find R=0
        do ll=1,h0_lat_num
            if ( h0_lat(1,ll)==0 .and. h0_lat(2,ll)==0 .and. h0_lat(3,ll)==0 ) exit
        enddo
        !
        do aa=1, natom
        do La=0, Lmax
        if (atom(aa)%nobtL(La)>0) then
            !
            do kk=1, lsize
                if (nint(soc_list(1,kk))==0) then
                    exit
                elseif ( nint(soc_list(1,kk))==atom(aa)%elmt .and. &
                         nint(soc_list(2,kk))==La ) then
                    !
                    do ii=1, atom(aa)%nobtL(La)
                    do jj=1, atom(aa)%nobtL(La)
                    do nn=1,nspin
                    do mm=1,nspin
                        !
                        ii_ = atom(aa)%idx(ii,La) + (nn-1)*nobt
                        jj_ = atom(aa)%idx(jj,La) + (mm-1)*nobt
                        !
                        ! add soc
                        h0_dat(ii_,jj_,ll) = h0_dat(ii_,jj_,ll) + 0.5*soc_list(3,kk) &
                               *socmat( atom(aa)%M(ii,La),nn, atom(aa)%M(jj,La),mm, La )
                        !
                    enddo
                    enddo
                    enddo
                    enddo
                    !
                endif
            enddo
            
            !
        endif
        enddo
        enddo
        !
    endif
endsubroutine get_hr

subroutine init_obtsym()
    !
    obtsym( 0,0) = 's'
    !
    obtsym(-1,1) = 'py'
    obtsym( 0,1) = 'pz'
    obtsym( 1,1) = 'px'
    !
    obtsym(-2,2) = 'dxy'
    obtsym(-1,2) = 'dyz'
    obtsym( 0,2) = 'dz2'
    obtsym( 1,2) = 'dzx'
    obtsym( 2,2) = 'dx2-y2'
    !
    obtsym(-3,3) = 'y(3x2-y2)'
    obtsym(-2,3) = 'xyz'
    obtsym(-1,3) = 'yz2'
    obtsym( 0,3) = 'z3'
    obtsym( 1,3) = 'xz2'
    obtsym( 2,3) = 'z(x2-y2)'
    obtsym( 3,3) = 'x(x2-3y2)'
endsubroutine init_obtsym

subroutine init_soc()
    !
    integer :: ii, jj
    !
    complex(kind=dp), parameter :: c_i3 = c_im*sqrt(3.0_dp), &
                                   c_i2 = c_im*2, &
                                   c_r3 = sqrt(3.0_dp), &
                                   c_r0 = cmplx(0.0_dp, 0.0_dp), &
                                   c_r1 = cmplx(1.0_dp, 0.0_dp)
    socmat(:,:, :,:, :) = 0.0_dp
    !
    ! soc_p up-up
    !                        pyu,  pzu,  pxu
    socmat(-1,1, -1:1,1, 1) = (/ c_r0, c_r0, c_im  /)   ! pyu
    socmat( 0,1, -1:1,1, 1) = (/ c_r0, c_r0, c_r0  /)   ! pzu
    socmat( 1,1, -1:1,1, 1) = (/-c_im, c_r0, c_r0  /)   ! pxu
    ! soc_p up-dw
    !                        pyd,  pzd,  pxd
    socmat(-1,1, -1:1,2, 1) = (/ c_r0,-c_im, c_r0  /)   ! pyu
    socmat( 0,1, -1:1,2, 1) = (/ c_im, c_r0,-c_r1  /)   ! pzu
    socmat( 1,1, -1:1,2, 1) = (/ c_r0, c_r1, c_r0  /)   ! pxu
    do ii=-1,1
    do jj=-1,1
        socmat(ii,2, jj,1, 1) = conjg( socmat(jj,1, ii,2, 1) )
    enddo
    enddo
    !
    socmat(:,2,:,2, 1) = -socmat(:,1,:,1, 1)
    !
    ! soc_d up-up
    !                        xyu,  yzu,  z2u,  zxu,  x2u
    socmat(-2,1, -2:2,1, 2) = (/ c_r0, c_r0, c_r0, c_r0, c_i2/)   ! xyu
    socmat(-1,1, -2:2,1, 2) = (/ c_r0, c_r0, c_r0, c_im, c_r0/)   ! yzu
    socmat( 0,1, -2:2,1, 2) = (/ c_r0, c_r0, c_r0, c_r0, c_r0/)   ! z2u
    socmat( 1,1, -2:2,1, 2) = (/ c_r0,-c_im, c_r0, c_r0, c_r0/)   ! zxu
    socmat( 2,1, -2:2,1, 2) = (/-c_i2, c_r0, c_r0, c_r0, c_r0/)   ! x2u
    ! soc_d up-dw
    !                        xyd,  yzd,  z2d,  zxd,  x2d
    socmat(-2,1, -2:2,2, 2) = (/ c_r0, c_r1, c_r0,-c_im, c_r0/)   ! xyu
    socmat(-1,1, -2:2,2, 2) = (/-c_r1, c_r0,-c_i3, c_r0,-c_im/)   ! yzu
    socmat( 0,1, -2:2,2, 2) = (/ c_r0, c_i3, c_r0,-c_r3, c_r0/)   ! z2u
    socmat( 1,1, -2:2,2, 2) = (/ c_im, c_r0, c_r3, c_r0,-c_r1/)   ! zxu
    socmat( 2,1, -2:2,2, 2) = (/ c_r0, c_im, c_r0, c_r1, c_r0/)   ! x2u
    !
    do ii=-2,2
    do jj=-2,2
        socmat(ii,2, jj,1, 2) = conjg( socmat(jj,1, ii,2, 2) )
    enddo
    enddo
    !
    socmat(:,2, :,2, 2) = -socmat(:,1, :,1, 2)
    !
endsubroutine init_soc

endmodule tghtbd
