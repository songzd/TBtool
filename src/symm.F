!
! This module offers the following subroutines/datas:
!   1. The space group operations.
!   2. Representations (reducible) on the wannier basis.
!   3. The subroutine to symmetrize input hamiltonian.
!   4. The subroutine to calculate little d matrices for integer angular 
!      momentum.
!
module symm
!
use typcnst, only : dp, LTmax, sym_op, c_pi, c_im, Lmax
use my_math, only : dot_prod, mix_prod, vect_len, othmat, unimat
use pntgrp,  only : TonA, make_group, pgrp_name
use reprsnt, only : dmat_rl, dmat_rl_op, dmat_spin
use tghtbd,  only : cell_a, cell_b, natom, atom, h0_dat, h0_lat, h0_lat_num, &
                    h0_LTidx, nobt, nobt2
!
implicit none
!
private
public :: get_sym, cstrct_rep, spgop, nspgop, time_rev, mtpltb, atmap, &
          atmapR, rep_at, reps_at

!
! 18 group generators
real(kind=dp), dimension(3,3,18)  :: op0    ! These generators are complete
character(len=30), dimension(18)  :: op0_n  !   for 32 point groups
!
! space group operation
type(sym_op), dimension(48)::  spgop        ! Space group operations
integer                    :: nspgop        ! Number of space group operations
integer                    :: nspgen        ! Number of generators
integer,dimension(48,48)   :: mtpltb        ! Multiplication table of the point group operations
character(len=30)          :: tau_conv      ! Coorindates convention for input of tau
logical                    :: time_rev      ! Time Reversal 
real(dp)                   :: symtol        ! tolerace of symmetry check
!
namelist /Symmetry/ spgop, tau_conv, time_rev
!
!  Variables to construct representation matrices
integer, allocatable, dimension(:,:)   :: atmap   ! Atoms mapping, gg*aa is the atmap(aa,gg)
integer, allocatable, dimension(:,:,:) :: atmapR  !    atom in the atmapR(:,aa,gg) lattice
real(kind=dp), allocatable, dimension(:,:,:,:,:)  :: rep_at 
complex(kind=dp), allocatable, dimension(:,:,:,:) :: reps_at
! rep_at(:,:,L,aa,gg) is the matrix D_{ap,alphp; a,alph}(gg) in my note, 
!   where, ap={g_R*a + g_T}
!
contains

subroutine get_sym()
    !
    integer :: ii
    character(len=30) :: gname
    !
    call read_sym()
    call make_group(nspgop, spgop, mtpltb, cell_a, cell_b)
    !
    ! Get Point Group Name
    !
    call pgrp_name( nspgop, spgop, gname )
    !
    ! Output ========================================================
    !
    !write(*,'(5X, "Number of generators:          ", I4)') nspgen
    !write(*,'(5X, "Number of symmetry operations: ", I4)') nspgop
    !
    write(*,'(5X, "Point Group:", A9)') trim(gname)
    write(*,'(5X, " iop", " det", "    alph", 2X, "axis (in cart)", &
              15X, "tau (in cart)")')
    do ii=1, nspgop
58      format(5X,2I4, F8.1, '  (',3F8.4, ')','   (',3F8.4,')')
        write(*, 58) ii, spgop(ii)%det, spgop(ii)%thet/c_pi*180, &
                     spgop(ii)%axis, spgop(ii)%tau
    enddo
    !
endsubroutine get_sym
!
!***** cstrct_rep() *************************************************
!
! Used to construct atmap, atmapR, rep_at and reps_at
!   g|RaLis> = \sum_j D(i',i,L,a,g) D(s',s,a,g) |R'a'Li's'>
!
! Here |Rais> is the i-th orbital with orbital momentum L and spinor s
!  at the a-th atom lattice R, a'=atmap, R'=atmapR are the rotated atom 
!  and lattice, respectively.
! rep_at(i',i,L,a,g) and reps_at(s',s,a,g) is the representation
!  matrice of real orbitals and spinor orbitals, respectively
!
!********************************************************************
subroutine cstrct_rep(nspin)
    integer,intent(in) :: nspin
    !
    integer :: aa, bb, aap, ii, jj, LL, gg, La , nn
    !
    real(kind=dp), dimension(3) :: aapos    ! used for atom mapping
    !
    ! D matrix of atoms
    real(kind=dp), dimension(-Lmax:Lmax,-Lmax:Lmax,0:Lmax,natom) :: dmat_at
    complex(kind=dp), dimension(1:nspin,1:nspin,natom)                   ::dmats_at
    !
    ! D matrix of operation
    real(kind=dp), dimension(-Lmax:Lmax,-Lmax:Lmax,0:Lmax) :: dmat_op
    complex(kind=dp), dimension(1:nspin,1:nspin) :: dmats_op
    !
    real(kind=dp), dimension(-Lmax:Lmax,-Lmax:Lmax) :: dmat
    complex(kind=dp), dimension(1:nspin,1:nspin) :: dmats
#if 0 
    ! check D matrix
    !
    real(kind=dp) :: alph, beta, gamm
    !
    alph=0.0_dp;  beta=0.0_dp;  gamm=c_pi/2
    LL=3
    write(*,*) 'L=3'
    call Dmat_rl(alph, beta, gamm, dmat(-LL:LL,-LL:LL), LL)
    write(*,'(7F10.5)') ( dmat(ii,-LL:LL), ii=-LL,LL )
    LL=2
    write(*,*) 'L=2'
    call Dmat_rl(alph, beta, gamm, dmat(-LL:LL,-LL:LL), LL)
    write(*,'(5F10.5)') ( dmat(ii,-LL:LL), ii=-LL,LL )
    LL=1
    write(*,*) 'L=1'
    call Dmat_rl(alph, beta, gamm, dmat(-LL:LL,-LL:LL), LL)
    write(*,'(3F10.5)') ( dmat(ii,-LL:LL), ii=-LL,LL )
    LL=0
    write(*,*) 'L=0'
    call Dmat_rl(alph, beta, gamm, dmat(-LL:LL,-LL:LL), LL)
    write(*,'(1F10.5)') ( dmat(ii,-LL:LL), ii=-LL,LL )
#endif
    !
    ! Local D matrices ==============================================
    !
    dmat_at(:,:,:,:) = 0.0_dp
    dmats_at(:,:,:)  = 0.0_dp
    !
    do aa=1,natom
        do LL=0, Lmax
            call Dmat_RL( atom(aa)%alph, atom(aa)%beta, atom(aa)%gamm, &
                          dmat_at(-LL:LL,-LL:LL, LL, aa), LL)
        enddo
        if (nspin==2) then
            call dmat_spin( atom(aa)%alph, atom(aa)%beta, atom(aa)%gamm, dmats_at(:,:,aa) )
        else
            dmats_at(1,1,aa)=1.0_dp
        endif
    enddo
    !
    ! Atoms mapping and rep_at ===============
    ! 
    !open(unit=100,file='out.rotatm',status='replace')
    !write(100,'(1X,"Rotation on atoms:")')
    !write(100,'(5X,"  iop", "    a", "   ap", "     Rp")')
    !
    open(unit=101,file='out.rep_at',status='replace')
    open(unit=102,file='out.reps_at',status='replace')
    !
    allocate(atmap(natom,nspgop) )
    allocate(atmapR(3,natom,nspgop))
    allocate(rep_at(1:2*Lmax+1, 1:2*Lmax+1, 0:Lmax, natom, nspgop))
    allocate(reps_at(1:nspin,1:nspin,natom,nspgop))
    !
    dmat_op(:,:,:) = 0.0_dp
    dmats_op(:,:)  = 0.0_dp
    atmap(:,:)=0
    atmapR(:,:,:)=0
    rep_at(:,:,:,:,:) = 0.0_dp
    reps_at(:,:,:,:) = 0.0_dp
    !
    do gg=1,nspgop
        !
        ! calculate D-matrix of gg
        do LL=0,Lmax
            call Dmat_RL_op( spgop(gg), dmat_op(-LL:LL,-LL:LL,LL), LL )
        enddo
        !
        if (nspin==2) then
            call dmat_spin(spgop(gg)%alph, spgop(gg)%beta, spgop(gg)%gamm, dmats_op(:,:) )
        else
            dmats_op(1,1)=1.0_dp
        endif
        !
        ! Atom Mapping
        do aa=1,natom
            call TonA( matmul(spgop(gg)%R,atom(aa)%pos) + spgop(gg)%tau, &
                       atmapR(:,aa,gg), aapos, cell_a, cell_b )
            do bb=1,natom
                if ( vect_len( aapos-atom(bb)%pos )<1.0e-7_dp .and. &
                     atom(bb)%elmt==atom(aa)%elmt ) exit
            enddo
            atmap(aa,gg)=bb
            !
            !write(100,'(5X,7I5)') gg, aa, atmap(aa,gg), atmapR(:,aa,gg)
            !
        enddo
        !
        ! rep_at
        do aa=1, natom
        do La=0, Lmax
            nn = atom(aa)%nobtL(La)
            aap = atmap(aa,gg)
            if (nn>0) then
                !
                dmat(-La:La,-La:La) = &
                    matmul(dmat_op(-La:La,-La:La,La),dmat_at(-La:La,-La:La,La,aa))
                dmat(-La:La,-La:La) = & 
                    matmul(transpose(dmat_at(-La:La,-La:La,La,aap)), dmat(-La:La,-La:La))
                do ii=1,nn
                do jj=1,nn
                    rep_at(ii,jj,La,aa,gg) = dmat(atom(aa)%M(ii,La), atom(aa)%M(jj,La))
                enddo
                enddo
                !
                if ( .not. othmat( rep_at(1:nn,1:nn,La,aa,gg), nn ) ) then
                    write(*,*) 'Error in cstrct_rep(): orbitals are not compatible with symmetry !!!'
                    stop
                endif
                !
                ! write to out.rep_at 
                write(101, '("iop=", I4, "  aa=", I4, "  aap=", I4, "  Rp=",3I4, "  ********************")') &
                    gg, aa, aap, atmapR(:,aa,gg)
                write(101, '(4X, "Orbitals on aa ", I4, "- ",I4, " : " )', advance='no') &
                    atom(aa)%idx(1,La), atom(aa)%idx(nn,La)
                do ii=1,nn
                    write(101,'( "(", I1, I2, ") ") ',advance='no') La, atom(aa)%M(ii,La)
                enddo
                write(101,*) ' '
                !
                write(101, '(4X, "Orbitals on aap", I4, "- ",I4, " : " )', advance='no') &
                    atom(aap)%idx(1,La), atom(aap)%idx(nn,La)
                do ii=1,nn
                    write(101,'( "(", I1, I2, ") ") ',advance='no') La, atom(aap)%M(ii,La)
                enddo
                write(101,*) ' '
                !
                do ii=1,nn
                    write(101,'(4X,  5F10.5 )') rep_at(ii,1:nn,La,aa,gg)
                enddo
            endif
            !
            ! reps_at
            !
            dmats(:,:) = matmul(dmats_op(:,:), dmats_at(:,:,aa))
            reps_at(:,:,aa,gg)=matmul(transpose(conjg(dmats_at(:,:,aap))), dmats)
            !
            if ( .not. unimat( reps_at(1:nspin,1:nspin,aa,gg),nspin ) ) then
                write(*,*) 'Error in cstrct_rep(): spins are not compatible with symmetry !!!'
                stop
            endif
            !
            write(102, '("iop=", I4, "  aa=", I4, "  aap=", I4, "  Rp=",3I4, "  ********************")') &
                gg, aa, aap, atmapR(:,aa,gg)
            do ii=1,nspin
                write(102,'(4X,  2("(",2F10.5,")")   )') reps_at(ii,1:nspin,aa,gg)
            enddo
            !
        enddo
        enddo
    enddo
    !
    !close(100)
    close(101)
    close(102)
    !
endsubroutine cstrct_rep



subroutine read_sym()
    !
    integer :: ii, jj, kk, cnt
    real(kind=dp), dimension(3, natom) :: atom_pos1
    integer,       dimension(3) :: ntmp
    !
    ! read_sym ======================================================
    !
    spgop(:)%str='null'
    spgop(:)%det=0
    time_rev = .false.
    !
    open(unit=100, file='system.in', status='old')
    read(100, Symmetry)
    close(100)
    !
    call op0ini()
    nspgop=0
    do
        if (spgop(nspgop+1)%str(1:4)=='null') exit
        !
        nspgop = nspgop+1
        do ii=1,18
            if ( spgop(nspgop)%str==op0_n(ii) ) then
                spgop(nspgop)%R=op0(:,:,ii)
                spgop(nspgop)%det = nint( mix_prod( op0(:,1,ii), op0(:,2,ii), &
                                                                 op0(:,3,ii) ) )
            endif
        enddo
        !
        if (abs(spgop(nspgop)%det)/=1) then
            write(*,*) 'Error in read_sym(): can not find generator !!!'
            stop
        endif
        !
        if (tau_conv(1:1)=='d') then
            spgop(nspgop)%tau = spgop(nspgop)%tau(1)*cell_a(:,1) &
                              + spgop(nspgop)%tau(2)*cell_a(:,2) &
                              + spgop(nspgop)%tau(3)*cell_a(:,3)
        endif
        !
    enddo
    !
    nspgen = nspgop
    !
    write(*,'(1X,"Symmetry:")')
    if (time_rev) then
        write(*,'(5X,"Time-reversal: Yes")')
    else
        write(*,'(5X,"Time-reversal: No")')
    endif
    !
    ! check symm ====================================================
    !
    do kk=1,nspgop
        !
        cnt=0   ! count the rotated atoms that can correspond to original atoms
        do ii=1,natom
            !
            ! calculate the new positions
            atom_pos1(:,ii)=matmul(spgop(kk)%R,atom(ii)%pos)+spgop(kk)%tau
            call TonA( atom_pos1(:,ii), ntmp, atom_pos1(:,ii), cell_a, cell_b)
            !
            ! check the new position
            do jj=1,natom
                if ( vect_len( atom(jj)%pos-atom_pos1(:,ii) )<1.0e-7_dp &
                     .and. atom(jj)%elmt==atom(ii)%elmt) cnt=cnt+1
            enddo
            !
        enddo
        !
        if ( cnt==natom ) then
111         format(5X,'checked generator: ', A6, ', tau=', 3F10.5, ' (in cart)')
            write(*,111) trim(spgop(kk)%str), spgop(kk)%tau 
        else
112         format(1X, 'Error in read_sym(): spgop(', I3, ') is not satisfied !!!')
            write(*,112) kk 
            stop
        endif
        !
    enddo
    !
endsubroutine read_sym
!
subroutine op0ini()
    op0_n(1)='E'
    op0(1,:,1) = (/ 1.0_dp ,   0.0_dp ,   0.0_dp  /)  
    op0(2,:,1) = (/ 0.0_dp ,   1.0_dp ,   0.0_dp  /)
    op0(3,:,1) = (/ 0.0_dp ,   0.0_dp ,   1.0_dp  /)
    !
    op0_n(2)='I'
    op0(1,:,2) = (/-1.0_dp ,   0.0_dp ,   0.0_dp  /)  
    op0(2,:,2) = (/ 0.0_dp ,  -1.0_dp ,   0.0_dp  /)
    op0(3,:,2) = (/ 0.0_dp ,   0.0_dp ,  -1.0_dp  /)
    !
    op0_n(3)='C4x'
    op0(1,:,3) = (/ 1.0_dp ,   0.0_dp ,   0.0_dp  /)
    op0(2,:,3) = (/ 0.0_dp ,   0.0_dp ,  -1.0_dp  /)
    op0(3,:,3) = (/ 0.0_dp ,   1.0_dp ,   0.0_dp  /)
    !
    op0_n(4)='C4y'
    op0(1,:,4) = (/ 0.0_dp ,   0.0_dp ,   1.0_dp  /)
    op0(2,:,4) = (/ 0.0_dp ,   1.0_dp ,   0.0_dp  /)
    op0(3,:,4) = (/-1.0_dp ,   0.0_dp ,   0.0_dp  /)
    !
    op0_n(5)='C4z'
    op0(1,:,5) = (/ 0.0_dp ,  -1.0_dp ,   0.0_dp  /)
    op0(2,:,5) = (/ 1.0_dp ,   0.0_dp ,   0.0_dp  /)
    op0(3,:,5) = (/ 0.0_dp ,   0.0_dp ,   1.0_dp  /)
    !
    op0_n(6)='C2x'
    op0(1,:,6) = (/ 1.0_dp ,   0.0_dp ,   0.0_dp  /)
    op0(2,:,6) = (/ 0.0_dp ,  -1.0_dp ,   0.0_dp  /)
    op0(3,:,6) = (/ 0.0_dp ,   0.0_dp ,  -1.0_dp  /)
    !
    op0_n(7)='C2y'
    op0(1,:,7) = (/-1.0_dp ,   0.0_dp ,   0.0_dp  /)
    op0(2,:,7) = (/ 0.0_dp ,   1.0_dp ,   0.0_dp  /)
    op0(3,:,7) = (/ 0.0_dp ,   0.0_dp ,  -1.0_dp  /)
    !
    op0_n(8)='C2z'
    op0(1,:,8) = (/-1.0_dp ,   0.0_dp ,   0.0_dp  /)
    op0(2,:,8) = (/ 0.0_dp ,  -1.0_dp ,   0.0_dp  /)
    op0(3,:,8) = (/ 0.0_dp ,   0.0_dp ,   1.0_dp  /)
    !
    op0_n(9)='C3xyz'
    op0(1,:,9) = (/ 0.0_dp ,   0.0_dp ,   1.0_dp  /)
    op0(2,:,9) = (/ 1.0_dp ,   0.0_dp ,   0.0_dp  /)
    op0(3,:,9) = (/ 0.0_dp ,   1.0_dp ,   0.0_dp  /)
    !
    op0_n(10)='Mx'
    op0(1,:,10) = (/-1.0_dp ,   0.0_dp ,   0.0_dp  /) 
    op0(2,:,10) = (/ 0.0_dp ,   1.0_dp ,   0.0_dp  /)
    op0(3,:,10) = (/ 0.0_dp ,   0.0_dp ,   1.0_dp  /)
    !
    op0_n(11)='My'
    op0(1,:,11) = (/ 1.0_dp ,   0.0_dp ,   0.0_dp  /) 
    op0(2,:,11) = (/ 0.0_dp ,  -1.0_dp ,   0.0_dp  /)
    op0(3,:,11) = (/ 0.0_dp ,   0.0_dp ,   1.0_dp  /)
    !
    op0_n(12)='Mz'
    op0(1,:,12) = (/ 1.0_dp ,   0.0_dp ,   0.0_dp  /) 
    op0(2,:,12) = (/ 0.0_dp ,   1.0_dp ,   0.0_dp  /)
    op0(3,:,12) = (/ 0.0_dp ,   0.0_dp ,  -1.0_dp  /)
    !
    op0_n(13)='Mxy'
    op0(1,:,13) = (/ 0.0_dp ,  -1.0_dp ,   0.0_dp  /) 
    op0(2,:,13) = (/-1.0_dp ,   0.0_dp ,   0.0_dp  /)
    op0(3,:,13) = (/ 0.0_dp ,   0.0_dp ,   1.0_dp  /)
    !
    op0_n(14)='S4z'
    op0(1,:,14) = (/ 0.0_dp ,  -1.0_dp ,   0.0_dp  /) 
    op0(2,:,14) = (/ 1.0_dp ,   0.0_dp ,   0.0_dp  /)
    op0(3,:,14) = (/ 0.0_dp ,   0.0_dp ,  -1.0_dp  /)
    !
    op0_n(15)='C2xy'
    op0(1,:,15) = (/ 0.0_dp ,   1.0_dp ,   0.0_dp  /) 
    op0(2,:,15) = (/ 1.0_dp ,   0.0_dp ,   0.0_dp  /)
    op0(3,:,15) = (/ 0.0_dp ,   0.0_dp ,  -1.0_dp  /)
    !
    op0_n(16)='C2_xy'
    op0(1,:,16) = (/ 0.0_dp ,  -1.0_dp ,   0.0_dp  /) 
    op0(2,:,16) = (/-1.0_dp ,   0.0_dp ,   0.0_dp  /)
    op0(3,:,16) = (/ 0.0_dp ,   0.0_dp ,  -1.0_dp  /)
    !
    op0_n(17)='C3z'
    op0(1,:,17)  = (/ cos(c_pi*2/3),-sin(c_pi*2/3), 0.0_dp  /)
    op0(2,:,17)  = (/ sin(c_pi*2/3), cos(c_pi*2/3), 0.0_dp  /)
    op0(3,:,17)  = (/ 0.0_dp ,       0.0_dp ,       1.0_dp  /)
    !
    op0_n(18)='C6z'
    op0(1,:,18) = (/ cos(c_pi/3), -sin(c_pi/3),   0.0_dp  /)
    op0(2,:,18) = (/ sin(c_pi/3),  cos(c_pi/3),   0.0_dp  /)
    op0(3,:,18) = (/ 0.0_dp ,      0.0_dp ,       1.0_dp  /)
    !
endsubroutine op0ini

endmodule symm
