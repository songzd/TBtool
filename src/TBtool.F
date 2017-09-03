program TBtool
    !
    use typcnst, only: DP, i12
    use tghtbd,  only: get_tb, nspin
    use symm,    only: get_sym, cstrct_rep
    !use reprsnt, only: dmat_check
    !
    implicit none
    !
    character(len=30)   :: task
    integer             :: lent, stat
    real(kind=dp),parameter :: c_pi=3.141592653589793_dp
    !
    ! Get task ======================================================
    !
    call GET_COMMAND_ARGUMENT(1,task,length=lent,status=stat)
    if (stat>0 .or. lent<2) then
        write(*,*) "Please set the task !!! For example: RPAmain.x band"
        stop
    endif
    !
    ! Get TB Model
    call get_tb()
    !
    ! Get symmtry
    call get_sym()
    call cstrct_rep(nspin)
    !
    select case( task(1:2) )
    case('ba')
        ! Juts plot bands
        call plot_band()
    case('sy')
        ! Impose Symmetry on TB
        call impose_sym_H0()
    case('ir')
        ! Band irreps
        call band_irreps()
    case default
        write(*,'(1X, "Nothing to be done")')
    endselect
    !
#if 0
    !call dmat_check(0.0_dp, 0.0_dp, c_pi/2)
    call irrep_test()
#endif
    !
contains

subroutine plot_band()
    use typcnst,    only: dp, Lmax, lsize, c_pi
    use my_math,    only: vect_len
    use band,       only: bnd_slv, bndspc_ini, bndspc_rls
    use tghtbd,     only: cell_a, cell_b, h0_dat, h0_lat, &
                          h0_lat_num, nobt, nobt2, nspin, natom, atom, nobt
    implicit none
    !
    ! parameters to plot bands
    character(len=30)                   :: k_conv
    integer                             :: nk_path, nklist
    real(kind=dp)                       :: mu_band
    real(kind=dp),dimension(3,lsize)    :: klist
    !
    ! band plot
    integer :: nkpt
    real(kind=dp),allocatable,dimension(:,:):: kpt
    real(kind=dp),allocatable,dimension(:)  :: klen
    real(kind=dp),allocatable,dimension(:,:):: keng
    !
    namelist  /BandPlot/    k_conv, nk_path, mu_band, klist
    !
    integer :: ii, jj
    complex(kind=dp), dimension(nobt2,nobt2) :: wave
    !
    ! initialize
    !
    klist(:,:) = -999
    !
    ! Read ==========================================================
    !
    open(unit=100, file='job.in', status='old')
    read(100, BandPlot)
    close(100)
    !
    nklist=0
    do ii=1,lsize
        if ( klist(1,ii)==-999 ) then
            exit
        else
            nklist = nklist + 1
        endif
    enddo
    if (nklist<2) then
        write(*,*) 'Error in read_band(): too short klist !!!'
        stop
    endif
    !
    if ( k_conv(1:1)=='d' ) then
        do ii=1,nklist
            klist(:,ii) = klist(1,ii)*cell_b(:,1) + klist(2,ii)*cell_b(:,2) &
                        + klist(3,ii)*cell_b(:,3)
        enddo
    endif
    !
    ! Generate kpt and klen =========================================
    !
    nkpt = nk_path*(nklist-1)
    allocate( kpt(3,nkpt), klen(nkpt) )
    !
    do ii=1, nklist-1
        do jj=1,nk_path
            kpt(:, (ii-1)*nk_path+jj) = klist(:,ii) + (klist(:,ii+1)-klist(:,ii)) &
                                                      *(jj-1)/(nk_path-1)
        enddo
    enddo
    !
    klen(:)=0.0_dp
    do ii=2, nkpt
        klen(ii) = klen(ii-1) + vect_len( kpt(:,ii)-kpt(:,ii-1) )
    enddo
    !
    write(*,'(1X, "Band:")')
    write(*,'(5X, "nk per path:", I5)') nk_path
    do ii=1,nklist
        write(*,'(8X, "(", 3F12.6, ")*2pi" )') klist(:,ii)/(2*c_pi)
    enddo
    !
#if 0 
    open(100,file='kpt.dat',status='replace')
    do ii=1, nkpt
        write(100,'(5X, 3F12.6)') kpt(:,ii)/(2*c_pi)
    enddo
    close(100)
#endif
    !
    !
    allocate(keng(nobt2,nkpt))
    !
    write(*,'(5X,"Plotting Bands ... ")')
    call bndspc_ini(nobt2)
    do ii=1, nkpt
        call bnd_slv( kpt(:,ii), wave, keng(:,ii), nobt, nobt2, nspin, h0_lat_num, h0_dat, &
                      h0_lat, natom, atom, cell_a )
    enddo
    call bndspc_rls()
    !
    ! output ========================================================
    !
    open(100,file='band.dat',status='replace')
    !
    do ii=1, nobt2
        do jj=1, nkpt
            write(100,'(5X, 2F20.10)') klen(jj), keng(ii,jj)-mu_band
        enddo
        write(100,*) ' '
    enddo
    !
    do ii=1,nklist-1
        write(100,'(5X, 2F20.10)') klen(ii*nk_path), minval(keng(1,:))-mu_band-0.5_dp
        write(100,'(5X, 2F20.10)') klen(ii*nk_path), maxval(keng(nobt2,:))-mu_band+0.5_dp
        write(100, *) ' '
    enddo
    close(100)
    !
    write(*,'(5X,"done")')
    !
endsubroutine plot_band
!
!
subroutine impose_sym_H0()
    use typcnst,    only: LTmax
    use band,       only: impose_sym
    use tghtbd,     only: h0_lat_num, nobt, nobt2, nspin, h0_dat, h0_lat, h0_LTidx, &
                          natom, atom, cell_a, cell_b
    use symm,       only: nspgop, spgop, time_rev, atmap, atmapR, rep_at, reps_at
    !
    integer :: RR,ii,jj
!    complex(dp),dimension(nobt2,nobt2,h0_lat_num) :: h0_dat_
!    real(dp) :: h0max
!    !
!    ! The input should be a symmetrized hamiltonian
!    h0_dat_(:,:,:)=h0_dat(:,:,:)
!
    !
    write(*,'(1X,"Impose symmetry on H0 ... ")' )
    call impose_sym(h0_lat_num, nobt, nobt2, nspin, h0_dat, h0_lat, LTmax, h0_LTidx, nspgop, spgop, &
                    time_rev, natom, atom, atmap, atmapR, rep_at, reps_at, &
                    cell_a, cell_b)
!
!    write(*,'("Debug imps_sym in impose_sym_H0() ")' )
!    !
!    h0max=0.0_dp
!    do RR=1,h0_lat_num
!    do ii=1,nobt2
!    do jj=1,nobt2
!        if( abs( h0_dat_(ii,jj,RR)-h0_dat(ii,jj,RR) )>h0max ) &
!            h0max = abs( h0_dat_(ii,jj,RR)-h0_dat(ii,jj,RR) )
!    enddo
!    enddo
!    enddo
!    !
!    write(*,'("h0_dat diffs by",E12.5)') h0max
!    !
!
    write(*,*) " done. Symmetrized Hamiltonian is written in hr.dat.sym"
    !
    ! Output ========================================================
    !
    open(100,file='hr.dat.sym',status='replace')
    do RR=1,h0_lat_num
        do jj=1,nobt2
        do ii=1,nobt2
273         format(5I5,2F20.10)
            if ( abs(h0_dat(ii,jj,RR))>1.0e-7_dp ) then
                write(100,273) h0_lat(:,RR), ii, jj, h0_dat(ii,jj,RR)
            endif
        enddo
        enddo
    enddo
    close(100)
    !
endsubroutine

subroutine band_irreps()
    use typcnst,    only: dp, Lmax, lsize, c_pi, sym_op
    use my_math,    only: vect_len, dot_prod
    use band,       only: bnd_slv, bndspc_ini, bndspc_rls
    use tghtbd,     only: cell_a, cell_b, h0_dat, h0_lat, &
                          h0_lat_num, nobt, nobt2, nspin, natom, atom, nobt
    use pntgrp,     only: find_ltgrp, pgrp_name
    use symm,       only: spgop, nspgop
    implicit none
    !
    ! controlling parameters 
    character(len=30)                   :: k_conv
    real(kind=dp),dimension(3,lsize)    :: klist
    !
    namelist  /BandIrrep/    k_conv, klist
    !   
    ! Little Group
    type(sym_op),dimension(48)  :: ltgop        ! Little group operations
    integer                     :: nltgop       ! Number of Little group operations
    integer,dimension(48)       :: ltg2spg      ! Index of little group operations in spgop()
    character(len=30)           :: ltgnm        ! Little group name
    real(dp)                    :: aax(3),cax(3)! Corrdinates for little group
    !
    integer :: ii, jj, kk, ss
    complex(kind=dp), dimension(nobt2,nobt2) :: wave
    real(kind=dp), dimension(nobt2) :: eng
    !
    ! initialize
    !
    klist(:,:) = -999
    !
    ! Read ==========================================================
    !
    open(unit=100, file='job.in', status='old')
    read(100, BandIrrep)
    close(100)
    !
    do kk=1,lsize
        if ( klist(1,kk)==-999 ) then
            exit
        endif
        !
        if ( k_conv(1:1)=='d' ) then
            klist(:,kk) = klist(1,kk)*cell_b(:,1) + klist(2,kk)*cell_b(:,2) &
                        + klist(3,kk)*cell_b(:,3)
        endif
        !
        ! find the little group
        !
        call find_ltgrp(klist(:,kk), nspgop, spgop, cell_a, cell_b, nltgop, ltgop, ltg2spg)
        call pgrp_name(nltgop, ltgop, ltgnm, cax, aax)
        if(dot_prod(cax,klist(:,kk))<0) cax=-cax
        !
        write(*,'(5X, "k = ")')
        write(*,'(5X,"Little group: ", A9)') trim(ltgnm)
        write(*,'(8X, "with caxis=",3F8.4, " aaxis=",3F8.4)') cax, aax
        ss=10
        do ii=1,nltgop/ss
            write(*,'(8X,"iop =", 16I5)') ltg2spg( (ii-1)*ss+1:min(ii*ss,nltgop) )
        enddo
        ii=nltgop/ss
        if ( nltgop>ss*ii ) then
            write(*,'(8X,"iop =", 16I5)') ltg2spg( ss*ii+1:nltgop )
        endif
        
    enddo


endsubroutine band_irreps

endprogram TBtool
