program TBtool
    !
    use typcnst,    only: DP, i12
    use tghtbd,     only: get_obt, get_hr, nspin
    use symm,       only: get_sym, cstrct_rep
    !
    implicit none
    !
    character(len=30)       :: task
    integer                 :: lent, stat
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
    call get_obt()
    call get_hr()
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
    case('sl')
        ! slice-wannier
        call build_swan()
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
    else
        klist = klist*2*c_pi
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
    use my_math,    only: vect_len
    use band,       only: impose_sym
    use tghtbd,     only: h0_lat_num, nobt, nobt2, nspin, obt, h0_dat, h0_lat, h0_LTidx, &
                          natom, atom, cell_a, cell_b
    use symm,       only: nspgop, spgop, time_rev, atmap, atmapR, rep_at, reps_at
    !
    ! controlling parameters
    real(kind=dp) ::  cut_dist, cut_eng
    namelist /ImposeSymm/    cut_dist, cut_eng
    !
    integer :: RR,ii,jj
    real(kind=dp), dimension(3) :: latt
    !
!    complex(dp),dimension(nobt2,nobt2,h0_lat_num) :: h0_dat_
!    real(dp) :: h0max
!    !
!    ! The input should be a symmetrized hamiltonian
!    h0_dat_(:,:,:)=h0_dat(:,:,:)
!
    write(*,'(1X,"Impose symmetry on H0 ... ")' )
    !
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
    open(unit=100, file='job.in', status='old')
    read(100, ImposeSymm)
    close(100)
    !
    open(100,file='hr.dat.sym',status='replace')
    do RR=1,h0_lat_num
        latt = cell_a(:,1)*h0_lat(1,RR) + cell_a(:,2)*h0_lat(2,RR) + cell_a(:,3)*h0_lat(3,RR)
        do jj=1,nobt2
        do ii=1,nobt2
273         format(5I5,2F20.10)
            if ( abs(h0_dat(ii,jj,RR))>=cut_eng .and. &
                 vect_len( latt+obt(jj)%pos-obt(ii)%pos )<cut_dist ) then
                write(100,273) h0_lat(:,RR), ii, jj, h0_dat(ii,jj,RR)
            endif
        enddo
        enddo
    enddo
    close(100)
    !
endsubroutine

subroutine band_irreps()
    use typcnst,    only: dp, Lmax, lsize, c_pi, sym_op, irep_type
    use my_math,    only: vect_len, dot_prod
    use band,       only: bnd_slv, bndspc_ini, bndspc_rls, get_sewing
    use tghtbd,     only: cell_a, cell_b, h0_dat, h0_lat, &
                          h0_lat_num, nobt, nobt2, nspin, natom, atom, nobt
    use pntgrp,     only: find_ltgrp_k, pgrp_name
    use symm,       only: spgop, nspgop, atmap, atmapR, rep_at, reps_at
    use reprsnt,    only: init_pntgrp_bst, set_pntgrp_irep
    !
    implicit none
    !
    ! controlling parameters 
    character(len=30)   :: k_conv
    real(dp)            :: deg_tol, klist(3,lsize)
    !
    namelist  /BandIrrep/    k_conv, deg_tol, klist
    !   
    ! Little Group
    type(sym_op),dimension(48)  :: ltgop        ! Little group operations
    integer                     :: nltgop       ! Number of Little group operations
    integer,dimension(48)       :: ltg2spg      ! Index of little group operations in spgop()
    character(len=30)           :: ltgnm        ! Little group name
    real(dp)                    :: aax(3),cax(3)! Corrdinates for little group
    !
    ! table of irreps
    type(irep_type), dimension(20) :: irep
    integer :: nirep
    !
    ! Bloch representations
    complex(dp)     :: sewing(nobt2,nobt2,48), Brep(nobt2,48), nxitmp
    integer         :: nxi, checkdim

    !
    integer :: ii, jj, jj_, kk, ss, xi, gg
    real(kind=dp), dimension(3) :: kvecd, kvecc
    complex(kind=dp), dimension(nobt2,nobt2) :: wave
    real(kind=dp), dimension(nobt2) :: eng
    !
    ! initialize
    !
    klist(:,:) = -999
    write(*,*) 'BandIrreps:'
    !
    ! Read ==========================================================
    !
    open(unit=100, file='job.in', status='old')
    read(100, BandIrrep)
    close(100)
    !
    call bndspc_ini(nobt2)
    do kk=1,lsize
        if ( klist(1,kk)==-999 ) exit
        !
        if ( k_conv(1:1)=='d' ) then
            kvecd = klist(:,kk)
            kvecc = kvecd(1)*cell_b(:,1) + kvecd(2)*cell_b(:,2) + kvecd(3)*cell_b(:,3)
        else
            kvecc = klist(:,kk)*2*c_pi
            kvecd = matmul(transpose(cell_a),kvecc)/(2*c_pi)
        endif
        !
        ! find the little group (relevant point group)
        !
        call find_ltgrp_k(kvecc, nspgop, spgop, cell_a, cell_b, nltgop, ltgop, ltg2spg)
        call pgrp_name(nltgop, ltgop, ltgnm, cax, aax)
        !if(dot_prod(cax,klist(:,kk))<0) cax=-cax
        !
        write(*,'(4X, "********************************************************&
                       *********************************")')
        write(*,'(6X, "k = (", 3F12.6, ")*2pi" )') kvecc/(2*c_pi)
        write(*,'(6X, "    (", 3F12.6, ") in cell_b" )') kvecd
        write(*,'(6X,"Little group: ", A9)') trim(ltgnm)
        write(*,'(8X, "with caxis=",3F8.4, " aaxis=",3F8.4)') cax, aax
        ss=10
        do ii=1,nltgop/ss
            write(*,'(8X,"iop =", 16I5)') ltg2spg( (ii-1)*ss+1:min(ii*ss,nltgop) )
        enddo
        ii=nltgop/ss
        if ( nltgop>ss*ii ) then
            write(*,'(8X,"iop =", 16I5)') ltg2spg( ss*ii+1:nltgop )
        endif
        !
        ! Set the point group irreps
        !
        call init_pntgrp_bst()
        call set_pntgrp_irep(nltgop,ltgop,ltgnm,nspin,cax,aax,nirep,irep)
        !
        write(*,'(6X, "Character table:")')
        do xi=1,nirep
            ss=4
            write(*,'(T9,A5,T1)',advance='no') irep(xi)%name
            do ii=1,nltgop/ss
                write(*,'(T15, 4("(", 2F8.4, ") ") )') irep(xi)%cha( (ii-1)*ss+1:min(ii*ss,nltgop) )
            enddo
            ii=nltgop/ss
            if ( nltgop>ss*ii ) then
                write(*,'(T15, 4("(", 2F8.4, ") ") )') irep(xi)%cha( ii*ss+1:nltgop )
                !write(*,'(8X,"iop =", 16I5)') ltg2spg( ss*ii+1:nltgop )
            endif
        enddo
        !
        ! calculate the representations of Bloch states
        !
        ! sewing matrice
        call get_sewing(nobt, nobt2, nspin, nltgop, ltgop, ltg2spg, natom, &
                        atom, atmap, atmapR, rep_at, reps_at, cell_a, kvecc, sewing)
        !
        ! diagonalize the hamiltonian
        call bnd_slv( kvecc, wave, eng, nobt, nobt2, nspin, h0_lat_num, h0_dat, &
                      h0_lat, natom, atom, cell_a )
        !
        ! Band representation
        Brep = 0.0_dp
        do gg=1,nltgop
            do ii=1,nobt2
                do jj =1,nobt2
                do jj_=1,nobt2
                    Brep(ii,gg) = Brep(ii,gg) + conjg(wave(jj_,ii))*sewing(jj_,jj,gg)*wave(jj,ii)
                enddo
                enddo
            enddo
        enddo
        !
        ! get irreps
        write(*,*) ' '
        write(*,'(8X, "Bands", 9X, "Energy", 10X, "Irreps")')
        ii=1
        do while (ii<=nobt2)
            !
            ! find the degenerate bands
            do jj=ii,nobt2
                if ( abs(eng(ii)-eng(jj))>deg_tol ) exit
            enddo
            write(*,'(8X, I4, " -", I4, 4X, F9.4, "  :")',advance='no') ii, jj-1, eng(ii)
            !
            ! get irreps
            checkdim = 0
            do xi=1,nirep
                nxitmp = 0.0_dp
                do gg=1,nltgop
                    nxitmp = nxitmp + conjg(sum(Brep(ii:jj-1,gg)))*irep(xi)%cha(gg)
                enddo
                nxi = nint(real(nxitmp)/nltgop)
                checkdim = checkdim + nxi*irep(xi)%D
                !
                if (nxi>0) then
                    write(*,'(2X, I4, "*", A5)',advance='no') nxi,  irep(xi)%name
                endif
            enddo
            write(*,*) ' '
            !
            ! check irreps dimension 
            if (checkdim/=jj-ii) then
                write(*,*) 'Error in band_irreps(): the dimensions of degenerate bands and &
                            irreps do not match !!!'
                stop
            endif
            !
            ! next band
            ii=jj
            !
        enddo
    enddo
    !
    call bndspc_rls()
    !
endsubroutine band_irreps

subroutine build_swan()
    !
    use typcnst,    only: dp, LTmax
    use tghtbd,     only: cell_a, cell_b, nat=>natom, at=>atom, nobt, nspin, nobt2, &
                          nlatt=>h0_lat_num, latt=>h0_lat, LTidx=>h0_LTidx
    use symm,       only: nop=>nspgop, op=>spgop, time_rev, atmap, atmapR, rep_at, reps_at
    use slice,      only: read_swan, cstrct_rep_slice, nsat, sat, maxrep_sat, nsproj, satmap, &
                          satmapR, rep_sat, Ut_sat, dist_cut
    !
    use slicewann,  only: random_ovlp, symm_ovlp
    !
    complex(dp)         :: ovlp(nsproj,nobt2,nlatt)
    !
    ! get slice projector information ===============================
    !
    call read_swan()
    call cstrct_rep_slice()
    !
    ! set random ovlp ===============================================
    !
    call random_ovlp( nsat, sat, maxrep_sat, nsproj, nat, at, nobt, nobt2, &
                      nspin, nlatt, latt, cell_a, dist_cut, ovlp)
    !



endsubroutine build_swan

endprogram TBtool
