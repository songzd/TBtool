!
module pntgrp
!
use typcnst,    only: dp, c_pi, sym_op
use my_math,    only: vect_len, mix_prod, dot_prod, x_prod, angle
!
implicit none
public
!
contains

subroutine euler_angle(c3,a3,alph,beta,gamm)
    real(dp),dimension(3),intent(in) :: c3, a3   ! a and c axis
    real(dp),intent(out) :: alph, beta, gamm     ! Euler's angles
    !
    real(dp), dimension(3) :: a3_, b3_, c3_, a2_, b2_
    real(dp)    ::  sc, sa
    !   
    sa=vect_len(a3)
    sc=vect_len(c3)
    !   
    if (sa<1e-9_dp .or. sc<1e-10_dp .or. dot_prod(a3,c3)>1e-4_dp) then
        write(*,*) 'Error in euler_angle !!!'
        write(*,'(3F12.6)') a3
        write(*,'(3F12.6)') c3
        stop
    endif
    !   
    a3_=a3/sa
    c3_=c3/sc
    b3_= x_prod(c3_,a3_)
    !   
    beta = acos(c3_(3))
    if (abs(beta)<1.0e-8_dp) then
    ! if c3 || z
        alph=0_dp
        gamm = angle( a3_(1), a3_(2) )
    elseif (abs(beta-c_pi)<1.0e-8_dp) then
    ! if c3 || -z
        alph = angle( b3_(2), -b3_(1) )
        gamm = 0.0_dp
    else
        alph = angle(c3_(1),c3_(2))
        a2_= (/ cos(beta)*cos(alph), cos(beta)*sin(alph), -sin(beta) /)
        b2_= (/-sin(alph), cos(alph), 0.0_dp /)
        gamm = angle( dot_prod(a3_,a2_), dot_prod(a3_,b2_) )
    endif
    !
endsubroutine
!
subroutine TonA( t, tint, frac, A, B ) 
    real(dp),dimension(3),intent(in)   :: t     ! a vector
    real(dp),dimension(3),intent(out)  :: frac  ! t's fractional coordinates on A
    integer,dimension(3),intent(out)   :: tint  ! t's integer coordinates on A
    real(dp),dimension(3,3),intent(in) :: A, B  ! A: lattice, B: reciprocal lattice
                                                ! A*B = 2pi
    !   
    integer :: ii
    real(dp),dimension(3) ::  frac_
    !   
    frac_(1)=dot_prod(t,B(:,1))/(2*c_pi)
    frac_(2)=dot_prod(t,B(:,2))/(2*c_pi)
    frac_(3)=dot_prod(t,B(:,3))/(2*c_pi)
    !   
    do ii=1,3
        if( abs(frac_(ii)-nint(frac_(ii)))<1.0e-4 ) then
            tint(ii)=nint(frac_(ii))
        else
            tint(ii)=floor(frac_(ii))
        endif
    enddo
    !   
    frac_ = frac_ - tint
    !   
    frac(:) = frac_(1)*A(:,1) + frac_(2)*A(:,2) + frac_(3)*A(:,3)
    !   
endsubroutine
!
subroutine get_rotation(R, angle, axis)
    real(dp),intent(in)    :: R(3,3)
    real(dp),intent(out)   :: angle, axis(3)
    ! local variables
    integer    :: i
    ! dgeev variables
    integer    :: info, iaxis
    real(dp)   :: det, mat(3,3), dvl(3,3), dvr(3,3), wi(3), dwork(12), wr(3), arg
    !----
    det = mix_prod(R(:,1),R(:,2),R(:,3)) 
    mat = R*det
    !
    arg=((mat(1,1)+mat(2,2)+mat(3,3)-1.0_dp)*0.5_dp)
    if(arg>1.0_dp)  arg=1.0_dp
    if(arg<-1.0_dp) arg=-1.0_dp
    angle=acos(arg)

    if(abs(abs(angle) - c_pi) .lt. 1e-4) then
       ! angle is 180 deg => can't find the axis the
       ! easy way. Diagonalize rotation matrix and 
       ! pick the eigenvector corresponding to 
       ! unity eigenvalue.
       call DGEEV( 'N', 'V', 3, mat, 3, wr, wi, dvl, 1, &
            dvr,  3,  dwork, 12,info)
       if(info .ne. 0) then
          write(*,*) 'error in dgeev. info=',info
          STOP
       end if
       ! find the axis...just pick the first one with e=1
       iaxis=0
       do i=1,3
          if(abs(wr(i)-1.0_dp) .lt. 1e-9 .and. abs(wi(i)) .lt. 1e-9) then
             iaxis=i
             exit
          end if
       end do
       if(iaxis .lt. 1) then
          write(*,*) 'Error in get_rotation(): could not find rotation axis'
          STOP
       end if
    else if(abs(angle) .gt. 1e-3) then
       ! standard case. See Altmann's book
       dvr(1,1)=mat(3,2)-mat(2,3)
       dvr(2,1)=mat(1,3)-mat(3,1)
       dvr(3,1)=mat(2,1)-mat(1,2)
       dvr=dvr/sin(angle)/2.0_dp
       iaxis=1
    else if(abs(angle) .lt. 1e-4) then
       dvr(1,1)=1.0_dp
       dvr(2,1)=0.0_dp
       dvr(3,1)=0.0_dp
       iaxis=1
    end if
    axis=dvr(:,iaxis)
endsubroutine get_rotation
!
function opmltp( op1, op2, cell_a, cell_b )
    type(sym_op) :: op1, op2, opmltp
    real(dp),dimension(3,3) :: cell_a, cell_b
    !
    integer, dimension(3) :: ntmp
    real(kind=dp), dimension(3) :: t
    !
    opmltp%det = op1%det*op2%det
    !
    t = matmul(op1%R, op2%tau) + op1%tau
    call TonA( t, ntmp, opmltp%tau, cell_a, cell_b )
    !
    opmltp%R = matmul(op1%R, op2%R)
    !
endfunction opmltp

function opeq(op1, op2)
    type(sym_op)  ::  op1, op2
    logical :: opeq
    !
    integer :: ii, jj
    !
    opeq = .true.
    !
    if ( op1%det /= op2%det ) opeq = .false.
    !
    do ii=1,3
        if ( abs(op1%tau(ii)-op2%tau(ii))>1.0e-3_dp ) then
            opeq = .false.
            exit
        endif
    enddo
    !
    do ii=1,3
        if ( .not. opeq) exit
        do jj=1,3
            if ( abs( op1%R(ii,jj)-op2%R(ii,jj) ) >1.0e-3_dp ) then
                opeq = .false.
                exit
            endif
        enddo
    enddo
    !
endfunction opeq
!
subroutine make_group(nop,op,mtb,cell_a,cell_b)
    integer, intent(inout) :: nop 
    type(sym_op), dimension(48), intent(inout) ::  op  
    integer,dimension(48,48), intent(out)      :: mtb 
    real(dp), dimension(3,3), intent(in)       :: cell_a, cell_b
    !
    logical :: new 
    type(sym_op) ::  op1 
    !   
    integer :: ii, jj, kk
    !   
    ii=1
    do while (ii<=nop)
        !   
        jj=1
        do while (jj<=nop)
            !   
            op1=opmltp(op(ii),op(jj),cell_a,cell_b)
            !   
            ! whether it is new
            new = .true.
            do kk=1,nop
                if ( opeq( op(kk),op1)  ) then
                    new = .false.
                    exit
                endif
            enddo
            !
            ! append it to op if new
            if (new) then
                nop = nop+1
                op(nop) = op1
            endif
            !
            jj=jj+1
        enddo
        !
        ii=ii+1
    enddo
    !
    ! Multiplication Table
    mtb(:,:)=0
    do ii=1,nop
    do jj=1,nop
        op1=opmltp(op(ii),op(jj),cell_a,cell_b)
        do kk=1,nop
            if ( opeq(op(kk),op1) ) then
                mtb(ii,jj)=kk
                exit
            endif
        enddo
    enddo
    enddo
    !   
    ! calculate Euler's angles
    !   
    do ii=1, nop
        call euler_angle( op(ii)%det*op(ii)%R(:,3), &
                          op(ii)%det*op(ii)%R(:,1), &
                          op(ii)%alph, op(ii)%beta, op(ii)%gamm)
    enddo
    !   
    ! Get rotation axes and angles
    !   
    do ii=1, nop
        call get_rotation( op(ii)%R, op(ii)%thet, op(ii)%axis )
        op(ii)%rotm = nint(2*c_pi/op(ii)%thet)
    enddo
endsubroutine make_group
!
!
subroutine pgrp_name( nop, op, gname, cx, ax)
    integer,intent(in)              :: nop          ! Number of operations
    type(sym_op),intent(in)         :: op(nop)      ! Operations
    character(len=*),intent(out)    :: gname        ! Name of the Group
    real(dp),optional,intent(out)   :: cx(3), ax(3) ! Principla axis
    !
    real(dp) :: S(3,3,48)
    integer :: ii, cc, aa, bb
    !
    ! Find group name ===============================================
    !
    do ii=1,nop
        S(:,:,ii) = op(ii)%R
    enddo
    call PGROUP_VASP(S,nop,ii,gname)
    !
    if( .not. present(cx) ) return
    if( .not. present(ax) ) return
    !
    ! Find Principal Axis ===========================================
    !
    select case(trim(gname))
    case('C1', 'Ci')
        cx=(/0.0, 0.0, 1.0/)
        ax=(/1.0, 0.0, 0.0/)
    case('T', 'Th', 'Td', 'D2h')
        ! If the point group is T, Th or Td
        !   cx is the 2-fold axis closest to (001)
        !   ax is the 2-fold axis closest to (100)
        cc=0
        do ii=1,nop
        if ( op(ii)%rotm==2 .and. op(ii)%det==1) then
            if( cc==0 ) then
                cc=ii
            elseif( abs(op(ii)%axis(3))>abs(op(cc)%axis(3)) ) then
                cc=ii
            endif
        endif
        enddo
        cx=op(cc)%axis*sign( 1.0_dp, op(cc)%axis(3) )
        !   
        aa=0
        do ii=1,nop
        if ( op(ii)%rotm==2 .and. op(ii)%det==1 .and. &
             abs(dot_prod( op(ii)%axis,cx ))<1.0e-3_dp) then
            if(aa==0) then
                aa=ii
            elseif( abs(op(ii)%axis(1))>abs(op(aa)%axis(1)) ) then
                aa=ii
            endif
        endif
        enddo
        ax=op(aa)%axis*sign( 1.0_dp, op(aa)%axis(1) )
        !
    !case('D2h')
    !    cc=0
    !    do ii=1,nop
    !        if( op(ii)%rotm==2 .and. op(ii)%det==-1) then
    !            cc=ii
    !            exit
    !        endif
    !    enddo
    !    cx = op(cc)%axis*sign(1.0_dp,op(cc)%axis(3))
    !    !
    !    aa=0
    !    do ii=1,nop
    !        if ( op(ii)%rotm==2 .and. &
    !             abs(dot_prod(op(ii)%axis,cx))<1.0e-3_dp ) then
    !            if (aa==0) then
    !                aa=ii
    !            elseif ( abs(op(ii)%axis(1))>abs(op(aa)%axis(1)) ) then
    !                aa=ii
    !            endif
    !            !
    !        endif
    !    enddo
    !    ax=op(aa)%axis*sign(1.0_dp,op(aa)%axis(1))
    !    !
    case('D2d')
        aa=0
        bb=0
        do ii=1,nop
            if ( op(ii)%rotm==2 .and. op(ii)%det==-1 ) then
                if (aa==0) then
                    aa=ii
                elseif (bb==0) then
                    bb=ii
                endif
            endif
        enddo
        !
        if( abs(op(bb)%axis(1))>abs(op(aa)%axis(1)) ) then
            cc=aa
            aa=bb
            bb=cc
        endif
        !
        ax=op(aa)%axis*sign(1.0_dp,op(aa)%axis(1))
        cx=x_prod(ax,op(bb)%axis)
        cx=cx*sign( 1.0_dp, cx(3) )
        !
    case default
        !------------------------------------------------------------
        ! The following rules work for: Oh, O, D6h, D4h, D3h, D3d, D6
        !   D4, D3, D2, S6, S4, C6h, C4h, C3h, C2h, Cs, C6v, C4v, C3v
        !   C2v, C6, C4, C3, C2
        !------------------------------------------------------------
        !
        ! cx:   1. higher symmetry
        !       2. rotation > pseudo rotation
        !       3. closest to 001
        cc=1
        do ii=1,nop
            if (op(ii)%rotm>op(cc)%rotm) then
                cc=ii
            elseif( op(ii)%rotm==op(cc)%rotm ) then
                if( op(ii)%det > op(cc)%det ) then
                    cc=ii
                elseif( op(ii)%det == op(cc)%det .and. &
                    abs(op(ii)%axis(3))>abs(op(cc)%axis(3)) ) then
                    cc=ii
                endif
            endif
        enddo
        cx = op(cc)%axis*sign(1.0_dp,op(cc)%axis(3))
        !
        ! ax: try all (pseudo) symmetry axis perpendicular to cx and consider:
        !       1. rotation > pseudo rotation
        !       2. higher symmetry
        !       3. closest to (100)
        aa=0
        do ii=1,nop
        if ( op(ii)%rotm>1 .and. &
             abs(dot_prod( op(ii)%axis,cx ))<1.0e-3_dp) then
            !
            if (aa==0) then
                ! first set
                aa=ii
            else
                if( op(ii)%det>op(aa)%det ) then
                    ! A rotation has a higher priority than pseudo rotation
                    aa=ii
                elseif(op(ii)%det==op(aa)%det) then
                    if (op(ii)%rotm>op(aa)%rotm) then
                        ! An axis has higher symmetry
                        aa=ii
                    elseif(op(ii)%rotm==op(aa)%rotm .and. &
                           abs(op(ii)%axis(1))>abs(op(aa)%axis(1)) ) then
                        aa=ii
                    endif
                endif
            endif
            !
        endif
        enddo
        !
        if (aa/=0) then
            ax=op(aa)%axis*sign(1.0_dp,op(aa)%axis(1))
        else
            ax=x_prod( (/0.0_dp, 1.0_dp, 0.0_dp/), cx )
            if ( vect_len(ax)>1.0e-2_dp ) then
                ax=ax/vect_len(ax)
            else
                ax=x_prod( (/1.0_dp, 0.0_dp, 0.0_dp/), cx )
                ax=ax/vect_len(ax)
            endif
        endif
        !
    endselect
    !
endsubroutine pgrp_name

!**** find_ltgrp_r() **************************************************
!
! Return the little group G(r),  which is defined as
!   (1)  \forall g=(R,tau) \in G(r), R*r+tau=r (module a reciprocal lattice)
!
!********************************************************************
subroutine find_ltgrp_r( rvec, nop, op, cell_a, cell_b, nlop, lop, lop2op)
    real(dp),intent(in)     :: rvec(3)      ! vector
    integer,intent(in)      :: nop          ! number of space group operations
    type(sym_op),intent(in) :: op(:)        ! space group operations
    real(dp),intent(in)     :: cell_a(3,3)  ! lattice
    real(dp),intent(in)     :: cell_b(3,3)  ! reciprocal lattice, A*B=2pi
    integer,intent(out)     :: nlop         ! number of little group operations
    type(sym_op),intent(out):: lop(:)       ! little group operations
    integer,intent(out)     :: lop2op(:)    ! lop's index in op 
    !
    integer :: ii, ntmp(3)
    real(dp) :: rlen, rtmp(3),ftmp(3)
    real(dp) :: rtol=1e-3_dp
    !
    ! if r=(000)
    rlen = vect_len(rvec)
    if ( rlen<1.0e-10_dp ) then
        nlop=nop
        lop(1:nop) =op(1:nop)
        do ii=1,nop
            lop2op(ii)=ii
        enddo
        return
    endif
    !
    nlop=0
    do ii=1,nop
        rtmp = matmul(op(ii)%R,rvec) + op(ii)%tau - rvec
        call TonA(rtmp,ntmp,ftmp,cell_a,cell_b)
        !
        if( vect_len(ftmp)/rlen<rtol ) then
            nlop = nlop+1
            lop(nlop) = op(ii)
            lop2op(nlop) = ii
        endif
    enddo
endsubroutine find_ltgrp_r
!
!**** find_ltgrp_k() **************************************************
!
! Return the relevant little group G(q),  which is defined as
!   (1)  \forall g=(R,tau) \in G(q), R*q=q (module a reciprocal lattice)
!   (2)  \forall g=(R,tau), g'=(R',tau') \in G(q)
!        there is  (Rq-q)*tau'=2*pi*n, n=0,\pm 1,\pm 2 ...
!
!********************************************************************
subroutine find_ltgrp_k( qvec, nop, op, cell_a, cell_b, nlop, lop, lop2op)
    real(dp),intent(in)     :: qvec(3)      ! q-vector
    integer,intent(in)      :: nop          ! number of space group operations
    type(sym_op),intent(in) :: op(:)        ! space group operations
    real(dp),intent(in)     :: cell_a(3,3)  ! lattice
    real(dp),intent(in)     :: cell_b(3,3)  ! reciprocal lattice, A*B=2pi
    integer,intent(out)     :: nlop         ! number of little group operations
    type(sym_op),intent(out):: lop(:)       ! little group operations
    integer,intent(out)     :: lop2op(:)    ! lop's index in op 
    !
    integer :: ii, jj, gint(3)
    logical :: acpt
    real(dp) :: qlen, rr, gtmp(3),gtmp2(3), ftmp(3)
    real(dp) :: qtol=1e-3_dp
    !
    ! if q=(000)
    qlen = vect_len(qvec)
    if ( qlen<1.0e-10_dp ) then
        nlop=nop
        lop(1:nop) =op(1:nop)
        do ii=1,nop
            lop2op(ii)=ii
        enddo
        return
    endif
    !
    ! Normal Operations
    nlop=0
    do ii=1,nop
        if( vect_len(matmul(op(ii)%R,qvec)-qvec)/qlen<qtol ) then
            nlop = nlop+1
            lop(nlop) = op(ii)
            lop2op(nlop) = ii
        endif
    enddo
    !
    ! Umklapp Operations 
    do ii=1,nop
        gtmp = matmul(op(ii)%R,qvec)-qvec
        if( vect_len(gtmp)/qlen<qtol ) cycle
        !
        ! Get an umklapp operation
        call TonA(gtmp,gint,ftmp,cell_b,cell_a)
        if( vect_len(ftmp)/qlen<qtol ) then
            !
            ! Whether it satisfy relevant irrep condition
            acpt=.true.
            jj=0
            do while(jj<nlop)
                jj=jj+1
                rr=dot_prod(gtmp,lop(jj)%tau)/(2.0_dp*c_pi)
                if ( abs(rr-nint(rr))>1.0e-3 ) then
                    acpt=.false.
                    exit
                endif
                !
                gtmp2=matmul(lop(jj)%R,qvec)-qvec
                rr=dot_prod(gtmp2,op(ii)%tau)/(2.0_dp*c_pi)
                if ( abs(rr-nint(rr))>1.0e-3 ) then
                    acpt=.false.
                    exit
                endif
            enddo
            !
            if (acpt) then
                nlop = nlop+1
                lop(nlop) = op(ii)
                lop2op(nlop) = ii
            endif
            !
        endif
    enddo
endsubroutine find_ltgrp_k
!***********************************************************************
!                                                                      *
      SUBROUTINE PGROUP_VASP(S,NROT,PGIND,SFNAME)
!                                                                      *
!***********************************************************************
!                                                                      *
!   Given some set of point symmetry operations routine PGROUP returns *
!   the name of the point group defined by these symmetry operations.  *
!                                                                      *
!                                                                      *
!   Input parameters:                                                  *
!   -----------------                                                  *
!                                                                      *
!      S(3,3,48) contains the INTEGER rotation matrices.               *
!      NROT contains the number of symmetry operations.                *
!                                                                      *
!                                                                      *
!   Output parameters:                                                 *
!   ------------------                                                 *
!                                                                      *
!      PGIND gives the "name" of the group defined by following key:   *
!       1 --> C_1       9 --> C_3      17 --> D_4      25 --> C_6v     *
!       2 --> S_2      10 --> S_6      18 --> C_4v     26 --> D_3h     *
!       3 --> C_2      11 --> D_3      19 --> D_2d     27 --> D_6h     *
!       4 --> C_1h     12 --> C_3v     20 --> D_4h     28 --> T        *
!       5 --> C_2h     13 --> D_3d     21 --> C_6      29 --> T_h      *
!       6 --> D_2      14 --> C_4      22 --> C_3h     30 --> O        *
!       7 --> C_2v     15 --> S_4      23 --> C_6h     31 --> T_d      *
!       8 --> D_2h     16 --> C_4h     24 --> D_6      32 --> O_h      *
!      SFNAME is the explicit name in form of a string (Schoenflies).  *
!                                                                      *
!                                                                      *
!***********************************************************************

      real(dp) :: S(3,3,48)
      INTEGER NROT,PGIND
      CHARACTER (LEN=*) SFNAME

      !
      CHARACTER (4) GNAME(32)
      SAVE GNAME
      DATA GNAME / 'C1', 'Ci', 'C2', 'Cs','C2h', 'D2','C2v', &
     &            'D2h', 'C3', 'S6', 'D3','C3v','D3d', 'C4', &
     &             'S4','C4h', 'D4','C4v','D2d','D4h', 'C6', &
     &            'C3h','C6h', 'D6','C6v','D3h','D6h',  'T', &
     &             'Th',  'O', 'Td', 'Oh'/
     !
     INTEGER :: TRACE,DET, INVERS, NC2,NC3,NC4,NC6,NS2,NS3,NS4,NS6,IR

! Need at least four characters to store the group name ... :
      IF (LEN(SFNAME)<4) then
        write(*,*)'Error in PGROUP(): Variable SFNAME declared too short in  &
                   the calling program!!!'
        stop    
      ENDIF
! Trivial case: Only one symmetry operation (can only be E --> PGIND=1):
      IF (NROT==1) THEN
         PGIND=1
         SFNAME=GNAME(PGIND)
         RETURN
      END IF
! There is an other trivial case: group O_h (PGIND=32), because it is
! the only group having 48 elements ... !
      IF (NROT==48) THEN
         PGIND=32
         SFNAME=GNAME(PGIND)
         RETURN
      END IF
! An other trivial case is group D_4h (PGIND=20), because it is the only
! group having 16 elements ... !
      IF (NROT==16) THEN
         PGIND=20
         SFNAME=GNAME(PGIND)
         RETURN
      END IF
! And finally there is a fourth trivial case: it is group C_3 (PGIND=9),
! because it is the only group having 3 elements ... !
      IF (NROT==3) THEN
         PGIND=9
         SFNAME=GNAME(PGIND)
         RETURN
      END IF
! All other groups need further investigations and detailed analysis ...
! First determine the type of elements and count them. Possible elements
! are E, I, C_2,3,4,6 and S_2,3,4,6 (S_2 = m), E is trivial and always
! present. The type of a symmetry operation can be identified simply by
! calculating the trace and the determinant of the rotation matrix. The
! combination of these two quantities is specific for specific elements:
!
! Element:         E    I  C_2  C_3  C_4  C_6  S_2  S_6  S_4  S_3
! Trace:          +3   -3   -1    0   +1   +2   +1    0   -1   -2
! Determinant:    +1   -1   +1   +1   +1   +1   -1   -1   -1   -1

      INVERS=0
      NC2=0
      NC3=0
      NC4=0
      NC6=0
      NS2=0
      NS6=0
      NS4=0
      NS3=0
      DO 1 IR=1,NROT
         TRACE=nint(S(1,1,IR)+S(2,2,IR)+S(3,3,IR))
! Found unity operator (trivial):
         IF (TRACE==3) GOTO 1
! Found inversion ...
         IF (TRACE==(-3)) THEN
            INVERS=1
            GOTO 1
         END IF
         DET=nint( S(1,1,IR)*(S(2,2,IR)*S(3,3,IR)-S(2,3,IR)*S(3,2,IR))+ &
     &             S(1,2,IR)*(S(2,3,IR)*S(3,1,IR)-S(2,1,IR)*S(3,3,IR))+ &
     &             S(1,3,IR)*(S(2,1,IR)*S(3,2,IR)-S(2,2,IR)*S(3,1,IR)) )
! Found C_2:
         IF ((TRACE==(-1)).AND.(DET==1)) NC2=NC2+1
! Found S_2:
         IF ((TRACE==1).AND.(DET==(-1))) NS2=NS2+1
! Found C_3:
         IF ((TRACE==0).AND.(DET==1)) NC3=NC3+1
! Found S_6:
         IF ((TRACE==0).AND.(DET==(-1))) NS6=NS6+1
! Found C_4:
         IF ((TRACE==1).AND.(DET==1)) NC4=NC4+1
! Found S_4:
         IF ((TRACE==(-1)).AND.(DET==(-1))) NS4=NS4+1
! Found C_6:
         IF ((TRACE==2).AND.(DET==1)) NC6=NC6+1
! Found S_3:
         IF ((TRACE==(-2)).AND.(DET==(-1))) NS3=NS3+1
    1 CONTINUE
! Now we know which elements we have and so we know the group ... :
      IF (NROT==2) THEN
! Groups with 2 elements:
         IF (INVERS==1) THEN
! Contains inversion --> S_2 (PGIND=2):
            PGIND=2
            SFNAME=GNAME(PGIND)
            RETURN
         END IF
         IF (NC2==1) THEN
! Contains twofold rotation --> C_2 (PGIND=3):
            PGIND=3
            SFNAME=GNAME(PGIND)
            RETURN
         END IF
         IF (NS2==1) THEN
! Contains mirror plane --> C_1h (PGIND=4):
            PGIND=4
            SFNAME=GNAME(PGIND)
            RETURN
         END IF
      END IF
      IF (NROT==4) THEN
! Groups with 4 elements:
         IF (INVERS==1) THEN
! Contains inversion --> C_2h (PGIND=5):
            PGIND=5
            SFNAME=GNAME(PGIND)
            RETURN
         END IF
         IF (NC2==3) THEN
! Contains three twofold rotations --> D_2 (PGIND=6):
            PGIND=6
            SFNAME=GNAME(PGIND)
            RETURN
         END IF
         IF (NS2==2) THEN
! Contains two mirror planes --> C_2v (PGIND=7):
            PGIND=7
            SFNAME=GNAME(PGIND)
            RETURN
         END IF
         IF (NC4==1) THEN
! Contains fourfold rotation --> C_4 (PGIND=14):
            PGIND=14
            SFNAME=GNAME(PGIND)
            RETURN
         END IF
         IF (NS4==2) THEN
! Contains fourfold improper rotation --> S_4 (PGIND=15):
            PGIND=15
            SFNAME=GNAME(PGIND)
            RETURN
         END IF
      END IF
      IF (NROT==6) THEN
! Groups with 6 elements:
         IF (INVERS==1) THEN
! Contains inversion --> S_6 (PGIND=10):
            PGIND=10
            SFNAME=GNAME(PGIND)
            RETURN
         END IF
         IF (NC2==3) THEN
! Contains three twofold rotations --> D_3 (PGIND=11):
            PGIND=11
            SFNAME=GNAME(PGIND)
            RETURN
         END IF
         IF (NS2==3) THEN
! Contains three mirror planes --> C_3v (PGIND=12):
            PGIND=12
            SFNAME=GNAME(PGIND)
            RETURN
         END IF
         IF (NC2==1) THEN
! Contains only one twofold rotations --> C_6 (PGIND=21):
            PGIND=21
            SFNAME=GNAME(PGIND)
            RETURN
         END IF
         IF (NS2==1) THEN
! Contains only one mirror plane --> C_3h (PGIND=22):
            PGIND=22
            SFNAME=GNAME(PGIND)
            RETURN
         END IF
      END IF
      IF (NROT==8) THEN
! Groups with 8 elements:
         IF (NS2==3) THEN
! Contains three mirror planes --> D_2h (PGIND=8):
            PGIND=8
            SFNAME=GNAME(PGIND)
            RETURN
         END IF
         IF (NS2==1) THEN
! Contains one mirror planes --> C_4h (PGIND=16):
            PGIND=16
            SFNAME=GNAME(PGIND)
            RETURN
         END IF
         IF (NS2==0) THEN
! Contains no mirror planes --> D_4 (PGIND=17):
            PGIND=17
            SFNAME=GNAME(PGIND)
            RETURN
         END IF
         IF (NS2==4) THEN
! Contains four mirror planes --> C_4v (PGIND=18):
            PGIND=18
            SFNAME=GNAME(PGIND)
            RETURN
         END IF
         IF (NS2==2) THEN
! Contains two mirror planes --> D_2d (PGIND=19):
            PGIND=19
            SFNAME=GNAME(PGIND)
            RETURN
         END IF
      END IF
      IF (NROT==12) THEN
! Groups with 12 elements:
         IF (NS2==3) THEN
! Contains three mirror planes --> D_3d (PGIND=13):
            PGIND=13
            SFNAME=GNAME(PGIND)
            RETURN
         END IF
         IF (NS2==1) THEN
! Contains one mirror planes --> C_6h (PGIND=23):
            PGIND=23
            SFNAME=GNAME(PGIND)
            RETURN
         END IF
         IF (NC2==7) THEN
! Contains seven twofold rotations --> D_6 (PGIND=24):
            PGIND=24
            SFNAME=GNAME(PGIND)
            RETURN
         END IF
         IF (NS2==6) THEN
! Contains six mirror planes --> C_6v (PGIND=25):
            PGIND=25
            SFNAME=GNAME(PGIND)
            RETURN
         END IF
         IF (NS2==4) THEN
! Contains four mirror planes --> D_3h (PGIND=26):
            PGIND=26
            SFNAME=GNAME(PGIND)
            RETURN
         END IF
         IF (NC3==8) THEN
! Contains eight threefold rotations --> T (PGIND=28):
            PGIND=28
            SFNAME=GNAME(PGIND)
            RETURN
         END IF
      END IF
      IF (NROT==24) THEN
! Groups with 24 elements:
         IF (NC6==2) THEN
! Contains two sixfold rotations --> D_6h (PGIND=27):
            PGIND=27
            SFNAME=GNAME(PGIND)
            RETURN
         END IF
         IF (INVERS==1) THEN
! Contains inversion --> T_h (PGIND=29):
            PGIND=29
            SFNAME=GNAME(PGIND)
            RETURN
         END IF
         IF (NC4==6) THEN
! Contains six fourfold rotations --> O (PGIND=30):
            PGIND=30
            SFNAME=GNAME(PGIND)
            RETURN
         END IF
         IF (NS4==6) THEN
! Contains six fourfold improper rotations --> T_d (PGIND=31):
            PGIND=31
            SFNAME=GNAME(PGIND)
            RETURN
         END IF
      END IF
! Ready!
      RETURN
      ENDSUBROUTINE

endmodule
