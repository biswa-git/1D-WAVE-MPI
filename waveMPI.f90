!!  MODULE par

module par
    implicit none
    
    integer (kind=4), parameter         :: nx      = 664
    double precision, parameter         :: pi      = 4.0d0*datan(1.0d0)
!!  GEOMETRY INFO
    double precision, parameter         :: x_start = 0.0d0
    double precision, parameter         :: x_end   = 10.0d0
    double precision, parameter         :: length  = x_end - x_start
    double precision, parameter         :: dx      = length/dble(nx)
    double precision, parameter         :: idx     = 1.0d0 /dx

!!  WAVE PACKET INFO
    double precision, parameter         :: kh      = 0.8382420d0
    double precision, parameter         :: k       = kh/dx

!!  WAVE PHSICAL PROPERTY
    double precision, parameter         :: c       = 0.50d0
    double precision, parameter         :: CFL     = 0.020d0
    double precision, parameter         :: dt      = CFL*dx/c
    double precision, parameter         :: t_max   = 1000.0d0
!!  RUNGE-KUTTA COEFFICIENTS
    double precision, parameter, dimension(1:4) ::  &
    &   aj =(/                                      &
    &            0.0d0        ,                     &
    &           -5.0d0 / 9.0d0,                     &
    &           -1.0d0        ,                     &
    &           -33.0d0/25.0d0                      &
    &       /)

    double precision, parameter, dimension(1:4) ::  &
    &   bj =(/                                      &
    &            1.0d0 / 9.0d0,                     &
    &            3.0d0 / 4.0d0,                     &
    &            2.0d0 / 5.0d0,                     &
    &            5.0d0 / 4.0d0                      &
    &       /)

end module par
!-----------------------------------------------------------------------------------------------!



!!  MODULE var : CONTAINS VARIABLES

module var
    use par
    implicit none

    double precision, allocatable, dimension(:) :: x, u, du, u_new, dudx
    double precision                            :: time
    integer (kind=4)                            :: iter

end module var
!-----------------------------------------------------------------------------------------------!



!!  MODULE oucs3_const : CONTAINS THE INFORMATION OF COEFFICIENTS OF OUCS3 STECILS

module oucs3_const
    use par, only: idx
    implicit none

    double precision, parameter                     ::  DD  =    0.37938949120d0
    double precision, parameter                     ::  EE  =    1.57557379000d0
    double precision, parameter                     ::  FF  =    0.18320519200d0
    double precision, parameter                     ::  eta =   -0.0d0, phi = 0.0d0

    double precision, parameter, dimension(-1:1)    ::  &
    &   p = (/                                          &
    &           DD - eta/60.0d0,                        &
    &           1.0d0          ,                        &
    &           DD + eta/60.0d0                         &
    &       /)

    double precision, parameter, dimension(-2:2)    ::  &
    &   q = (/                                          &
    &          -FF *0.250d0 + eta/300.0d0,              &
    &          -EE *0.500d0 + eta/030.0d0,              &
    &          -eta*(11.0d0/150.0d0)     ,              &
    &           EE *0.500d0 + eta/030.0d0,              &
    &           FF *0.250d0 + eta/300.0d0               &
    &       /)

    double precision, parameter, dimension(-4:4)    ::  &
    &   r = (/                                          &
    &          idx*( 1.0d0/280.0d0 + phi*1.00d0),       &
    &          idx*(-4.0d0/105.0d0 - phi*8.00d0),       &
    &          idx*( 1.0d0/5.0d0   + phi*28.0d0),       &
    &          idx*(-4.0d0/5.0d0   - phi*56.0d0),       &
    &          idx*( 0.0d0         + phi*70.0d0),       &
    &          idx*( 4.0d0/5.0d0   - phi*56.0d0),       &
    &          idx*(-1.0d0/5.0d0   + phi*28.0d0),       &
    &          idx*( 4.0d0/105.0d0 - phi*8.00d0),       &
    &          idx*(-1.0d0/280.0d0 + phi*1.00d0)        &
    &       /)

end module oucs3_const
!-----------------------------------------------------------------------------------------------!



!! MPI VARIABLES

module MPI_VAR
    use MPI
    implicit none

!!  MPI        
!---------------------------------------------------------------!
    integer (kind=4), parameter     :: ROOT = 0
    integer (kind=4)                :: MPI_ERR, MPI_COMM_CART
    integer (kind=4)                :: NUMPROC
    integer (kind=4)                :: RANK
    integer (kind=4)                :: STATUS(MPI_STATUS_SIZE)
!---------------------------------------------------------------!

!! DATA RELATED TO TOPOLOGY
!---------------------------------------------------------------!
    integer (kind=4), parameter     :: NDIM = 1, CART_NX = 4
    logical                         :: PERIOD(1:NDIM) = .true.
    integer (kind=4)                :: CART  (1:NDIM)
    integer (kind=4)                :: COORDS(1:NDIM),COORDS_QUERY(1:NDIM)

    integer (kind=4)                :: NEIGHBOUR_EAST 
    integer (kind=4)                :: NEIGHBOUR_WEST

    logical                         :: L_BC_EAST  = .false.
    logical                         :: L_BC_WEST  = .false.
!---------------------------------------------------------------!

    integer (kind=4), parameter     :: L_SN = 1                 !----> SHARED NODE NUMBER
    integer (kind=4)                :: DDT_SN                   !----> DERIVED DATA TYPE FOR SHARED NODE
    integer (kind=4)                :: L_NX                     !----> LOCAL NX
    integer (kind=4)                :: L_START_X                !----> LOCAL STARTIONG INDEX


!!  MY-IDEA PARAMETERS
    integer (kind=4), parameter                 :: stencil_nx = 31, stencil_mid = stencil_nx/2+1
    double precision, dimension(1:stencil_nx)   :: stencil
    double precision                            :: EAST_RIGHT, EAST_LEFT
    double precision                            :: WEST_RIGHT, WEST_LEFT

end module MPI_VAR




!! PROGRAM wave : MAIN PROGRAM BLOCK

program wave
    use par
    use var
    use MPI_VAR
    implicit none

    call prog_initialize
    call tag_process

    call initial_condition

!! START OF MAIN TIME LOOP
!-------------------------------------------------------!
    do while(time .le. t_max)

        call Runge_Kutta

        iter = iter + 1
        time = dt*dble(iter)

        if(mod(iter,500) .eq. 0) then
            call write_IO
        end if

    end do
!-------------------------------------------------------!

    call write_IO


    call prog_finalize
end program wave






subroutine prog_initialize
    use MPI
    use MPI_VAR, only: MPI_ERR, NUMPROC, RANK, stencil_nx, stencil
    implicit none
    
    integer (kind=4)    :: i,j,ios, iounit = 45

    call MPI_INIT(MPI_ERR)                                                  !----> INITIALIZE MPI COMM 
    call MPI_COMM_SIZE(MPI_COMM_WORLD,NUMPROC,MPI_ERR)                      !----> DETERMINE NUMBER OF PROCESSOR
 !  call MPI_COMM_RANK(MPI_COMM_WORLD,RANK,MPI_ERR)                         !----> DETERMINE RANK

    call oucs3_stencil(stencil_nx,stencil)

    return
end subroutine prog_initialize




!!  SUBROUTINE TAGPROCESS
!---------------------------------------------------------------------------!

subroutine tag_process
    use MPI
    use MPI_VAR
    use par
    implicit none

!!  LOCAL VARIABLE
!---------------------------------------------------------------------------!
    integer (kind=4)    :: MODX
!---------------------------------------------------------------------------!

    CART(1:NDIM) = (/CART_NX/)

    call MPI_CART_CREATE(MPI_COMM_WORLD,NDIM,CART,PERIOD,.false.,MPI_COMM_CART,MPI_ERR)
    call MPI_COMM_RANK  (MPI_COMM_CART ,RANK,MPI_ERR)
    call MPI_CART_COORDS(MPI_COMM_CART ,RANK,NDIM,COORDS,MPI_ERR)

!!  ASSIGNING NEIGHBOUR
!---------------------------------------------------------------------------!
    call MPI_CART_SHIFT(MPI_COMM_CART,0,1,NEIGHBOUR_WEST ,NEIGHBOUR_EAST ,MPI_ERR)
!---------------------------------------------------------------------------!

!!  LOCATE BOUNDARY (ENQUIRE WHETHER THERE IS ANY BOUNDARY OR NOT)
!---------------------------------------------------------------------------!
    call boundary_locator
!---------------------------------------------------------------------------!

!!  OBTAIN LOCAL NX AND NY
!---------------------------------------------------------------------------!

    L_NX      = nx/CART_NX + 1

    MODX      = mod(nx,CART_NX)

    L_START_X = COORDS(1)*(L_NX-1)+1            ! +1 ACCOUNTS FOR THE OFFSET FROM ZERO

    if(COORDS(1) .ge. CART_NX-MODX) L_START_X = L_START_X + (COORDS(1)-CART_NX+MODX)
    if(COORDS(1) .ge. CART_NX-MODX) L_NX      = L_NX      + 1

!---------------------------------------------------------------------------!

    return
end subroutine tag_process
!---------------------------------------------------------------------------!




!! HELPS TO LOCATE THE BOUNDARY
!---------------------------------------------------------------------------!

subroutine boundary_locator
    use MPI    , only: MPI_PROC_NULL
    use MPI_VAR, only: NEIGHBOUR_EAST , &
    &                  NEIGHBOUR_WEST , &
    &                  L_BC_EAST      , &
    &                  L_BC_WEST      
    implicit none

    if(NEIGHBOUR_EAST  .eq. MPI_PROC_NULL) L_BC_EAST  = .true.
    if(NEIGHBOUR_WEST  .eq. MPI_PROC_NULL) L_BC_WEST  = .true.

    return
end subroutine boundary_locator
!---------------------------------------------------------------------------!




subroutine initial_condition
    use MPI
    use MPI_VAR
!   use par
    use var
    implicit none

!-----------------------local variable--------------------------!
    integer                         :: i, err
    double precision                :: func
!-----------------------local variable--------------------------!

!!  ALLOCATING 1D ARRAY OF VARIABLES 
!!--------------------------------------------------------------!  
    allocate(x(1:L_NX), stat=err)
    if (err /= 0) print *, "x: Allocation request denied"

    allocate(u(1-L_SN :L_NX+L_SN), stat=err)
    if (err /= 0) print *, "u: Allocation request denied"

    allocate(dudx(1:L_NX), stat=err)
    if (err /= 0) print *, "dudx: Allocation request denied"

!!--------------------------------------------------------------!

!!  GRID GENERATION
    forall(i=1:L_NX) x(i) = dx*dble(L_START_X+i-2)  !! UNIFORM GRID

!!  INITIAL FUNCTION
    forall(i=1:L_NX) u(i) = func(x(i))

!! INITIAL TIME AND ITERATION 
    time = 0.00d0
    iter = 0

    return
end subroutine initial_condition




subroutine prog_finalize
    use MPI
    use MPI_VAR, only: MPI_ERR, RANK, L_NX
    use var,     only: x,u,dudx
    implicit none

!-----------------------local variable--------------------------!
    integer                         :: i,err
!-----------------------local variable--------------------------!

    if (allocated(x)) deallocate(x, stat=err)
    if (err /= 0) print *, "x: Deallocation request denied"

    if (allocated(u)) deallocate(u, stat=err)
    if (err /= 0) print *, "u: Deallocation request denied"

    if (allocated(dudx)) deallocate(dudx, stat=err)
    if (err /= 0) print *, "dudx: Deallocation request denied"

    call MPI_FINALIZE(MPI_ERR)

    return
end subroutine prog_finalize




subroutine meet_neighbour
    use MPI
    use MPI_VAR
    use var
    implicit none


    call MPI_TYPE_CONTIGUOUS(L_SN,MPI_DOUBLE_PRECISION,DDT_SN,MPI_ERR)
    call MPI_TYPE_COMMIT(DDT_SN,MPI_ERR)

!!  INFORMATION IS PASSED FROM RIGHT SIDE OF DONOR TO LEFT SIDE OF RECEIVER

    call MPI_SENDRECV(u(L_NX-L_SN),1,DDT_SN,NEIGHBOUR_EAST,1, &
    &                 u(1   -L_SN),1,DDT_SN,NEIGHBOUR_WEST,1,MPI_COMM_WORLD,STATUS,MPI_ERR)


    call MPI_SENDRECV(u(2        ),1,DDT_SN,NEIGHBOUR_WEST,2, &
    &                 u(L_NX+1   ),1,DDT_SN,NEIGHBOUR_EAST,2,MPI_COMM_WORLD,STATUS,MPI_ERR)


    call MPI_TYPE_FREE(DDT_SN,MPI_ERR)

    call MPI_BARRIER(MPI_COMM_WORLD, MPI_ERR)

    return
end subroutine meet_neighbour




!! OBTAINING DERIVATIVE WITH COMPACT SCHEME (OUCS3)

subroutine spatial_devative
    use oucs3_const,  only: p,q,r
    use var
    use MPI_VAR
    use MPI
    implicit none

    interface tdma
        function tdma(a) result(b)
            double precision, intent(in)    :: a(:,:)
            double precision                :: b(1:size(a,1))
        end function tdma
    end interface ! tdma


    integer (kind=4)                        :: i,j
    double precision, dimension(1:L_NX,1:4) :: a
    double precision, dimension(1:L_NX    ) :: b
    double precision :: tempdat

    a = 0.0d0

!!  DEFINING A AND B MATRIX MATRIX
    do i=2,L_NX-1

        a(i,1) = p(-1)
        a(i,2) = 1.0d0
        a(i,3) = p( 1)

        do j=-2,2
            a(i,4) = a(i,4) + u(i+j)*q(j)*idx
        end do

    end do


    EAST_LEFT  = 0.0d0
    WEST_RIGHT = 0.0d0
    do j=1,stencil_nx/2
        EAST_LEFT  = EAST_LEFT  + u(L_NX-j)*stencil(stencil_mid-j)
        WEST_RIGHT = WEST_RIGHT + u(1   +j)*stencil(stencil_mid+j)
    end do


!!  INFORMATION IS PASSED FROM RIGHT SIDE OF DONOR TO LEFT SIDE OF RECEIVER AND VICEVERSA

    call MPI_SENDRECV(EAST_LEFT ,1,MPI_DOUBLE_PRECISION,NEIGHBOUR_EAST,1, &
    &                 WEST_LEFT ,1,MPI_DOUBLE_PRECISION,NEIGHBOUR_WEST,1,MPI_COMM_WORLD,STATUS,MPI_ERR)

    call MPI_SENDRECV(WEST_RIGHT,1,MPI_DOUBLE_PRECISION,NEIGHBOUR_WEST,1, &
    &                 EAST_RIGHT,1,MPI_DOUBLE_PRECISION,NEIGHBOUR_EAST,1,MPI_COMM_WORLD,STATUS,MPI_ERR)

    call MPI_BARRIER(MPI_COMM_WORLD, MPI_ERR)

!!  EAST BOUNDARY CLOSURE : SUGGESTED BY BISWAJIT
!------------------------------------------------------!
    if (L_BC_EAST) then
        print*, "Non-Periodic conditon not available!"
    else
!!  PERODIC BOUNDARY CONDITION AT EAST
        a(L_NX,2) = 1.0d0
        a(L_NX,4) = idx*(EAST_RIGHT + EAST_LEFT)
    end if
!------------------------------------------------------!


!!  WEST BOUNDARY CLOSURE : SUGGESTED BY BISWAJIT
!------------------------------------------------------!
    if (L_BC_WEST) then
        print*, "Non-Periodic conditon not available!"
    else
!!  PERODIC BOUNDARY CONDITION AT WEST
        a(1   ,2) = 1.0d0
        a(1   ,4) = idx*(WEST_RIGHT + WEST_LEFT)
    end if
!------------------------------------------------------!

!!  OBTAINING 1ST DERIVATIVE
    !-------------------!
        dudx = tdma(a)  !
    !-------------------!

    return
end subroutine spatial_devative




!!  PERFORMING TIME INTEGRATION WITH RUNGE KUTTA METHOD

subroutine Runge_Kutta
    use par, only: aj,bj
    use var
    use MPI_VAR
    implicit none

    integer (kind=4)                    :: i,irk,err

    allocate(du(1:L_NX), stat=err)
    if (err /= 0) print *, "du: Allocation request denied"

    allocate(u_new(1:L_NX), stat=err)
    if (err /= 0) print *, "u_new: Allocation request denied"
    
!------------------rk-------------------!
    du(:)  = 0.0d0

    do irk = 1,4

        call meet_neighbour
        call spatial_devative


        do i   = 1,L_NX
            du(i) = aj(irk)*du(i) - c*dt*dudx(i)
        end do

        do i = 1, L_NX
            u_new(i) = u(i) + bj(irk)*du(i)
        end do

        u(1:L_NX) = u_new(1:L_NX)

    end do
!------------------rk-------------------!

    if (allocated(du)) deallocate(du, stat=err)
    if (err /= 0) print *, "du: Deallocation request denied"

    if (allocated(u_new)) deallocate(u_new, stat=err)
    if (err /= 0) print *, "u_new: Deallocation request denied"

    return
end subroutine Runge_Kutta





subroutine write_IO
    use par
    use var
    use MPI_VAR
    use MPI
    implicit none

    integer (KIND=4), dimension(1:1)                    :: G_SIZE, L_SIZE, L_ORG
    integer (KIND=4)                                    :: DDT_CHAR_STR, DDT_ARRAY
    integer (kind=mpi_offset_kind)                      :: OFFSET
    integer (KIND=4)                                    :: FH,INT_LENG
    integer (KIND=4), parameter                         :: CHAR_STR_LEN = 60
    character(len = CHAR_STR_LEN), dimension(1:L_NX-1)  :: DATASET
    character(len=20)                                   :: FILENAME, ITER_STMP, filen, FORMAT_STR
    character(len=30)                                   :: SOLUTIONTIME, TITLE, VARIABLES, ZONE

    integer (KIND=4)                                    :: i


    if(RANK .eq. ROOT)  print*, "writing data at time =", time

    write(ITER_STMP,"(I10)") iter  
    INT_LENG = len(trim(adjustl(ITER_STMP))) 

    write(FORMAT_STR,"(A5,I1,A4)") "(A7,A",INT_LENG,",A4)"
    write(FILENAME,FORMAT_STR) "output_",adjustl(ITER_STMP),".dat"

    do i = 1,L_NX-1
        write(DATASET(i),*),  NEW_LINE('A'), x(i),u(i)
    end do


    write(TITLE,"(A17)") "TITLE = WAVE PLOT"
    write(VARIABLES,"(A15)") "VARIABLES = X,U"
    write(ZONE,"(A7,I4)")"ZONE I=",nx
    write(SOLUTIONTIME,"(A15,F10.5)") "SOLUTIONTIME = ",time

    open(80,file=trim(FILENAME))
        write(80,"(A30)") TITLE
        write(80,"(A30)") VARIABLES
        write(80,"(A30)") ZONE
        write(80,"(A30)") SOLUTIONTIME
    close(80)




    OFFSET = 120
    G_SIZE = (/nx/)
    L_SIZE = (/L_NX-1/)
    L_ORG  = (/L_START_X-1/)


    call MPI_TYPE_CONTIGUOUS(CHAR_STR_LEN, MPI_CHARACTER, DDT_CHAR_STR, MPI_ERR)
    call MPI_TYPE_COMMIT(DDT_CHAR_STR, MPI_ERR)


    call MPI_TYPE_CREATE_SUBARRAY(1,G_SIZE,L_SIZE,L_ORG,MPI_ORDER_FORTRAN,DDT_CHAR_STR,DDT_ARRAY,MPI_ERR)
    call MPI_TYPE_COMMIT(DDT_ARRAY, MPI_ERR)

    call MPI_BARRIER(MPI_COMM_WORLD, MPI_ERR)


!-----------------------------------------------------------------------------------!
    call MPI_FILE_OPEN(MPI_COMM_WORLD,trim(FILENAME),MPI_MODE_CREATE+MPI_MODE_WRONLY,MPI_INFO_NULL,FH,MPI_ERR)
    call MPI_FILE_SET_VIEW(FH,OFFSET,DDT_CHAR_STR,DDT_ARRAY,"NATIVE",MPI_INFO_NULL,MPI_ERR)

    call MPI_FILE_WRITE_ALL(FH,DATASET,L_NX-1,DDT_CHAR_STR,STATUS,MPI_STATUS_IGNORE)
    call MPI_FILE_CLOSE(FH, MPI_ERR)
!-----------------------------------------------------------------------------------!

    call mpi_type_free(ddt_array, mpi_err)
    call mpi_type_free(ddt_char_str, mpi_err)

    return
end subroutine write_IO




function func(x) result(f)
    use par, only: k, x_start, x_end
    implicit none

    double precision, intent(in)    :: x
    double precision                :: f
    double precision, parameter     :: alpha = 50.0d0
    double precision, parameter     :: xm    = 1.50d0
!
!----------------------wave packet----------------------!
!                                                       !
    f = dexp( -alpha*(x-xm)**2 ) * dsin(k*x)            !
!                                                       !
!----------------------wave packet----------------------!

    return
end function func




!----------------------------FUNCTION tdma------------------------------!

function tdma(a) result(b)
    implicit none

    double precision, intent(in)    :: a(:,:)
    double precision                :: b(1:size(a,1))

!---------------------------local variable------------------------------!
    integer (kind=4)                :: i,n,m
    double precision                :: c(1:size(a,1)),d(1:size(a,1))
!---------------------------local variable------------------------------!

    n = size(a,1)
    m = size(a,2)
    if(m/=4) stop "ERROR: Mismatch in 'tdma' array dimension!"

    c(1)=a(1,3)/a(1,2)
    d(1)=a(1,4)/a(1,2)

    do i=2,n-1
        c(i)=a(i,3)/(a(i,2)-a(i,1)*c(i-1))
        d(i)=(a(i,4)-a(i,1)*d(i-1))/(a(i,2)-a(i,1)*c(i-1))
    end do
        
    d(n)=(a(n,4)-a(n,1)*d(n-1))/(a(n,2)-a(n,1)*c(n-1))

    b(n)=d(n)

    do i=n-1,1,-1
        b(i)=d(i)-c(i)*b(i+1)
    end do

    return
end function tdma

!----------------------------END FUNCTION tdma--------------------------!




subroutine oucs3_stencil(n,stencil)
    implicit none

    interface inverse
        function inverse(a) result(c)
            double precision, intent(in )   :: a(:,:)
            double precision                :: c(size(a,dim=1),size(a,dim=2))
        end function inverse
    end interface ! inverse

    integer (kind=4), intent(in )                   ::  n
    double precision, dimension(1:n),intent(out)    ::  stencil

    double precision, dimension(1:n,1:n)            ::  A, B, A_inverse, C
    integer (kind=4)                                ::  i,j
    double precision , parameter                    ::  E = 0.7877868950d0, &
    &                                                   F = 0.0458012980d0, &
    &                                                   D = 0.3793894912d0

    if(mod(n,2) .eq. 0) stop "ERROR: Please enter ODD value of 'n'"

    A = 0.0d0
    B = 0.0d0

!! DEFINATION OF A
    forall (j=2:n  ) A(j  ,j-1) = D
    forall (j=1:n  ) A(j  ,j  ) = 1.0d0
    forall (j=1:n-1) A(j  ,j+1) = D

    A(1,n) = D
    A(n,1) = D

!! DEFINATION OF B
    forall (j=3:n  ) B(j  ,j-2) = -F
    forall (j=2:n  ) B(j  ,j-1) = -E
    forall (j=1:n-1) B(j  ,j+1) =  E
    forall (j=1:n-2) B(j  ,j+2) =  F

    B(1,n  ) = -E
    B(1,n-1) = -F
    B(2,n  ) = -F

    B(n  ,1) = E
    B(n  ,2) = F
    B(n-1,1) = F


    A_inverse = inverse(A)

    C = matmul(A_inverse,B)

    do i=1,n 
    stencil(i) = C(n/2+1,i)
    end do 

end subroutine oucs3_stencil




function inverse(a) result(c)
    implicit none

    double precision, intent(inout)                             :: a(:,:)
    double precision                                            :: c(size(a,dim=1),size(a,dim=2))

    double precision, dimension(size(a,dim=1),size(a,dim=2))    :: l,u
    double precision, dimension(size(a,dim=1))                  :: b,d,x

    double precision                                            :: coeff
    integer (kind=4)                                            :: i,j,k,n

    if(size(a,dim=1) .ne. size(a,dim=2)) stop "ERROR: please enter a square matrix!"
    n = size(a,dim=1)

    l = 0.0d0;  u = 0.0d0;  b = 0.0d0

!! STEP 1: FORWARD ELIMINATION
    do k=1, n-1
        do i=k+1,n
            coeff  = a(i,k)/a(k,k)
            l(i,k) = coeff
            do j=k+1,n
                a(i,j) = a(i,j)-coeff*a(k,j)
            end do
        end do
    end do

!!  STEP 2: PREPARE l AND u MATRICES 
    do i=1,n
        l(i,i) = 1.0d0
    end do

!!  U MATRIX IS THE UPPER TRIANGULAR PART OF a
    do j=1,n
        do i=1,j
            u(i,j) = a(i,j)
        end do
    end do

!!  STEP 3: COMPUTE COLUMNS OF THE INVERSE MATRIX c
    do k=1,n
        b(k) = 1.0d0
        d(1) = b(1)
!!  STEP 3A: SOLVE ld=b USING THE FORWARD SUBSTITUTION
        do i=2,n
            d(i)=b(i)
            do j=1,i-1
                d(i) = d(i) - l(i,j)*d(j)
            end do
        end do
!!  STEP 3B: SOLVE ux=d USING THE BACK SUBSTITUTION
        x(n)=d(n)/u(n,n)
        do i = n-1,1,-1
            x(i) = d(i)
            do j=n,i+1,-1
                x(i)=x(i)-u(i,j)*x(j)
            end do
            x(i) = x(i)/u(i,i)
        end do
!!  STEP 3C: FILL THE SOLUTIONS x(n) INTO COLUMN k OF c
        do i=1,n
            c(i,k) = x(i)
        end do
        b(k) = 0.0d0
    end do

    return
end function inverse