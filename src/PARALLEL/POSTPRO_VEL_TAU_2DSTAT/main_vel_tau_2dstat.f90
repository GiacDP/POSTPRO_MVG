PROGRAM postpro_vel_tau_serialstats
! 
! This program reads ../stat.bin and write sub1 vtr.
! 
! ATTENTION: All the variables are normalised by U_infty and U_infty^2
!
! INPUT FILES
! ./input_subvolume_serial.dat
!   nxmax, nymax, nzmax      Streamwise/Wall-normal/Spanwise total number of points
!   nprocx, nprocy, nprocz        "          "         "     number of blocks
!   re_delta                 Reynolds number from simulation "flow_params.dat"
!   mach                     Mach number       "       "            " 
!   ig_start                 Min global streamwise index considered
!   jg_end                   Max   "    wall-normal  "      " 
!   nzsmax                   Spanwise number of points considered in symmetric average
! ./dxg, ./dyg, ./dzg        Coords file from simulation
! ../stat.bin                Time- and span-wise averaged 2D stats for uncontrolled case
!
! OUTPUT FILES
!  sub1_0001.vtr (rho_bar, (u_i)_tilde, p_bar)
!
! 11/01/2023 G. Della Posta
!
 USE modpostpro
 IMPLICIT NONE
!
!------------------------------------------------------------------------------
! Subvolume

 INTEGER :: i,j,k,l,ks,kk,m,jj,ii
 INTEGER :: ig_start, ig_end, num_blocks, nv_sub, jg_start, jg_end
 INTEGER :: kproc, ksym_start, ksym_end
 INTEGER :: nprocx, nprocy, nprocz, nxserial, nyserial, nzserial
 INTEGER :: nxsmax, nysmax, nzsmax, nxs, nys, nzs
 INTEGER :: k_start, i_l, i_r, j_l, j_r, k_l, k_r
 INTEGER :: idestin, isource, ii1, ii2
 INTEGER :: mp_cart_xz

 INTEGER, DIMENSION(3,3) :: it 

 REAL(rkind) :: duiidx, duiidy, duiidz, pp, dxm1, dym1, dzm1

 REAL(rkind), DIMENSION(:,:,:), ALLOCATABLE :: w_stat2

 REAL(rkind), DIMENSION(:,:,:,:), ALLOCATABLE :: w_sym
 REAL(rkind), DIMENSION(:,:,:,:), ALLOCATABLE :: production
 
 CHARACTER(len = 4) :: chx, chz
 CHARACTER(len = 12), DIMENSION(5) :: names1 = &
         & ["Reyn_Density", "FaVelocity_x", "FaVelocity_y", & 
         & "FaVelocity_z", "ReynPressure"]

!------------------------------------------------------------------------------
! MPI

 LOGICAL :: reord
 
 INTEGER, DIMENSION(3) :: sizes    ! Dimensions of the total grid
 INTEGER, DIMENSION(3) :: subsizes ! Dimensions of grid local to a procs
 INTEGER, DIMENSION(3) :: starts   ! Starting coordinates
 INTEGER :: ntot3d
 INTEGER :: mpi_io_file
 INTEGER :: filetype
 INTEGER :: size_real
 INTEGER (kind=mpi_offset_kind) :: offset

!------------------------------------------------------------------------------

 call mpi_init(iermpi)
 call mpi_comm_rank(mpi_comm_world,nrank,iermpi)
 call mpi_comm_size(mpi_comm_world,nproc,iermpi)

 masterproc = .false.
 if (nrank==0) masterproc = .true.

 IF (masterproc) WRITE(*,*) '!--------START PROGRAM POSTPRO_VEL_TAU-------'

 OPEN(12, FILE='./input_subvolume_serial.dat', ACTION='read', STATUS='old')
 READ(12,*) 
 READ(12,*) nxmax, nymax, nzmax, nprocx, nprocy, nprocz
 READ(12,*) 
 READ(12,*) 
 READ(12,*) re_delta, mach, ig_start, jg_end, nzsmax
 CLOSE(12)

 nv_stat_3d = 14
 nv_sub     = 12

 u0         = SQRT(gamma)*mach

 nblocks(1)  = nproc
 nblocks(2)  = 1
 nblocks(3)  = 1

 pbc(1) = .false.
 pbc(2) = .false.
 pbc(3) = .true.

! Create 3D topology

 reord = .false.
 call mpi_cart_create(mpi_comm_world,ndims,nblocks,pbc,reord,mp_cart,iermpi)
 call mpi_cart_coords(mp_cart,nrank,ndims,ncoords,iermpi)

! Sizes of stats files and procs data

 nxserial = nxmax/nprocx
 nyserial = nymax/nprocy
 nzserial = nzmax/nprocz

 nx       = nxserial ! Parallel in x
 ny       = nymax/1 
 nz       = nzmax/1

 i_block_start = ig_start/nxserial
 i_block_end   = nprocx - 1
 num_blocks    = i_block_end - i_block_start + 1

 IF (num_blocks /= nproc) THEN
  IF (masterproc) WRITE(*,*) 'num_blocks, nproc = ', num_blocks, nproc
  IF (masterproc) WRITE(*,*) 'i_block_start, i_block_end = ', i_block_start, i_block_end
  RETURN
 ENDIF

 ig_start = nxserial*i_block_start + 1
 ig_end   = nxmax 
 jg_start = 1

!------------------------------------------------------------------------------
! Subvolume inputs
!

 IF (MOD(nzsmax,2) /= 0) THEN
  nzsmax = CEILING(REAL(nzsmax)/2)*2
  IF (masterproc) WRITE(*,*) 'nzsmax (multiple of 2!) = ', nzsmax
 ENDIF

 k_start = (nzmax - nzsmax)/2 ! First point before central subvolume

 nxsmax  = ig_end - ig_start + 1
 nysmax  = jg_end - jg_start + 1
 nzsmax  = nzsmax/2

 nxs     = nx
 nys     = nysmax
 nzs     = nzsmax

 allocate(xg(nxmax))
 allocate(yg(nymax), dyg(nymax-1))
 allocate(zg(nzmax))

 open(18,file='./dxg.dat',action='read',status='old')
 do i=1,nxmax
  read(18,*) xg(i)
 enddo
 close(18)
 open(18,file='./dyg.dat',action='read',status='old')
 do j=1,nymax
  read(18,*) yg(j)
 enddo
 close(18)
 open(18,file='./dzg.dat',action='read',status='old')
 do k=1,nzmax
  read(18,*) zg(k)
 enddo
 close(18)
 
 call mpi_barrier(mpi_comm_world,iermpi)

!------------------------------------------------------------------------------
! Import time and spanwise averaged solution without ramp (stat.bin)
!------------------------------------------------------------------------------

 IF (masterproc) WRITE(*,*) 'Start read clean stats.'

 ALLOCATE(w_stat2(nxmax,nymax,70), w_sym(5,nx,ny,1))

 OPEN(11, FILE='../stat.bin', FORM='unformatted', STATUS='old', ACTION='read', & 
 &    ACCESS='stream')
 READ(11) w_stat2
 ii1 = ig_start + nrank*nx
 ii2 = ii1 + nx - 1
 w_sym(1,:,:,1) = w_stat2(ii1:ii2,:, 1)                   ! rho/rho_infty
 w_sym(2,:,:,1) = w_stat2(ii1:ii2,:,13)/w_sym(1,:,:,1)/u0 ! Favre u/u_infty
 w_sym(3,:,:,1) = w_stat2(ii1:ii2,:,14)/w_sym(1,:,:,1)/u0 ! Favre v/u_infty
 w_sym(4,:,:,1) = w_stat2(ii1:ii2,:,15)/w_sym(1,:,:,1)/u0 ! Favre w/u_infty
 w_sym(5,:,:,1) = w_stat2(ii1:ii2,:, 5)                   ! p/p_infty
 DEALLOCATE(w_stat2)
 CLOSE(11)

 IF (masterproc) WRITE(*,*) 'End read clean stats.'

 CALL MPI_BARRIER(mpi_comm_world,iermpi)

!------------------------------------------------------------------------------
! Write output
!------------------------------------------------------------------------------

 ! Redefine grid => same formal operations
 nxmax       = nxsmax
 nymax       = nysmax
 nzmax       = 1
 nx          = nxs
 ny          = nys
 nz          = 1
 i_l         = ig_start
 i_r         = ig_end
 j_l         = 1
 j_r         = nymax
 k_l         = 1
 k_r         = 1
 xg(1:nxmax) = xg(i_l:i_r)
 yg(1:nymax) = yg(j_l:j_r)
 zg(1:nzmax) = zg(k_l:k_r)
 nblocks(1)  = num_blocks
 
 IF (masterproc) THEN
   WRITE(*,*) 'nxsmax, nysmax, nzsmax, ig_start, ig_end, jg_start, jg_end, kg_start, kg_end'
   WRITE(*,*)  nxmax, nymax, nzmax, i_l, i_r, j_l, j_r, k_l, k_r
   WRITE(*,*) 'nblocks = ', nblocks
 ENDIF

 IF (masterproc) WRITE(*,*) 'Writing primitive variables.' 
 CALL write_vtk_reduced_extended_x(5, w_sym(     1:5,:,:,:), 'sub1', names1)

 IF (masterproc) WRITE(*,*) '!-------------END WRITING OUTPUT-------------'

 IF (masterproc) WRITE(*,*) '!---------END PROGRAM POSTPRO_VEL_TAU--------'

 call mpi_finalize(iermpi)

END PROGRAM postpro_vel_tau_serialstats
