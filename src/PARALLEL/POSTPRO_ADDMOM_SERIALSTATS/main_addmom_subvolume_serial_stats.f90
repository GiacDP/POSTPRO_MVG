PROGRAM postpro_addmom_subvolume_serial_stats
!
! This program evaluates the added momentum on the xz plane.
!
! E_add is already normalised by U_inf**2 (rho_inf*U_inf**2) AND yg_edge (h/delta_0).
! ==> e_add = E_add/h
! 
! OSS. 
!   h = height of the ramp
!   yg_edge = h/delta_0
!   h_delta = h/delta_99 with delta_99 blay thickness at shock impingment
!   yg_edge/h_delta = delta_99/delta_0
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
!   mu0                      Dynamic viscosity at reference temperature from "flow_params.dat"
!   yg_edge                  h/delta_0
!   h_delta                  h/delta_99
! ./dxg, ./dyg, ./dzg        Coords file from simulation
! ./stat3d0_'//chx//'_0000_'//chz//'.bin  Time-averaged 3D stats for controlled case (block chx, chz)
! ../stat.bin                Time- and span-wise averaged 2D stats for uncontrolled case
! 
! OUTPUT FILES
! ./Epln_0001.vtr            Added momentum on xz plane .vtr
!   E_adh_ufavr2             E_add based on Favre-averaged velocity
!   E_adh_Rrhou2               "     "    " rho*u**2 average
!
! 23/01/2023 G. Della Posta
!

USE modpostpro
IMPLICIT NONE

!------------------------------------------------------------------------------

INTEGER :: i, j, k, l, ll, kk, m, ks, jj, ii, ii1, ii2, j99
INTEGER :: ig_start, ig_end, num_blocks, nv_sub, jg_start, jg_end
INTEGER :: kproc, ksym_start, ksym_end
INTEGER :: nprocx, nprocy, nprocz, nxserial, nyserial, nzserial
INTEGER :: nxsmax, nysmax, nzsmax, nxs, nys, nzs, k_start
INTEGER :: i_block_start, i_block_end
INTEGER :: j_giepman

INTEGER, DIMENSION(1) :: jtmp

REAL(rkind) :: cf_c, dy, ygg, int1, int2, um, up, u2m, u2p, delta99, tmp
REAL(rkind) :: uu, dely, unum, uden
REAL(rkind) :: u_infty, mu0
REAL(rkind) :: dd

REAL(rkind), DIMENSION(:), ALLOCATABLE :: dudy_w_bl
REAL(rkind), DIMENSION(:), ALLOCATABLE :: delta_inc_bl, theta_inc_bl, delta_99_bl

REAL(rkind), DIMENSION(:,:), ALLOCATABLE :: w_sym_tmp
REAL(rkind), DIMENSION(:,:), ALLOCATABLE :: e_add_tmp

REAL(rkind), DIMENSION(:,:,:), ALLOCATABLE :: ww_stat, w_stat2

REAL(rkind), DIMENSION(:,:,:,:), ALLOCATABLE :: w_sym
REAL(rkind), DIMENSION(:,:,:,:), ALLOCATABLE :: diff_mom, e_add

CHARACTER(len = 4) :: chx, chz
CHARACTER(len = 12), DIMENSION(2) :: names = ["E_adh_ufavr2", "E_adh_Rrhou2"]

!------------------------------------------------------------------------------

LOGICAL :: reord

INTEGER (kind=mpi_offset_kind) :: offset

INTEGER :: proc_last, ny_last, nymax_paral
INTEGER :: ntot3d, ntotxy
INTEGER :: mpi_io_file
INTEGER :: filetype
INTEGER :: size_real

INTEGER, DIMENSION(3) :: sizes    ! Dimensions of the total grid
INTEGER, DIMENSION(3) :: subsizes ! Dimensions of grid local to a procs
INTEGER, DIMENSION(3) :: starts   ! Starting coordinates

!------------------------------------------------------------------------------

 call mpi_init(iermpi)
 call mpi_comm_rank(mpi_comm_world,nrank,iermpi)
 call mpi_comm_size(mpi_comm_world,nproc,iermpi)

 masterproc = .false.
 if (nrank==0) masterproc = .true.

 OPEN(12, FILE='./input_subvolume_serial.dat', ACTION='read', STATUS='old')
 READ(12,*)
 READ(12,*) nxmax, nymax, nzmax, nprocx, nprocy, nprocz
 READ(12,*)
 READ(12,*)
 READ(12,*) re_delta, mach, ig_start, jg_end, nzsmax, mu0, yg_edge, h_delta
 CLOSE(12)

 u_infty    = SQRT(gamma)*mach

 nblocks(1) = nproc
 nblocks(2) = 1
 nblocks(3) = 1

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

 nxsmax  = ig_end - ig_start + 1
 nysmax  = jg_end - jg_start + 1
 nzsmax  = nzsmax/2

 nxs     = nx
 nys     = nysmax
 nzs     = nzsmax

 k_start    = 0

 nv_stat_3d = 14 ! Reynolds rho, Favre u, Favre w
 nv_sub     = 3

!------------------------------------------------------------------------------

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

 dxg     = xg( 2) - xg( 1)
 dzg     = zg( 2) - zg( 1)
 dyg     = yg(2:nymax) - yg(1:nymax-1)
 
 call mpi_barrier(mpi_comm_world,iermpi)

 pbc(1) = .false.
 pbc(2) = .false.
 pbc(3) = .true.

! Create 3D topology

 reord = .false.
 call mpi_cart_create(mpi_comm_world,ndims,nblocks,pbc,reord,mp_cart,iermpi)
 call mpi_cart_coords(mp_cart,nrank,ndims,ncoords,iermpi)

!------------------------------------------------------------------------------
! Read stat3d.bin
!------------------------------------------------------------------------------

 IF (masterproc) WRITE(*,*) 'Start read stat3d.bin.'

 ! Read stat3d.bin

 IF (masterproc) WRITE(*,*) 'Start read stat3d and average subvolume.'
 IF (masterproc) WRITE(*,*) nxs, nxserial, nzs 

 ALLOCATE(w_stat_3d(1:nv_stat_3d,nxserial,nyserial,nzserial))
 ALLOCATE(w_sym(nv_sub, nxs, nys, nzs+1))

 DO kproc = 0, nprocz-1

  WRITE(chx,'(I4.4)') nrank + i_block_start
  WRITE(chz,'(I4.4)') kproc

  IF (masterproc) WRITE(*,*) 'Reading kproc = ', kproc

  w_stat_3d = 0._rkind

  OPEN(11,file='stat3d0_'//chx//'_0000_'//chz//'.bin',form='unformatted', &
      & STATUS='old', ACTION='read')
  READ(11) w_stat_3d(1:nv_stat_3d,1:nxserial,1:nyserial,1:nzserial)
  CLOSE(11)

  IF (kproc < nprocz/2) THEN
    ksym_start = kproc*nzserial + 1
    ksym_end   = ksym_start + nzserial - 1
    w_sym(1,:,:,ksym_start:ksym_end) = w_stat_3d(1,1:nxs,1:nys,1:nzserial)
    w_sym(2,:,:,ksym_start:ksym_end) = w_stat_3d(2,1:nxs,1:nys,1:nzserial)
    w_sym(3,:,:,ksym_start:ksym_end) = w_stat_3d(7,1:nxs,1:nys,1:nzserial)
  ELSE
    ksym_start = (nprocz - kproc - 1)*nzserial + 1
    ksym_end   = ksym_start + nzserial - 1
    IF (kproc == nprocz/2) THEN
      w_sym(1,:,:,ksym_end+1) = w_stat_3d(1, 1:nxs, 1:nys, 1)
      w_sym(2,:,:,ksym_end+1) = w_stat_3d(2, 1:nxs, 1:nys, 1)
      w_sym(3,:,:,ksym_end+1) = w_stat_3d(7, 1:nxs, 1:nys, 1)
    ENDIF
    ! rho
    w_sym(1,:,:,ksym_start+1:ksym_end) = &
       & 0.5_rkind * (w_sym(1,:,:,ksym_start+1:ksym_end) + w_stat_3d(1,1:nxs, 1:nys, nzserial:2:-1))
    w_sym(1,:,:,ksym_end+1) = 0.5_rkind * (w_sym(1,:,:,ksym_end+1) + w_stat_3d(1, 1:nxs, 1:nys, 1))
    ! rho u
    w_sym(2,:,:,ksym_start+1:ksym_end) = &
       & 0.5_rkind * (w_sym(2,:,:,ksym_start+1:ksym_end) + w_stat_3d(2,1:nxs, 1:nys, nzserial:2:-1))
    w_sym(2,:,:,ksym_end+1) = 0.5_rkind * (w_sym(2,:,:,ksym_end+1) + w_stat_3d(2, 1:nxs, 1:nys, 1))
    ! rho u^2
    w_sym(3,:,:,ksym_start+1:ksym_end) = &
       & 0.5_rkind * (w_sym(3,:,:,ksym_start+1:ksym_end) + w_stat_3d(7,1:nxs, 1:nys, nzserial:2:-1))
    w_sym(3,:,:,ksym_end+1) = 0.5_rkind * (w_sym(3,:,:,ksym_end+1) + w_stat_3d(7, 1:nxs, 1:nys, 1))
  ENDIF

 ENDDO
 DEALLOCATE(w_stat_3d)

 w_sym(2,:,:,:) = w_sym(2,:,:,:)/w_sym(1,:,:,:)
 w_sym(1,:,:,:) = w_sym(2,:,:,:)**2   ! u_tilde^2
 w_sym(2,:,:,:) = w_sym(3,:,:,:)      ! (rho u^2)_bar

 IF (masterproc) WRITE(*,*) 'End read stat3d.bin and average.'

!------------------------------------------------------------------------------
! Import time and spanwise averaged solution without ramp (stat.bin)
 !------------------------------------------------------------------------------

 IF (masterproc) WRITE(*,*) 'Start read clean stats.'

 ALLOCATE(w_stat2(nxmax,nymax,70), ww_stat(2,nx,ny))

 OPEN(11, FILE='../stat.bin', FORM='unformatted', STATUS='old', ACTION='read', ACCESS='stream')
 READ(11) w_stat2
 CLOSE(11)
 ii1 = ig_start + nrank*nx
 ii2 = ii1 + nx - 1
 ww_stat(1,:,:) = w_stat2(ii1:ii2,:,13)/w_stat2(ii1:ii2,:,1)
 ww_stat(2,:,:) = w_stat2(ii1:ii2,:,16)                ! (rho u^2)_bar
 ww_stat(1,:,:) = ww_stat(1,:,:)**2                    ! u_tilde^2
 DEALLOCATE(w_stat2)

 IF (masterproc) WRITE(*,*) 'End read clean stats.'

 !------------------------------------------------------------------------------
 ! Added momentum 
 !------------------------------------------------------------------------------

 IF (masterproc) WRITE(*,*) 'E_add 0.43 delta evaluation.'

 ALLOCATE(diff_mom(  2, nxs,   nys, nzs+1))
 DO k = 1, nzs+1
   DO j = 1, nys
     diff_mom(1,:,j,k) = w_sym(1,:,j,k) - ww_stat(1,:,j)
     diff_mom(2,:,j,k) = w_sym(2,:,j,k) - ww_stat(2,:,j)
   ENDDO
 ENDDO
 DEALLOCATE(w_sym, ww_stat)

 IF (masterproc) WRITE(*,*) 'Diff mom defined.'

 ALLOCATE(e_add_tmp(2, nzs+1), e_add(2, nx, 1, nzs+1))
 dd = 0.43_rkind*yg_edge/h_delta
 jtmp = MINLOC(ABS(yg-dd))
 j_giepman = jtmp(1)
 DO i = 1, nx
   e_add_tmp = 0._rkind
   DO j = 2, j_giepman
     e_add_tmp = e_add_tmp + 0.5_rkind*dyg(j-1)*(diff_mom(:,i,j,:)+diff_mom(:,i,j-1,:))
   ENDDO
   e_add(:,i,1,:) = e_add_tmp
 ENDDO
 e_add = e_add/u_infty**2/yg_edge
 DEALLOCATE(diff_mom, e_add_tmp)

 IF (masterproc) WRITE(*,*) 'End E_add 0.43 delta evaluation.'

 CALL MPI_BARRIER(mpi_comm_world,iermpi)

 !------------------------------------------------------------------------------
 ! Finalize and deallocate
 !------------------------------------------------------------------------------
 ! WRITE OUTPUT

 IF (masterproc) WRITE(*,*) '!------------START WRITING OUTPUT------------'

 IF (masterproc) WRITE(*,*) 'Writing added momentum plane. ALREADY NORMALISED by (rho_inf)*U_inf^2*h.'

 !OPEN(12,file='./e_rho_add_tecplot_'//chx//'.dat')
 !WRITE(12,*) 'zone i=',nxs,', j=',nzs+1
 !DO k = 1, nzs+1
 ! DO i = 1, nx
 !  WRITE(12,'(3ES20.10)') xg(i+(i_block_start+nrank)*nx), zg(k), e_add(2,i,1,k)
 ! ENDDO
 !ENDDO
 !CLOSE(12)

 CALL write_vtk_general(ig_start, k_start, 2, nx, 1, nzs+1, &
        & e_add(1:2, 1:nx, 1, 1:nzs+1), 'Epln', names, &
        & nblocks, ncoords, mp_cart)

 IF (masterproc) WRITE(*,*) '!-------------END WRITING OUTPUT-------------'

 IF (masterproc) WRITE(*,*) '!---------END PROGRAM POSTPRO_ADDMOM---------'

 call mpi_finalize(iermpi)

END PROGRAM postpro_addmom_subvolume_serial_stats
