PROGRAM postpro_blay_serial_stats
!
! This program evaluates boundary layer props 
! (cf_dudy, cf_dwdy, delta/theta_inc, delta99, rhow) for blay and ramp cases.
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
! ./dxg, ./dyg, ./dzg        Coords file from simulation
! ./stat3d0_'//chx//'_0000_'//chz//'.bin  Time-averaged 3D stats for controlled case (block chx, chz)
! ../stat.bin                Time- and span-wise averaged 2D stats for uncontrolled case
! 
! OUTPUT FILES
! ./blpr_0001.vtr            Boundary layer properties on xz plane .vtr
!    *_bla/*_ram              referred to 2D uncontrolled case/3D ramp-controlled case
! fricdudy*                  Streamwise skin friction coefficient
! fricdwdy* (only ram)       Spanwise    "       "       " 
! deltainc*                  Incompressible displacement boundary layer thickness/delta_0
! thetainc*                         "       momentum        "      "      "   
! delta_99*                  Boundary layer thickness based on 99% velocity/delta_0
! rhow_inf*                  Density at the wall/rho_infty
! pwallinf*                  Mean wall pressure/p_infty
! psdevinf*                  Wall-pressure standard deviation/p_infty
!
!  delta_0 is the reference length = delta_99 at inlet section
!
! 11/01/2023 G. Della Posta
!

USE modpostpro
IMPLICIT NONE

!------------------------------------------------------------------------------

INTEGER :: i, j, k, l, ll, kk, m, ks, jj, ii, ii1, ii2, j99
INTEGER :: ig_start, ig_end, num_blocks, nv_sub, jg_start, jg_end
INTEGER :: kproc, ksym_start, ksym_end
INTEGER :: nprocx, nprocy, nprocz, nxserial, nyserial, nzserial
INTEGER :: nxsmax, nysmax, nzsmax, nxs, nys, nzs, k_start

REAL(rkind) :: cf_c, dy, ygg, int1, int2, um, up, u2m, u2p, delta99, tmp
REAL(rkind) :: uu, dely, unum, uden
REAL(rkind) :: f_t_wall, mu0, s_t, t_wall, u_infty

REAL(rkind), DIMENSION(:), ALLOCATABLE :: dudy_w_bl
REAL(rkind), DIMENSION(:), ALLOCATABLE :: delta_inc_bl, theta_inc_bl, delta_99_bl

REAL(rkind), DIMENSION(:,:), ALLOCATABLE :: w_sym_tmp !, w_stat2
REAL(rkind), DIMENSION(:,:), ALLOCATABLE :: dudy_w, dwdy_w
REAL(rkind), DIMENSION(:,:), ALLOCATABLE :: delta_inc, theta_inc, delta_99

REAL(rkind), DIMENSION(:,:,:), ALLOCATABLE :: ww_stat, w_stat2

REAL(rkind), DIMENSION(:,:,:,:), ALLOCATABLE :: w_sym
REAL(rkind), DIMENSION(:,:,:,:), ALLOCATABLE :: output_array

CHARACTER(len = 4) :: chx, chz
CHARACTER(len = 12), DIMENSION(15) :: names = &
        & ["fricdudy_bla",                 "deltainc_bla", &
        &  "thetainc_bla", "delta_99_bla", "rhow_inf_bla", & 
        &  "pwallinf_bla", "psdevinf_bla",                  & 
        &  "fricdudy_ram", "fricdwdy_ram", "deltainc_ram", & 
        &  "thetainc_ram", "delta_99_ram", "rhow_inf_ram", &
        &  "pwallinf_ram", "psdevinf_ram"]

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
 READ(12,*) re_delta, mach, ig_start, jg_end, nzsmax, mu0
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
 nv_sub     = 5

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

 ALLOCATE(w_stat_3d(nv_stat_3d,nxserial,nyserial,nzserial))
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
  w_stat_3d(3,:,:,:) = w_stat_3d(4,:,:,:) ! rho*w overwrites rho*v
  w_stat_3d(4,:,:,:) = w_stat_3d(14,:,:,:) - w_stat_3d(5,:,:,:)**2 ! 4: sigma_p/p_infty, 5: p_w/p_infty
  w_stat_3d(4,:,:,:) = SQRT(w_stat_3d(4,:,:,:))

  IF (kproc < nprocz/2) THEN
    ksym_start = kproc*nzserial + 1
    ksym_end   = ksym_start + nzserial - 1
    w_sym(1:3,:,:,ksym_start:ksym_end) = w_stat_3d(1:3,1:nxs,1:nys,1:nzserial)
    w_sym(4,:,:,ksym_start:ksym_end) = w_stat_3d(5,1:nxs,1:nys,:) ! 5: sigma_p/p_infty, 4: p_w/p_infty
    w_sym(5,:,:,ksym_start:ksym_end) = w_stat_3d(4,1:nxs,1:nys,:)
  ELSE
    ksym_start = (nprocz - kproc - 1)*nzserial + 1
    ksym_end   = ksym_start + nzserial - 1
    IF (kproc == nprocz/2) THEN
      w_sym(1,:,:,ksym_end+1) = w_stat_3d(1, 1:nxs, 1:nys, 1)
      w_sym(2,:,:,ksym_end+1) = w_stat_3d(2, 1:nxs, 1:nys, 1)
      w_sym(3,:,:,ksym_end+1) = w_stat_3d(3, 1:nxs, 1:nys, 1)
      w_sym(4,:,:,ksym_end+1) = w_stat_3d(5, 1:nxs, 1:nys, 1)
      w_sym(5,:,:,ksym_end+1) = w_stat_3d(4, 1:nxs, 1:nys, 1)
    ENDIF
    ! rho
    w_sym(1,:,:,ksym_start+1:ksym_end) = &
      & 0.5_rkind * (w_sym(1,:,:,ksym_start+1:ksym_end) + w_stat_3d(1,1:nxs, 1:nys, nzserial:2:-1))
    w_sym(1,:,:,ksym_end+1) = 0.5_rkind * (w_sym(1,:,:,ksym_end+1) + w_stat_3d(1, 1:nxs, 1:nys, 1))
    ! u
    w_sym(2,:,:,ksym_start+1:ksym_end) = &
      & 0.5_rkind * (w_sym(2,:,:,ksym_start+1:ksym_end) + w_stat_3d(2,1:nxs, 1:nys, nzserial:2:-1))
    w_sym(2,:,:,ksym_end+1) = 0.5_rkind * (w_sym(2,:,:,ksym_end+1) + w_stat_3d(2, 1:nxs, 1:nys, 1))
    ! w
    w_sym(3,:,:,ksym_start+1:ksym_end) = &
      & 0.5_rkind * (w_sym(3,:,:,ksym_start+1:ksym_end) - w_stat_3d(3,1:nxs, 1:nys, nzserial:2:-1))
    w_sym(3,:,:,ksym_end+1) = 0.5_rkind * (w_sym(3,:,:,ksym_end+1) - w_stat_3d(3, 1:nxs, 1:nys, 1))
    ! p_w/p_infty
    w_sym(4,:,:,ksym_start+1:ksym_end) = &
      & 0.5_rkind * (w_sym(4,:,:,ksym_start+1:ksym_end) + w_stat_3d(5,1:nxs, 1:nys, nzserial:2:-1))
    w_sym(4,:,:,ksym_end+1) = 0.5_rkind * (w_sym(4,:,:,ksym_end+1) + w_stat_3d(5, 1:nxs, 1:nys, 1))
    ! sigma_pw/p_infty
    w_sym(5,:,:,ksym_start+1:ksym_end) = &
      & 0.5_rkind * (w_sym(5,:,:,ksym_start+1:ksym_end) + w_stat_3d(4,1:nxs, 1:nys, nzserial:2:-1))
    w_sym(5,:,:,ksym_end+1) = 0.5_rkind * (w_sym(5,:,:,ksym_end+1) + w_stat_3d(4, 1:nxs, 1:nys, 1))
  ENDIF

 ENDDO
 DEALLOCATE(w_stat_3d)

 w_sym(2,:,:,:) = w_sym(2,:,:,:)/w_sym(1,:,:,:)
 w_sym(3,:,:,:) = w_sym(3,:,:,:)/w_sym(1,:,:,:)

 ! w_sym(1,:,:,:): rho_bar, w_sym(2,:,:,:): u_tilde, w_sym(3,:,:,:): w_tilde
 ! w_sym(4,:,:,:): pw_bar/p_infty, w_sym(5,:,:,:): sigma_pw/p_infty

 IF (masterproc) WRITE(*,*) 'End read stat3d.bin and average.'

!------------------------------------------------------------------------------
! Import time and spanwise averaged solution without ramp (stat.bin)
 !------------------------------------------------------------------------------

 IF (masterproc) WRITE(*,*) 'Start read clean stats.'

 ALLOCATE(w_stat2(nxmax,nymax,70), ww_stat(4,nx,ny))

 OPEN(11, FILE='../stat.bin', FORM='unformatted', STATUS='old', ACTION='read', ACCESS='stream')
 READ(11) w_stat2
 ii1 = ig_start + nrank*nx
 ii2 = ii1 + nx - 1
 ww_stat(1,:,:) = w_stat2(ii1:ii2,:,1)
 ww_stat(2,:,:) = w_stat2(ii1:ii2,:,13)/ww_stat(1,:,:) ! Favre u_clean
 ww_stat(3,:,:) = w_stat2(ii1:ii2,:, 5)                ! p/p_infty
 ww_stat(4,:,:) = w_stat2(ii1:ii2,:,11) - w_stat2(ii1:ii2,:,5)**2 ! (sigma_p/p_infty)^2
 ww_stat(4,:,:) = SQRT(ww_stat(4,:,:))
 DEALLOCATE(w_stat2)
 CLOSE(11)

 IF (masterproc) WRITE(*,*) 'End read clean stats.'

 !------------------------------------------------------------------------------
 ! Skin friction
 !------------------------------------------------------------------------------
 
 IF (masterproc) WRITE(*,*) 'Start friction.' 

 t_wall    = 1._rkind + 0.5_rkind * gm1 * (pr)**(1./3.) * mach**2 ! Adim temp = recovery temp
 s_t       = s0_dim/t_ref_dim
 f_t_wall  = t_wall**1.5_rkind * (1._rkind+s_t)/(t_wall+s_t)
 cf_c      = mu0 * f_t_wall / (0.5_rkind * u_infty**2)
 dy        = (-22._rkind*yg(1)+36._rkind*yg(2) &
         &    -18._rkind*yg(3)+ 4._rkind*yg(4))/12._rkind

 ALLOCATE(dudy_w_bl(nx), dudy_w(nx,nzs+1), dwdy_w(nx,nzs+1))     

 ! Clean wall streamwise vel
 dudy_w_bl  = (-22._rkind*ww_stat(2,:,1)+36._rkind*ww_stat(2,:,2) &
         &     -18._rkind*ww_stat(2,:,3)+ 4._rkind*ww_stat(2,:,4))/12._rkind
 dudy_w_bl  = cf_c*dudy_w_bl/dy 
 IF (masterproc) WRITE(*,*) 'Done cf_dudy clean.' 

 ! Ramp streamwise vel 
 dudy_w = (-22._rkind*w_sym(2,:,1,:)+36._rkind*w_sym(2,:,2,:) &
         & -18._rkind*w_sym(2,:,3,:)+ 4._rkind*w_sym(2,:,4,:))/12._rkind
 dudy_w = cf_c*dudy_w/dy 
 IF (masterproc) WRITE(*,*) 'Done cf_dudy ramp.' 

 ! Ramp spanwise vel
 dwdy_w = (-22._rkind*w_sym(3,:,1,:)+36._rkind*w_sym(3,:,2,:) &
         & -18._rkind*w_sym(3,:,3,:)+ 4._rkind*w_sym(3,:,4,:))/12._rkind
 dwdy_w = cf_c*dwdy_w/dy 
 IF (masterproc) WRITE(*,*) 'Done cf_dwdy ramp.' 

 IF (masterproc) WRITE(*,*) 'End friction.' 

 !------------------------------------------------------------------------------
 ! Boundary layer thickness
 !------------------------------------------------------------------------------

 IF (masterproc) WRITE(*,*) 'Start boundary layer thicknesses.'

 ygg = yg(nymax)

 ! Clean boundary layer
 ALLOCATE(delta_inc_bl(nx), theta_inc_bl(nx), delta_99_bl(nx))
 DO i = 1, nx
   int1  = 0._rkind
   int2  = 0._rkind
   DO j = 2, ny
     dy      = dyg(j-1)
     up      = ww_stat(2,i,  j)/u_infty
     um      = ww_stat(2,i,j-1)/u_infty
     u2p     = up*up
     u2m     = um*um
     int1    = int1 + 0.5_rkind*( up+ um)*dy
     int2    = int2 + 0.5_rkind*(u2p+u2m)*dy
   ENDDO
   delta_inc_bl(i) = ygg - int1
   theta_inc_bl(i) = int1  - int2

   j99    = 1
   do j=1,ny-1
    uu = ww_stat(2,i,j)/u_infty
    if (uu>0.99_rkind) then
     j99 = j-1
     exit
    endif
   enddo
   dely     = yg(j99+1)-yg(j99)
   unum     = 0.99_rkind*u_infty-ww_stat(2,i,j99)
   uden     = ww_stat(2,i,j99+1)-ww_stat(2,i,j99)
   delta_99_bl(i) = yg(j99)+dely*(unum/uden)

 ENDDO
 IF (masterproc) WRITE(*,*) 'End clean blay.'

 ! Ramp
 ALLOCATE(delta_inc(nx,nzs+1), theta_inc(nx, nzs+1), delta_99(nx, nzs+1))
 DO i = 1, nx
   DO k = 1, nzs+1
     int1  = 0._rkind
     int2  = 0._rkind
     DO j = 2, ny
       dy      = dyg(j-1)
       up      = w_sym(2,i,  j,k)/u_infty
       um      = w_sym(2,i,j-1,k)/u_infty
       u2p     = up*up
       u2m     = um*um
       int1    = int1 + 0.5_rkind*( up+ um)*dy
       int2    = int2 + 0.5_rkind*(u2p+u2m)*dy
     ENDDO
     delta_inc(i,k) = ygg - int1
     theta_inc(i,k) = int1  - int2
      
     j99    = 1
     do j=1,ny-1
      uu = w_sym(2,i,j,k)/u_infty
      if (uu>0.99_rkind) then
       j99 = j-1
       exit
      endif
     enddo
     dely     = yg(j99+1)-yg(j99)
     unum     = 0.99_rkind*u_infty-w_sym(2,i,j99,k)
     uden     = w_sym(2,i,j99+1,k)-w_sym(2,i,j99,k)
     delta_99(i,k) = yg(j99)+dely*(unum/uden)

   ENDDO
 ENDDO
 IF (masterproc) WRITE(*,*) 'End ramp.'

 IF (masterproc) WRITE(*,*) 'End boundary layer thicknesses.'

 CALL MPI_BARRIER(mpi_comm_world,iermpi)
 
 !------------------------------------------------------------------------------
 ! Gather output in single array
 !------------------------------------------------------------------------------

 ALLOCATE(output_array(15, nx, 1, nzs+1))
 DO k = 1, nzs+1
  output_array(1,:,1,k) = dudy_w_bl
  output_array(2,:,1,k) = delta_inc_bl
  output_array(3,:,1,k) = theta_inc_bl
  output_array(4,:,1,k) = delta_99_bl
  output_array(5,:,1,k) = ww_stat(1,:,1)
  output_array(6,:,1,k) = ww_stat(3,:,1)
  output_array(7,:,1,k) = ww_stat(4,:,1)
 ENDDO
 output_array( 8,:,1,:) = dudy_w
 output_array( 9,:,1,:) = dwdy_w
 output_array(10,:,1,:) = delta_inc
 output_array(11,:,1,:) = theta_inc
 output_array(12,:,1,:) = delta_99
 output_array(13,:,1,:) = w_sym(1,:,1,:)
 output_array(14,:,1,:) = w_sym(4,:,1,:)
 output_array(15,:,1,:) = w_sym(5,:,1,:)
 DEALLOCATE(w_sym, ww_stat)
 DEALLOCATE(dudy_w_bl,    dudy_w,       dwdy_w     )
 DEALLOCATE(delta_inc,    theta_inc,    delta_99   )
 DEALLOCATE(delta_inc_bl, theta_inc_bl, delta_99_bl)

 IF (masterproc) WRITE(*,*) 'End friction.'

 CALL MPI_BARRIER(mpi_comm_world,iermpi)

 !------------------------------------------------------------------------------
 ! Finalize and deallocate
 !------------------------------------------------------------------------------
 ! WRITE OUTPUT
 
 IF (masterproc) WRITE(*,*) '!------------START WRITING OUTPUT------------'
 
 IF (masterproc) WRITE(*,*) 'Writing blay properties on xz plane.'
 CALL write_vtk_general(ig_start, k_start, 15, nx, 1, nzs+1, &
         & output_array(:, 1:nx, 1, 1:nzs+1), 'blpr', names, &
         & nblocks, ncoords, mp_cart)

 IF (masterproc) WRITE(*,*) '!-------------END WRITING OUTPUT-------------'

 call mpi_finalize(iermpi)

END PROGRAM postpro_blay_serial_stats
