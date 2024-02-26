PROGRAM postpro_timespanav_slice_bl
!
! This program reads the stat.bin file (averaged in time and span) 
! and writes a binary file containing the xy grid and 
! (1) rho_bar/rho_infty
! (2) u_tilde/u_infty
! (3) v_tilde/u_infty, 
! (4) w_tilde/u_infty,
! (5) p_bar/p_infty
! 
! bar = Reynolds average, tilde = Favre average
!
! INPUT FILES:
! ./input_subvolume_serial.dat
!     nxmax, nymax, nzmax     Streamwise/wall-normal/spanwise number of points
!     nprocx, nprocy, nprocz      "           "          "    number of blocks
!     re_delta                Reynolds number from simulation "flow_params.dat"
!     mach                    Mach number      "        "           " 
!     ig_start                Minimum global streamwise  index of subvolume of interest
!     jg_end                  Maximum    "    wall-normal  "   "       "     "    "
! ./dxg.dat, ./dyg.dat        Grid files from simulation
! ./stat.bin                  Stats file   "     " 
!
! OUTPUT FILES: 
! time_span_av_plane_infty_norm_bl.bin
!
! 24/01/2024 G. Della Posta
!

USE, intrinsic :: iso_fortran_env!, only : error_unit
USE iso_c_binding

IMPLICIT NONE

!------------------------------------------------------------------------------

INTEGER, PARAMETER :: singtype = selected_real_kind(6,37)    ! single precision
INTEGER, PARAMETER :: doubtype = selected_real_kind(15,307)  ! double precision
INTEGER, PARAMETER :: rkind    = doubtype
INTEGER, PARAMETER :: ikind    = INT32
!
INTEGER, PARAMETER :: nv_stat  = 70 
INTEGER :: nxmax, nymax, nzmax, ig_start, jg_end
INTEGER :: nprocx, nprocy, nprocz, nxserial
INTEGER :: i_block_start
INTEGER :: nxsmax, nysmax
INTEGER :: i, j

REAL(rkind), PARAMETER :: gamma = 1.4_rkind
REAL(rkind) :: re_delta, u_infty, mach

REAL(rkind), ALLOCATABLE, DIMENSION(:) :: xg, yg

REAL(rkind), ALLOCATABLE, DIMENSION(:,:,:) :: w_stat

!------------------------------------------------------------------------------

 WRITE(*,*) '!---------END PROGRAM POSTPRO_TIMESPANAV_SLICE_BL--------'

 OPEN(12, FILE='./input_subvolume_serial.dat', ACTION='read', STATUS='old')
 READ(12,*)
 READ(12,*) nxmax, nymax, nzmax, nprocx, nprocy, nprocz
 READ(12,*)
 READ(12,*)
 READ(12,*) re_delta, mach, ig_start, jg_end
 CLOSE(12)

 u_infty = sqrt(gamma)*mach

 nxserial      = nxmax/nprocx
 i_block_start = ig_start/nxserial
 ig_start      = nxserial*i_block_start + 1

 nxsmax  = nxmax  - ig_start + 1
 nysmax  = jg_end 

!------------------------------------------------------------------------------
! Read grid
!------------------------------------------------------------------------------

 allocate(xg(nxmax))
 allocate(yg(nymax))

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

!------------------------------------------------------------------------------
! Import time and spanwise averaged solution (stat.bin)
!------------------------------------------------------------------------------

 WRITE(*,*) 'Start read clean stats.'

 ALLOCATE(w_stat(nxmax,nymax,nv_stat))

 OPEN(11, FILE='./stat.bin', FORM='unformatted', STATUS='old', ACTION='read', ACCESS='stream')
 READ(11) w_stat
 CLOSE(11)

 !w_stat(:,:,1)                                       ! rho_bar/rho_infty
 w_stat(:,:,2) = w_stat(:,:,13)/w_stat(:,:,1)/u_infty ! u_tilde/u_infty 
 w_stat(:,:,3) = w_stat(:,:,14)/w_stat(:,:,1)/u_infty ! v
 w_stat(:,:,4) = w_stat(:,:,15)/w_stat(:,:,1)/u_infty ! w
 w_stat(:,:,5) = w_stat(:,:,5)                        ! p_bar/p_infty

 WRITE(*,*) 'End read clean stats.'
 WRITE(*,*) 'Writing primitive variables.'

 OPEN(666,FILE='time_span_av_plane_infty_norm_bl.bin', &
     & ACTION='write', STATUS='new', FORM='unformatted')
 WRITE(666) nxsmax, nysmax, 5
 WRITE(666) xg(ig_start:ig_start+nxsmax-1)
 WRITE(666) yg(1:nysmax)
 DO j = 1, 5
   WRITE(666) w_stat(ig_start:ig_start+nxsmax-1,1:nysmax, j)
 ENDDO
 CLOSE(666)
 DEALLOCATE(w_stat)

 WRITE(*,*) 'End write clean stats.'

 WRITE(*,*) '!---------END PROGRAM POSTPRO_TIMESPANAV_SLICE_BL--------'

END PROGRAM postpro_timespanav_slice_bl
