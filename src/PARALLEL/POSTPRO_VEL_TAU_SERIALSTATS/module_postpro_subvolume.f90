module modpostpro

 use mpi
 use, intrinsic :: iso_fortran_env!, only : error_unit
 use iso_c_binding
 implicit none
 save
!
 integer, parameter :: singtype = selected_real_kind(6,37)    ! single precision
 integer, parameter :: doubtype = selected_real_kind(15,307)  ! double precision
 integer, parameter :: rkind    = doubtype
 integer, parameter :: ikind    = INT32
 integer, parameter :: mpi_prec = mpi_real8
!
 integer, parameter :: ndims            = 3
 integer, parameter :: i_stat_3d_vtk    = 0    
 integer, parameter :: i_stat_3d_PLOT3D = 0    

 integer :: nv_stat  = 19 ! 70 
 integer :: nxmax, nymax, nzmax
 integer :: nx,ny,nz,nv,nv_stat_3d, j_block_start, j_block_end, i_block_start, i_block_end
!
 logical :: masterproc
!
! MPI related parameters 
!
 integer :: iermpi
 integer :: nrank,nproc
!
 integer, dimension(2) :: dims
 integer, dimension(3) :: nblocks
 logical, dimension(3) :: pbc
 integer, dimension(mpi_status_size) :: istatus
!
 integer, dimension(3) :: ncoords
 integer :: mp_cart
!
 real(rkind), parameter :: gamma     = 1.4_rkind
 real(rkind), parameter :: s0_dim    = 111.0_rkind
 real(rkind), parameter :: t_ref_dim = 160.0_rkind
 real(rkind), parameter :: pr        = 0.72_rkind
 real(rkind), parameter :: gm1       = gamma-1._rkind
 real(rkind), parameter :: gm        = 1._rkind/gm1
 real(rkind), parameter :: pi2       = 2.*ATAN(1._rkind)
 !real(rkind), parameter :: theta_r   = 8.64_rkind * pi2 / 90._rkind
 !real(rkind), parameter :: theta_s   = 24.0_rkind * pi2 / 90._rkind
 !real(rkind), parameter :: domain_semiwidth = 5._rkind
 !real(rkind), parameter :: mach      = 2.0_rkind
 !real(rkind), parameter :: u0        = SQRT(gamma)*mach
 !real(rkind), parameter :: xg_edge   = 60._rkind

 real(rkind) :: theta_r, theta_s, domain_semiwidth
 real(rkind) :: mach, u0, xg_edge
 real(rkind) :: re_delta, yg_edge, h_delta
 real(rkind) :: dxg,dzg
 real(rkind) :: ramp_length, ramp_semiwidth
 real(rkind) :: ramp_x_start, ramp_z_start
!
 real(rkind), allocatable, dimension(:)       :: xg,yg,zg,dyg

 real(rkind), allocatable, dimension(:,:,:)   :: wmean, corrz
 real(rkind), allocatable, dimension(:,:,:)   :: w_stat

 real(rkind), allocatable, dimension(:,:,:,:) :: w, w_stat_3d
!
CONTAINS 

  function int2str(int_num)
       implicit none
       integer :: int_num
       character(len=16) :: int2str, ret_value
       write(ret_value, "(I0)") int_num
       int2str = ret_value
  endfunction int2str

  function int2str_o(int_num)
       use mpi
       implicit none
       integer(KIND=MPI_OFFSET_KIND) :: int_num
       character(len=32) :: int2str_o, ret_value
       write(ret_value, "(I0)") int_num
       int2str_o = ret_value
  endfunction int2str_o
  !
  !-------------------------------------------------------------------------------------------------
  !-------------------------------------------------------------------------------------------------
  !
  SUBROUTINE write_vtk_reduced_extended_x(nv_io, w_aux_io, file_str, names)
  
   implicit none
  
   integer, parameter :: int64_kind = selected_int_kind(2*range(1))
   
   INTEGER, INTENT(in) :: nv_io
   integer :: istore = 1
   integer :: l
   integer :: mpi_io_file
   integer :: filetype
   integer :: size_real, size_integer, size_integer8
   integer :: ntot
   integer(kind=mpi_offset_kind)  :: offset
   integer(kind=mpi_offset_kind)  :: offset_x,offset_y,offset_z,delta_offset_w
   integer(kind=mpi_address_kind) :: stride
   integer(int64_kind) :: gridsize_64
   
   INTEGER :: comm_print, nrank_print, nproc_print, color, group
  
   integer, dimension(3) :: sizes     ! Dimensions of the total grid
   integer, dimension(3) :: subsizes  ! Dimensions of grid local to a procs
   integer, dimension(3) :: starts    ! Starting coordinates
   
   real(rkind), dimension(1:nx,1:ny,1:nz) :: w_tmp
  
   REAL(rkind), DIMENSION(1:nv_io,1:nx,1:ny,1:nz), INTENT(in) :: w_aux_io
  
   CHARACTER(len=    4), INTENT(in) :: file_str
   character(len=    4) :: chstore
   character(len=    7) :: vtk_float
   character(len=   32) :: file_prefix_
   character(len=65536) :: xml_part
   
   CHARACTER(len = 12), DIMENSION(nv_io), INTENT(in) :: names
  
  !-------------------------------------------------------------------------------------------------
  
   size_real = storage_size(1._rkind)/8
   size_integer = storage_size(nxmax)/8
   size_integer8 = storage_size(gridsize_64)/8
   
   write(chstore(1:4), '(I4.4)') istore
   
   if (masterproc) print *, 'Storing VTK stat.'
   
   call MPI_TYPE_SIZE(mpi_prec,size_real,iermpi)
   if(size_real == 4) then
    vtk_float = "Float32"
   elseif(size_real == 8) then
    vtk_float  = "Float64"
   else
    if(masterproc) write(*,*) "Error on VTK write! size_real must be either 4 or 8"
    call MPI_ABORT(MPI_COMM_WORLD,iermpi,iermpi)
   endif
   gridsize_64 = int(size_real,int64_kind)*int(nxmax,int64_kind)*int(nymax,int64_kind)*&
                 int(nzmax,int64_kind)
   
   if(storage_size(gridsize_64) /= 64) then
    if(masterproc) write(*,*) "Error on VTK write! Size of int64_kind integers is not 8 bytes!"
    call MPI_ABORT(MPI_COMM_WORLD, iermpi, iermpi)
   endif
   
   file_prefix_ = file_str//"_"
   
   if (masterproc) then
    offset_x = 0
    offset_y = size_real*nxmax + storage_size(gridsize_64)/8
    offset_z = offset_y + size_real*nymax + storage_size(gridsize_64)/8
    delta_offset_w = gridsize_64 + storage_size(gridsize_64)/8 ! the second part is because of the header of bytes before data
   
    open(unit=123, file=trim(file_prefix_)//chstore//'.vtr', access="stream", form="unformatted", status="replace")
   
    xml_part = ' &
     & <?xml version="1.0"?> &
     & <VTKFile type="RectilinearGrid" version="1.0" byte_order="LittleEndian" header_type="UInt64"> &
     &  <RectilinearGrid WholeExtent="+1 +'&
     &   //int2str(nxmax)//' +1 +'//int2str(nymax)//' +1 +'//int2str(nzmax)//'"> &
     &   <Piece Extent="+1 +'//int2str(nxmax)//&
     &    ' +1 +'//int2str(nymax)//' +1 +'//int2str(nzmax)//'"> &
     &    <Coordinates> &
     &     <DataArray type="'//vtk_float//'" NumberOfComponents="1" Name="X" format="appended" offset="'//&
     &       int2str_o(offset_x)//'"/> &
     &     <DataArray type="'//vtk_float//'" NumberOfComponents="1" Name="Y" format="appended" offset="'//&
     &       int2str_o(offset_y)//'"/> &
     &     <DataArray type="'//vtk_float//'" NumberOfComponents="1" Name="Z" format="appended" offset="'//&
     &       int2str_o(offset_z)//'"/> &
     &    </Coordinates> &
     &   <PointData> '
   
    offset = offset_z + size_real*nzmax + storage_size(gridsize_64)/8
    do l=1,nv_io
       xml_part = trim(adjustl(xml_part)) // ' &
      & <DataArray type="'//vtk_float//'" NumberOfComponents="1" Name="'//trim(names(l))//'" format="appended" &
      &  offset="'//int2str_o(offset)//'"/>'
       offset = offset + delta_offset_w
    enddo
   
    xml_part = trim(adjustl(xml_part)) // ' &
    &       </PointData> &
    &     </Piece> &
    &   </RectilinearGrid> &
    &   <AppendedData encoding="raw"> '
   
    write(123) trim(adjustl(xml_part))
   
    write(123) "_"
    write(123) size_real*int(nxmax,int64_kind) , xg(1:nxmax)
    write(123) size_real*int(nymax,int64_kind) , yg(1:nymax)
    write(123) size_real*int(nzmax,int64_kind) , zg(1:nzmax)
    flush(123)
    close(123)
   endif
   
   sizes(1) = nblocks(1)*nx
   sizes(2) = nblocks(2)*ny
   sizes(3) = nblocks(3)*nz

   !subsizes(1) = 0
   !IF ((nrank >= i_block_start) .AND. (nrank <= i_block_end)) THEN
   !  subsizes(1) = nx
   !ENDIF
   subsizes(1) = nx
   subsizes(2) = ny
   subsizes(3) = nz

   !IF ((nrank >= i_block_start) .AND. (nrank <= i_block_end)) THEN
   ! starts(1) = 0 + (ncoords(1)-i_block_start)*subsizes(1) 
   !ELSE
   ! starts(1) = 0 !+ i_block_start*nx
   !ENDIF
   starts(1) = 0 + ncoords(1)*subsizes(1)
   starts(2) = 0 + ncoords(2)*subsizes(2)
   starts(3) = 0 + ncoords(3)*subsizes(3)

   ntot = subsizes(1)*subsizes(2)*subsizes(3)
   !WRITE(*,*) 'i_block_start, i_block_end, nrank, sizes, subsizes, starts = ', &
   !  & i_block_start, i_block_end, nrank, sizes, subsizes, starts
  
   !stride = sizes(1)*8_mpi_address_kind
   !CALL MPI_TYPE_CREATE_HVECTOR(ny*nz, subsizes(1), stride, mpi_prec, filetype, iermpi)
   call MPI_TYPE_CREATE_SUBARRAY(3,sizes,subsizes,starts,MPI_ORDER_FORTRAN,mpi_prec,filetype,iermpi)
   call MPI_TYPE_COMMIT(filetype,iermpi)
   
   do l=1,nv_io
    call MPI_BARRIER(mp_cart,iermpi)
    if (masterproc) then
     open(unit=123, file=trim(file_prefix_)//chstore//'.vtr', access="stream", form="unformatted", position="append")
     write(123) gridsize_64
     flush(123)
     close(123)
    endif
    call MPI_BARRIER(mp_cart, iermpi)
    call MPI_FILE_OPEN(mp_cart,trim(file_prefix_)//chstore//'.vtr',MPI_MODE_RDWR,MPI_INFO_NULL,mpi_io_file,iermpi)
    call MPI_FILE_OPEN(mp_cart,trim(file_prefix_)//chstore//'.vtr',MPI_MODE_RDWR,MPI_INFO_NULL,mpi_io_file,iermpi)
    call MPI_FILE_GET_SIZE(mpi_io_file, offset, iermpi)
    call MPI_BARRIER(mp_cart, iermpi)
    call MPI_FILE_SET_VIEW(mpi_io_file,offset,mpi_prec,filetype,"native",MPI_INFO_NULL,iermpi)
    w_tmp = w_aux_io(l,1:nx,1:ny,1:nz)
    call MPI_FILE_WRITE_ALL(mpi_io_file,w_tmp,ntot,mpi_prec,istatus,iermpi)
    call MPI_FILE_CLOSE(mpi_io_file,iermpi)
   enddo
   
   call MPI_TYPE_FREE(filetype,iermpi)
   if (masterproc) then
    open(unit=123, file=trim(file_prefix_)//chstore//'.vtr', access="stream", position="append", form="unformatted")
     write(123) ' &
      &    </AppendedData> &
      &  </VTKFile>'
    close(123)
   endif
  
  END SUBROUTINE write_vtk_reduced_extended_x

end module modpostpro
!
!-------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------
