program vanadis2netcdf
   use, intrinsic :: iso_fortran_env, only : real32, real64
   use netcdf
   implicit none

   integer :: t_start, t_end, t, first_t, nfiles, it
   integer :: npts, nx, ny, nz
   integer :: ncid, dim_x, dim_y, dim_z, dim_t
   integer :: var_x, var_y, var_z, var_t, var_c
   integer :: dimids4(4)
   logical :: exists
   character(len=512) :: input_dir, output_file, first_file, fname, arg
   real(real64), allocatable :: xall(:), yall(:), zall(:)
   real(real64), allocatable :: xaxis(:), yaxis(:), zaxis(:)
   real(real32), allocatable :: field4(:,:,:,:)
   real(real64) :: time_one(1)

   ! Defaults: Vanadis OUTPUT directory, 100 s ... 2000 s.
   input_dir  = 'OUTPUT'
   output_file = 'vanadis_00100_02000.nc'
   t_start = 100
   t_end   = 2000

   ! Usage:
   !   vanadis2netcdf [input_dir] [t_start] [t_end] [output.nc]
   ! Example:
   !   vanadis2netcdf OUTPUT 100 2000 vanadis.nc
   call get_command_argument(1, arg)
   if (len_trim(arg) > 0) input_dir = trim(arg)

   call get_command_argument(2, arg)
   if (len_trim(arg) > 0) then
      read(arg, *, err=900) t_start
   end if

   call get_command_argument(3, arg)
   if (len_trim(arg) > 0) then
      read(arg, *, err=900) t_end
   end if

   call get_command_argument(4, arg)
   if (len_trim(arg) > 0) output_file = trim(arg)

   if (t_end < t_start) then
      write(*,'(A)') 'ERROR: t_end must be >= t_start.'
      stop 1
   end if

   ! Find all Vanadis files T_XXXXX.TXT in the requested integer-second range.
   first_t = -1
   nfiles = 0
   do t = t_start, t_end
      call make_filename(input_dir, t, fname)
      inquire(file=trim(fname), exist=exists)
      if (exists) then
         nfiles = nfiles + 1
         if (first_t < 0) then
            first_t = t
            first_file = fname
         end if
      end if
   end do

   if (nfiles == 0) then
      write(*,'(A,I0,A,I0,A,A)') 'ERROR: no T_XXXXX.TXT files found for ', &
           t_start, ' ... ', t_end, ' s in ', trim(input_dir)
      stop 1
   end if

   write(*,'(A,I0)') 'Found time files: ', nfiles
   write(*,'(A,A)')  'First file:       ', trim(first_file)

   ! Read coordinates from the first available file and infer the structured grid.
   call read_first_grid(trim(first_file), npts, xall, yall, zall)
   call infer_grid(npts, xall, yall, zall, nx, ny, nz, xaxis, yaxis, zaxis)

   write(*,'(A,I0)') 'Nodes:             ', npts
   write(*,'(A,I0,A,I0,A,I0)') 'Detected grid:     nx=', nx, ', ny=', ny, ', nz=', nz

   deallocate(xall, yall, zall)
   allocate(field4(nx, ny, nz, 1))

   ! Create NetCDF-4 output.  In the Fortran API x is the fastest varying
   ! in-memory dimension.  NetCDF readers such as xarray normally expose the
   ! variable in conventional time/z/y/x order.
   call nc_check(nf90_create(trim(output_file), ior(nf90_clobber, nf90_netcdf4), ncid), &
                 'nf90_create')

   call nc_check(nf90_def_dim(ncid, 'x',    nx,             dim_x), 'def_dim x')
   call nc_check(nf90_def_dim(ncid, 'y',    ny,             dim_y), 'def_dim y')
   call nc_check(nf90_def_dim(ncid, 'z',    nz,             dim_z), 'def_dim z')
   call nc_check(nf90_def_dim(ncid, 'time', nf90_unlimited, dim_t), 'def_dim time')

   call nc_check(nf90_def_var(ncid, 'x', nf90_double, dim_x, var_x), 'def_var x')
   call nc_check(nf90_def_var(ncid, 'y', nf90_double, dim_y, var_y), 'def_var y')
   call nc_check(nf90_def_var(ncid, 'z', nf90_double, dim_z, var_z), 'def_var z')
   call nc_check(nf90_def_var(ncid, 'time', nf90_double, dim_t, var_t), 'def_var time')

   dimids4 = (/ dim_x, dim_y, dim_z, dim_t /)
   call nc_check(nf90_def_var(ncid, 'concentration', nf90_float, dimids4, var_c, &
                              deflate_level=2, shuffle=.true.), 'def_var concentration')

   call nc_check(nf90_put_att(ncid, var_x, 'long_name', 'Vanadis X coordinate'), 'att x long_name')
   call nc_check(nf90_put_att(ncid, var_x, 'units', 'm'), 'att x units')
   call nc_check(nf90_put_att(ncid, var_x, 'axis', 'X'), 'att x axis')

   call nc_check(nf90_put_att(ncid, var_y, 'long_name', 'Vanadis Y coordinate'), 'att y long_name')
   call nc_check(nf90_put_att(ncid, var_y, 'units', 'm'), 'att y units')
   call nc_check(nf90_put_att(ncid, var_y, 'axis', 'Y'), 'att y axis')

   call nc_check(nf90_put_att(ncid, var_z, 'long_name', 'Vanadis Z coordinate'), 'att z long_name')
   call nc_check(nf90_put_att(ncid, var_z, 'units', 'm'), 'att z units')
   call nc_check(nf90_put_att(ncid, var_z, 'axis', 'Z'), 'att z axis')

   call nc_check(nf90_put_att(ncid, var_t, 'long_name', 'simulation time'), 'att time long_name')
   call nc_check(nf90_put_att(ncid, var_t, 'units', 's'), 'att time units')
   call nc_check(nf90_put_att(ncid, var_t, 'axis', 'T'), 'att time axis')

   call nc_check(nf90_put_att(ncid, var_c, 'long_name', 'pollutant concentration'), &
                 'att concentration long_name')
   call nc_check(nf90_put_att(ncid, var_c, 'units', 'mg m-3'), 'att concentration units')
   call nc_check(nf90_put_att(ncid, var_c, 'coordinates', 'x y z'), 'att concentration coordinates')

   call nc_check(nf90_put_att(ncid, nf90_global, 'title', &
                 'Vanadis 3D atmospheric transport output'), 'global title')
   call nc_check(nf90_put_att(ncid, nf90_global, 'source', &
                 'Vanadis 3D FEM atmospheric transport solver'), 'global source')
   call nc_check(nf90_put_att(ncid, nf90_global, 'converter', &
                 'vanadis2netcdf'), 'global converter')
   call nc_check(nf90_put_att(ncid, nf90_global, 'coordinate_note', &
                 'Coordinates are Vanadis model coordinates; no geodetic CRS is encoded.'), &
                 'global coordinate_note')

   call nc_check(nf90_enddef(ncid), 'nf90_enddef')

   call nc_check(nf90_put_var(ncid, var_x, xaxis), 'write x')
   call nc_check(nf90_put_var(ncid, var_y, yaxis), 'write y')
   call nc_check(nf90_put_var(ncid, var_z, zaxis), 'write z')

   ! Read and write one timestep at a time.  No time-step interval needs to be
   ! supplied: every existing T_XXXXX.TXT file in the requested range is used.
   it = 0
   do t = t_start, t_end
      call make_filename(input_dir, t, fname)
      inquire(file=trim(fname), exist=exists)
      if (.not. exists) cycle

      it = it + 1
      call read_field_checked(trim(fname), nx, ny, nz, xaxis, yaxis, zaxis, &
                              field4(:,:,:,1))

      time_one(1) = real(t, real64)
      call nc_check(nf90_put_var(ncid, var_t, time_one, start=(/it/), count=(/1/)), &
                    'write time')
      call nc_check(nf90_put_var(ncid, var_c, field4, &
                    start=(/1,1,1,it/), count=(/nx,ny,nz,1/)), &
                    'write concentration')

      write(*,'(A,I0,A,A)') 'Written t=', t, ' s  <-  ', trim(fname)
   end do

   call nc_check(nf90_close(ncid), 'nf90_close')

   write(*,'(A)') 'Done.'
   write(*,'(A,A)') 'NetCDF file: ', trim(output_file)
   write(*,'(A,I0)') 'Time levels: ', it
   write(*,'(A,I0,A,I0,A,I0)') 'Grid: nx=', nx, ', ny=', ny, ', nz=', nz
   stop

900 continue
   write(*,'(A)') 'ERROR: invalid command-line integer argument.'
   write(*,'(A)') 'Usage: vanadis2netcdf [input_dir] [t_start] [t_end] [output.nc]'
   stop 1

contains

   subroutine nc_check(status, where)
      integer, intent(in) :: status
      character(len=*), intent(in) :: where
      if (status /= nf90_noerr) then
         write(*,'(A,A)') 'NetCDF ERROR at ', trim(where)
         write(*,'(A)') trim(nf90_strerror(status))
         stop 2
      end if
   end subroutine nc_check


   subroutine make_filename(dir, time_s, name)
      character(len=*), intent(in)  :: dir
      integer,          intent(in)  :: time_s
      character(len=*), intent(out) :: name
      character(len=32) :: base

      if (time_s < 0 .or. time_s > 99999) then
         write(*,'(A,I0)') 'ERROR: filename convention T_XXXXX.TXT cannot represent time ', time_s
         stop 1
      end if

      write(base,'("T_",I5.5,".TXT")') time_s
      if (len_trim(dir) == 0 .or. trim(dir) == '.') then
         name = trim(base)
      else
         name = trim(dir)//'/'//trim(base)
      end if
   end subroutine make_filename


   subroutine read_first_grid(name, n, xa, ya, za)
      character(len=*), intent(in) :: name
      integer, intent(out) :: n
      real(real64), allocatable, intent(out) :: xa(:), ya(:), za(:)

      integer :: u, ios, i
      real(real64) :: xx, yy, zz, ss

      open(newunit=u, file=trim(name), status='old', action='read', iostat=ios)
      if (ios /= 0) then
         write(*,'(A,A)') 'ERROR: cannot open ', trim(name)
         stop 1
      end if

      n = 0
      do
         read(u, *, iostat=ios) xx, yy, zz, ss
         if (ios < 0) exit
         if (ios > 0) then
            write(*,'(A,A)') 'ERROR: malformed record in ', trim(name)
            close(u)
            stop 1
         end if
         n = n + 1
      end do

      if (n <= 0) then
         write(*,'(A,A)') 'ERROR: empty file ', trim(name)
         close(u)
         stop 1
      end if

      rewind(u)
      allocate(xa(n), ya(n), za(n))

      do i = 1, n
         read(u, *, iostat=ios) xa(i), ya(i), za(i), ss
         if (ios /= 0) then
            write(*,'(A,A)') 'ERROR: cannot reread ', trim(name)
            close(u)
            stop 1
         end if
      end do
      close(u)
   end subroutine read_first_grid


   subroutine infer_grid(n, xa, ya, za, nx_o, ny_o, nz_o, xv, yv, zv)
      integer, intent(in) :: n
      real(real64), intent(in) :: xa(n), ya(n), za(n)
      integer, intent(out) :: nx_o, ny_o, nz_o
      real(real64), allocatable, intent(out) :: xv(:), yv(:), zv(:)

      integer :: nxy, i, j, k, p

      ! Vanadis test output is ordered with X varying fastest, then Y, then Z.
      nx_o = 1
      do while (nx_o < n)
         if (.not. same_coord(ya(nx_o+1), ya(1))) exit
         if (.not. same_coord(za(nx_o+1), za(1))) exit
         nx_o = nx_o + 1
      end do

      nxy = 1
      do while (nxy < n)
         if (.not. same_coord(za(nxy+1), za(1))) exit
         nxy = nxy + 1
      end do

      if (mod(nxy, nx_o) /= 0) then
         write(*,'(A)') 'ERROR: first Z-plane is not an integer number of X rows.'
         stop 1
      end if
      ny_o = nxy / nx_o

      if (mod(n, nxy) /= 0) then
         write(*,'(A)') 'ERROR: file is not a complete tensor-product X/Y/Z grid.'
         stop 1
      end if
      nz_o = n / nxy

      if (nx_o * ny_o * nz_o /= n) then
         write(*,'(A)') 'ERROR: inconsistent inferred grid dimensions.'
         stop 1
      end if

      allocate(xv(nx_o), yv(ny_o), zv(nz_o))
      xv = xa(1:nx_o)
      do j = 1, ny_o
         yv(j) = ya(1 + (j-1)*nx_o)
      end do
      do k = 1, nz_o
         zv(k) = za(1 + (k-1)*nx_o*ny_o)
      end do

      ! Full validation: every record must map to x(i), y(j), z(k).
      p = 0
      do k = 1, nz_o
         do j = 1, ny_o
            do i = 1, nx_o
               p = p + 1
               if (.not. same_coord(xa(p), xv(i)) .or. &
                   .not. same_coord(ya(p), yv(j)) .or. &
                   .not. same_coord(za(p), zv(k))) then
                  write(*,'(A,I0)') 'ERROR: non-structured or unexpected node ordering at record ', p
                  write(*,'(A)') 'Expected Vanadis tensor grid with X fastest, then Y, then Z.'
                  stop 1
               end if
            end do
         end do
      end do
   end subroutine infer_grid


   subroutine read_field_checked(name, nx_i, ny_i, nz_i, xv, yv, zv, field)
      character(len=*), intent(in) :: name
      integer, intent(in) :: nx_i, ny_i, nz_i
      real(real64), intent(in) :: xv(nx_i), yv(ny_i), zv(nz_i)
      real(real32), intent(out) :: field(nx_i,ny_i,nz_i)

      integer :: u, ios, i, j, k, p
      real(real64) :: xx, yy, zz, ss

      open(newunit=u, file=trim(name), status='old', action='read', iostat=ios)
      if (ios /= 0) then
         write(*,'(A,A)') 'ERROR: cannot open ', trim(name)
         stop 1
      end if

      p = 0
      do k = 1, nz_i
         do j = 1, ny_i
            do i = 1, nx_i
               read(u, *, iostat=ios) xx, yy, zz, ss
               p = p + 1
               if (ios /= 0) then
                  write(*,'(A,A,A,I0)') 'ERROR: too few or malformed records in ', trim(name), &
                                        ' at record ', p
                  close(u)
                  stop 1
               end if

               if (.not. same_coord(xx, xv(i)) .or. &
                   .not. same_coord(yy, yv(j)) .or. &
                   .not. same_coord(zz, zv(k))) then
                  write(*,'(A,A,A,I0)') 'ERROR: grid in ', trim(name), &
                                        ' differs from first file at record ', p
                  close(u)
                  stop 1
               end if

               ! Preserve Vanadis output exactly in sign.  Do not clip negative
               ! concentrations: the converter is a format converter, not a filter.
               field(i,j,k) = real(ss, real32)
            end do
         end do
      end do

      ! There must not be additional numeric records.
      read(u, *, iostat=ios) xx, yy, zz, ss
      if (ios == 0) then
         write(*,'(A,A)') 'ERROR: extra records found in ', trim(name)
         close(u)
         stop 1
      end if

      close(u)
   end subroutine read_field_checked


   logical function same_coord(a, b)
      real(real64), intent(in) :: a, b
      real(real64), parameter :: atol = 1.0e-8_real64
      real(real64), parameter :: rtol = 1.0e-10_real64
      same_coord = abs(a-b) <= max(atol, rtol*max(abs(a),abs(b)))
   end function same_coord

end program vanadis2netcdf
