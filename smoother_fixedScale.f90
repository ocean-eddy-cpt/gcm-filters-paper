! compile this program with gfortran using
! gfortran -o Smooth.x smooth.f90 -O3 -fopenmp -mcmodel=medium -I/usr/include -L/usr/lib -lnetcdff -Wl,-Bsymbolic-functions -Wl,-z,relro -lnetcdf
! Everything after mcmodel=medium can be replaced by the output of nf-config --flibs and nf-config --fflags
! On cheyenne you'll need to load modules gnu and netcdf
! Before running the compiled executable, put the following at your bash prompt (or bash script): export OMP_NUM_THREADS=16
! export OMP_STACKSIZE=1536M
! ulimit -s unlimited
! Change 16 to however many cores you have available on the node.
! To read PDUSmooth3D.dat in matlab: fid = fopen('PDUSmooth3D.dat','r');tmp = fread(fid,1/0,'double');fclose(fid); clear ans fid; PDUSmooth=reshape(tmp,[2400 3600 62]);clear tmp
program main
  use netcdf
  implicit none

  ! We are reading 3D data, a 2400 x 3600 x 62 grid. 
  integer, parameter :: NX = 3600, NY = 2400, NZ = 62
  integer, parameter :: rec_length = 8*NX*NY ! This only works for gfortran, not ifort
  double precision, allocatable, dimension(:,:) :: smoothed, field, fieldT
  double precision, allocatable, dimension(:,:,:) :: TEMP
  double precision, dimension(NX,NY) :: mask
  double precision, dimension(NY,NX) :: maskT
  double precision, dimension(NX,NY) :: TAREA, TLON, TLAT
  double precision :: weightsX(NX), weightsY(NY), dLONX(NX), dLONY(NY)
  double precision :: d2r, r2d, pi, TMPX(NX), TMPY(NY)
  integer :: i, j ! Loop indices
  integer :: KMT(NX,NY)

  ! This will be the netCDF ID for the file and data variable.
  integer :: ncid, varid

  ! Some timing variables
  integer :: tCount0, tCount1, tCountRate

! Read smaller auxiliary fields
  ! Open the file. NF90_NOWRITE tells netCDF we want read-only access to
  ! the file.
  call check( nf90_open("current_data", NF90_NOWRITE, ncid) )
  ! Get the varid of the data variable, based on its name.
  ! Read TAREA
  call check( nf90_inq_varid(ncid, "TAREA", varid) )
  call check( nf90_get_var(ncid, varid, TAREA) )

  ! Read Longitude (degrees)
  call check( nf90_inq_varid(ncid, "TLONG", varid) )
  call check( nf90_get_var(ncid, varid, TLON) )

  ! Read Latitude (degrees)
  call check( nf90_inq_varid(ncid, "TLAT", varid) )
  call check( nf90_get_var(ncid, varid, TLAT) )

  ! Read KMT
  call check( nf90_inq_varid(ncid, "KMT", varid) )
  call check( nf90_get_var(ncid, varid, KMT) )

  call system_clock(tCount0,tCountRate)
! Read TEMP
  allocate(TEMP(NX,NY,NZ)) ! All the alocate/deallocate stuff is more
                           ! useful when dealing with multiple 3D fields; not
                           ! really important here
  call check( nf90_inq_varid(ncid, "TEMP", varid) )
  call check( nf90_get_var(ncid, varid, TEMP) )
  ! Close the file, freeing all resources.
  call check( nf90_close(ncid) )
  call system_clock(tCount1)
  print *, 'Reading TEMP: ', real(tCount1-tCount0)/real(tCountRate)
  allocate(field(-44:NX+45,NY))
  field(1:NX,:) = TEMP(:,:,1)
  deallocate(TEMP)

! Time the whole computation
  call system_clock(tCount0)
  ! degrees to radians
  pi = acos(-1.)
  d2r = pi/180.
  r2d = 180./pi

  TAREA = TAREA/MAXVAL(TAREA)
  TLON = TLON*d2r
  TLAT = TLAT*d2r

  where( KMT == -1 )
    KMT = 0
  end where

! Set land mask
  mask = 1.
  do j=1,NY
    do i=1,NX
      if ( KMT(i,j) < 1 ) then
        mask(i,j) = 0.
      end if
    end do
  end do
  maskT = transpose(mask)

! Set values on land to 0; these values shouldn't be used
  do j=1,NY
    do i=1,NX
      if ( KMT(i,j) < 1 ) then
        field(i,j) = 0.
      end if
    end do
  end do

! Apply filter over longitude
  allocate(smoothed(NX,NY))
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(weightsX,dLONX,TMPX,i,j)
  do j=1,NY
    do i=1,NX
      ! Only over ocean
      if( KMT(i,j) >= 1 ) then
        ! Compute the weights
        dLONX = abs(TLON(i,j)-TLON(:,j))
        dLONX = min(dLONX,2.*pi-dLONX)
        TMPX = sin(TLAT(i,j))*sin(TLAT(:,j)) + cos(TLAT(i,j))*cos(TLAT(:,j))*cos(dLONX)
        where(TMPX > 1.) ! if i=j roundoff puts TMPX>1, which leads to acos(>1)=NaN
          TMPX = 1.
        end where
        weightsX = mask(:,j)*TAREA(:,j)*exp(-6.*(r2d*acos(TMPX))**2)
        smoothed(i,j) = sum(weightsX*field(:,j))/sum(weightsX)
      end if
    end do
  end do
!$OMP END PARALLEL DO
  deallocate(field)

! Aapply filter over latitude
  allocate(fieldT(NY,NX))
  fieldT = transpose(smoothed)
  deallocate(smoothed)
  allocate(smoothed(NY,NX))

!$OMP PARALLEL DO PRIVATE(weightsY,dLONY,TMPY,i,j)
  do i=1,NX
    do j=1,NY
      ! Only over ocean 
      if( KMT(i,j) >= 1  ) then
        ! Compute the weights
        dLONY = abs(TLON(i,j)-TLON(i,:))
        dLONY = min(dLONY,2.*pi-dLONY)
        TMPY = sin(TLAT(i,j))*sin(TLAT(i,:)) + cos(TLAT(i,j))*cos(TLAT(i,:))*cos(dLONY)
        where( TMPY > 1. )
          TMPY = 1.
        end where
        weightsY = maskT(:,i)*TAREA(i,:)*exp(-6.*(r2d*acos(TMPY))**2 )
        smoothed(j,i) = sum(weightsY*fieldT(:,i))/sum(weightsY)
        end if
    end do
  end do
!$OMP END PARALLEL DO
  call system_clock(tCount1)
  print *, 'Total computing time: ', real(tCount1-tCount0)/real(tCountRate)
  deallocate(fieldT)

! Write smoothed TEMP
!  open(unit=21,file='TEMPSmooth3D.dat',access='DIRECT',form='UNFORMATTED',status='UNKNOWN',RECL=rec_length)
!  do k=1,1
!    write(21,REC=k) transpose(smoothed(:,:,k))
!  end do
!  close(21)
!  deallocate(smoothed)

contains
  subroutine check(status)
    integer, intent ( in) :: status
    
    if(status /= nf90_noerr) then 
      print *, trim(nf90_strerror(status))
      stop "Stopped"
    end if
  end subroutine check  

end program main
