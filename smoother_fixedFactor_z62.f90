! compile this program with gfortran using
! gfortran -o Smooth.x smooth.f90 -O3 -fopenmp -mcmodel=medium -I/usr/include -L/usr/lib -lnetcdff -Wl,-Bsymbolic-functions -Wl,-z,relro -lnetcdf
! Everything after mcmodel=medium can be replaced by the output of nf-config --flibs and nf-config --fflags
! On cheyenne you'll need to load modules gnu and netcdf
! Before running the compiled executable, put the following at your bash prompt (or bash script): export OMP_NUM_THREADS=16
! export OMP_STACKSIZE=1024M
! ulimit -s unlimited
! Change 16 to however many cores you have available on the node.
program main 
  use netcdf
  implicit none

  ! We are reading 3D data, a 2400 x 3600 x 62 grid. 
  integer, parameter :: NX = 3600, NY = 2400, NZ = 62
  integer, parameter :: rec_length = 8*NX*NY ! This only works for gfortran, not ifort
  double precision, allocatable, dimension(:,:,:) :: smoothed, field, fieldT
  double precision, allocatable, dimension(:,:,:) :: TEMP
  double precision, dimension(-8:NX+9,NY,NZ) :: mask
  double precision, dimension(-8:NY+9,NX,NZ) :: maskT
  double precision :: weights(19)
  integer :: i, j, k ! Loop indices
  integer :: KMT(NX,NY)

  ! This will be the netCDF ID for the file and data variable.
  integer :: ncid, varid

  ! Some timing variables
  integer :: tCount0, tCount1, tCountRate

! Read smaller auxiliary fields
  ! Open the file. NF90_NOWRITE tells netCDF we want read-only access to
  ! the file.
  call check( nf90_open("current_data", NF90_NOWRITE, ncid) )

  ! Read KMT
  call check( nf90_inq_varid(ncid, "KMT", varid) )
  call check( nf90_get_var(ncid, varid, KMT) )

  call system_clock(tCount0,tCountRate)
! Read TEMP
  allocate(TEMP(NX,NY,NZ)) 
  call check( nf90_inq_varid(ncid, "TEMP", varid) )
  call check( nf90_get_var(ncid, varid, TEMP) )
  ! Close the file, freeing all resources.
  call check( nf90_close(ncid) )
  call system_clock(tCount1)
  print *, 'Reading TEMP: ', real(tCount1-tCount0)/real(tCountRate)
  allocate(field(-8:NX+9,NY,NZ))
  do k=1,NZ
    field(1:NX,:,k) = TEMP(:,:,k)
  end do
  deallocate(TEMP)

! Time the whole computation
  call system_clock(tCount0)
  ! set weights
  do i=-9,9
    weights(i+10) = exp( -(3./50.)*(i**2) )
  end do
  weights = weights/sum(weights)

  where( KMT == -1 )
    KMT = 0
  end where
! Set land mask, and values on land to 0
  mask = 1.
  do k=1,NZ
    do j=1,NY
      do i=1,NX
        if ( KMT(i,j) < k ) then
          mask(i,j,k) = 0.
          field(i,j,k) = 0. ! These values not used
        end if
      end do
    end do
    call per_ext(mask(:,:,k))
    call per_ext(field(:,:,k))
    maskT(1:NY,:,k) = transpose(mask(1:NX,:,k))
    call tri_ext(maskT(:,:,k))
  end do

! Apply filter over longitude
  allocate(smoothed(NX,NY,NZ))
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,k)
  do k=1,NZ
    do j=1,NY
      do i=1,NX
        ! Only over ocean 
        if( KMT(i,j) >= k ) then
          smoothed(i,j,k) = sum(weights*mask(i-9:i+9,j,k)*field(i-9:i+9,j,k))/sum(weights*mask(i-9:i+9,j,k))
        end if
      end do
    end do
  end do
!$OMP END PARALLEL DO
  deallocate(field)

! Apply filter over latitude
  allocate(fieldT(-8:NY+9,NX,NZ))
  do k=1,NZ
    fieldT(1:NY,:,k) = transpose(smoothed(:,:,k))
    call tri_ext(fieldT(:,:,k))
  end do
  deallocate(smoothed)
  allocate(smoothed(NY,NX,NZ))

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,k)
  do k=1,NZ
    do i=1,NX
      do j=1,NY
        ! Only over ocean 
        if( KMT(i,j) >= 1  ) then
          smoothed(j,i,k) = sum(weights*maskT(j-9:j+9,i,k)*fieldT(j-9:j+9,i,k))/sum(weights*maskT(j-9:j+9,i,k))
        end if
      end do
    end do
  end do
!$OMP END PARALLEL DO
  call system_clock(tCount1)
  print *, 'Total computing time: ', real(tCount1-tCount0)/real(tCountRate)
  deallocate(fieldT)

! Write smoothed TEMP
!  open(unit=21,file='TEMPSmooth3D.dat',access='DIRECT',form='UNFORMATTED',status='UNKNOWN',RECL=rec_length)
!  do k=1,1
!    write(21,REC=k) transpose(smoothed(1:NY,:,k))
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

  subroutine per_ext(dummy)
    double precision, dimension(-8:NX+9,NY), intent(inout) :: dummy
    dummy(-8:0,:) = dummy(NX-8:NX,:)
    dummy(NX+1:NX+9,:) = dummy(1:9,:)
  end subroutine per_ext

  subroutine tri_ext(dummy)
    double precision, dimension(-8:NY+9,NX), intent(inout) :: dummy
    integer :: i, j
    do j=1,9
      do i=1,NX
        dummy(1-j,i) = dummy(j,NX+1-i)
        dummy(NY+j,i) = dummy(NY+1-j,NX+1-i)
      end do
    end do
  end subroutine tri_ext

end program main 
