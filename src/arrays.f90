! -*- F90 -*-

module arrays
  implicit none

  ! fortran array pointer
  real*8, pointer, save :: cord(:,:,:), temp(:,:), vel(:,:,:), stress0(:,:,:,:), &
       force(:,:,:), balance(:,:,:), amass(:,:), rmass(:,:), &
       area(:,:,:), dvol(:,:,:), strain(:,:,:), bc(:,:,:)

  integer, pointer, save :: ncod(:,:,:)

  ! theromchronology
  real*8, pointer, save :: chron_if(:,:,:), chron_temp(:,:,:), &
      chron_time(:,:,:), max_temp(:,:)

  integer, pointer, save :: nmark(:,:)

  ! temporary array
  real*8, pointer, save :: junk2(:,:)


contains

  subroutine allocate_arrays(nz, nx)
    implicit none

    integer, intent(in) :: nz, nx

    allocate(cord(nz, nx, 2))
    allocate(temp(nz, nx))
    allocate(vel(nz, nx, 2))
    allocate(stress0(nz-1, nx-1, 4, 4))
    allocate(force(nz, nx, 2))
    allocate(balance(nz, nx, 2))
    allocate(amass(nz, nx))
    allocate(rmass(nz, nx))
    allocate(area(nz-1, nx-1, 4))
    allocate(dvol(nz-1, nx-1, 4))
    allocate(strain(nz-1, nx-1, 3))
    allocate(bc(nz, nx, 2))

    allocate(ncod(nz, nx, 2))

    allocate(max_temp(nz,nx))
    allocate(nmark(nz-1,nx-1))
    allocate(chron_if(3,nz-1,nx-1))
    allocate(chron_temp(3,nz-1,nx-1))
    allocate(chron_time(3,nz-1,nx-1))

    allocate(junk2(nz, nx))

  end subroutine allocate_arrays

end module arrays
