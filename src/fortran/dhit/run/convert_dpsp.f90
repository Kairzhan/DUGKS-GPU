program main
implicit none
real(8), allocatable, dimension(:,:,:,:) :: f
real(8), allocatable, dimension(:,:,:) :: u, v, w, rho
real(4), allocatable, dimension(:,:,:,:) :: f1
real(4), allocatable, dimension(:,:,:) :: u1, v1, w1, rho1
integer, parameter :: nx=256
integer, parameter :: ny=256
integer, parameter :: nz=256
integer, parameter :: npop=19
integer :: i, j, k, ip
character(16) :: filename

allocate(f(0:npop-1, nx, ny, nz))
allocate(u(nx, ny, nz))
allocate(v(nx, ny, nz))
allocate(w(nx, ny, nz))
allocate(rho(nx, ny, nz))

allocate(f1(0:npop-1, nx, ny, nz))
allocate(u1(nx, ny, nz))
allocate(v1(nx, ny, nz))
allocate(w1(nx, ny, nz))
allocate(rho1(nx, ny, nz))

write(filename, "('duv',i9.9,'.dat')") 20050
open(unit=100, file=filename, status="old", access="stream")
do k=1,nz
   do j=1,ny
      do i=1,nx
         read(100) rho(i,j,k), u(i,j,k), v(i,j,k), w(i,j,k)
      enddo
   enddo
enddo
close(100)

write(filename, "('ddf',i9.9,'.dat')") 20050
open(unit=100, file=filename, status="old", access="stream")
read(100) f(ip, i,j,k)
!do k=1,nz
!   do j=1,ny
!      do i=1,nx
!         do ip=0, npop-1
!           read(100) f(ip, i,j,k)
!         enddo
!      enddo
!   enddo
!enddo
close(100)

f1=f
u1=u
v1=v
w1=w
rho1=rho

write(filename, "('suv',i9.9,'.dat')") 20050
open(unit=100, file=filename, status="replace", access="stream")
do k=1,nz
   print *, "k=", k
   do j=1,ny
      do i=1,nx
         write(100) rho1(i,j,k), u1(i,j,k), v1(i,j,k), w1(i,j,k)
      enddo
   enddo
enddo
close(100)

write(filename, "('sdf',i9.9,'.dat')") 20050
open(unit=100, file=filename, status="replace", access="stream")
write(100) f1
!do k=1,nz
!   do j=1,ny
!      do i=1,nx
!         do ip=0, npop-1
!           write(100) f1(ip, i,j,k)
!         enddo
!      enddo
!   enddo
!enddo
close(100)

endprogram main
