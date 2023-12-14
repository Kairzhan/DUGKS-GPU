module stat
  !===========================================================================
contains

  subroutine stat_te(step, dt, xwidth, visc, nx, ny, nz, ux_global, uy_global, uz_global)
    use, intrinsic :: iso_c_binding 
    use dft
    implicit none
    integer, intent(in) :: step
    real, intent(in) :: dt, xwidth, visc
    real :: dx, cwidth, cwidthh
    integer :: nx, ny, nz
    integer, parameter :: prec=4
    real(prec), dimension(nx, ny, nz) :: ux_global
    real(prec), dimension(nx, ny, nz) :: uy_global
    real(prec), dimension(nx, ny, nz) :: uz_global
    real(prec), allocatable, dimension(:,:,:) :: utemp

    integer :: lx, ly, lz
    character(len=16) :: filename
    integer :: i, j, k

    real, parameter :: pi=3.14159, pi2=pi*2

    ! First
    real(prec) :: te

    ! Second
    real(prec), allocatable, dimension(:, :, :) :: et, ets
    real(prec), allocatable, dimension(:) :: ett, etts
    real(prec) :: uave, vave, wave

    real(prec) :: top, bottom, skew, flat

    real(prec), allocatable, dimension(:, :, :) :: kx, ky, kz, ksqr

    integer :: wn
    real :: knorm

    ! Third
    real(prec) :: etotal, ediss_total, ediss_total0
    integer :: iwn
    real(prec) :: eta, vk, tk, uprm, tmse, Retms, xl, et2, kmxeta

    real(prec) :: L

    lx=nx
    ly=ny
    lz=nz

    cwidth=xwidth
    cwidthh=cwidth/2
    dx=xwidth/nx

    allocate(utemp(lx, ly, lz))

    ! Compute KE in a classical way
    te=0.5*sum(ux_global**2+uy_global**2+uz_global**2)/(nx*ny*nz)
    open(unit=100, file="stat_te.dat", status="unknown", form="formatted", position="append")
    write(100, "(I7,1X,F16.8,1X,1pe15.6)") step, step*dt, te
    close(100)

    ! Now let's compute statistics using FFT
    L=xwidth

    uave=sum(ux_global)/(nx*ny*nz)
    vave=sum(uy_global)/(nx*ny*nz)
    wave=sum(uz_global)/(nx*ny*nz)

    ux_global=ux_global-uave
    uy_global=uy_global-vave
    uz_global=uz_global-wave

    print *, "average vel.:", uave, vave, wave

    knorm=pi2/xwidth

    allocate(et(nx, ny, nz))
    allocate(ets(nx, ny, nz))
    allocate(ett(0:nx*2))
    allocate(etts(0:nx*2))

    allocate(kx(nx, ny, nz))
    allocate(ky(nx, ny, nz))
    allocate(kz(nx, ny, nz))
    allocate(ksqr(nx, ny, nz))

    ! wavenumbers
    do i=1, nx
       if (i<=nx/2) then
          kx(i, :, :)=2*pi/L*(i-1)
       else
          kx(i, :, :)=2*pi/L*(-nx+i-1)
       endif
    enddo
    do i=1, ny
       if (i<=ny/2) then
          ky(:, i, :)=2*pi/L*(i-1)
       else
          ky(:, i, :)=2*pi/L*(-ny+i-1)
       endif
    enddo
    do i=1, nz
       if (i<=nz/2) then
          kz(:, :, i)=2*pi/L*(i-1)
       else
          kz(:, :, i)=2*pi/L*(-nz+i-1)
       endif
    enddo

    ksqr=kx*kx+ky*ky+kz*kz

    do k=1, nz
       do j=1, ny
          do i=1, nx
             uin(i,j,k)=cmplx(ux_global(i,j,k), 0)
             vin(i,j,k)=cmplx(uy_global(i,j,k), 0)
             win(i,j,k)=cmplx(uz_global(i,j,k), 0)
          enddo
       enddo
    enddo

    call fftw_execute_dft(planu, uin, uout)
    call fftw_execute_dft(planv, vin, vout)
    call fftw_execute_dft(planw, win, wout)

    uout=uout/(nx*ny*nz)
    vout=vout/(nx*ny*nz)
    wout=wout/(nx*ny*nz)

    et=0.5*real(uout*conjg(uout)+vout*conjg(vout)+wout*conjg(wout))
    ets=2*visc*et*ksqr

    et=cshift(et, shift=-nx/2, dim=1)
    et=cshift(et, shift=-ny/2, dim=2)
    et=cshift(et, shift=-nz/2, dim=3)

    ets=cshift(ets, shift=-nx/2, dim=1)
    ets=cshift(ets, shift=-ny/2, dim=2)
    ets=cshift(ets, shift=-nz/2, dim=3)

    ett=0
    etts=0
    do k=1, nz
       do j=1, ny
          do i=1, nx
             wn=nint(sqrt((i-nx/2-1)**2+(j-ny/2-1)**2+(k-nz/2-1)**2+0.0))
             ett(wn)=ett(wn)+et(i,j,k)!/(knorm+1e-20)
             etts(wn)=etts(wn)+ets(i,j,k)!*knorm**2
          enddo
       enddo
    enddo

    ! energy from spectral
    etotal=sum(ett)
    ediss_total0=sum(etts)
    ediss_total=0

    eta   = (visc**3/ediss_total0)**0.25       ! Kolmogorov length scale
    vk    = (visc*ediss_total0)**0.25          ! Kolmogorov velocity scale
    tk    = sqrt(visc/ediss_total0)            ! Kolmogorov time scale
    uprm  = sqrt(2./3.*etotal)                 ! u_prime = u_rms
    tmse  = sqrt(15*visc*uprm**2/ediss_total0) ! Taylor microscale
    Retms = uprm*tmse/visc                     ! Taylor scale Reynolds number
    xl    = uprm**3/ediss_total0               ! Large eddy length scale
    kmxeta= (cwidthh-1.5*dx)*pi2/cwidth*eta    ! kmax*Kolmogorov length scale

    open(unit=100, file="stat_te_spectral.dat", status="unknown", form="formatted", position="append")
    write(100, "(I7,1X,F16.8,1X,3(1pe15.6,1X),8(1pe15.6,1X))") step, step*dt, etotal, ediss_total0, ediss_total, &
         eta, vk, tk, uprm, tmse, Retms, xl, kmxeta
    close(100)

    write(filename, "('spe',i9.9,'.dat')") step
    open(100, file=filename, status="replace")
    do i=0, nx
       write(100, "(I4,1X,1pe15.6)") i, ett(i)
    enddo
    close(100)

    ! compute first derivative
    uout=kx*(cmplx(0,1)*uout)
    vout=ky*(cmplx(0,1)*vout)
    wout=kz*(cmplx(0,1)*wout)

    call fftw_execute_dft(planudx, uout, uout2)
    call fftw_execute_dft(planvdy, vout, vout2)
    call fftw_execute_dft(planwdz, wout, wout2)

    ! compute skewness and flatness
    bottom=sum(real(uout2)**2+real(vout2)**2+real(wout2)**2)/(3*nx*ny*nz)
    top=sum(real(uout2)**3+real(vout2)**3+real(wout2)**3)/(3*nx*ny*nz)
    skew=top/bottom**(3./2.)
    top=sum(real(uout2)**4+real(vout2)**4+real(wout2)**4)/(3*nx*ny*nz)
    flat=top/bottom**2

    open(unit=100, file="skew-flat.dat", status="unknown", position="append")
    write(100, "(I7,1X,1pe15.6,1X,1pe15.6,1X,1pe15.6)") step, step*dt, skew, flat
    close(100)

    deallocate(kx, ky, kz, ksqr)
    deallocate(etts)
    deallocate(ets)
    deallocate(ett)
    deallocate(et)
  endsubroutine stat_te
endmodule stat
