module dft
  use, intrinsic :: iso_c_binding 
  include 'fftw3.f03'

  type(C_PTR) :: planu, planv, planw
  type(C_PTR) :: planudx, planvdy, planwdz

  ! complex(C_FLOAT_COMPLEX), allocatable, dimension(:, :, :) :: uin, uout
  ! complex(C_FLOAT_COMPLEX), allocatable, dimension(:, :, :) :: vin, vout
  ! complex(C_FLOAT_COMPLEX), allocatable, dimension(:, :, :) :: win, wout

  ! complex(C_FLOAT_COMPLEX), allocatable, dimension(:, :, :) :: uout2, vout2, wout2

  complex(C_DOUBLE_COMPLEX), allocatable, dimension(:, :, :) :: uin, uout
  complex(C_DOUBLE_COMPLEX), allocatable, dimension(:, :, :) :: vin, vout
  complex(C_DOUBLE_COMPLEX), allocatable, dimension(:, :, :) :: win, wout

  complex(C_DOUBLE_COMPLEX), allocatable, dimension(:, :, :) :: uout2, vout2, wout2

contains
  
  subroutine init_dft(nx, ny, nz)
    implicit none
    integer :: nx, ny, nz

    allocate(uin(nx, ny, nz), uout(nx, ny, nz))
    allocate(vin(nx, ny, nz), vout(nx, ny, nz))
    allocate(win(nx, ny, nz), wout(nx, ny, nz))

    allocate(uout2(nx, ny, nz))
    allocate(vout2(nx, ny, nz))
    allocate(wout2(nx, ny, nz))

    planu = fftw_plan_dft_3d(nz, ny, nx, uin, uout, FFTW_FORWARD,FFTW_MEASURE)
    planv = fftw_plan_dft_3d(nz, ny, nx, vin, vout, FFTW_FORWARD,FFTW_MEASURE)
    planw = fftw_plan_dft_3d(nz, ny, nx, win, wout, FFTW_FORWARD,FFTW_MEASURE)

    planudx = fftw_plan_dft_3d(nz, ny, nx, uout, uout2, FFTW_BACKWARD,FFTW_MEASURE)
    planvdy = fftw_plan_dft_3d(nz, ny, nx, vout, vout2, FFTW_BACKWARD,FFTW_MEASURE)
    planwdz = fftw_plan_dft_3d(nz, ny, nx, wout, wout2, FFTW_BACKWARD,FFTW_MEASURE)
  endsubroutine init_dft

  subroutine free_dft
    implicit none
    
    call fftw_destroy_plan(planu)
    call fftw_destroy_plan(planv)
    call fftw_destroy_plan(planw)
    call fftw_destroy_plan(planudx)
    call fftw_destroy_plan(planvdy)
    call fftw_destroy_plan(planwdz)

    deallocate(uin, uout)
    deallocate(vin, vout)
    deallocate(win, wout)
    deallocate(uout2, vout2, wout2)
  endsubroutine free_dft
endmodule dft
