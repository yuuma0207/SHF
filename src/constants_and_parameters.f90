module constants_and_parameters
    use,intrinsic :: iso_fortran_env, only: int32, real64
    implicit none
    
    !-- math and physical constants --
    real(real64),    parameter :: PI                 = 4.0d0*atan(1.0d0)
    real(real64),    parameter :: HBARC              = 197.32705d0          ! [MeV.fm]
    real(real64),    parameter :: mass_of_proton     = 938.272088d0         ! [MeV]
    real(real64),    parameter :: mass_of_neutron    = 939.565413d0         ! [MeV]
    real(real64),    parameter :: fine_structure     = 1d0/137.04d0         ! Fine Structure Constant
    real(real64),    parameter :: elementary_charge2 = HBARC*fine_structure ! [MeV.fm]
    real(real64),    parameter :: HBAR2_over_2m      = 20.75250d0           ! [MeV.fm^2]
    complex(real64), parameter :: imaginary_unit     = (0.0d0,1.0d0)
    
    
    !-- Parameters of Skyrme force --
    real(real64), parameter :: t_0   = -2645.0d0
    real(real64), parameter :: t_1   =  410.0d0
    real(real64), parameter :: t_2   = -135.0d0
    real(real64), parameter :: t_3   =  15595.0d0
    real(real64), parameter :: t_4   =  130.0d0
    real(real64), parameter :: x_0   =  0.090d0
    real(real64), parameter :: x_1   =  0.0d0
    real(real64), parameter :: x_2   =  0.0d0
    real(real64), parameter :: x_3   =  0.0d0
    real(real64), parameter :: param_gamma = 1.0d0/6.0d0
    
    real(real64),parameter :: b_0  = t_0*(1d0 + x_0/2d0)
    real(real64),parameter :: b_0p = t_0*(1d0/2d0 + x_0)
    real(real64),parameter :: b_1  = 1d0/4d0*(t_1*(1d0 + x_1/2d0) + t_2*(1d0 + x_2/2d0))
    real(real64),parameter :: b_1p = 1d0/4d0*(t_1*(1d0/2d0 + x_1) - t_2*(1d0/2d0 + x_2))
    real(real64),parameter :: b_2  = 1d0/8d0*(3d0*t_1*(1d0 + x_1/2d0) - t_2*(1d0 + x_2/2d0))
    real(real64),parameter :: b_2p = 1d0/8d0*(3d0*t_1*(1d0/2d0 + x_1) + t_2*(1d0/2d0 + x_2))
    real(real64),parameter :: b_3  = 1d0/4d0*t_3*(1d0 + x_3/2d0)
    real(real64),parameter :: b_3p = 1d0/4d0*t_3*(1d0/2d0 + x_3)
    real(real64),parameter :: b_4  = t_4/2d0
    real(real64),parameter :: b_4p = t_4/2d0
    real(real64),parameter :: eta = 1d0
    real(real64),parameter :: c_1  = eta*1d0/16d0*(t_1*x_1 + t_2*x_2)
    real(real64),parameter :: c_1p = eta*1d0/16d0*(t_1 - t_2)
    
    !-- Numerical parameters (can override by inputfile) --
    integer(int32)            :: Nr=100                           ! The number of grid points in r-direction (r^2=x^2+y^2)
    integer(int32)            :: Nz=100                           ! The number of grid points in z-direction
    integer(int32)            :: num_p=2                          ! The number of protons
    integer(int32)            :: num_n=2                          ! The number of neutrons
    real(real64)  ,parameter  :: r_max=20d0                       ! [fm] The maximum value of r (0<=r<=r_max)
    real(real64)  ,parameter  :: z_max=20d0                       ! [fm] The maximum value of z (-z_max<=z<=z_max)
    real(real64)  ,parameter  :: deformation_degree = 0.0d0       ! deformation degree of nucleus ! 0.1 がのーまる
    integer(int32)            :: prep_max_m = 3                   ! 方位量子数の最大値
    integer(int32)            :: prep_max_r = 3                   ! 主量子数の最大値
    integer(int32)            :: prep_max_z = 3                   ! z成分の最大値
    integer(int32)            :: prep_number = 2                  ! 用意する軌道の数

    real(real64)              :: dr                               ! [fm] The grid spacing in r-direction
    real(real64)              :: dz                               ! [fm] The grid spacing in z-direction
    real(real64)              :: z_center = z_max/2d0             ! [fm] The center of z-direction
    character(len=8)          :: str_z                            ! center of z-direction
    real(real64)              :: imaginary_time_step = 1d-2       ! imaginary time step E=50MeVとしてE*ct/hbarc <1になるように設定
    real(real64)              :: cou_eps = 1d0                    ! coulobm interactionを考慮するしきい値
    real(real64)              :: SD_eps = 1d-6
    real(real64)              :: prev_eps = 1d-12                 ! 計算の収束判定に用いる値
    real(real64)              :: lagrange_multiplier = 1d0        ! ラグランジュの未定乗数法におけるラグランジュ乗数    
                  
    character(len=8)          :: nuc_name
    integer(int32) :: max_iter  = 5000 ! 最大反復回数
    integer(int32) :: plot_iter = 500
    logical        :: is_output = .true.
    logical        :: is_use_dd = .true.
    logical        :: not_divJ  = .true.
    logical        :: is_use_B  = .true.
    logical        :: is_use_W  = .true.
    logical        :: is_TimeR  = .true.
    logical        :: is_converged = .false.
    logical        :: is_minus_var = .false.
    logical        :: is_coulomb   = .false.
    logical        :: is_update_coulomb    = .true.  ! 最初にクーロンを計算する
    logical        :: coulomb_updated_flag = .false. ! クーロンを更新したかどうかのフラグ
    logical        :: is_use_second_imag   = .false.

    !-- use output file(for python) --
    character(len=8)           :: date
    character(len=10)          :: time
    character(len=12)          :: MMDDHHMMSS
    character(len=9),parameter :: resultdir = "./../res/"   ! outputdir
    character(len=21)          :: foldername
end module constants_and_parameters