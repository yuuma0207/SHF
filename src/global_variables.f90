module global_variables
    use,intrinsic :: iso_fortran_env
    use constants_and_parameters
    implicit none

!-- Density and potential (caluculated by Occupied Orbitals) --
    real(real64)   , allocatable :: rho_n(:,:)                 ! [fm^-3]  (r,z)
    real(real64)   , allocatable :: rho_p(:,:)                 ! [fm^-3]  (r,z)
    real(real64)   , allocatable :: dr_rho_n(:,:)
    real(real64)   , allocatable :: dz_rho_n(:,:)
    real(real64)   , allocatable :: dr_rho_p(:,:)
    real(real64)   , allocatable :: dz_rho_p(:,:)
    real(real64)   , allocatable :: ddr_rho_n(:,:)
    real(real64)   , allocatable :: ddz_rho_n(:,:)
    real(real64)   , allocatable :: ddr_rho_p(:,:)
    real(real64)   , allocatable :: ddz_rho_p(:,:)
    real(real64)   , allocatable :: lap_rho_n(:,:)             ! [fm^-5]  (r,z)
    real(real64)   , allocatable :: lap_rho_p(:,:)             ! [fm^-5]  (r,z)
    real(real64)   , allocatable :: kinetic_dens_n(:,:)        ! [fm^-3]  (r,z)
    real(real64)   , allocatable :: kinetic_dens_p(:,:)        ! [fm^-3]  (r,z)
    real(real64)   , allocatable :: thomas_fermi_n(:,:)        ! [fm^-4]  (r,z)
    real(real64)   , allocatable :: thomas_fermi_p(:,:)        ! [fm^-4]  (r,z)
    real(real64)   , allocatable :: spin_orbit_dens_n(:,:,:)   ! [fm^-3]  (r,z,vec component)
    real(real64)   , allocatable :: spin_orbit_dens_p(:,:,:)   ! [fm^-3]  (r,z,vec component)
    real(real64)   , allocatable :: div_spin_orbit_dens_n(:,:) ! [fm^-4]  (r,z)
    real(real64)   , allocatable :: div_spin_orbit_dens_p(:,:) ! [fm^-4]  (r,z)
    real(real64)   , allocatable :: U_n(:,:)                   ! [MeV] (r,z)
    real(real64)   , allocatable :: U_p(:,:)                   ! [MeV] (r,z)
    real(real64)   , allocatable :: B_n(:,:)                   ! [MeV] (r,z)
    real(real64)   , allocatable :: B_p(:,:)                   ! [MeV] (r,z)
    real(real64)   , allocatable :: W_n(:,:,:)                 ! [MeV] (r,z,component)
    real(real64)   , allocatable :: W_p(:,:,:)                 ! [MeV] (r,z,component)
    real(real64)   , allocatable :: direct_Coulomb_pot(:,:)    ! [MeV] (r,z)
    real(real64)   , allocatable :: exchange_Coulomb_pot(:,:)  ! [MeV] (r,z)

! -- All Orbitals(Occupied state / All state : prefix = prep(last component is n or p)) --
    integer(int32) , allocatable :: prep_idx_array(:,:)        ! nucleon index (prepared)
    integer(int32) , allocatable :: magnetic_q_num(:)          ! magnetic quantum number
    integer(int32) , allocatable :: prep_magnetic_q_num(:,:)   ! magnetic quantum number (prepared)
    character(4)   , allocatable :: init_rzm_to_write(:)       ! write quantum number (only to write)
    real(real64)   , allocatable :: prep_one_particle_E(:,:)   ! [MeV] 1粒子エネルギー (prepared)
    real(real64)   , allocatable :: prep_one_particle_E2(:,:)  ! [MeV2] <α|h^2|α> (prepared) prep, n/p, u/d
    real(real64)   , allocatable :: prep_one_E_Kin(:,:,:)      ! [MeV] 1粒子エネルギー (prepared) prep, n/p, u/d
    real(real64)   , allocatable :: prep_one_E_U(:,:,:)        ! [MeV] 1粒子エネルギー (prepared)
    real(real64)   , allocatable :: prep_one_E_B(:,:,:)        ! [MeV] 1粒子エネルギー (prepared)
    real(real64)   , allocatable :: prep_one_E_W(:,:,:)        ! [MeV] 1粒子エネルギー (prepared)
    real(real64)   , allocatable :: wf_plus(:,:,:)             ! [fm^-3/2]  (r,z,nucleon)
    real(real64)   , allocatable :: wf_minus(:,:,:)            ! [fm^-3/2]  (r,z,nucleon)
    real(real64)   , allocatable :: prep_wf_plus(:,:,:,:)      ! [fm^-3/2]  (r,z,prepsize, n/p)
    real(real64)   , allocatable :: prep_wf_minus(:,:,:,:)     ! [fm^-3/2]  (r,z,prepsize, n/p)
    !----- derivative of wave function -----
    real(real64)   , allocatable :: dr_wf_plus(:,:,:)
    real(real64)   , allocatable :: dz_wf_plus(:,:,:)
    real(real64)   , allocatable :: ddr_wf_plus(:,:,:)
    real(real64)   , allocatable :: ddz_wf_plus(:,:,:)
    real(real64)   , allocatable :: lap_wf_plus(:,:,:)
    real(real64)   , allocatable :: prep_dr_wf_plus(:,:,:,:)
    real(real64)   , allocatable :: prep_dz_wf_plus(:,:,:,:)
    real(real64)   , allocatable :: prep_ddr_wf_plus(:,:,:,:)
    real(real64)   , allocatable :: prep_ddz_wf_plus(:,:,:,:)
    real(real64)   , allocatable :: prep_lap_wf_plus(:,:,:,:)
    real(real64)   , allocatable :: dr_wf_minus(:,:,:)
    real(real64)   , allocatable :: dz_wf_minus(:,:,:)
    real(real64)   , allocatable :: ddr_wf_minus(:,:,:)
    real(real64)   , allocatable :: ddz_wf_minus(:,:,:)
    real(real64)   , allocatable :: lap_wf_minus(:,:,:)
    real(real64)   , allocatable :: prep_dr_wf_minus(:,:,:,:)
    real(real64)   , allocatable :: prep_dz_wf_minus(:,:,:,:)
    real(real64)   , allocatable :: prep_ddr_wf_minus(:,:,:,:)
    real(real64)   , allocatable :: prep_ddz_wf_minus(:,:,:,:)
    real(real64)   , allocatable :: prep_lap_wf_minus(:,:,:,:)
    real(real64)   , allocatable :: prep_Upsi_plus(:,:,:,:)
    real(real64)   , allocatable :: prep_Upsi_minus(:,:,:,:)
    real(real64)   , allocatable :: prep_Bpsi_plus(:,:,:,:)
    real(real64)   , allocatable :: prep_Bpsi_minus(:,:,:,:)
    real(real64)   , allocatable :: prep_Wpsi_plus(:,:,:,:)
    real(real64)   , allocatable :: prep_Wpsi_minus(:,:,:,:)
    !----- use for calculation -----
    ! SDNWP/M = sigma dot nabla plus/minus
    real(real64)   , allocatable ::  SDNWP(:,:,:)               ! [fm^-3/2]  (r,z,nucleon)
    real(real64)   , allocatable ::  SDNWM(:,:,:)               ! [fm^-3/2]  (r,z,nucleon)
    real(real64)   , allocatable ::  prep_SDNWP(:,:,:,:)        ! [fm^-3/2]  (r,z,nucleon)
    real(real64)   , allocatable ::  prep_SDNWM(:,:,:,:)        ! [fm^-3/2]  (r,z,nucleon)
    ! Hpsi
    real(real64)   , allocatable ::  prep_Hpsi2_plus(:,:,:,:)   ! [MeV fm^-3/2]  (r,z,nucleon)
    real(real64)   , allocatable ::  prep_Hpsi2_minus(:,:,:,:)  ! [MeV fm^-3/2]  (r,z,nucleon)
    real(real64)   , allocatable ::  prep_Hpsi_plus(:,:,:,:)    ! [MeV fm^-3/2]  (r,z,nucleon)
    real(real64)   , allocatable ::  prep_Hpsi_minus(:,:,:,:)   ! [MeV fm^-3/2]  (r,z,nucleon)

    !-- Grid variables --
    real(real64)   , allocatable :: r_vec(:)                    ! [fm] r position
    real(real64)   , allocatable :: z_vec(:)                    ! [fm] z position

    !-- Observables --
    real(real64) :: Skyrme_Energy
    real(real64) :: Kinetic_Energy
    real(real64) :: Coulomb_Energy
    real(real64) :: CM_Energy
    real(real64) :: Total_Energy
    real(real64) :: Avarage_variance
    real(real64) :: b_0_term,b_3_term,b_2_term,b_1_term,b_4_term,c_1_term
    
    
end module global_variables