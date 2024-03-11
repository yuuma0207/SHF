module determine_meanfield_potential
    use :: constants_and_parameters,  only: is_coulomb
    implicit none
    contains
    subroutine calc_pot_U
        implicit none
        !call calc_Woods_Saxon

        ! pについてはクーロン項を付け足す必要がある
        call calc_pot_U_np("n",.false.)
        call calc_pot_U_np("p",is_coulomb)
        
    end subroutine calc_pot_U

    subroutine calc_Woods_Saxon
        use,intrinsic :: iso_fortran_env, only: int32, real64
        use :: constants_and_parameters,  only: Nr,Nz,z_center
        use :: global_variables,          only: U_n,U_p,r_vec,z_vec
        implicit none

        real(real64) :: r0
        real(real64) :: a
        real(real64) :: VC
        real(real64) :: radius
        integer(int32) :: i,j

        r0 = 1.27d0 ! fm
        a  = 0.67d0 ! fm
        VC = -51d0  !+ 33d0*(2d0-2d0)/4d0  ! MeV
        radius = r0*(4d0)**(1d0/3d0) ! fm

        do j=1,Nz
            do i=1,Nr
                U_n(i,j) = VC/(1d0+exp((sqrt(r_vec(i)**2 + (z_vec(j)-z_center)**2)-radius)/a))
                U_p(i,j) = VC/(1d0+exp((sqrt(r_vec(i)**2 + (z_vec(j)-z_center)**2)-radius)/a))
            end do
        end do

    end subroutine calc_Woods_Saxon

    subroutine calc_pot_U_np(np_str,is_coulomb)
        use,intrinsic :: iso_fortran_env, only: int32, real64
        use :: write_to_file,             only: write_U_term
        use :: constants_and_parameters,  only: e2 => elementary_charge2,PI &
                                               ,param_gamma,b_0,b_0p,b_1,b_1p,b_2,b_2p,b_3,b_3p,b_4,b_4p,Nr,Nz &
                                               ,is_update_coulomb, coulomb_updated_flag
        use :: global_variables,          only: tau_n => kinetic_dens_n &
                                              , tau_p => kinetic_dens_p &
                                              , j_n => spin_orbit_dens_n &
                                              , j_p => spin_orbit_dens_p &
                                              , div_j_n => div_spin_orbit_dens_n &
                                              , div_j_p => div_spin_orbit_dens_p &
                                              , rho_n,rho_p,lap_rho_n,lap_rho_p &
                                              , U_n,U_p,direct_Coulomb_pot,exchange_Coulomb_pot
        implicit none
        character(1),intent(in)   :: np_str
        logical,intent(in),optional :: is_coulomb
        real(real64), allocatable :: U_np(:,:)                  ! [MeV] (r,z)
        real(real64), allocatable :: rho_np(:,:)                ! [fm^-3] (r,z)
        real(real64), allocatable :: tau_np(:,:)                ! [fm^-3] (r,z)
        real(real64), allocatable :: j_np(:,:,:)                ! [fm^-3] (r,z,component)
        real(real64), allocatable :: laplacian_rho_np(:,:)      ! [fm^-5] (r,z)
        real(real64), allocatable :: div_j_np(:,:)              ! [fm^-4] (r,z,component)
        real(real64), allocatable :: rho(:,:)                   ! [fm^-3] (r,z)
        real(real64), allocatable :: tau(:,:)                   ! [fm^-3] (r,z)
        real(real64), allocatable :: laplacian_rho(:,:)         ! [fm^-5] (r,z)
        real(real64), allocatable :: div_j(:,:)                 ! [fm^-4] (r,z,component)

        if (np_str /= "n" .and. np_str /= "p") then
            write(*,*) "Error: np_str must be 'n' or 'p'."
            stop
        end if
        
        allocate(U_np(Nr,Nz))
        allocate(rho_np(Nr,Nz))
        allocate(tau_np(Nr,Nz))
        allocate(j_np(Nr,Nz,3))
        allocate(laplacian_rho_np(Nr,Nz))
        allocate(div_j_np(Nr,Nz))
        allocate(rho(Nr,Nz))
        allocate(tau(Nr,Nz))
        allocate(laplacian_rho(Nr,Nz))
        allocate(div_j(Nr,Nz))

        if(np_str == "n")then
            rho_np(:,:) = rho_n(:,:)
            tau_np(:,:) = tau_n(:,:)
            j_np(:,:,:) = j_n(:,:,:)
            laplacian_rho_np(:,:) = lap_rho_n(:,:)
            div_j_np(:,:) = div_j_n(:,:)
        else if(np_str == "p")then
            rho_np(:,:) = rho_p(:,:)
            tau_np(:,:) = tau_p(:,:)
            j_np(:,:,:) = j_p(:,:,:)
            laplacian_rho_np(:,:) = lap_rho_p(:,:)
            div_j_np(:,:) = div_j_p(:,:)
        end if

        laplacian_rho(:,:)   = lap_rho_n(:,:) + lap_rho_p(:,:)
        rho(:,:)             = rho_n(:,:) + rho_p(:,:)
        tau(:,:)             = tau_n(:,:) + tau_p(:,:)
        div_j(:,:)           = div_j_n(:,:) + div_j_p(:,:)
        

        U_np(:,:) = b_0*rho(:,:) - b_0p*rho_np(:,:) &
        + b_1*tau(:,:)           - b_1p*tau_np(:,:) & 
        - b_2*laplacian_rho(:,:) + b_2p*laplacian_rho_np(:,:) &
        +b_3*(param_gamma+2d0)/3d0*rho(:,:)**(param_gamma+1d0) &
        -b_3p*2d0/3d0*rho(:,:)**(param_gamma)*rho_np(:,:) &
        -b_3p*param_gamma/3d0*rho(:,:)**(param_gamma-1d0)*(rho_n(:,:)**2 + rho_p(:,:)**2) &
        -b_4*div_j(:,:)          - b_4p*div_j_np(:,:) ! 9bのdelta_{qp}の前までの式を計算している。

        
        if(np_str == "n")then
            U_n(:,:) = U_np(:,:)
        else if(np_str == "p")then ! クーロンポテンシャルを追加する必要あり
            if(is_coulomb)then
                if(is_update_coulomb)then
                    call calc_coulomb_pot ! update direct_Coulomb_pot
                    is_update_coulomb = .false.
                    coulomb_updated_flag = .true.
                end if
                exchange_Coulomb_pot(:,:) = e2*(3d0/PI)**(1d0/3d0)*rho_p(:,:)**(1d0/3d0)
                ! direct_Coulomb_pot   = e**2 * \int dV' rho(r')/|r-r'|
                ! exchange_Coulomb_pot = e**2 * (3/PI)**(1/3) * rho_p(r')**(1/3)
            end if
            U_p(:,:) = U_np(:,:) + direct_Coulomb_pot(:,:) - exchange_Coulomb_pot(:,:)
        end if

        block 
            real(real64),allocatable :: b0t(:,:),b0pt(:,:),b1t(:,:),b1pt(:,:),b2t(:,:),b2pt(:,:)&
                                        ,b3t(:,:),b3pt(:,:),b3pt2(:,:),b4t(:,:),b4pt(:,:)
            allocate(b0t(Nr,Nz))
            allocate(b0pt(Nr,Nz))
            allocate(b1t(Nr,Nz))
            allocate(b1pt(Nr,Nz))
            allocate(b2t(Nr,Nz))
            allocate(b2pt(Nr,Nz))
            allocate(b3t(Nr,Nz))
            allocate(b3pt(Nr,Nz))
            allocate(b3pt2(Nr,Nz))
            allocate(b4t(Nr,Nz))
            allocate(b4pt(Nr,Nz))
            b0t(:,:)   = b_0*rho(:,:)
            b0pt(:,:)  = b_0p*rho_np(:,:)
            b1t(:,:)   = b_1*tau(:,:)
            b1pt(:,:)  = b_1p*tau_np(:,:)
            b2t(:,:)   = b_2*laplacian_rho(:,:)
            b2pt(:,:)  = b_2p*laplacian_rho_np(:,:)
            b3t(:,:)   = b_3*(param_gamma+2d0)/3d0*rho(:,:)**(param_gamma+1d0)
            b3pt(:,:)  = b_3p*2d0/3d0*rho(:,:)**(param_gamma)*rho_np(:,:)
            b3pt2(:,:) = b_3p*param_gamma/3d0*rho(:,:)**(param_gamma-1d0)*(rho_n(:,:)**2 + rho_p(:,:)**2)
            b4t(:,:)   = b_4*div_j(:,:)
            b4pt(:,:)  = b_4p*div_j_np(:,:)

            call write_U_term(b0t(:,:),b0pt(:,:),b1t(:,:),b1pt(:,:),b2t(:,:)&
                             ,b2pt(:,:),b3t(:,:),b3pt(:,:),b3pt2(:,:),b4t(:,:),b4pt(:,:)&
                             ,direct_Coulomb_pot(:,:),exchange_Coulomb_pot(:,:),np_str)

            deallocate(b0t,b0pt,b1t,b1pt,b2t,b2pt,b3t,b3pt,b3pt2,b4t,b4pt)
        end block
        deallocate(U_np)
        deallocate(rho_np)
        deallocate(tau_np)
        deallocate(j_np)
        deallocate(laplacian_rho_np)
        deallocate(div_j_np)
        deallocate(rho)
        deallocate(tau)
        deallocate(laplacian_rho)
        deallocate(div_j)
    end subroutine calc_pot_U_np

    subroutine calc_coulomb_pot  ! 後でモーメントの足し算に変更するかも．今はlapackで解いている
        use,intrinsic :: iso_fortran_env, only: int32
        use global_variables,only: rho_p,direct_Coulomb_pot
        use poisson_problem, only: nine_point_poisson_by_lapack
        implicit none
        ! 時間計測を行う
        !integer(int32) :: ti,tf,tr
        !call system_clock(ti)
        call nine_point_poisson_by_lapack(rho_p(:,:),direct_Coulomb_pot(:,:))
        !call system_clock(tf,tr)
        !write(6,'(f10.3,A)')(tf-ti)/dble(tr),'[s]'
    end subroutine calc_coulomb_pot
    
    subroutine calc_pot_B
        use :: constants_and_parameters,only: HBAR2_over_2m,b_1,b_1p
        use :: global_variables,only:rho_n,rho_p,B_n,B_p
        implicit none
        B_n(:,:) = b_1*(rho_n(:,:)+rho_p(:,:)) - b_1p*rho_n(:,:)
        B_p(:,:) = b_1*(rho_n(:,:)+rho_p(:,:)) - b_1p*rho_p(:,:)
    end subroutine calc_pot_B

    subroutine calc_pot_W
        use,intrinsic :: iso_fortran_env, only: int32, real64
        use :: constants_and_parameters,  only: Nr,Nz,b_4,b_4p,c_1,c_1p
        use :: global_variables,          only: j_n => spin_orbit_dens_n, j_p => spin_orbit_dens_p &
                                              , rho_n,rho_p,dr_rho_n,dr_rho_p,dz_rho_n,dz_rho_p &
                                              , W_n,W_p
        use :: math_derivation ,          only: math_diff_r_9point,math_diff_z_9point
        implicit none
        real(real64), allocatable :: rho(:,:) ! [fm^-3] (r,z)
        real(real64), allocatable :: spin_orbit_dens(:,:,:)   ! [fm^-3] (r,z,component)
        real(real64), allocatable :: dr_rho(:,:)
        real(real64), allocatable :: dz_rho(:,:)
        integer(int32) :: i,j

        allocate(rho(Nr,Nz))
        allocate(spin_orbit_dens(Nr,Nz,3))
        allocate(dr_rho,mold=rho_n)
        allocate(dz_rho,mold=rho_n)
        
        rho(:,:) = rho_n(:,:) + rho_p(:,:)
        spin_orbit_dens(:,:,:) = j_n(:,:,:) + j_p(:,:,:)
        dr_rho(:,:) = dr_rho_n(:,:) + dr_rho_p(:,:)
        dz_rho(:,:) = dz_rho_n(:,:) + dz_rho_p(:,:)
        
        do j=1,Nz
            do i=1,Nr
                W_n(i,j,1) =  b_4*dr_rho(i,j) + b_4p*dr_rho_n(i,j) &
                             -2d0*c_1*spin_orbit_dens(i,j,1) + 2d0*c_1p*j_n(i,j,1)
                            
                W_n(i,j,2) = 0d0
                
                W_n(i,j,3) =  b_4*dz_rho(i,j) + b_4p*dz_rho_n(i,j) &
                             -2d0*c_1*spin_orbit_dens(i,j,3) + 2d0*c_1p*j_n(i,j,3)
                            
                W_p(i,j,1) =  b_4*dr_rho(i,j) + b_4p*dr_rho_p(i,j) &
                             -2d0*c_1*spin_orbit_dens(i,j,1) + 2d0*c_1p*j_p(i,j,1)
                             
                W_p(i,j,2) = 0d0
                
                W_p(i,j,3) =  b_4*dz_rho(i,j) + b_4p*dz_rho_p(i,j) &
                             -2d0*c_1*spin_orbit_dens(i,j,3) + 2d0*c_1p*j_p(i,j,3)
            end do
        end do
    end subroutine calc_pot_W
end module determine_meanfield_potential