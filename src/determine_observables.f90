module determine_observables
    implicit none

    contains
    
    
    subroutine calc_Skyrme_Energy_by_integrate
        use,intrinsic :: iso_fortran_env,only:int32,real64
        use constants_and_parameters,only: b_0,b_0p,b_3,b_3p,b_2,b_2p,b_1,b_1p,b_4,b_4p,c_1,c_1p,param_gamma &
                                          ,Nr,Nz
        use global_variables,only:Skyrme_Energy,b_0_term,b_3_term,b_2_term,b_1_term,b_4_term,c_1_term
        use math_integrate,only: volume_integrate
        use global_variables,only: rho_n, rho_p, lap_rho_n, lap_rho_p &
                                  ,tau_n => kinetic_dens_n, tau_p => kinetic_dens_p &
                                  ,j_n => spin_orbit_dens_n, j_p => spin_orbit_dens_p &
                                  ,div_j_n => div_spin_orbit_dens_n, div_j_p => div_spin_orbit_dens_p
        implicit none
        real(real64),allocatable :: rho(:,:),lap_rho(:,:), tau(:,:), j(:,:,:), div_j(:,:)
                                  
        allocate(rho(Nr,Nz))
        allocate(lap_rho(Nr,Nz))
        allocate(tau(Nr,Nz))
        allocate(j(Nr,Nz,3))
        allocate(div_j(Nr,Nz))
        
        rho(:,:) = rho_n(:,:) + rho_p(:,:)
        lap_rho(:,:) = lap_rho_n(:,:) + lap_rho_p(:,:)
        tau(:,:) = tau_n(:,:) + tau_p(:,:)
        j(:,:,:) = j_n(:,:,:) + j_p(:,:,:)
        div_j(:,:) = div_j_n(:,:) + div_j_p(:,:)
        
        b_0_term = volume_integrate(b_0/2.0d0*(rho(:,:))**2 -b_0p/2.0d0*(rho_n(:,:)**2 + rho_p(:,:)**2))
        
        b_3_term = volume_integrate(b_3/3.0d0*(rho(:,:))**(param_gamma+2.0d0) &
                                    - b_3p/3.0d0*rho(:,:)**(param_gamma)*(rho_n(:,:)**2 + rho_p(:,:)**2))
                                    
        b_2_term = volume_integrate(-b_2/2.0d0*rho(:,:)*lap_rho(:,:) &
                                    + b_2p/2.0d0*(rho_n(:,:)*lap_rho_n(:,:) + rho_p(:,:)*lap_rho_p(:,:)))
                                    
        b_1_term = volume_integrate(b_1*rho(:,:)*tau(:,:) - b_1p*(rho_n(:,:)*tau_n(:,:) + rho_p(:,:)*tau_p(:,:)))
        
        b_4_term = volume_integrate(b_4*rho(:,:)*div_j(:,:) + b_4p*(rho_n(:,:)*div_j_n(:,:) + rho_p(:,:)*div_j_p(:,:)))
        
        c_1_term = volume_integrate(c_1*(j(:,:,1)**2 + j(:,:,3)**2) &
                                    - c_1p*(j_n(:,:,1)**2 + j_n(:,:,3)**2) &
                                    - c_1p*(j_p(:,:,1)**2 + j_p(:,:,3)**2))
        
        Skyrme_Energy = b_0_term + b_3_term + b_2_term + b_1_term - b_4_term - c_1_term
        
        deallocate(rho)
        deallocate(lap_rho)
        deallocate(tau)
        deallocate(j)
        deallocate(div_j)
    end subroutine calc_Skyrme_Energy_by_integrate
    
    
    subroutine calc_Kinetic_Energy_by_integrate
        use constants_and_parameters,only:HBAR2_over_2m
        use global_variables,only:Kinetic_Energy
        use math_integrate,only: volume_integrate
        use global_variables,only: tau_n => kinetic_dens_n, tau_p => kinetic_dens_p
        implicit none
        
        Kinetic_Energy = HBAR2_over_2m*volume_integrate(tau_n(:,:) + tau_p(:,:))
        
    end subroutine calc_Kinetic_Energy_by_integrate

    subroutine calc_Coulomb_Energy
        use global_variables,only:direct_Coulomb_pot,exchange_Coulomb_pot,Coulomb_Energy,rho_p
        use math_integrate,only: volume_integrate
        implicit none

        Coulomb_Energy = volume_integrate(rho_p(:,:)*(direct_Coulomb_pot(:,:)/2d0 - 3d0/4d0*exchange_Coulomb_pot(:,:)))
    end subroutine calc_Coulomb_Energy

    subroutine calc_CM_Energy
        use,intrinsic :: iso_fortran_env, only : int32, real64
        use constants_and_parameters, only : mass_of_neutron, mass_of_proton,num_n,num_p
        use global_variables, only : CM_Energy
        implicit none
        integer(int32) :: k,kp
        real(real64) :: MC2_tot

        CM_Energy = 0d0

        do k=1,num_n
            CM_Energy = CM_Energy + calc_PC2_ExV(k)
        end do
        do k=1,num_p
            kp = k + num_n
            CM_Energy = CM_Energy + calc_PC2_ExV(kp)
        end do
        MC2_tot = num_n*mass_of_neutron + num_p*mass_of_proton
        CM_Energy = CM_Energy/(2.0d0*MC2_tot)
    end subroutine calc_CM_Energy

    function calc_PC2_ExV(nuc_idx) result(expected_value) ! calculate (PC)^2 expectation value
        use,intrinsic :: iso_fortran_env, only : int32, real64
        use constants_and_parameters, only : HBARC
        use global_variables, only : wf_plus, wf_minus, lap_wf_plus, lap_wf_minus
        use math_integrate, only : volume_integrate
        implicit none
        integer(int32),intent(in) :: nuc_idx
        real(real64) :: expected_value

        expected_value = -(HBARC)**2*volume_integrate(wf_plus(:,:,nuc_idx)*lap_wf_plus(:,:,nuc_idx) &
                                                     + wf_minus(:,:,nuc_idx)*lap_wf_minus(:,:,nuc_idx))

    end function calc_PC2_ExV

    
    subroutine calc_Total_Energy
        use global_variables,only:Total_Energy, Skyrme_Energy, Kinetic_Energy, Coulomb_Energy, CM_Energy
        
        Total_Energy = Skyrme_Energy + Kinetic_Energy + Coulomb_Energy - CM_Energy
        
        
    end subroutine calc_Total_Energy
    
end module determine_observables