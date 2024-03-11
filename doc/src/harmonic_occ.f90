module harmonic_occ_wf
    use,intrinsic :: iso_fortran_env
    use           :: constants_and_parameters
    implicit none
    
    real(real64) :: hboma
    real(real64) :: xpi=sqrt(4.0d0/(5.0d0*pi))
    
    contains
    function set_harmonic_occ(node_r,node_z,mag,r,z,betain) result(normilized_wf)
        use global_variables, only : PI,HBAR2_over_2m
        use spacial_function
        integer(int32), intent(in) :: node_z,node_r,mag
        real(real64), intent(in) :: r,z
        real(real64), intent(in) :: betain
        real(real64) :: alpha_z,alpha_r
        real(real64) :: normilized_wf
        real(real64) :: z_func
        real(real64) :: rad_func
        real(real64) :: rad_norm
        real(real64) :: z_norm
        real(real64) :: hbomz,hbomr
        hboma  = 41.0d0/dble(num_n+num_p)**(1d0/3d0) 
        hbomz = hboma*exp(-xpi*betain)
        hbomr=hboma*exp(xpi/2.0d0*betain)
        alpha_z = sqrt(2.0d0*HBAR2_over_2m/hbomz) 
        alpha_r = sqrt(2.0d0*HBAR2_over_2m/hbomr) 
        
        z_func = Hermite_n(node_z,z/alpha_z)*exp(-z**2/(2d0*alpha_z**2))
        rad_func = Laguerre_nl(node_r,abs(mag),(r/alpha_r)**2)*exp(-r**2/(2d0*alpha_r**2))*(r/alpha_r)**abs(mag)
        z_norm   = sqrt(1d0/alpha_z)*sqrt(1d0/(2d0**node_z*gamma(dble(node_z+1))))*sqrt(1d0/sqrt(PI))
        rad_norm = sqrt(1d0/alpha_r**2)*sqrt(gamma(dble(node_r+1))/gamma(dble(node_r+abs(mag)+1)))*sqrt(1d0/PI)
        
        normilized_wf = z_norm*rad_norm*z_func*rad_func
        
        if(normilized_wf>1)then
            print*, node_r,node_z,abs(mag),r,z,betain
            print*, z_func,rad_func
            print*, z_norm,rad_norm
            print*, "Error in harmonic_occ_wf: normilized_wf>1"
            stop
        end if
    end function set_harmonic_occ

    function calc_occ_energy(node_r,node_z,mag,betain) result(occ_energy)
        integer(int32), intent(in) :: node_z,node_r,mag
        real(real64), intent(in) :: betain
        real(real64) :: occ_energy
        real(real64) :: hbomz,hbomr

        hbomz = hboma*exp(-xpi*betain)
        hbomr = hboma*exp(xpi/2.0d0*betain)
        occ_energy = hbomz*(0.5d0 + node_z) + hbomr*dble((2*node_r + abs(mag) + 1))
    end function calc_occ_energy
    
end module harmonic_occ_wf