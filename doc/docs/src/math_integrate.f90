module math_integrate
    use,intrinsic :: iso_fortran_env
    implicit none

    
    contains
    function volume_integrate(f) result (result)
        use constants_and_parameters, only : Nr, Nz, dr, dz, PI
        use global_variables, only: r_vec
        real(real64), intent(in) :: f(:,:)
        real(real64)             :: result
        integer(int32)           :: r_size, z_size
        !real(real64)             :: Wr(Nr), Wz(Nz)
        
        r_size = size(f,1)
        z_size = size(f,2)
        if(r_size /= Nr) write(*,*) 'Error: r_size /= Nr'
        if(z_size /= Nz) write(*,*) 'Error: z_size /= Nz'

        
        !Wr(:) = r_vec(:)*dr;  !論文ではr=0でdr/2とかいているが、間違っている。円筒座標系ではr=0の体積要素が0なのでこのままでよい
        !Wz(:) = dz    ! ; Wz(1) = dz/2d0; Wz(Nz) = dz/2d0

        result = 0d0
        
        block
            integer(int32) :: i,j
            do j=1,Nz
                do i=1,Nr
                    result = result + f(i,j)*r_vec(i)*dr*dz
                end do 
            end do
        end block
        result = result*2d0*PI
        return
    end function volume_integrate

    function sympson_volume_integrate(f) result (result)
        use constants_and_parameters, only: Nr,Nz,dr,dz,PI
        use global_variables, only: r_vec
        real(real64), intent(in) :: f(:,:)
        real(real64)             :: result
        integer(int32)           :: r_size, z_size
        real(real64)             :: Wr(Nr), Wz(Nz)

        r_size = size(f,1)
        z_size = size(f,2)
        if(r_size /= Nr) write(*,*) 'Error: r_size /= Nr'
        if(z_size /= Nz) write(*,*) 'Error: z_size /= Nz'




    end function sympson_volume_integrate
    
end module math_integrate