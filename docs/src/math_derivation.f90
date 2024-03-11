module math_derivation
    use,intrinsic :: iso_fortran_env,only: int32,real64
    implicit none
    real(real64) :: eps = 0.03d0

    interface math_diff_r_9point
        module procedure math_diff_r_9point_realfunc, math_diff_r_9point_complexfunc
    end interface math_diff_r_9point

    interface math_diff_z_9point
        module procedure math_diff_z_9point_realfunc, math_diff_z_9point_complexfunc
    end interface math_diff_z_9point
    
    interface math_diff_2r_9point
        module procedure math_diff_2r_9point_realfunc
    end interface math_diff_2r_9point
    
    interface math_diff_2z_9point
        module procedure math_diff_2z_9point_realfunc
    end interface math_diff_2z_9point


contains
    ! r方向の微分はr=0で偶関数であることを利用している
    ! z方向の微分とr遠方は関数値がゼロになることを利用している
    function math_diff_2r_9point_realfunc(f) result(f_2r)
        use constants_and_parameters, only : dr,z_center
        implicit none
        real(real64),intent(in) :: f(:,:)
        real(real64),allocatable :: f_2r(:,:)

        real(real64) :: weight(-4:4) = [-9d0,128d0,-1008d0,8064d0,-14350d0,8064d0,-1008d0,128d0,-9d0]
        real(real64) :: f_tmp(-4:4) = 0d0
        integer(int32) :: i,j,k
        integer(int32) :: Nr,Nz

        allocate(f_2r,mold=f)

        f_2r(:,:) = 0d0
        Nr = size(f,1)
        Nz = size(f,2)
        do i=1,Nr
        do j=1,Nz
            do k=-4,-1 ! 負の方向について
                if(i+k<1) then ! f(0) -> f(1), f(-1) -> f(2),...
                    if(abs(f(1,int(z_center/dr))) < eps)then
                        f_tmp(k) = -f(-k+1-i,j)
                    else
                        f_tmp(k) = f(-k+1-i,j)
                    end if
                else
                    f_tmp(k) = f(i+k,j)
                end if
                f_2r(i,j) = f_2r(i,j) + weight(k)*f_tmp(k)
            end do
            f_2r(i,j) = f_2r(i,j) + weight(0)*f(i,j)
            do k=1,4 ! 正の方向について
                if(i+k>Nr) then
                    f_tmp(k) = 0d0 ! 境界の外側はゼロになる
                else
                    f_tmp(k) = f(i+k,j)
                end if
                f_2r(i,j) = f_2r(i,j) + weight(k)*f_tmp(k)
            end do
            f_2r(i,j) = f_2r(i,j)/(5040d0*dr*dr)
        end do 
        end do
        
    end function math_diff_2r_9point_realfunc
    
    function math_diff_2z_9point_realfunc(f) result(f_2z)
        use constants_and_parameters, only : dz
        implicit none
        real(real64),intent(in) :: f(:,:)
        real(real64),allocatable :: f_2z(:,:)
    
        real(real64) :: weight(-4:4) = [-9d0,128d0,-1008d0,8064d0,-14350d0,8064d0,-1008d0,128d0,-9d0]
        real(real64) :: f_tmp(-4:4) = 0d0
        integer(int32) :: i,j,k
        integer(int32) :: Nr,Nz
    
        allocate(f_2z,mold=f)
    
        f_2z(:,:) = 0d0
        Nr = size(f,1)
        Nz = size(f,2)
    
        do i=1,Nr
        do j=1,Nz
            do k=-4,-1
                if(j+k<1) then
                    f_tmp(k) = 0d0
                else
                    f_tmp(k) = f(i,j+k)
                end if
                f_2z(i,j) = f_2z(i,j) + weight(k)*f_tmp(k)
            end do
            f_2z(i,j) = f_2z(i,j) + weight(0)*f(i,j)
            do k=1,4
                if(j+k>Nz) then
                    f_tmp(k) = 0d0
                else
                    f_tmp(k) = f(i,j+k)
                end if
                f_2z(i,j) = f_2z(i,j) + weight(k)*f_tmp(k)
            end do
            f_2z(i,j) = f_2z(i,j)/(5040d0*dz*dz)
        end do 
        end do
      end function math_diff_2z_9point_realfunc
    
    function math_diff_r_9point_realfunc(f) result(f_r)
        use constants_and_parameters, only : dr,z_center
        implicit none
        real(real64),intent(in) :: f(:,:)
        real(real64),allocatable :: f_r(:,:)

        real(real64) :: weight(-4:4) = [3d0,-32d0,168d0,-672d0,0d0,672d0,-168d0,32d0,-3d0]
        real(real64) :: f_tmp(-4:4) = 0d0
        integer(int32) :: i,j,k
        integer(int32) :: Nr,Nz

    allocate(f_r,mold=f)

    f_r(:,:) = 0d0
    Nr = size(f,1)
    Nz = size(f,2)
    do i=1,Nr
    do j=1,Nz
        do k=-4,-1
            if(i+k<1) then ! f(-1)を参照するならf(1)を参照する
                if(abs(f(1,int(z_center/dr))) < eps)then
                    f_tmp(k) = -f(-k+1-i,j)
                else
                    f_tmp(k) = f(-k+1-i,j)
                end if
            else
                f_tmp(k) = f(i+k,j)
            end if            
            f_r(i,j) = f_r(i,j) + weight(k)*f_tmp(k)
        end do
        do k=1,4
            if(i+k>Nr) then
                f_tmp(k) = 0d0
            else
                f_tmp(k) = f(i+k,j)
            end if
            f_r(i,j) = f_r(i,j) + weight(k)*f_tmp(k)
        end do
        f_r(i,j) = f_r(i,j)/(840d0*dr)
    end do 
    end do
    end function math_diff_r_9point_realfunc

    function math_diff_r_9point_complexfunc(f) result(f_r)
        use constants_and_parameters, only : dr
        implicit none
        complex(real64),intent(in) :: f(:,:)
        complex(real64),allocatable :: f_r(:,:)

        real(real64) :: weight(-4:4) = [3d0,-32d0,168d0,-672d0,0d0,672d0,-168d0,32d0,-3d0]
        complex(real64) :: f_tmp(-4:4) = 0d0
        integer(int32) :: i,j,k
        integer(int32) :: Nr,Nz

    allocate(f_r,mold=f)

    f_r(:,:) = 0d0
    Nr = size(f,1)
    Nz = size(f,2)

    do i=1,Nr
    do j=1,Nz
        do k=-4,-1
            if(i+k<1) then
                f_tmp(k) = f(-k+1-i,j)
            else
                f_tmp(k) = f(i+k,j)
            end if
            f_r(i,j) = f_r(i,j) + weight(k)*f_tmp(k)
        end do
        do k=1,4
            if(i+k>Nr) then
                f_tmp(k) = 0d0
            else
                f_tmp(k) = f(i+k,j)
            end if
            f_r(i,j) = f_r(i,j) + weight(k)*f_tmp(k)
        end do
        f_r(i,j) = f_r(i,j)/(840d0*dr)
    end do 
    end do
    end function math_diff_r_9point_complexfunc 

  function math_diff_z_9point_realfunc(f) result(f_z)
    use constants_and_parameters, only : dz
    implicit none
    real(real64),intent(in) :: f(:,:)
    real(real64),allocatable :: f_z(:,:)

    real(real64) :: weight(-4:4) = [3d0,-32d0,168d0,-672d0,0d0,672d0,-168d0,32d0,-3d0]
    real(real64) :: f_tmp(-4:4) = 0d0
    integer(int32) :: i,j,k
    integer(int32) :: Nr,Nz

    allocate(f_z,mold=f)

    f_z(:,:) = 0d0
    Nr = size(f,1)
    Nz = size(f,2)

    do i=1,Nr
    do j=1,Nz
        do k=-4,-1
            if(j+k<1) then
                f_tmp(k) = f(i,-k+1-j)
            else
                f_tmp(k) = f(i,j+k)
            end if
            f_z(i,j) = f_z(i,j) + weight(k)*f_tmp(k)
        end do
        do k=1,4
            if(j+k>Nz) then
                f_tmp(k) = 0d0
            else
                f_tmp(k) = f(i,j+k)
            end if
            f_z(i,j) = f_z(i,j) + weight(k)*f_tmp(k)
        end do
        f_z(i,j) = f_z(i,j)/(840d0*dz)
    end do 
    end do

  end function math_diff_z_9point_realfunc

  function math_diff_z_9point_complexfunc(f) result(f_z)
    use constants_and_parameters, only : dz
    implicit none
    complex(real64),intent(in) :: f(:,:)
    complex(real64),allocatable :: f_z(:,:)

    real(real64) :: weight(-4:4) = [3d0,-32d0,168d0,-672d0,0d0,672d0,-168d0,32d0,-3d0]
    complex(real64) :: f_tmp(-4:4) = 0d0
    integer(int32) :: i,j,k
    integer(int32) :: Nr,Nz

    allocate(f_z,mold=f)

    f_z(:,:) = 0d0
    Nr = size(f,1)
    Nz = size(f,2)

    do i=1,Nr
    do j=1,Nz
        do k=-4,-1
            if(j+k<1) then
                f_tmp(k) = f(i,-k+1-j)
            else
                f_tmp(k) = f(i,j+k)
            end if
            f_z(i,j) = f_z(i,j) + weight(k)*f_tmp(k)
        end do
        do k=1,4
            if(j+k>Nz) then
                f_tmp(k) = 0d0
            else
                f_tmp(k) = f(i,j+k)
            end if
            f_z(i,j) = f_z(i,j) + weight(k)*f_tmp(k)
        end do
        f_z(i,j) = f_z(i,j)/(840d0*dz)
    end do 
    end do
  end function math_diff_z_9point_complexfunc
  

  
end module math_derivation