module spacial_function
    use,intrinsic :: iso_fortran_env
    implicit none
    contains
    function Laguerre_nl(n,l,x) result(res)
        implicit none
        real(real64) :: res,tes
        integer(int32),intent(in) :: n,l
        real(real64),intent(in) :: x
        integer(int32) :: i
        
        res = 0d0
        do i=0,n
            res = res + (-x)**i/(gamma(dble(i+l+1))*gamma(dble(n-i+1))*gamma(dble(i+1)))
        end do
        tes = (gamma(dble(n+l+1)))
        res = res*(gamma(dble(n+l+1))) ! 足立先生のプリントとは違う定義
    end function Laguerre_nl
      
    function Hermite_n(n,x)  result(res)
        implicit none
        real(real64) :: res
        integer(int32),intent(in) :: n
        real(real64),intent(in) :: x
        integer(int32) :: i
        
        res = 0d0
        do i=0,floor(dble(n)/2d0)
            res = res + (-1)**i*gamma(dble(n+1))/(gamma(dble(i+1))*gamma(dble(n-2*i+1)))*(2d0*x)**(n-2*i)
        end do
    end function Hermite_n 
end module spacial_function