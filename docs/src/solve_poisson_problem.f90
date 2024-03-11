module poisson_problem
    contains

    subroutine use_DGESV(N,A,b,res_x)
        use,intrinsic :: iso_fortran_env
        implicit none
        integer,intent(in) :: N
        real(kind=8),intent(in) :: A(:,:)
        real(kind=8),intent(in) :: b(:)
        real(kind=8),intent(out) :: res_x(:)
        
        integer(int32) :: lda, ldb, NRHS, INFO
        integer(int32),allocatable :: IPIV(:)
        real(real64),allocatable :: A_lapack(:,:), b_lapack(:)
        
        allocate(IPIV(N))
        allocate(A_lapack(N,N))
        allocate(b_lapack(N))
        
        A_lapack(:,:) = A(:,:)
        b_lapack(:) = b(:)
        
        lda = N
        ldb = N
        NRHS = 1
        call DGESV(N, NRHS, A_lapack, lda, IPIV, b_lapack, ldb, INFO)
        if(INFO/=0)then
            write(*,*) "Error in DGESV"
            write(*,*) INFO
            stop
        end if
        res_x(:) = b_lapack(:)
        write(*,*) "DGESV finished"
        !write(*,*) b_lapack
        deallocate(IPIV)
        deallocate(A_lapack)
        deallocate(b_lapack)
    end subroutine use_DGESV

    ! r=0付近は偶関数の条件→行列の生成時に取り入れる
    ! r=Nr,z=0,z=Nzはディリクレ境界条件→解くべき行列の大きさは(Nr-1)*(Nz-2)

    subroutine nine_point_poisson_by_lapack(rho_p,coulomb)
        use,intrinsic :: iso_fortran_env
        use constants_and_parameters, e2 => elementary_charge2
        use global_variables,only:r_vec
        real(real64),intent(in)  :: rho_p(:,:)
        real(real64),intent(out) :: coulomb(:,:)
        real(real64),allocatable :: A(:,:),b(:),res_x(:) ! solve Ax=b. res_x is the solution
        real(real64) :: dd_weight(-4:4) = [-9d0,128d0,-1008d0,8064d0,-14350d0,8064d0,-1008d0,128d0,-9d0]
        real(real64) :: d_weight(-4:4) = [3d0,-32d0,168d0,-672d0,0d0,672d0,-168d0,32d0,-3d0]
        real(real64) :: ddr_coef(-4:4), ddz_coef(-4:4)
        real(real64) :: dr_coef(-4:4), dz_coef(-4:4)
        real(real64) :: f_tmp(-4:4)
        integer(int32) :: size_A
        integer(int32) :: i,j,k
        integer(int32) :: r_idx,z_idx,ir_idx,iz_idx
        logical        :: is_outside_r,is_outside_z
        size_A = (Nr-1)*(Nz-2)
        allocate(A(size_A,size_A))
        allocate(b(size_A))
        allocate(res_x(size_A))
        f_tmp(:) = 0d0
        ddr_coef(:) = dd_weight(:)/(5040d0*dr*dr)
        ddz_coef(:) = dd_weight(:)/(5040d0*dz*dz)
        dr_coef(:)  = d_weight(:)/(840d0*dr)
        dz_coef(:)  = d_weight(:)/(840d0*dz)

        A(:,:) = 0d0

        do k=1,size_A
            r_idx = modulo(k-1,Nr-1)+1
            z_idx = (k-1)/(Nr-1)+2

            do i=-4,4
                is_outside_r = .false.
                ir_idx = r_idx+i
                if(ir_idx>=1.and.ir_idx<=Nr-1)then ! rは1からNr-1までを考えている
                    ir_idx = ir_idx
                elseif(ir_idx<1)then ! 偶関数の条件
                    ir_idx = -(ir_idx-1)
                elseif(ir_idx>Nr-1)then ! rの外に出るとき（加わるものをゼロにする）
                    is_outside_r = .true.
                end if
                if(is_outside_r)then
                    ! 何もしない
                else 
                    A(k,get_k(ir_idx,z_idx)) = A(k,get_k(ir_idx,z_idx)) + ddr_coef(i) + dr_coef(i)*r_vec(r_idx)
                end if
            end do
            do j=-4,4
                is_outside_z = .false.
                iz_idx = z_idx+j
                if(iz_idx>=2.and.iz_idx<=Nz-1)then ! zは2からNz-1までを考えている  
                    iz_idx = iz_idx
                elseif(iz_idx<2)then 
                    is_outside_z = .true.
                elseif(iz_idx>Nr-1)then
                    is_outside_z = .true.
                end if
                if(is_outside_z)then
                    ! 何もしない
                else 
                    A(k,get_k(r_idx,iz_idx)) = A(k,get_k(r_idx,iz_idx)) + ddz_coef(j)
                end if
            end do
        end do

        ! ディリクレ境界上では値を0とするので，b内で値を引き算する必要はなし．
        ! r=1の条件は行列に反映されている．
        b(:) = 0.0d0 
        ! ここでbを定義する
        do k=1, size_A
            r_idx = modulo(k-1,Nr-1)+1
            z_idx = (k-1)/(Nr-1)+2
            b(k) = -4d0*PI*e2*rho_p(r_idx,z_idx) ! \laplacian \phi = -4\pi e^2 \rho_p
        end do

        call use_DGESV(size_A,A(:,:),b(:),res_x(:))

        ! ここでres_xをcoulombに代入する
        do k=1, size_A
            r_idx = modulo(k-1,Nr-1)+1
            z_idx = (k-1)/(Nr-1)+2
            coulomb(r_idx,z_idx) = res_x(k)
        end do

        deallocate(A)
        deallocate(b)
        deallocate(res_x)
    end subroutine nine_point_poisson_by_lapack



    
    subroutine three_point_poisson_by_lapack(rho_p,A,b)
        use constants_and_parameters,e2 => elementary_charge2
        use,intrinsic :: iso_fortran_env
        implicit none
        real(real64),intent(in)  :: rho_p(:,:)
        real(real64),intent(out) :: A(:,:)
        real(real64),intent(out) :: b(:)
        integer(int32) :: N ! (Nr-1)*(Nz-2)
        real(real64)   :: r
        real(real64)   :: A0, Ap1, Am1, A2 ! 差分方程式の係数
        integer(int32) :: r_idx,z_idx ! 1次元配列のインデックスからi,jを復元するための変数
        integer(int32) :: k
        
        N = (Nr-1)*(Nz-2)
        if(size(A,1)/=N.or.size(A,2)/=N)then
            write(*,*) "Error in make_matrix_of_poisson"
            write(*,*) "size of A is not correct"
            stop
        end if
        if(size(b)/=N)then
            write(*,*) "Error in make_matrix_of_poisson"
            write(*,*) "size of b is not correct"
            stop
        end if        

        A0  = -(2.0d0/(dr*dr) + 2.0d0/(dz*dz))
        Ap1 =  0.0d0 !ここだけはループ内で変更
        Am1 =  0.0d0 !ここだけはループ内で変更
        A2  =  1.0d0/(dz*dz)
        
        
        ! 行列Aの生成
        ! A(k,k)がf(i,j)に対応する。jを固定して先にiを変えていく。
        ! r_idx = modulo(k-1,Nr-1)+1, z_idx = (k-1)/(Nr-1)+2
        ! (z_idx-2)*(Nr-1) + r_idx で一次元配列の位置を復元できる
        A(:,:) = 0.0d0
        
        do k=1,N
            r_idx = modulo(k-1,Nr-1)+1
            z_idx = (k-1)/(Nr-1)+2
            r = dble(r_idx)*dr
            Ap1 = 1d0/(dr)**2 + 1d0/(2d0*dr*r)
            Am1 = 1d0/(dr)**2 - 1d0/(2d0*dr*r)
            
            if(r_idx==1)then
                if(z_idx==2)then
                A(k,                   k) = A(k,k)                    + A0
                A(k,get_k(r_idx+1,z_idx)) = A(k,get_k(r_idx+1,z_idx)) + Ap1 ! f(i+1,j)の係数
                A(k,get_k(1      ,z_idx))   = A(k,get_k(1    ,z_idx))   + Am1 ! f(i-1,j)の係数
                A(k,get_k(r_idx,z_idx+1)) = A(k,get_k(r_idx,z_idx+1)) + A2  ! f(i,j+1)の係数
                !A(k,get_k(r_idx,z_idx-1)) = A(k,get_k(r_idx,z_idx-1)) + A2  ! f(i,j-1)の係数
                else if(z_idx==Nz-1)then
                A(k,                   k) = A(k,k)                    + A0
                A(k,get_k(r_idx+1,z_idx)) = A(k,get_k(r_idx+1,z_idx)) + Ap1 ! f(i+1,j)の係数
                A(k,get_k(1      ,z_idx))   = A(k,get_k(1    ,z_idx))   + Am1 ! f(i-1,j)の係数
                !A(k,get_k(r_idx,z_idx+1)) = A(k,get_k(r_idx,z_idx+1)) + A2  ! f(i,j+1)の係数
                A(k,get_k(r_idx,z_idx-1)) = A(k,get_k(r_idx,z_idx-1)) + A2  ! f(i,j-1)の係数
                else
                A(k,                   k) = A(k,k)                    + A0
                A(k,get_k(r_idx+1,z_idx)) = A(k,get_k(r_idx+1,z_idx)) + Ap1 ! f(i+1,j)の係数
                A(k,get_k(1      ,z_idx))   = A(k,get_k(1    ,z_idx))   + Am1 ! f(i-1,j)の係数
                A(k,get_k(r_idx,z_idx+1)) = A(k,get_k(r_idx,z_idx+1)) + A2  ! f(i,j+1)の係数
                A(k,get_k(r_idx,z_idx-1)) = A(k,get_k(r_idx,z_idx-1)) + A2  ! f(i,j-1)の係数
                end if
            else if(r_idx==Nr-1)then
                if(z_idx==2)then
                A(k,                   k) = A(k,k)                    + A0
                !A(k,get_k(r_idx+1,z_idx)) = A(k,get_k(r_idx+1,z_idx)) + Ap1 ! f(i+1,j)の係数
                A(k,get_k(r_idx-1,z_idx))   = A(k,get_k(r_idx-1,z_idx))   + Am1 ! f(i-1,j)の係数
                A(k,get_k(r_idx,z_idx+1)) = A(k,get_k(r_idx,z_idx+1)) + A2  ! f(i,j+1)の係数
                !A(k,get_k(r_idx,z_idx-1)) = A(k,get_k(r_idx,z_idx-1)) + A2  ! f(i,j-1)の係数
                else if(z_idx==Nz-1)then
                A(k,                   k) = A(k,k)                    + A0
                !A(k,get_k(r_idx+1,z_idx)) = A(k,get_k(r_idx+1,z_idx)) + Ap1 ! f(i+1,j)の係数
                A(k,get_k(r_idx-1,z_idx))   = A(k,get_k(r_idx-1,z_idx))   + Am1 ! f(i-1,j)の係数
                !A(k,get_k(r_idx,z_idx+1)) = A(k,get_k(r_idx,z_idx+1)) + A2  ! f(i,j+1)の係数
                A(k,get_k(r_idx,z_idx-1)) = A(k,get_k(r_idx,z_idx-1)) + A2  ! f(i,j-1)の係数
                else 
                A(k,                   k) = A(k,k)                    + A0
                !A(k,get_k(r_idx+1,z_idx)) = A(k,get_k(r_idx+1,z_idx)) + Ap1 ! f(i+1,j)の係数
                A(k,get_k(r_idx-1,z_idx))   = A(k,get_k(r_idx-1,z_idx))   + Am1 ! f(i-1,j)の係数
                A(k,get_k(r_idx,z_idx+1)) = A(k,get_k(r_idx,z_idx+1)) + A2  ! f(i,j+1)の係数
                A(k,get_k(r_idx,z_idx-1)) = A(k,get_k(r_idx,z_idx-1)) + A2  ! f(i,j-1)の係数
                end if
            else
                if(z_idx==2)then
                A(k,                   k) = A(k,k)                    + A0
                A(k,get_k(r_idx+1,z_idx)) = A(k,get_k(r_idx+1,z_idx)) + Ap1 ! f(i+1,j)の係数
                A(k,get_k(r_idx-1,z_idx))   = A(k,get_k(r_idx-1,z_idx))   + Am1 ! f(i-1,j)の係数
                A(k,get_k(r_idx,z_idx+1)) = A(k,get_k(r_idx,z_idx+1)) + A2  ! f(i,j+1)の係数
                !A(k,get_k(r_idx,z_idx-1)) = A(k,get_k(r_idx,z_idx-1)) + A2  ! f(i,j-1)の係数
                else if(z_idx==Nz-1)then
                A(k,                   k) = A(k,k)                    + A0
                A(k,get_k(r_idx+1,z_idx)) = A(k,get_k(r_idx+1,z_idx)) + Ap1 ! f(i+1,j)の係数
                A(k,get_k(r_idx-1,z_idx))   = A(k,get_k(r_idx-1,z_idx))   + Am1 ! f(i-1,j)の係数
                !A(k,get_k(r_idx,z_idx+1)) = A(k,get_k(r_idx,z_idx+1)) + A2  ! f(i,j+1)の係数
                A(k,get_k(r_idx,z_idx-1)) = A(k,get_k(r_idx,z_idx-1)) + A2  ! f(i,j-1)の係数
                else 
                A(k,                   k) = A(k,k)                    + A0
                A(k,get_k(r_idx+1,z_idx)) = A(k,get_k(r_idx+1,z_idx)) + Ap1 ! f(i+1,j)の係数
                A(k,get_k(r_idx-1,z_idx))   = A(k,get_k(r_idx-1,z_idx))   + Am1 ! f(i-1,j)の係数
                A(k,get_k(r_idx,z_idx+1)) = A(k,get_k(r_idx,z_idx+1)) + A2  ! f(i,j+1)の係数
                A(k,get_k(r_idx,z_idx-1)) = A(k,get_k(r_idx,z_idx-1)) + A2  ! f(i,j-1)の係数
                end if
            end if
        end do
        
        !call write_mat(A(:,:),N)
        
        ! ディリクレ境界条件は、値が0なのでそのままでよく、r=1の条件は行列に反映されている
        b(:) = 0.0d0 
        ! ここでbを定義する
        do k=1, N
            r_idx = modulo(k-1,Nr-1)+1
            z_idx = (k-1)/(Nr-1)+2
            r = dble(r_idx)*dr
            b(k) = -4d0*PI*e2*rho_p(r_idx,z_idx)
        end do
        
        ! TODO: ここでAとbをDGESVに渡して解く
    end subroutine three_point_poisson_by_lapack
    
    subroutine write_mat(A,N)
        use,intrinsic :: iso_fortran_env
        real(real64),intent(in) :: A(:,:)
        integer,intent(in) :: N
        integer(int32) :: i,j
        
        write(*,fmt='(a)') "Matrix A:"
        do j=1, N
           do i=1, N
              write(*,fmt='(x,f8.3)',advance='no') A(j,i)
           end do
           write(*,*)""
        end do
    end subroutine write_mat
    
    function get_k(r_idx,z_idx) result(k_idx)
        use constants_and_parameters,only: Nr
        implicit none
        integer,intent(in) :: r_idx,z_idx
        integer :: k_idx
        k_idx = (z_idx-2)*(Nr-1) + r_idx
    end function get_k
    
    
end module poisson_problem
