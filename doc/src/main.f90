program main            
    use,intrinsic :: iso_fortran_env,only:int32
    use,intrinsic :: ieee_arithmetic,only:ieee_is_nan
    use :: write_to_file,only: outputsettings   !フォルダ作成
    use :: constants_and_parameters,only:foldername,max_iter,str_z
    implicit none
    integer(int32) :: ti1,tf1,tr1
    character(40) :: str
    integer(int32) :: num_iter
    call system_clock(ti1)
    call outputsettings
    call read_input_parameters
    call alloc_fields
    call initial_settings   ! 波動関数の初期化，調和f振動子に基づく1粒子エネルギーを設定．ソートも行っている
    call calc_ground_state
    call dealloc_fields
    call system_clock(tf1,tr1)
    write(str,'(i0)') int((tf1-ti1)/dble(tr1))
    print*, "実行時間は"//trim(str)//"秒"
    print*, "python3 ./../src/plotall.py "//foldername//" -z "//trim(str_z)
    print*, "endrecord"
    contains

    subroutine is_converge_for_one_E(is_converged)
        use,intrinsic :: iso_fortran_env,only:int32,real64
        use constants_and_parameters,only:SD_eps, prev_eps, num_n,num_p, is_minus_var, is_coulomb &
                                        , is_update_coulomb, coulomb_updated_flag, cou_eps
        use math_integrate,only:volume_integrate
        use global_variables,only:one_E => prep_one_particle_E, one_E2 => prep_one_particle_E2 &
                                 ,idx=>prep_idx_array, Avarage_variance
        implicit none
        logical,intent(out) :: is_converged
        integer(int32) :: k
        real(real64) :: avg, E_Variance, tmp_var, one_e_diff, prev_e_sum, one_e_sum, diff_SD
        character(6) :: str_num
        write(str_num,'(i6)') num_iter

        E_Variance = 0d0
        one_e_diff = 0d0
        do k=1,num_n
            one_e_sum = one_e_sum + one_e(k,1)
        end do
        do k=1,num_p
            one_e_sum = one_e_sum + one_e(k,2)
        end do
        one_e_diff = (one_e_sum - prev_e_sum)/dble(num_n+num_p)
        prev_e_sum = one_e_sum/dble(num_n+num_p)

        do k=1,num_n
            tmp_var = one_E2(idx(k,1),1) - one_E(idx(k,1),1)**2
            E_Variance = E_Variance + tmp_var
            if(tmp_var<0) then
                is_minus_var = .true.
                print*, tmp_var
            end if
        end do

        do k=1,num_p
            tmp_var = one_E2(idx(k,2),2) - one_E(idx(k,2),2)**2
            E_Variance = E_Variance + tmp_var
            if(tmp_var<0) then
                is_minus_var = .true.
                print*, tmp_var
            end if
        end do

        avg = E_Variance/dble(num_n+num_p)
        diff_SD = abs(sqrt(Avarage_variance) - sqrt(avg))
        Avarage_variance = avg
        ! クーロンを更新するかどうかの判定
        if( (is_coulomb .eqv. .true.) .and. (modulo(num_iter,500)==0) )then
            is_update_coulomb = .true.
            print*, trim(str_num)//"回目にクーロン力を更新．"
        end if
        
        ! 収束判定 coulomb入れる場合は，アップデート後の分散がepsを下回っているかどうか
        if(is_coulomb .eqv. .false.)then
            if(diff_SD < SD_eps)then
                is_converged = .true.
            end if
        else
            if(diff_SD < SD_eps .and. one_e_diff < prev_eps)then
                is_converged = .true.
            end if
        end if
        ! アップデートフラグをfalseにもどしておく．trueにするのはcoulombを更新したときのみ
        if(coulomb_updated_flag)then
            coulomb_updated_flag = .false.
        end if
    end subroutine is_converge_for_one_E
    
    subroutine calc_ground_state
        use,intrinsic :: iso_fortran_env,only:int32,real64
        use constants_and_parameters,only : max_iter, plot_iter, is_converged, is_minus_var
        use global_variables,only:Total_Energy, Avarage_variance
        use write_to_file,only:write_prep_sort &
                              ,write_wavefunction,write_d_wavefunction &
                              ,write_density_and_derivative,write_pot,write_meanfield_Energy &
                              ,write_iter_Energy
        implicit none
        character(40)  :: str_iter, str_ave, str_time
        integer(int32) :: num_write  = 100
        integer(int32) :: write_iter = 1
        integer(int32) :: ti,tf,tr
        real(real64)   :: tmp_energy = 100000d0
        real(real64)   :: prev_var = 100000d0
        real(real64)   :: two_prev_var = 100000d0
        if(num_write > max_iter)then
            write_iter = max_iter/num_write
        end if
        call system_clock(ti)
        ! このdo文に入るときには既にidxはソート済み．
        ! 初期波動関数はinitial_settingsで調和振動子のエネルギーに基づいたソートを行っている
        ! IMAG_Evolutionのサブルーチンで「虚時間発展→前回の密度を使ってエネルギー計算→ソート→正規直交化」を行っている．
        do num_iter=0,max_iter
            if(is_minus_var)then
                print*, "分散がゼロ以下になった．"
                exit
            end if
            if(ieee_is_nan(Total_Energy))then
                write(*,*) "エネルギーがNan．"
                exit
            end if
            if(is_converged)then
                write(str,"(f20.15)") Total_Energy
                write(str_iter,'(i0)') num_iter
                write(*,*) trim(str_iter)//"回目でマイナス"//trim(str)//"MeVに収束した．"
                exit
            end if
            if(modulo(num_iter,500)==0 .and. num_iter /=0)then
                call system_clock(tf,tr)
                write(str_time,'(i0)') int((tf-ti)/dble(tr))
                call system_clock(ti)
                write(str_iter,'(i0)') num_iter
                write(str_ave,"(f20.15)") Avarage_variance
                ! char型変数の長さも考慮しないといけない
                write(*,*) trim(str_iter)//"回目の計算で分散は"//trim(str_ave)//"MeV&
                            前回から"//trim(str_time)//"秒かかった．"
                !if(prev_var < Avarage_variance .and. two_prev_var < prev_var)then
                !    write(*,*) "分散が増えた．"
                !    exit
                !end if
                two_prev_var = prev_var
                prev_var = Avarage_variance
            end if
            tmp_energy = Total_Energy
            call calc_derivation_of_wavefunction       ! prep_**配列の計算
            call copy_to_occupied_wf                   ! 更新されているidxを使って占有軌道の配列にコピー
            call calc_densities                        ! 密度の計算
            call calc_meanfield                        ! 占有軌道を使って平均場を計算
            call calc_observables                      ! 密度で決定される物理量を計算（Skyrmeエネルギーや四重極モーメントなど）
            if(modulo(num_iter,write_iter) == 0)then
                call write_prep_sort(num_iter)
                call write_iter_Energy(num_iter)
                call write_meanfield_Energy(num_iter)                ! write meanfield energy
            end if
            if(modulo(num_iter,plot_iter) == 0 .or. num_iter == 0 .or. num_iter == max_iter)then
                call write_wavefunction(num_iter)          ! 波動関数の書き出し
                call write_d_wavefunction(num_iter)        ! 波動関数の微分の書き出し
                call write_density_and_derivative(num_iter) ! 密度とその微分の書き出し
                call write_pot(num_iter)                   ! 平均場の書き出し
                print*, "．重心の位置は"//trim(str_z)//"fm"
                print*, "python3 ./../src/plotall.py "//foldername//" -z "//trim(str_z)
            end if
            call calc_hf_by_imaginary_time_evolution   ! 虚時間発展，波動関数の更新，エネルギーの計算，idxソート，正規直交化
        end do
        write(str,"(f20.15)") Total_Energy
        if(is_converged .eqv. .false.)then
            if(Total_Energy>0)then
                print*, "収束しなかった．エネルギーは"//trim(str)//"MeV"
            else
                print*, "収束しなかった．エネルギーは"//trim(str)//"MeV．"
            end if
        end if
        call write_wavefunction(num_iter)          ! 波動関数の書き出し
        call write_d_wavefunction(num_iter)        ! 波動関数の微分の書き出し
        call write_density_and_derivative(num_iter) ! 密度とその微分の書き出し
        call write_pot(num_iter)                   ! 平均場の書き出し
    end subroutine calc_ground_state

    subroutine calc_densities
        implicit none
        call calc_density_and_lap               ! ρとΔρの計算 
        call calc_kinetic_density               ! τの計算
        call calc_spin_orbit_density_and_div    ! JとdivJの計算

    end subroutine calc_densities
    
    subroutine calc_observables
        use determine_observables
        implicit none
        call calc_Skyrme_Energy_by_integrate
        call calc_Kinetic_Energy_by_integrate
        call calc_CM_Energy
        call calc_Coulomb_Energy
        call calc_Total_Energy
        
    end subroutine calc_observables
    
    
    subroutine read_input_parameters
        use,intrinsic :: iso_fortran_env,only:int32
        use global_variables
        use constants_and_parameters
        implicit none
        integer(int32) :: fi,is
        character(4)  :: str
        namelist /nuclear_param/ num_n, num_p, prep_max_r, prep_max_z, prep_max_m, nuc_name
        namelist /calc_param/ max_iter, lagrange_multiplier, imaginary_time_step, cou_eps, SD_eps, prev_eps &
                             ,is_use_dd, is_output, not_divJ, is_use_B, is_use_W &
                             ,is_TimeR, is_coulomb, is_use_second_imag, Nr, Nz, plot_iter
        open(unit=fi,file="./../src/input_parameters.txt",status="old",iostat=is)
        read(fi,nuclear_param,iostat=is)
        if (is /= 0) stop "Error reading nuclear_param."
            
        read(fi,calc_param,iostat=is)
        if (is /= 0) stop "Error reading calc_param."
        close(fi)
        prep_number = (prep_max_r+1)*(prep_max_z+1)*(2*prep_max_m+1)
        write(str,'(i4)') prep_number
        print*, "計算する軌道の数は"//trim(str)//"個"
        if(prep_number < max(num_n,num_p))then
            write(*,*) "Error: prep_number must be larger than max(num_n,num_p)."
            stop
        end if
        dr = r_max/dble(Nr)        ! [fm] The grid spacing in r-direction
        dz = z_max/dble(Nz)        ! [fm] The grid spacing in z-direction
        ! is_memoを用意して，readしてターミナルから文字列を受け取るようにしてメモを付ける

    end subroutine read_input_parameters

    subroutine initial_settings
        use,intrinsic :: iso_fortran_env,only:int32
        use global_variables,only : r_vec, z_vec, dr, dz, Nr, Nz, z_max, is_TimeR
        implicit none
        integer(int32) :: i

        ! 格子点の位置を決める．特異点を回避するためにrはdr/2だけずらしている．
        do i=1,Nr
            r_vec(i) = dr*dble(i) - dr/2.0d0 
        end do
        do i=1,Nz
            z_vec(i) = dz*dble(i)
        end do
        if(is_TimeR)then
            call prep_time_reversal_wavefunction
        else
            call prep_initial_wavefunction
        end if

    end subroutine initial_settings

    subroutine prep_time_reversal_wavefunction ! 規格化はできている（overlap(k,k)=1）
        use,intrinsic :: iso_fortran_env,only:int32,real64
        use global_variables,only: prep_wf_plus, prep_wf_minus, r_vec,z_vec, z_max &
                                  ,prep_max_r, prep_max_z, prep_max_m &
                                  ,prep_mag => prep_magnetic_q_num, prep_idx_array &
                                  ,prep_E => prep_one_particle_E

        use constants_and_parameters,only: betain => deformation_degree, z_center, prep_number,Nr,Nz, nuc_name &
                                          ,prep_max_m, prep_max_r, prep_max_z
        use harmonic_occ_wf,only : set_harmonic_occ, calc_occ_energy
        use math_integrate,only: volume_integrate
        implicit none
        
        integer(int32) :: i,j,k,np,tmp
        integer(int32)   :: prep_r,prep_z,prep_m
        integer(int32),allocatable :: rzm_ar(:,:)
        integer(int32) :: prep
        real(real64)   :: z0  ! 原子核の中心（0, z_max/2）
        allocate(rzm_ar(prep_number,3))
        z0 = z_center

        ! 取りうる量子数の設定
        prep=1
        do k=-prep_max_m,prep_max_m
        do i=0,prep_max_r
        do j=0,prep_max_z
            rzm_ar(prep,1) = i
            rzm_ar(prep,2) = j
            rzm_ar(prep,3) = k
            prep_idx_array(prep,1) = prep
            prep_idx_array(prep,2) = prep
            prep = prep + 1
        end do
        end do
        end do

        ! Time-reversalを課して方位量子数を設定，配列も直す
        do np=1,2
            do prep=1,prep_number,2
                prep_mag(prep,np) = rzm_ar(prep,3)
                prep_mag(prep+1,np) = -prep_mag(prep,np) - 1 ! time-reversal
                rzm_ar(prep+1,3) = -prep_mag(prep,np) - 1 ! time-reversal
            end do
        end do

        do np=1,2
        do prep=1,prep_number,2
            do j=1,Nz
            do i=1,Nr
            prep_wf_plus(i,j,prep,np) &
            = set_harmonic_occ(rzm_ar(prep,1),rzm_ar(prep,2),prep_mag(prep,np),r_vec(i),z_vec(j)-z0,betain)/sqrt(2d0)
            prep_wf_minus(i,j,prep,np) &
            = set_harmonic_occ(rzm_ar(prep,1),rzm_ar(prep,2),prep_mag(prep,np)+1,r_vec(i),z_vec(j)-z0,betain)/sqrt(2d0)

            prep_wf_plus(i,j,prep+1,np)  = -prep_wf_minus(i,j,prep,np) ! time-reversal
            prep_wf_minus(i,j,prep+1,np) =  prep_wf_plus(i,j,prep,np) ! time-reversal
            end do
            end do
        end do
        end do
        do np=1,2
            do prep=1,prep_number
                prep_E(prep,np) = calc_occ_energy(rzm_ar(prep,1),rzm_ar(prep,2),prep_mag(prep,np),betain)
            end do
        end do  


        ! エネルギーの小さい順にソート
        call sort_energy_and_idx(prep_E(:,:),prep_idx_array(:,:))
    end subroutine prep_time_reversal_wavefunction
    

    ! 波動関数の初期化，調和振動子に基づく1粒子エネルギーを設定
    subroutine prep_initial_wavefunction ! 規格化はできている（overlap(k,k)=1）
        use,intrinsic :: iso_fortran_env,only:int32,real64
        use global_variables,only: prep_wf_plus, prep_wf_minus, r_vec,z_vec, z_max &
                                    ,prep_max_r, prep_max_z, prep_max_m &
                                    ,prep_mag => prep_magnetic_q_num, prep_idx_array &
                                    ,prep_E => prep_one_particle_E, init_rzm_to_write

        use constants_and_parameters,only: betain => deformation_degree, z_center, prep_number,Nr,Nz &
                                            ,prep_max_m, prep_max_r, prep_max_z
        use harmonic_occ_wf,only : set_harmonic_occ, calc_occ_energy
        use math_integrate,only: volume_integrate
        implicit none
        
        integer(int32) :: i,j,k,np
        integer(int32),allocatable :: rzm_ar(:,:)
        integer(int32) :: prep
        real(real64)   :: z0  ! 原子核の中心（0, z_max/2）
        character(1)   :: str_r
        character(1)   :: str_z
        character(2)   :: str_m
        allocate(rzm_ar(prep_number,3))
        z0 = z_center
        !!!!!!--------- prep_numberも変更してます -----------！！！！！！！！！！！

        ! 取りうる量子数の設定
        prep=1
        do k=-prep_max_m,prep_max_m
        do i=0,prep_max_r
        do j=0,prep_max_z
            rzm_ar(prep,1) = i
            rzm_ar(prep,2) = j
            rzm_ar(prep,3) = k
            prep_idx_array(prep,1) = prep
            prep_idx_array(prep,2) = prep
            write(str_r,'(i1)') i
            write(str_z,"(i1)") j
            write(str_m,"(i2)") k
            init_rzm_to_write(prep) = str_r//str_z//str_m
            prep = prep + 1
        end do
        end do
        end do

        do np=1,2
            do prep=1,prep_number
                prep_mag(prep,np) = rzm_ar(prep,3)
            end do
        end do

        do np=1,2
        do prep=1,prep_number
            do j=1,Nz
            do i=1,Nr
            prep_wf_plus(i,j,prep,np) &
            = set_harmonic_occ(rzm_ar(prep,1),rzm_ar(prep,2),prep_mag(prep,np),r_vec(i),z_vec(j)-z0,betain)/sqrt(2d0)
            prep_wf_minus(i,j,prep,np) &
            = set_harmonic_occ(rzm_ar(prep,1),rzm_ar(prep,2),prep_mag(prep,np)+1,r_vec(i),z_vec(j)-z0,betain)/sqrt(2d0)
            end do
            end do
        end do
        end do
        do np=1,2
            do prep=1,prep_number
                prep_E(prep,np) = calc_occ_energy(rzm_ar(prep,1),rzm_ar(prep,2),prep_mag(prep,np),betain)
            end do
        end do

        do np=1,2
            do prep=1,prep_number
                do j=1,Nz
                do i=1,Nr
                prep_wf_plus(i,j,prep,np) &
                = set_harmonic_occ(rzm_ar(prep,1),rzm_ar(prep,2),prep_mag(prep,np),r_vec(i),z_vec(j)-z0,betain)/sqrt(2d0)
                prep_wf_minus(i,j,prep,np) &
                = set_harmonic_occ(rzm_ar(prep,1),rzm_ar(prep,2),prep_mag(prep,np)+1,r_vec(i),z_vec(j)-z0,betain)/sqrt(2d0)
                end do
                end do
            end do
            end do
            do np=1,2
                do prep=1,prep_number
                    prep_E(prep,np) = calc_occ_energy(rzm_ar(prep,1),rzm_ar(prep,2),prep_mag(prep,np),betain)
                end do
            end do

        ! エネルギーの小さい順にソート
        call sort_energy_and_idx(prep_E(:,:),prep_idx_array(:,:))
    end subroutine prep_initial_wavefunction

    subroutine calc_derivation_of_wavefunction
        use :: constants_and_parameters
        use :: math_derivation , only: math_diff_r_9point,math_diff_z_9point,math_diff_2r_9point,math_diff_2z_9point
        use :: global_variables, only:prep_wf_plus,prep_wf_minus &
                                     ,prep_dr_wf_plus,prep_dr_wf_minus,prep_dz_wf_plus,prep_dz_wf_minus &
                                     ,prep_ddr_wf_plus,prep_ddr_wf_minus,prep_ddz_wf_plus,prep_ddz_wf_minus &
                                     ,prep_lap_wf_plus,prep_lap_wf_minus &
                                     ,r_vec , prep_mag => prep_magnetic_q_num &
                                     ,prep_SDNWP, prep_SDNWM

        implicit none
        integer(int32) :: prep, i, j, l
        do l=1,2
        do prep=1,prep_number
            prep_dr_wf_plus(:,:,prep,l)   = math_diff_r_9point(f=prep_wf_plus(:,:,prep,l))
            prep_dz_wf_plus(:,:,prep,l)   = math_diff_z_9point(f=prep_wf_plus(:,:,prep,l))
            if(is_use_dd)then
                prep_ddr_wf_plus(:,:,prep,l)  = math_diff_2r_9point(f=prep_wf_plus(:,:,prep,l))
                prep_ddz_wf_plus(:,:,prep,l)  = math_diff_2z_9point(f=prep_wf_plus(:,:,prep,l))
            else
                prep_ddr_wf_plus(:,:,prep,l)  = math_diff_r_9point(f=prep_dr_wf_plus(:,:,prep,l))
                prep_ddz_wf_plus(:,:,prep,l)  = math_diff_z_9point(f=prep_dz_wf_plus(:,:,prep,l))
            end if
            prep_dr_wf_minus(:,:,prep,l)  = math_diff_r_9point(f=prep_wf_minus(:,:,prep,l))
            prep_dz_wf_minus(:,:,prep,l)  = math_diff_z_9point(f=prep_wf_minus(:,:,prep,l))
            if(is_use_dd)then
                prep_ddr_wf_minus(:,:,prep,l) = math_diff_2r_9point(f=prep_wf_minus(:,:,prep,l))
                prep_ddz_wf_minus(:,:,prep,l) = math_diff_2z_9point(f=prep_wf_minus(:,:,prep,l))
            else
                prep_ddr_wf_minus(:,:,prep,l) = math_diff_r_9point(f=prep_dr_wf_minus(:,:,prep,l))
                prep_ddz_wf_minus(:,:,prep,l) = math_diff_z_9point(f=prep_dz_wf_minus(:,:,prep,l))
            end if
        end do
        end do
        ! calculation laplacian of wavefunction
        ! ddf + df/r + ddz - (m/r)^2 f
        do l=1,2
        do prep=1,prep_number
            do j=1,Nz
                do i=1,Nr
                    prep_lap_wf_plus(i,j,prep,l)  = prep_ddr_wf_plus(i,j,prep,l) + prep_dr_wf_plus(i,j,prep,l)/r_vec(i) &
                    + prep_ddz_wf_plus(i,j,prep,l) - (dble(prep_mag(prep,l))/r_vec(i))**2*prep_wf_plus(i,j,prep,l)

                    prep_lap_wf_minus(i,j,prep,l) = prep_ddr_wf_minus(i,j,prep,l) + prep_dr_wf_minus(i,j,prep,l)/r_vec(i)&
                    + prep_ddz_wf_minus(i,j,prep,l) - (dble(prep_mag(prep,l)+1)/r_vec(i))**2*prep_wf_minus(i,j,prep,l)
                end do
            end do
        end do
        end do
        ! caluculation sigma dot nabla of wavefunction
        do l=1,2
        do prep=1,prep_number
            do j=1,Nz
                do i=1,Nr
                    prep_SDNWP(i,j,prep,l) = prep_dz_wf_plus(i,j,prep,l) + prep_dr_wf_minus(i,j,prep,l) &
                                         + dble(prep_mag(prep,l)+1)/r_vec(i)*prep_wf_minus(i,j,prep,l)
                    prep_SDNWM(i,j,prep,l) = -prep_dz_wf_minus(i,j,prep,l) + prep_dr_wf_plus(i,j,prep,l) &
                                         - dble(prep_mag(prep,l))/r_vec(i)*prep_wf_plus(i,j,prep,l)

                end do
            end do
        end do 
        end do
        !print*, "prep_SDNWP(i,j,prep,l) = ",maxval(prep_dz_wf_plus), maxval(prep_dr_wf_minus)!, &
        !!maxval(dble(prep_mag(prep,l)+1)/r_vec(i)*prep_wf_minus)
        !print*, "prep_SDNWM(i,j,prep,l) = ",maxval(prep_dz_wf_minus), maxval(prep_dr_wf_plus)!, &
        !!maxval(-dble(prep_mag(prep,l))/r_vec(i)*prep_wf_plus)
        
    end subroutine calc_derivation_of_wavefunction

    subroutine copy_to_occupied_wf
        use :: constants_and_parameters, betain => deformation_degree
        use :: math_derivation , only: math_diff_r_9point,math_diff_z_9point,math_diff_2r_9point,math_diff_2z_9point
        use :: global_variables, only:wf_plus,wf_minus,dr_wf_plus &
                                     ,dr_wf_minus,dz_wf_plus,dz_wf_minus &
                                     ,ddr_wf_plus,ddr_wf_minus,ddz_wf_plus,ddz_wf_minus &
                                     ,lap_wf_plus,lap_wf_minus &
                                     ,prep_wf_plus,prep_wf_minus &
                                     ,prep_dr_wf_plus,prep_dr_wf_minus,prep_dz_wf_plus,prep_dz_wf_minus &
                                     ,prep_ddr_wf_plus,prep_ddr_wf_minus,prep_ddz_wf_plus,prep_ddz_wf_minus &
                                     ,prep_lap_wf_plus,prep_lap_wf_minus &
                                     ,prep_mag => prep_magnetic_q_num, mag => magnetic_q_num, prep_idx_array &
                                     ,prep_SDNWP, prep_SDNWM, SDNWP, SDNWM 

        implicit none
        integer(int32) :: k, kp



         !波動関数のコピー
        do k = 1,num_n
            mag(k) = prep_mag(prep_idx_array(k,1),1)
            wf_plus(:,:,k)  = prep_wf_plus(:,:,prep_idx_array(k,1),1)
            wf_minus(:,:,k) = prep_wf_minus(:,:,prep_idx_array(k,1),1)
        end do
        do k = 1,num_p
            kp = k + num_n
            mag(kp) = prep_mag(prep_idx_array(k,2),2)
            wf_plus(:,:,kp)  = prep_wf_plus(:,:,prep_idx_array(k,2),2)
            wf_minus(:,:,kp) = prep_wf_minus(:,:,prep_idx_array(k,2),2)
        end do
        ! 導関数とかのコピー
        do k=1,num_n
            dr_wf_plus(:,:,k)   = prep_dr_wf_plus(:,:,prep_idx_array(k,1),1)
            dz_wf_plus(:,:,k)   = prep_dz_wf_plus(:,:,prep_idx_array(k,1),1)
            dr_wf_minus(:,:,k)  = prep_dr_wf_minus(:,:,prep_idx_array(k,1),1)
            dz_wf_minus(:,:,k)  = prep_dz_wf_minus(:,:,prep_idx_array(k,1),1)
            ddr_wf_plus(:,:,k)  = prep_ddr_wf_plus(:,:,prep_idx_array(k,1),1)
            ddz_wf_plus(:,:,k)  = prep_ddz_wf_plus(:,:,prep_idx_array(k,1),1)
            ddr_wf_minus(:,:,k) = prep_ddr_wf_minus(:,:,prep_idx_array(k,1),1)
            ddz_wf_minus(:,:,k) = prep_ddz_wf_minus(:,:,prep_idx_array(k,1),1)

            lap_wf_plus(:,:,k)  = prep_lap_wf_plus(:,:, prep_idx_array(k,1),1)
            lap_wf_minus(:,:,k) = prep_lap_wf_minus(:,:,prep_idx_array(k,1),1)
            SDNWP(:,:,k) = prep_SDNWP(:,:,prep_idx_array(k,1),1)
            SDNWM(:,:,k) = prep_SDNWM(:,:,prep_idx_array(k,1),1)
        end do
        do k=1,num_p
            kp = k + num_n
            dr_wf_plus(:,:,kp)   = prep_dr_wf_plus(:,:,  prep_idx_array(k,2),2)
            dz_wf_plus(:,:,kp)   = prep_dz_wf_plus(:,:,  prep_idx_array(k,2),2)
            dr_wf_minus(:,:,kp)  = prep_dr_wf_minus(:,:, prep_idx_array(k,2),2)
            dz_wf_minus(:,:,kp)  = prep_dz_wf_minus(:,:, prep_idx_array(k,2),2)
            ddr_wf_plus(:,:,kp)  = prep_ddr_wf_plus(:,:, prep_idx_array(k,2),2)
            ddz_wf_plus(:,:,kp)  = prep_ddz_wf_plus(:,:, prep_idx_array(k,2),2)
            ddr_wf_minus(:,:,kp) = prep_ddr_wf_minus(:,:,prep_idx_array(k,2),2)
            ddz_wf_minus(:,:,kp) = prep_ddz_wf_minus(:,:,prep_idx_array(k,2),2)

            lap_wf_plus(:,:,kp)  = prep_lap_wf_plus(:,:, prep_idx_array(k,2),2)
            lap_wf_minus(:,:,kp) = prep_lap_wf_minus(:,:,prep_idx_array(k,2),2)
            SDNWP(:,:,kp) = prep_SDNWP(:,:,prep_idx_array(k,2),2)
            SDNWM(:,:,kp) = prep_SDNWM(:,:,prep_idx_array(k,2),2)
        end do
    end subroutine copy_to_occupied_wf

    subroutine sort_energy_and_idx(energy_array,idx)
        use,intrinsic :: iso_fortran_env,only:int32,real64
        implicit none
        real(real64),intent(inout)    :: energy_array(:,:)
        integer(int32), intent(inout) :: idx(:,:)
        integer(int32) :: i,j,N,l
        real(real64)   :: tmp_energy
        real(real64) :: eps
        integer(int32) :: tmp_idx
        logical,save :: is_first = .true.
        if(size(energy_array,2) /= 2 .or. size(idx,2) /=2)then
            write(*,*) "Error: array (dim2) must be 2."
            stop
        end if
        !stop "adawdas"
        eps = 0d0
        N=size(energy_array,1)
        if(size(idx,1) /= N)then
            write(*,*) "Error: array (dim1) must be same."
            stop
        end if
        ! インデックスの順番を整える（前回の結果が保存されているため初期化する）
        idx(:,1) = [(i,i=1,N)]
        idx(:,2) = [(i,i=1,N)]
        do l=1,2
        do i=1,N-1
            do j=i+1,N
                if(energy_array(i,l) > energy_array(j,l) + eps)then
                    !print*, tmp_energy
                    !stop
                    tmp_energy = energy_array(i,l)
                    tmp_idx = idx(i,l)
                    energy_array(i,l) = energy_array(j,l)
                    idx(i,l) = idx(j,l)
                    energy_array(j,l) = tmp_energy
                    idx(j,l) = tmp_idx
                end if
            end do
        end do
        end do
        !if(is_first)then
        !    is_first = .false.
        !else
        !print*, "energy_array = ",energy_array(:,1)
        !print*, "idx = ",idx(:,1)
        !stop "return sort energy_and_idx"
        !end if
    end subroutine sort_energy_and_idx
    

    subroutine calc_meanfield
        use determine_meanfield_potential, only: calc_pot_U,calc_pot_B,calc_pot_W
        implicit none

        call calc_pot_U ! 今はWood-Saxon型
        call calc_pot_B
        call calc_pot_W
                    
    end subroutine calc_meanfield
    

    subroutine calc_hf_by_imaginary_time_evolution
        use constants_and_parameters,only:is_converged, is_minus_var
        use global_variables,only : prep_E => prep_one_particle_E,prep_idx_array
        implicit none


        call calc_Hpsi
        call calc_one_particle_energy
        call calc_Hpsi2 ! calc  h^2|wf>
        call is_converge_for_one_E(is_converged)
        if(is_minus_var .neqv. .true.)then
            call calc_IMAG_Evolution(1)
            call calc_IMAG_Evolution(2)
            call calc_Hpsi
            call calc_one_particle_energy
            call sort_energy_and_idx(prep_E(:,:),prep_idx_array(:,:)) 
            call Schmidt_orthogonalization(1)
            call Schmidt_orthogonalization(2)
        end if
    end subroutine calc_hf_by_imaginary_time_evolution

    subroutine calc_Hpsi2
        use,intrinsic :: iso_fortran_env,only:int32,real64
        use constants_and_parameters,only:HBAR2_over_2m,Nr,Nz,is_use_B,is_use_W
        use math_derivation,only: math_diff_r_9point,math_diff_z_9point,math_diff_2r_9point,math_diff_2z_9point
        use global_variables,only:prep_number&
                                 ,prep_Hpsi2_minus,prep_Hpsi2_plus &
                                 ,prep_Hpsi_plus, prep_Hpsi_minus, U_n, B_n, W_n, U_p, B_p, W_p &
                                 ,r_vec,prep_mag=>prep_magnetic_q_num 
        integer(int32) :: prep,l,i,j
        real(real64),allocatable :: tmp_UHpsi_plus(:,:,:,:), tmp_UHpsi_minus(:,:,:,:)
        real(real64),allocatable :: tmp_BHpsi_plus(:,:,:,:), tmp_BHpsi_minus(:,:,:,:)
        real(real64),allocatable :: tmp_WHpsi_plus(:,:,:,:), tmp_WHpsi_minus(:,:,:,:)
        real(real64),allocatable :: tmp_lap_Hpsi_plus(:,:,:,:), tmp_lap_Hpsi_minus(:,:,:,:)
        real(real64),allocatable :: tmp_dr_Hpsi_plus(:,:,:,:), tmp_dr_Hpsi_minus(:,:,:,:)
        real(real64),allocatable :: tmp_dz_Hpsi_plus(:,:,:,:), tmp_dz_Hpsi_minus(:,:,:,:)
        real(real64),allocatable :: tmp_ddr_Hpsi_plus(:,:,:,:), tmp_ddr_Hpsi_minus(:,:,:,:)
        real(real64),allocatable :: tmp_ddz_Hpsi_plus(:,:,:,:), tmp_ddz_Hpsi_minus(:,:,:,:)
        real(real64),allocatable :: tmp_SDNHWP(:,:,:,:), tmp_SDNHWM(:,:,:,:)

        real(real64),allocatable  :: dr_B(:,:), dz_B(:,:)
        allocate(tmp_UHpsi_plus,mold=prep_Hpsi_plus)
        allocate(tmp_UHpsi_minus,mold=prep_Hpsi_minus)
        allocate(tmp_BHpsi_plus,mold=prep_Hpsi_plus)
        allocate(tmp_BHpsi_minus,mold=prep_Hpsi_minus)
        allocate(tmp_WHpsi_plus,mold=prep_Hpsi_plus)
        allocate(tmp_WHpsi_minus,mold=prep_Hpsi_minus)
        allocate(tmp_lap_Hpsi_plus,mold=prep_Hpsi_plus)
        allocate(tmp_lap_Hpsi_minus,mold=prep_Hpsi_minus)
        allocate(tmp_dr_Hpsi_plus,mold=prep_Hpsi_plus)
        allocate(tmp_dr_Hpsi_minus,mold=prep_Hpsi_minus)
        allocate(tmp_dz_Hpsi_plus,mold=prep_Hpsi_plus)
        allocate(tmp_dz_Hpsi_minus,mold=prep_Hpsi_minus)
        allocate(tmp_ddr_Hpsi_plus,mold=prep_Hpsi_plus)
        allocate(tmp_ddr_Hpsi_minus,mold=prep_Hpsi_minus)
        allocate(tmp_ddz_Hpsi_plus,mold=prep_Hpsi_plus)
        allocate(tmp_ddz_Hpsi_minus,mold=prep_Hpsi_minus)
        allocate(tmp_SDNHWP(Nr,Nz,prep_number,2)       )  ! [fm^-3/2]  (r,z,nucleon)
        allocate(tmp_SDNHWM(Nr,Nz,prep_number,2)       )  ! [fm^-3/2]  (r,z,nucleon)

        allocate(dr_B(Nr,Nz))
        allocate(dz_B(Nr,Nz))

        do prep=1,prep_number
            tmp_UHpsi_plus(:,:,prep,1)  = U_n(:,:)*prep_Hpsi_plus(:,:,prep,1)
            tmp_UHpsi_minus(:,:,prep,1) = U_n(:,:)*prep_Hpsi_minus(:,:,prep,1)
            tmp_UHpsi_plus(:,:,prep,2)  = U_p(:,:)*prep_Hpsi_plus(:,:,prep,2)
            tmp_UHpsi_minus(:,:,prep,2) = U_p(:,:)*prep_Hpsi_minus(:,:,prep,2)
        end do

        dr_B(:,:) = math_diff_r_9point(B_n(:,:))
        dz_B(:,:) = math_diff_z_9point(B_n(:,:))

        do l=1,2
            do prep=1,prep_number
                tmp_dr_Hpsi_plus(:,:,prep,l)   = math_diff_r_9point(f=prep_Hpsi_plus(:,:,prep,l))
                tmp_dz_Hpsi_plus(:,:,prep,l)   = math_diff_z_9point(f=prep_Hpsi_plus(:,:,prep,l))

                tmp_ddr_Hpsi_plus(:,:,prep,l)  = math_diff_2r_9point(f=prep_Hpsi_plus(:,:,prep,l))
                tmp_ddz_Hpsi_plus(:,:,prep,l)  = math_diff_2z_9point(f=prep_Hpsi_plus(:,:,prep,l))

                tmp_dr_Hpsi_minus(:,:,prep,l)  = math_diff_r_9point(f=prep_Hpsi_minus(:,:,prep,l))
                tmp_dz_Hpsi_minus(:,:,prep,l)  = math_diff_z_9point(f=prep_Hpsi_minus(:,:,prep,l))

                tmp_ddr_Hpsi_minus(:,:,prep,l) = math_diff_2r_9point(f=prep_Hpsi_minus(:,:,prep,l))
                tmp_ddz_Hpsi_minus(:,:,prep,l) = math_diff_2z_9point(f=prep_Hpsi_minus(:,:,prep,l))
            end do
        end do

        do l=1,2
            do prep=1,prep_number
                do j=1,Nz
                    do i=1,Nr
                        tmp_lap_Hpsi_plus(i,j,prep,l)  = tmp_ddr_Hpsi_plus(i,j,prep,l) + tmp_dr_Hpsi_plus(i,j,prep,l)/r_vec(i) &
                        + tmp_ddz_Hpsi_plus(i,j,prep,l) - (dble(prep_mag(prep,l))/r_vec(i))**2*prep_Hpsi_plus(i,j,prep,l)
    
                        tmp_lap_Hpsi_minus(i,j,prep,l) = tmp_ddr_Hpsi_minus(i,j,prep,l) + tmp_dr_Hpsi_minus(i,j,prep,l)/r_vec(i)&
                        + tmp_ddz_Hpsi_minus(i,j,prep,l) - (dble(prep_mag(prep,l)+1)/r_vec(i))**2*prep_Hpsi_minus(i,j,prep,l)
                    end do
                end do
            end do
            end do
        do prep=1,prep_number
            tmp_BHpsi_plus(:,:,prep,1) = B_n(:,:)*tmp_lap_Hpsi_plus(:,:,prep,1) &
                                      + dr_B(:,:)*tmp_dr_Hpsi_plus(:,:,prep,1) &
                                      + dz_B(:,:)*tmp_dz_Hpsi_plus(:,:,prep,1)

            tmp_BHpsi_minus(:,:,prep,1) = B_n(:,:)*tmp_lap_Hpsi_minus(:,:,prep,1) &
                                        + dr_B(:,:)*tmp_dr_Hpsi_minus(:,:,prep,1) &
                                        + dz_B(:,:)*tmp_dz_Hpsi_minus(:,:,prep,1)
        end do

        dr_B(:,:) = math_diff_r_9point(B_p(:,:))
        dz_B(:,:) = math_diff_z_9point(B_p(:,:))
        do prep=1,prep_number
            tmp_BHpsi_plus(:,:,prep,2) =  B_p(:,:)*tmp_lap_Hpsi_plus(:,:,prep,2) &
                                        + dr_B(:,:)*tmp_dr_Hpsi_plus(:,:,prep,2) &
                                        + dz_B(:,:)*tmp_dz_Hpsi_plus(:,:,prep,2)

            tmp_BHpsi_minus(:,:,prep,2) = B_p(:,:)*tmp_lap_Hpsi_minus(:,:,prep,2) &
                                        + dr_B(:,:)*tmp_dr_Hpsi_minus(:,:,prep,2) &
                                        + dz_B(:,:)*tmp_dz_Hpsi_minus(:,:,prep,2)
        end do
        do l=1,2
            do prep=1,prep_number
                do j=1,Nz
                    do i=1,Nr
                        tmp_SDNHWP(i,j,prep,l) = tmp_dz_Hpsi_plus(i,j,prep,l) + tmp_dr_Hpsi_minus(i,j,prep,l) &
                                             + dble(prep_mag(prep,l)+1)/r_vec(i)*prep_Hpsi_minus(i,j,prep,l)
                        tmp_SDNHWM(i,j,prep,l) = -tmp_dz_Hpsi_minus(i,j,prep,l) + tmp_dr_Hpsi_plus(i,j,prep,l) &
                                             - dble(prep_mag(prep,l))/r_vec(i)*prep_Hpsi_plus(i,j,prep,l)
                    end do
                end do
            end do 
        end do

        do prep=1,prep_number
            tmp_WHpsi_plus(:,:,prep,1)  = W_n(:,:,3)*(tmp_SDNHWP(:,:,prep,1)-tmp_dz_Hpsi_plus(:,:,prep,1)) &
                                        + W_n(:,:,1)*(tmp_SDNHWM(:,:,prep,1)-tmp_dr_Hpsi_plus(:,:,prep,1))
            tmp_WHpsi_minus(:,:,prep,1) = W_n(:,:,1)*(tmp_SDNHWP(:,:,prep,1)-tmp_dr_Hpsi_minus(:,:,prep,1)) &
                                        - W_n(:,:,3)*(tmp_SDNHWM(:,:,prep,1)+tmp_dz_Hpsi_minus(:,:,prep,1))

            tmp_WHpsi_plus(:,:,prep,2)  = W_p(:,:,3)*(tmp_SDNHWP(:,:,prep,2)-tmp_dz_Hpsi_plus(:,:,prep,2)) &
                                        + W_p(:,:,1)*(tmp_SDNHWM(:,:,prep,2)-tmp_dr_Hpsi_plus(:,:,prep,2))
            tmp_WHpsi_minus(:,:,prep,2) = W_p(:,:,1)*(tmp_SDNHWP(:,:,prep,2)-tmp_dr_Hpsi_minus(:,:,prep,2)) &
                                        - W_p(:,:,3)*(tmp_SDNHWM(:,:,prep,2)+tmp_dz_Hpsi_minus(:,:,prep,2))
        end do

        prep_Hpsi2_plus(:,:,:,:)  = tmp_UHpsi_plus(:,:,:,:)  - HBAR2_over_2m*tmp_lap_Hpsi_plus(:,:,:,:)
        prep_Hpsi2_minus(:,:,:,:) = tmp_UHpsi_minus(:,:,:,:) - HBAR2_over_2m*tmp_lap_Hpsi_minus(:,:,:,:)
        ! Bpsiには，密度依存部分のみを含めている．
        if(is_use_B)then
            prep_Hpsi2_plus(:,:,:,:)  = prep_Hpsi2_plus(:,:,:,:)  -  tmp_BHpsi_plus(:,:,:,:)
            prep_Hpsi2_minus(:,:,:,:) = prep_Hpsi2_minus(:,:,:,:) - tmp_BHpsi_minus(:,:,:,:)
        end if
        ! BHpsiはマイナスでたさないといけない
        if(is_use_W)then
            prep_Hpsi2_plus(:,:,:,:)  = prep_Hpsi2_plus(:,:,:,:) +  tmp_WHpsi_plus(:,:,:,:)
            prep_Hpsi2_minus(:,:,:,:) = prep_Hpsi2_minus(:,:,:,:) + tmp_WHpsi_minus(:,:,:,:)
        end if

    end subroutine calc_Hpsi2

    subroutine calc_one_particle_energy
        use,intrinsic :: iso_fortran_env,only:int32
        use constants_and_parameters,only:prep_number,HBAR2_over_2m
        use math_integrate,only: volume_integrate
        use global_variables,only: prep_wf_plus,prep_wf_minus &
                                  ,prep_Upsi_plus,prep_Upsi_minus,prep_Bpsi_plus,prep_Bpsi_minus,prep_Wpsi_plus, prep_Wpsi_minus &
                                  ,prep_E => prep_one_particle_E, prep_E_U => prep_one_E_U, prep_E_Kin => prep_one_E_Kin &
                                  ,prep_E_B => prep_one_E_B, prep_E_W => prep_one_E_W &
                                  ,prep_lap_wf_plus,prep_lap_wf_minus, prep_Hpsi_plus, prep_Hpsi_minus, prep_one_particle_E2
        implicit none
        integer(int32) :: prep,np

        ! 引数は  prep_E_U(prep,np,ud)である．
        do np=1,2
            do prep=1,prep_number
                prep_E_U(prep,np,1)   = volume_integrate(f=prep_wf_plus(:,:,prep,np)*prep_Upsi_plus(:,:,prep,np)) 
                prep_E_U(prep,np,2)   = volume_integrate(f=prep_wf_minus(:,:,prep,np)*prep_Upsi_minus(:,:,prep,np))
                prep_E_Kin(prep,np,1) = volume_integrate(f=prep_wf_plus(:,:,prep,np)&
                                                         *(-1d0*HBAR2_over_2m*prep_lap_wf_plus(:,:,prep,np)))
                prep_E_Kin(prep,np,2) = volume_integrate(f=prep_wf_minus(:,:,prep,np)&
                                                         *(-1d0*HBAR2_over_2m*prep_lap_wf_minus(:,:,prep,np)))
                prep_E_B(prep,np,1)   = volume_integrate(f=prep_wf_plus(:,:,prep,np)*prep_Bpsi_plus(:,:,prep,np))
                prep_E_B(prep,np,2)   = volume_integrate(f=prep_wf_minus(:,:,prep,np)*prep_Bpsi_minus(:,:,prep,np))
                prep_E_W(prep,np,1)   = volume_integrate(f=prep_wf_plus(:,:,prep,np)*prep_Wpsi_plus(:,:,prep,np))
                prep_E_W(prep,np,2)   = volume_integrate(f=prep_wf_minus(:,:,prep,np)*prep_Wpsi_minus(:,:,prep,np))

                !prep_E(prep,np) = prep_E_U(prep,np,1) + prep_E_U(prep,np,2) + prep_E_Kin(prep,np,1) + prep_E_Kin(prep,np,2) &
                !                 -prep_E_B(prep,np,1) - prep_E_B(prep,np,2) + prep_E_W(prep,np,1) + prep_E_W(prep,np,2)
                prep_E(prep,np) = volume_integrate(f=prep_wf_plus(:,:,prep,np)*prep_Hpsi_plus(:,:,prep,np)) &
                                 +volume_integrate(f=prep_wf_minus(:,:,prep,np)*prep_Hpsi_minus(:,:,prep,np))
            end do
        end do

        do np=1,2
            do prep=1,prep_number
                prep_one_particle_E2(prep,np) = volume_integrate(f=prep_Hpsi_plus(:,:,prep,np)**2) &
                                              +volume_integrate(f=prep_Hpsi_minus(:,:,prep,np)**2)
            end do
        end do

    end subroutine calc_one_particle_energy

    subroutine calc_Hpsi
        use,intrinsic :: iso_fortran_env,only:int32
        use constants_and_parameters,only:HBAR2_over_2m,is_use_B,is_use_W,lam => lagrange_multiplier
        use global_variables,only: prep_Hpsi_plus, prep_Hpsi_minus,prep_lap_wf_plus,prep_lap_wf_minus &
                                  ,prep_Upsi_plus,prep_Upsi_minus,prep_Bpsi_plus,prep_Bpsi_minus,prep_Wpsi_plus, prep_Wpsi_minus &
                                  ,Nr,Nz,z_vec,dz, prep_wf_plus, prep_wf_minus
        implicit none
        integer(int32) :: j
        call calc_Upsi
        call calc_Bpsi
        call calc_Wpsi
        
        prep_Hpsi_plus(:,:,:,:)  = prep_Upsi_plus(:,:,:,:)  - HBAR2_over_2m*prep_lap_wf_plus(:,:,:,:)
        prep_Hpsi_minus(:,:,:,:) = prep_Upsi_minus(:,:,:,:) - HBAR2_over_2m*prep_lap_wf_minus(:,:,:,:)
        ! Bpsiには，密度依存部分のみを含めている．
        if(is_use_B)then
            prep_Hpsi_plus(:,:,:,:)  = prep_Hpsi_plus(:,:,:,:)  - prep_Bpsi_plus(:,:,:,:)
            prep_Hpsi_minus(:,:,:,:) = prep_Hpsi_minus(:,:,:,:) - prep_Bpsi_minus(:,:,:,:)
        end if
        ! BHpsiはマイナスでたさないといけない
        if(is_use_W)then
            prep_Hpsi_plus(:,:,:,:)  = prep_Hpsi_plus(:,:,:,:) + prep_Wpsi_plus(:,:,:,:)
            prep_Hpsi_minus(:,:,:,:) = prep_Hpsi_minus(:,:,:,:) + prep_Wpsi_minus(:,:,:,:)
        end if

        ! 未定定数を引く
        !do j=1,Nz
        !    prep_Hpsi_plus(:,j,:,:) = prep_Hpsi_plus(:,j,:,:) &
        !                            - lam*(z_vec(j)**2*prep_wf_plus(:,j,:,:) &
        !                            - 2d0*z_vec(j)*(Nz/2d0*dz)*prep_wf_plus(:,j,:,:) &
        !                            + (Nz/2d0*dz)**2*prep_wf_plus(:,j,:,:))
        !    prep_Hpsi_minus(:,j,:,:) = prep_Hpsi_minus(:,j,:,:) &
        !                            - lam*(z_vec(j)**2*prep_wf_minus(:,j,:,:) &
        !                            - 2d0*z_vec(j)*(Nz/2d0*dz)*prep_wf_minus(:,j,:,:) &
        !                            + (Nz/2d0*dz)**2*prep_wf_minus(:,j,:,:))
        !end do
    end subroutine calc_Hpsi



    subroutine calc_IMAG_Evolution(np_int)
        use constants_and_parameters, cdt => imaginary_time_step
        use global_variables,only: prep_wf_plus, prep_wf_minus, prep_Hpsi_plus, prep_Hpsi_minus,prep_Hpsi2_plus,prep_Hpsi2_minus
        implicit none
        integer(int32),intent(in) :: np_int

        prep_wf_plus(:,:,:,np_int)  = prep_wf_plus(:,:,:,np_int)  - cdt/HBARC*prep_Hpsi_plus(:,:,:,np_int)
        prep_wf_minus(:,:,:,np_int) = prep_wf_minus(:,:,:,np_int) - cdt/HBARC*prep_Hpsi_minus(:,:,:,np_int)

        if(is_use_second_imag)then
            prep_wf_plus(:,:,:,np_int)  = prep_wf_plus(:,:,:,np_int)  + cdt**2/HBARC**2/2d0*prep_Hpsi2_plus(:,:,:,np_int)
            prep_wf_minus(:,:,:,np_int) = prep_wf_minus(:,:,:,np_int) + cdt**2/HBARC**2/2d0*prep_Hpsi2_minus(:,:,:,np_int)
        end if
    end subroutine calc_IMAG_Evolution
    
    subroutine Schmidt_orthogonalization(np_int) ! 陽子と中性子別々で正規直交化
        use constants_and_parameters
        use global_variables,only:prep_wf_plus,prep_wf_minus,idx=>prep_idx_array
        implicit none
        integer(int32),intent(in) :: np_int
        integer(int32) :: k,prep
        integer(int32) :: tar,sub ! 直交化の対象と引き算をするインデックス
        real(real64)   :: ovlp,norm
        real(real64)   :: coef

        do prep=1,prep_number
            tar = idx(prep,np_int) ! 直交化の対象インデックス. idxはエネルギーの昇順にソート済み
            ! 直交化パート
            do k=1,prep-1
                sub = idx(k,np_int) ! 引き算をするインデックス
                ovlp = prep_wf_overlap(tar,sub,np_int)
                norm = prep_wf_overlap(sub,sub,np_int)
                coef = ovlp/norm
                prep_wf_plus(:,:,tar,np_int) = prep_wf_plus(:,:,tar,np_int) - coef*prep_wf_plus(:,:,sub,np_int)
                prep_wf_minus(:,:,tar,np_int) = prep_wf_minus(:,:,tar,np_int) - coef*prep_wf_minus(:,:,sub,np_int)
            end do
            ! 規格化パート
            norm = prep_wf_overlap(tar,tar,np_int)
            prep_wf_plus(:,:,tar,np_int) = prep_wf_plus(:,:,tar,np_int)/sqrt(norm)
            prep_wf_minus(:,:,tar,np_int) = prep_wf_minus(:,:,tar,np_int)/sqrt(norm)
        end do
    end subroutine Schmidt_orthogonalization


    function prep_wf_overlap(idx_1,idx_2,np_int) result(ovlp)
        use,intrinsic :: iso_fortran_env,only:int32,real64
        use math_integrate,only: volume_integrate
        use global_variables,only:prep_wf_plus,prep_wf_minus
        implicit none
        integer(int32),intent(in) :: idx_1
        integer(int32),intent(in) :: idx_2
        integer(int32),intent(in) :: np_int
        real(real64)              :: ovlp
        real(real64)              :: term_plus,term_minus
        
        term_plus = volume_integrate(f=prep_wf_plus(:,:,idx_1,np_int)*prep_wf_plus(:,:,idx_2,np_int))
        term_minus= volume_integrate(f=prep_wf_minus(:,:,idx_1,np_int)*prep_wf_minus(:,:,idx_2,np_int))
        
        ovlp = term_plus + term_minus
    end function prep_wf_overlap
    
    
    
    subroutine calc_Upsi
        use,intrinsic :: iso_fortran_env,only:int32
        use constants_and_parameters,only:prep_number
        use global_variables,only:prep_wf_plus,prep_wf_minus,U_n,U_p,prep_Upsi_plus,prep_Upsi_minus
        implicit none
        integer(int32) :: prep

        do prep=1,prep_number
            prep_Upsi_plus(:,:,prep,1)  = U_n(:,:)*prep_wf_plus(:,:,prep,1)
            prep_Upsi_minus(:,:,prep,1) = U_n(:,:)*prep_wf_minus(:,:,prep,1)
            prep_Upsi_plus(:,:,prep,2)  = U_p(:,:)*prep_wf_plus(:,:,prep,2)
            prep_Upsi_minus(:,:,prep,2) = U_p(:,:)*prep_wf_minus(:,:,prep,2)
            !print*, maxval(U_n(:,:))
            !stop
        end do
    end subroutine calc_Upsi
    
    subroutine calc_Wpsi
        use,intrinsic :: iso_fortran_env,only:int32,real64
        use constants_and_parameters,only: zI => imaginary_unit, prep_number
        use global_variables,only:W_n, W_p &
                                 ,prep_dz_wf_minus, prep_dz_wf_plus, prep_dr_wf_minus, prep_dr_wf_plus&
                                 ,prep_SDNWP, prep_SDNWM, prep_Wpsi_plus, prep_Wpsi_minus
        implicit none
        integer(int32) :: prep

        do prep=1,prep_number
            prep_Wpsi_plus(:,:,prep,1)  = W_n(:,:,3)*(prep_SDNWP(:,:,prep,1)-prep_dz_wf_plus(:,:,prep,1)) &
                                        + W_n(:,:,1)*(prep_SDNWM(:,:,prep,1)-prep_dr_wf_plus(:,:,prep,1))
            prep_Wpsi_minus(:,:,prep,1) = W_n(:,:,1)*(prep_SDNWP(:,:,prep,1)-prep_dr_wf_minus(:,:,prep,1)) &
                                        - W_n(:,:,3)*(prep_SDNWM(:,:,prep,1)+prep_dz_wf_minus(:,:,prep,1))

            prep_Wpsi_plus(:,:,prep,2)  = W_p(:,:,3)*(prep_SDNWP(:,:,prep,2)-prep_dz_wf_plus(:,:,prep,2)) &
                                        + W_p(:,:,1)*(prep_SDNWM(:,:,prep,2)-prep_dr_wf_plus(:,:,prep,2))
            prep_Wpsi_minus(:,:,prep,2) = W_p(:,:,1)*(prep_SDNWP(:,:,prep,2)-prep_dr_wf_minus(:,:,prep,2)) &
                                        - W_p(:,:,3)*(prep_SDNWM(:,:,prep,2)+prep_dz_wf_minus(:,:,prep,2))
        end do
    end subroutine calc_Wpsi
    
    subroutine calc_Bpsi
        use,intrinsic :: iso_fortran_env,only:int32,real64
        use constants_and_parameters,    only: Nr,Nz,prep_number
        use math_derivation,             only: math_diff_r_9point, math_diff_z_9point
        use global_variables,            only: B_n, B_p &
                                             , prep_dr_wf_plus, prep_dr_wf_minus, prep_dz_wf_plus, prep_dz_wf_minus &
                                             , prep_lap_wf_plus,prep_lap_wf_minus,prep_Bpsi_plus,prep_Bpsi_minus
        implicit none
 
        real(real64),allocatable  :: dr_B(:,:), dz_B(:,:)
        integer(int32)            :: prep,np_int
        allocate(dr_B(Nr,Nz))
        allocate(dz_B(Nr,Nz))
        

        dr_B(:,:) = math_diff_r_9point(B_n(:,:))
        dz_B(:,:) = math_diff_z_9point(B_n(:,:))
        do prep=1,prep_number
            prep_Bpsi_plus(:,:,prep,1) = B_n(:,:)*prep_lap_wf_plus(:,:,prep,1) &
                                        + dr_B(:,:)*prep_dr_wf_plus(:,:,prep,1) &
                                        + dz_B(:,:)*prep_dz_wf_plus(:,:,prep,1)

            prep_Bpsi_minus(:,:,prep,1) = B_n(:,:)*prep_lap_wf_minus(:,:,prep,1) &
                                        + dr_B(:,:)*prep_dr_wf_minus(:,:,prep,1) &
                                        + dz_B(:,:)*prep_dz_wf_minus(:,:,prep,1)
        end do

        dr_B(:,:) = math_diff_r_9point(B_p(:,:))
        dz_B(:,:) = math_diff_z_9point(B_p(:,:))
        do prep=1,prep_number
            prep_Bpsi_plus(:,:,prep,2) =  B_p(:,:)*prep_lap_wf_plus(:,:,prep,2) &
                                        + dr_B(:,:)*prep_dr_wf_plus(:,:,prep,2) &
                                        + dz_B(:,:)*prep_dz_wf_plus(:,:,prep,2)

            prep_Bpsi_minus(:,:,prep,2) = B_p(:,:)*prep_lap_wf_minus(:,:,prep,2) &
                                        + dr_B(:,:)*prep_dr_wf_minus(:,:,prep,2) &
                                        + dz_B(:,:)*prep_dz_wf_minus(:,:,prep,2)
        end do
    end subroutine calc_Bpsi
    
    subroutine calc_density_and_lap
        use,intrinsic :: ieee_arithmetic
        use,intrinsic :: iso_fortran_env,only:int32,real64
        use constants_and_parameters,only: Nr,Nz,num_n,num_p,is_use_dd,z_center,str_z
        use global_variables,only:rho_n,rho_p,wf_plus,wf_minus, lap_rho_n,lap_rho_p, r_vec&
                                  , dr_rho_n, dr_rho_p, dz_rho_n, dz_rho_p,ddr_rho_n,ddr_rho_p,ddz_rho_n,ddz_rho_p,z_vec
        use math_integrate,only: volume_integrate
        use math_derivation,only: math_diff_r_9point,math_diff_z_9point,math_diff_2r_9point,math_diff_2z_9point
        implicit none
        integer(int32) :: k,i,j
        real(real64),allocatable   :: tmp(:,:)
        allocate(tmp(Nr,Nz))

        rho_n(:,:) = 0.0d0
        rho_p(:,:) = 0.0d0
        
        do k=1,num_n
            rho_n(:,:) = rho_n(:,:) + abs(wf_plus(:,:,k))**2 + abs(wf_minus(:,:,k))**2
        end do
        
        do k=1,num_p
            rho_p(:,:) = rho_p(:,:) + abs(wf_plus(:,:,num_n+k))**2 + abs(wf_minus(:,:,num_n+k))**2
        end do

        ! Correction center of the nucleus
        do j=1,Nz
            do i=1,Nr
                tmp(i,j) = (rho_n(i,j) + rho_p(i,j))*z_vec(j)
            end do
        end do
        z_center = volume_integrate(f=tmp(:,:))/volume_integrate(f=rho_n(:,:) + rho_p(:,:))
        !print*, 'z_center = ', z_center
        write(str_z,'(f0.2)') z_center


        dr_rho_n(:,:) = math_diff_r_9point(f=rho_n(:,:))
        dr_rho_p(:,:) = math_diff_r_9point(f=rho_p(:,:))
        dz_rho_n(:,:) = math_diff_z_9point(f=rho_n(:,:))
        dz_rho_p(:,:) = math_diff_z_9point(f=rho_p(:,:))
        if(is_use_dd)then
            ddr_rho_n(:,:) = math_diff_2r_9point(f=rho_n(:,:))
            ddr_rho_p(:,:) = math_diff_2r_9point(f=rho_p(:,:))
            ddz_rho_n(:,:) = math_diff_2z_9point(f=rho_n(:,:))
            ddz_rho_p(:,:) = math_diff_2z_9point(f=rho_p(:,:))
        else
            ddr_rho_n(:,:) = math_diff_r_9point(f=dr_rho_n(:,:))
            ddr_rho_p(:,:) = math_diff_r_9point(f=dr_rho_p(:,:))
            ddz_rho_n(:,:) = math_diff_z_9point(f=dz_rho_n(:,:))
            ddz_rho_p(:,:) = math_diff_z_9point(f=dz_rho_p(:,:))
        end if
        do j=1,Nz
            do i=1,Nr
                lap_rho_n(i,j) = ddr_rho_n(i,j) + dr_rho_n(i,j)/r_vec(i) + ddz_rho_n(i,j)
                lap_rho_p(i,j) = ddr_rho_p(i,j) + dr_rho_p(i,j)/r_vec(i) + ddz_rho_p(i,j)
            end do
        end do
        deallocate(tmp)
    end subroutine calc_density_and_lap
    

    subroutine calc_spin_orbit_density_and_div ! J13=J31=0
        use,intrinsic :: iso_fortran_env,only:int32,real64
        use,intrinsic :: ieee_arithmetic
        use constants_and_parameters
        use math_derivation,only: math_diff_r_9point,math_diff_z_9point
        use global_variables, only: j_n=>spin_orbit_dens_n, j_p=>spin_orbit_dens_p, &
                              div_j_n=>div_spin_orbit_dens_n, div_j_p=>div_spin_orbit_dens_p, &
                              wf_plus, wf_minus, m=>magnetic_q_num, dz_wf_plus, dz_wf_minus, dr_wf_plus, dr_wf_minus, &
                              r_vec, tau_n => kinetic_dens_n, tau_p => kinetic_dens_p &
                                ,SDNWP, SDNWM, r_vec
        implicit none
        integer(int32) :: i,j,k,kp
        real(real64),allocatable :: J12(:,:), J21(:,:), J23(:,:), J32(:,:)
        real(real64),allocatable :: dr_j(:,:), dz_j(:,:)
        allocate(J12(Nr,Nz))
        allocate(J21(Nr,Nz))
        allocate(J23(Nr,Nz))
        allocate(J32(Nr,Nz))
        allocate(dr_j(Nr,Nz))
        allocate(dz_j(Nr,Nz))


        J12=0d0;J21=0d0;J23=0d0;J32=0d0
        do k=1,num_n
            do j=1,Nz
                do i=1,Nr
                    J12(i,j) = J12(i,j) + wf_minus(i,j,k)*dr_wf_plus(i,j,k) - wf_plus(i,j,k)*dr_wf_minus(i,j,k)
                    J21(i,j) = J21(i,j) + dble(2*m(k)+1)/r_vec(i)*wf_plus(i,j,k)*wf_minus(i,j,k)
                    J23(i,j) = J23(i,j) + (dble(m(k))*wf_plus(i,j,k)**2 + dble(m(k)+1)*wf_minus(i,j,k)**2)/r_vec(i)
                    J32(i,j) = J32(i,j) + wf_minus(i,j,k)*dz_wf_plus(i,j,k) - wf_plus(i,j,k)*dz_wf_minus(i,j,k)
                end do
            end do
        end do
        j_n(:,:,1) = J23(:,:) - J32(:,:)
        j_n(:,:,2) = 0d0
        j_n(:,:,3) = J12(:,:) - J21(:,:)

        if(not_divJ)then
            div_j_n(:,:) = 0d0
        else
            !do k=1,num_n
            !    div_j_n(:,:) = div_j_n(:,:) + SDNWP(:,:,k)**2 + SDNWM(:,:,k)**2
            !end do
            !div_j_n(:,:) = div_j_n(:,:) - tau_n(:,:)
            dr_j(:,:) = math_diff_r_9point(f=j_n(:,:,1))
            dz_j(:,:) = math_diff_z_9point(f=j_n(:,:,3))
            do j=1,Nz
                do i=1,Nr
                div_j_n(i,j) = j_n(i,j,1)/r_vec(i) + dr_j(i,j) + dz_j(i,j)
                end do
            end do
        end if

        J12=0d0;J21=0d0;J23=0d0;J32=0d0
        do k=1,num_p
            kp=num_n+k
            do j=1,Nz
                do i=1,Nr
                    J12(i,j) = J12(i,j) + wf_minus(i,j,kp)*dr_wf_plus(i,j,kp) - wf_plus(i,j,kp)*dr_wf_minus(i,j,kp)
                    J21(i,j) = J21(i,j) + dble(2*m(kp)+1)/r_vec(i)*wf_plus(i,j,kp)*wf_minus(i,j,kp)
                    J23(i,j) = J23(i,j) + (dble(m(kp))*wf_plus(i,j,kp)**2 + dble(m(kp)+1)*wf_minus(i,j,kp)**2)/r_vec(i)
                    J32(i,j) = J32(i,j) + wf_minus(i,j,kp)*dz_wf_plus(i,j,kp) - wf_plus(i,j,kp)*dz_wf_minus(i,j,kp)
                end do
            end do
        end do
        j_p(:,:,1) = J23(:,:) - J32(:,:)
        j_p(:,:,2) = 0d0
        j_p(:,:,3) = J12(:,:) - J21(:,:)

        if(not_divJ)then
            div_j_p(:,:) = 0d0
        else
            !do k=1,num_p
            !    kp = num_n + k
            !    div_j_p(:,:) = div_j_p(:,:) + SDNWP(:,:,kp)**2 + SDNWM(:,:,kp)**2
            !end do
            !div_j_p(:,:) = div_j_p(:,:) - tau_p(:,:)
            dr_j(:,:) = math_diff_r_9point(f=j_p(:,:,1))
            dz_j(:,:) = math_diff_z_9point(f=j_p(:,:,3))
            do j=1,Nz
                do i=1,Nr
                div_j_p(i,j) = j_p(i,j,1)/r_vec(i) + dr_j(i,j) + dz_j(i,j)
                end do
            end do
        end if
        deallocate(J12)
        deallocate(J21)
        deallocate(J23)
        deallocate(J32)
    end subroutine calc_spin_orbit_density_and_div

    subroutine calc_kinetic_density
        use,intrinsic :: iso_fortran_env,only:int32,real64
        use global_variables,only: tau_n => kinetic_dens_n, tau_p => kinetic_dens_p, m => magnetic_q_num &
                                  ,rho_n,rho_p, dr_wf_plus, dr_wf_minus, dz_wf_plus, dz_wf_minus, wf_plus, wf_minus, r_vec &
                                  ,thomas_fermi_n,thomas_fermi_p
        use constants_and_parameters,only:Nr,Nz,num_n,num_p,pi
        implicit none
        
        integer(int32) :: i,j,k
        real(real64), allocatable :: tau_per_nuc(:,:) ! [fm^-3] (r,z)
        allocate(tau_per_nuc(Nr,Nz))
        tau_n(:,:) = 0.0d0
        tau_p(:,:) = 0.0d0
        do k=1,num_n
            ! 32b式を参照して、各核子のtauを計算している。各項が正しいか検算
            do j=1,Nz
                do i=1,Nr
                    tau_per_nuc(i,j) = abs(dr_wf_plus(i,j,k))**2  + abs(dz_wf_plus(i,j,k))**2 &
                                     + abs(dr_wf_minus(i,j,k))**2 + abs(dz_wf_minus(i,j,k))**2 &
                                     + (dble(m(k))/r_vec(i))**2*abs(wf_plus(i,j,k))**2 &
                                     + (dble(m(k)+1)/r_vec(i))**2*abs(wf_minus(i,j,k))**2
                end do
            end do
            tau_n(:,:) = tau_n(:,:) + tau_per_nuc(:,:) ! 1つの核子のtauを足し合わせている。
        end do
        
        do k=1,num_p
            do j=1,Nz
                do i=1,Nr
                    tau_per_nuc(i,j) = abs(dr_wf_plus(i,j,num_n+k))**2 + abs(dz_wf_plus(i,j,num_n+k))**2 &
                                     + abs(dr_wf_minus(i,j,num_n+k))**2 + abs(dz_wf_minus(i,j,num_n+k))**2 &
                                     + (dble(m(num_n+k))/r_vec(i))**2*abs(wf_plus(i,j,num_n+k))**2 &
                                     + (dble(m(num_n+k)+1)/r_vec(i))**2*abs(wf_minus(i,j,num_n+k))**2
                end do
            end do
            tau_p(:,:) = tau_p(:,:) + tau_per_nuc(:,:)
        end do


        thomas_fermi_n(:,:) = 3d0/5d0*(3d0*pi**2)**(2d0/3d0)*rho_n(:,:)**(5d0/3d0) ! ちゃんとした密度だったら計算しているのと同じような形になるはず
        thomas_fermi_p(:,:) = 3d0/5d0*(3d0*pi**2)**(2d0/3d0)*rho_p(:,:)**(5d0/3d0)

    end subroutine calc_kinetic_density
    
    subroutine alloc_fields
        use,intrinsic :: iso_fortran_env
        use :: global_variables
        use :: constants_and_parameters
        implicit none
        integer(int32) :: nuc  
        integer(int32) :: prep 
        nuc = num_n + num_p
        prep = prep_number
!-- Density and potential (caluculated by Occupied Orbitals) --
        allocate(rho_n(Nr,Nz)                 ) ! [fm^-3]  (r,z)
        allocate(rho_p(Nr,Nz)                 ) ! [fm^-3]  (r,z)
        allocate(dr_rho_n(Nr,Nz)              ) ! [fm^-4]  (r,z)
        allocate(dr_rho_p(Nr,Nz)              ) ! [fm^-4]  (r,z)
        allocate(dz_rho_n(Nr,Nz)              ) ! [fm^-4]  (r,z)
        allocate(dz_rho_p(Nr,Nz)              ) ! [fm^-4]  (r,z)
        allocate(ddr_rho_n(Nr,Nz)             ) ! [fm^-5]  (r,z)
        allocate(ddr_rho_p(Nr,Nz)             ) ! [fm^-5]  (r,z)
        allocate(ddz_rho_n(Nr,Nz)             ) ! [fm^-5]  (r,z)
        allocate(ddz_rho_p(Nr,Nz)             ) ! [fm^-5]  (r,z)
        allocate(lap_rho_n(Nr,Nz)             ) ! [fm^-5]  (r,z)
        allocate(lap_rho_p(Nr,Nz)             ) ! [fm^-5]  (r,z)
        allocate(kinetic_dens_n(Nr,Nz)        ) ! [fm^-3]  (r,z)
        allocate(kinetic_dens_p(Nr,Nz)        ) ! [fm^-3]  (r,z)
        allocate(thomas_fermi_n(Nr,Nz)        ) ! [fm^-3]  (r,z)
        allocate(thomas_fermi_p(Nr,Nz)        ) ! [fm^-3]  (r,z)
        allocate(spin_orbit_dens_n(Nr,Nz,3)   ) ! [fm^-3]  (r,z,vec component)
        allocate(spin_orbit_dens_p(Nr,Nz,3)   ) ! [fm^-3]  (r,z,vec component)
        allocate(div_spin_orbit_dens_n(Nr,Nz) ) ! [fm^-4]  (r,z)
        allocate(div_spin_orbit_dens_p(Nr,Nz) ) ! [fm^-4]  (r,z)
        allocate(U_n(Nr,Nz)                   ) ! [MeV] (r,z)
        allocate(U_p(Nr,Nz)                   ) ! [MeV] (r,z)
        allocate(B_n(Nr,Nz)                   ) ! [MeV] (r,z)
        allocate(B_p(Nr,Nz)                   ) ! [MeV] (r,z)
        allocate(W_n(Nr,Nz,3)                 ) ! [MeV] (r,z,component)
        allocate(W_p(Nr,Nz,3)                 ) ! [MeV] (r,z,component)
        allocate(direct_Coulomb_pot(Nr,Nz),source=0d0) ! [MeV] (r,z)
        allocate(exchange_Coulomb_pot(Nr,Nz),source=0d0  ) ! [MeV] (r,z

! -- All Orbitals(Occupied state / All state : prefix = prep) --
        allocate(prep_idx_array(prep,2)              ) ! index of nucleon
        allocate(magnetic_q_num(nuc)          ) ! magnetic quantum number
        allocate(prep_magnetic_q_num(prep,2)    ) ! magnetic quantum number (prepared)
        allocate(init_rzm_to_write(prep)    ) ! (idx,nlm,p/n) (prepared)
        allocate(prep_one_particle_E(prep,2)    ) ! [MeV] 1粒子エネルギー (prepared)
        allocate(prep_one_particle_E2(prep,2))
        allocate(prep_one_E_Kin(prep,2,2))
        allocate(prep_one_E_U(prep,2,2))
        allocate(prep_one_E_B(prep,2,2))
        allocate(prep_one_E_W(prep,2,2))
        allocate(wf_plus(Nr,Nz,nuc)           ) ! [fm^-3/2]  (r,z,nucleon)
        allocate(wf_minus(Nr,Nz,nuc)          ) ! [fm^-3/2]  (r,z,nucleon)
        allocate(prep_wf_plus(Nr,Nz,prep ,2)     ) ! [fm^-3/2]  (r,z,prepsize)
        allocate(prep_wf_minus(Nr,Nz,prep,2)    ) ! [fm^-3/2]  (r,z,prepsize)
        !----- derivative of wave function -----
        allocate(dr_wf_plus(Nr,Nz,nuc)        )
        allocate(dz_wf_plus(Nr,Nz,nuc)        )
        allocate(ddr_wf_plus(Nr,Nz,nuc)       )
        allocate(ddz_wf_plus(Nr,Nz,nuc)       )
        allocate(lap_wf_plus(Nr,Nz,nuc)       )
        allocate(prep_dr_wf_plus(Nr,Nz,prep ,2)  )
        allocate(prep_dz_wf_plus(Nr,Nz,prep ,2)  )
        allocate(prep_ddr_wf_plus(Nr,Nz,prep,2) )
        allocate(prep_ddz_wf_plus(Nr,Nz,prep,2) )
        allocate(prep_lap_wf_plus(Nr,Nz,prep,2) )
        allocate(dr_wf_minus(Nr,Nz,nuc)       )
        allocate(dz_wf_minus(Nr,Nz,nuc)       )
        allocate(ddr_wf_minus(Nr,Nz,nuc)      )
        allocate(ddz_wf_minus(Nr,Nz,nuc)      )
        allocate(lap_wf_minus(Nr,Nz,nuc)      )
        allocate(prep_dr_wf_minus(Nr,Nz,prep ,2) )
        allocate(prep_dz_wf_minus(Nr,Nz,prep ,2) )
        allocate(prep_ddr_wf_minus(Nr,Nz,prep,2))
        allocate(prep_ddz_wf_minus(Nr,Nz,prep,2))
        allocate(prep_lap_wf_minus(Nr,Nz,prep,2))
        allocate(prep_Upsi_plus(Nr,Nz,prep,2)      )
        allocate(prep_Upsi_minus(Nr,Nz,prep,2)     )
        allocate(prep_Wpsi_plus(Nr,Nz,prep,2)      )
        allocate(prep_Wpsi_minus(Nr,Nz,prep,2)     )
        allocate(prep_Bpsi_plus(Nr,Nz,prep,2)      )
        allocate(prep_Bpsi_minus(Nr,Nz,prep,2)     )
        !----- use for calculation -----
        ! SDNWP/M = sigma dot nabla plus/minus
        allocate(SDNWP(Nr,Nz,nuc)             )  ! [fm^-3/2]  (r,z,nucleon)
        allocate(SDNWM(Nr,Nz,nuc)             )  ! [fm^-3/2]  (r,z,nucleon)
        allocate(prep_SDNWP(Nr,Nz,prep,2)       )  ! [fm^-3/2]  (r,z,nucleon)
        allocate(prep_SDNWM(Nr,Nz,prep,2)       )  ! [fm^-3/2]  (r,z,nucleon)
        ! Hpsi
        allocate(prep_Hpsi2_plus(Nr,Nz,prep,2)         ) ! [MeV fm^-3/2]  (r,z,nucleon)
        allocate(prep_Hpsi2_minus(Nr,Nz,prep,2)        ) ! [MeV fm^-3/2]  (r,z,nucleon)
        allocate(prep_Hpsi_plus(Nr,Nz,prep ,2)   ) ! [MeV fm^-3/2]  (r,z,nucleon)
        allocate(prep_Hpsi_minus(Nr,Nz,prep,2)  ) ! [MeV fm^-3/2]  (r,z,nucleon)
    
        !-- Grid variables --
        allocate(r_vec(Nr)                    ) ! [fm] r position
        allocate(z_vec(Nz)                    ) ! [fm] z position
    end subroutine alloc_fields

    
    subroutine dealloc_fields
        use :: global_variables
        implicit none

!-- Density and potential (caluculated by Occupied Orbitals) --
        deallocate(rho_n                 ) ! [fm^-3]  (r,z)
        deallocate(rho_p                 ) ! [fm^-3]  (r,z)
        deallocate(dr_rho_n              ) ! [fm^-4]  (r,z)
        deallocate(dr_rho_p              ) ! [fm^-4]  (r,z)
        deallocate(dz_rho_n              ) ! [fm^-4]  (r,z)
        deallocate(dz_rho_p              ) ! [fm^-4]  (r,z)
        deallocate(ddr_rho_n             ) ! [fm^-5]  (r,z)
        deallocate(ddr_rho_p             ) ! [fm^-5]  (r,z)
        deallocate(ddz_rho_n             ) ! [fm^-5]  (r,z)
        deallocate(ddz_rho_p             ) ! [fm^-5]  (r,z)
        deallocate(lap_rho_n             ) ! [fm^-5]  (r,z)
        deallocate(lap_rho_p             ) ! [fm^-5]  (r,z)
        deallocate(kinetic_dens_n        ) ! [fm^-3]  (r,z)
        deallocate(kinetic_dens_p        ) ! [fm^-3]  (r,z)
        deallocate(thomas_fermi_n)
        deallocate(thomas_fermi_p)
        deallocate(spin_orbit_dens_n   ) ! [fm^-3]  (r,z,vec component)
        deallocate(spin_orbit_dens_p   ) ! [fm^-3]  (r,z,vec component)
        deallocate(div_spin_orbit_dens_n ) ! [fm^-4]  (r,z)
        deallocate(div_spin_orbit_dens_p ) ! [fm^-4]  (r,z)
        deallocate(U_n                   ) ! [MeV] (r,z)
        deallocate(U_p                  ) ! [MeV] (r,z)
        deallocate(B_n                  ) ! [MeV] (r,z)
        deallocate(B_p                   ) ! [MeV] (r,z)
        deallocate(W_n                 ) ! [MeV] (r,z,component)
        deallocate(W_p                 ) ! [MeV] (r,z,component)
        deallocate(direct_Coulomb_pot    ) ! [MeV] (r,z)
        deallocate(exchange_Coulomb_pot  ) ! [MeV] (r,z)

! -- All Orbitals(Occupied state / All state : prefix = prep) --
        deallocate(prep_idx_array              ) ! index of nucleon
        deallocate(magnetic_q_num          ) ! magnetic quantum number
        deallocate(prep_magnetic_q_num    ) ! magnetic quantum number (prepared)
        deallocate(init_rzm_to_write    ) ! (idx,nlm,p/n) (prepared)
        deallocate(prep_one_particle_E    ) ! [MeV] 1粒子エネルギー (prepared)
        deallocate(prep_one_particle_E2)
        deallocate(prep_one_E_Kin)
        deallocate(prep_one_E_U)
        deallocate(prep_one_E_B)
        deallocate(prep_one_E_W)
        deallocate(wf_plus          ) ! [fm^-3/2]  (r,z,nucleon)
        deallocate(wf_minus          ) ! [f^-3/2]  (r,z,nucleon)
        deallocate(prep_wf_plus    ) ! [fm^-3/2]  (r,z,prepsize)
        deallocate(prep_wf_minus    ) ! [fm^-3/2]  (r,z,prepsize)
        !----- derivative of wave function -----
        deallocate(dr_wf_plus     )
        deallocate(dz_wf_plus     )
        deallocate(ddr_wf_plus     )
        deallocate(ddz_wf_plus     )
        deallocate(lap_wf_plus     )
        deallocate(prep_dr_wf_plus )
        deallocate(prep_dz_wf_plus )
        deallocate(prep_ddr_wf_plus )
        deallocate(prep_ddz_wf_plus )
        deallocate(prep_lap_wf_plus )
        deallocate(dr_wf_minus      )
        deallocate(dz_wf_minus      )
        deallocate(ddr_wf_minus      )
        deallocate(ddz_wf_minus      )
        deallocate(lap_wf_minus      )
        deallocate(prep_dr_wf_minus)
        deallocate(prep_dz_wf_minus)
        deallocate(prep_ddr_wf_minus)
        deallocate(prep_ddz_wf_minus)
        deallocate(prep_lap_wf_minus)
        deallocate(prep_Upsi_plus      )
        deallocate(prep_Wpsi_plus      )
        deallocate(prep_Bpsi_plus      )
        deallocate(prep_Upsi_minus     )
        deallocate(prep_Wpsi_minus     )
        deallocate(prep_Bpsi_minus     )
        !----- use for calculation -----
        ! SDNWP/M = sigma dot nabla plus/minus
        deallocate(SDNWP      )  ! [fm^-3/2]  (r,z,nucleon)
        deallocate(SDNWM      )  ! [fm^-3/2]  (r,z,nucleon)
        deallocate(prep_SDNWP      )  ! [fm^-3/2]  (r,z,nucleon)
        deallocate(prep_SDNWM      )  ! [fm^-3/2]  (r,z,nucleon)
        ! Hpsi
        deallocate(prep_Hpsi2_plus ) ! [MeV fm^-3/2]  (r,z,nucleon)
        deallocate(prep_Hpsi2_minus ) ! [MeV fm^-3/2]  (r,z,nucleon)
        deallocate(prep_Hpsi_plus ) ! [MeV fm^-3/2]  (r,z,nucleon)
        deallocate(prep_Hpsi_minus ) ! [MeV fm^-3/2]  (r,z,nucleon)
    
        !-- Grid variables --
        deallocate(r_vec          ) ! [fm] r position
        deallocate(z_vec          ) ! [fm] z position
    end subroutine dealloc_fields

end program main