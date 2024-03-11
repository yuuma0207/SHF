module write_to_file
    use,intrinsic :: iso_fortran_env,only:int32,real64
    use :: constants_and_parameters,only:resultdir,foldername
    implicit none
    contains
    
    subroutine outputsettings
        use constants_and_parameters,only:date,time,MMDDHHMMSS
        implicit none
        ! 現在の日付と時刻を取得
        call date_and_time(date, time)

        ! フォルダ名を MMDD/HHMMSS 形式で設定
        MMDDHHMMSS = date(5:8) // "/" // time(1:6) // "/"
        foldername = resultdir // MMDDHHMMSS

        ! フォルダ作成 (システム依存)
        call system("mkdir -p " // foldername)
        call system("mkdir -p " // foldername // "orbital_wavefunction/")
        call system("mkdir -p " // foldername // "densities/")
        call system("mkdir -p " // foldername // "potentials/")


        print *, "フォルダ名は"//MMDDHHMMSS//"です"
        ! make output files(use append mode)
        call make_one_file("sort_information.txt")
        call make_one_file("meanfield_Energy.txt")
        call make_one_file("mag_log.txt")
        call make_one_file("change_orbital.txt")
        call make_one_file("iter_Energy.txt")
    end subroutine outputsettings
    
    subroutine make_one_file(filename)
        integer(int32) :: fi,is=1
        character(*)              :: filename
        character(:), allocatable :: alloc_filename
        allocate(alloc_filename, source=filename)
        open(newunit=fi,file=foldername//filename,action="write",status="replace",iostat=is)
        if(is/=0) then
            print*, "Error: cannot make file."
            stop
        else if(is==0) then
            !print*, "File created: ", filename
        end if
        close(fi)
    end subroutine make_one_file
    
    subroutine write_meanfield_Energy(iter)
        use global_variables,only: Skyrme_Energy, Kinetic_Energy, Coulomb_Energy, CM_Energy, Total_Energy
        implicit none
        integer(int32) :: fi,is
        integer(int32),intent(in) :: iter
        logical,save :: make_index = .true.
        character(:), allocatable :: filename

        allocate(filename, source="meanfield_Energy.txt")
        open(newunit=fi,file=foldername//filename,action="write",position="append",status="old",iostat=is)
        if(is/=0) print*, "Error: cannot open file.",filename
        if(make_index) then
            write(fi,*) "iter Kinetic Skyrme Coulomb CM Total"
            make_index = .false.
        end if
        !write(fi, '(I5, G20.10, G20.10, G20.10, G20.10, G20.10)') &
        write(fi,*) &
                   iter, Kinetic_Energy, Skyrme_Energy, Coulomb_Energy, CM_Energy, Total_Energy
        close(fi)
    end subroutine write_meanfield_Energy

    subroutine write_prep_sort(iter)
        use global_variables, only: one_E => prep_one_particle_E, one_Kin => prep_one_E_Kin, one_U => prep_one_E_U &
                                   ,one_B => prep_one_E_B, one_W => prep_one_E_W &
                                   ,m => prep_magnetic_q_num, idx => prep_idx_array, prep_number, init_rzm_to_write
        use constants_and_parameters, only: num_n,num_p
        implicit none
        integer(int32) :: fi,is,prep
        integer(int32),intent(in) :: iter
        logical,save :: make_index = .true.
        logical,save :: make_index2 = .true.
        character(:), allocatable :: filename,filename2,filename3

        allocate(filename, source="sort_information.txt")
        allocate(filename2, source="change_orbital.txt")
        allocate(filename3, source="mag_log.txt")
        open(newunit=fi,file=foldername//filename,action="write",position="append",status="old",iostat=is)
        if(is/=0) then 
            print*, "Error: cannot open file.",filename
            stop
        end if
        if(make_index) then
            write(fi,*) "one_E_n Kin(u) Kin(d) U(u) U(d) B(u) B(d) W(u) W(d) m_n idx_n  |  &
            one_E_p Kin(u) Kin(d) U(u) U(d) B(u) B(d) W(u) W(d) m_p idx_p"
            make_index = .false.
        else
            do prep=1,2 ! 最初の2軌道だけを出力 もしかしたらidx(prep,1) -> idx(prep) が正しい
                !write(fi,*) &
                !            one_E(idx(prep,1),1), one_Kin(idx(prep,1),1,1), one_Kin(idx(prep,1),1,2) &
                !           ,one_U(idx(prep,1),1,1), one_U(idx(prep,1),1,2) &
                !           ,one_B(idx(prep,1),1,1), one_B(idx(prep,1),1,2) &
                !           ,one_W(idx(prep,1),1,1), one_W(idx(prep,1),1,2) &
                !           ,m(idx(prep,1),1), idx(prep,1), "|"&
                !           ,one_E(idx(prep,2),2), one_Kin(idx(prep,2),2,1), one_Kin(idx(prep,2),2,2) &
                !           ,one_U(idx(prep,2),2,1), one_U(idx(prep,2),2,2) &
                !           ,one_B(idx(prep,2),2,1), one_B(idx(prep,2),2,2) &
                !           ,one_W(idx(prep,2),2,1), one_W(idx(prep,2),2,2) &
                !           ,m(idx(prep,2),2), idx(prep,2)
                write(fi,*) &
                            one_E(prep,1), one_Kin(idx(prep,1),1,1), one_Kin(idx(prep,1),1,2) &
                           ,one_U(idx(prep,1),1,1), one_U(idx(prep,1),1,2) &
                           ,one_B(idx(prep,1),1,1), one_B(idx(prep,1),1,2) &
                           ,one_W(idx(prep,1),1,1), one_W(idx(prep,1),1,2) &
                           ,m(idx(prep,1),1), idx(prep,1), "|"&
                           ,one_E(prep,2), one_Kin(idx(prep,2),2,1), one_Kin(idx(prep,2),2,2) &
                           ,one_U(idx(prep,2),2,1), one_U(idx(prep,2),2,2) &
                           ,one_B(idx(prep,2),2,1), one_B(idx(prep,2),2,2) &
                           ,one_W(idx(prep,2),2,1), one_W(idx(prep,2),2,2) &
                           ,m(idx(prep,2),2), idx(prep,2)
            end do
        end if
        write(fi,*) "---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----"
        close(fi)
        open(newunit=fi,file=foldername//filename2,action="write",position="append",status="old",iostat=is)
        if(is/=0) then 
            print*, "Error: cannot open file.",filename2
            stop
        end if
        if(make_index2) then
            write(fi,*) "one_E_n m_n idx_n one_E_p m_p idx_p"
            make_index2 = .false.
        else
            do prep=1,prep_number
                write(fi,*) one_E(prep,1), m(idx(prep,1),1), idx(prep,1)&
                , one_E(prep,2), m(idx(prep,2),2), idx(prep,2)
            end do
        end if
        write(fi,*) "---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----"
        close(fi)

        open(newunit=fi,file=foldername//filename3,action="write",position="append",status="old",iostat=is)
        if(is/=0) then 
            print*, "Error: cannot open file.",filename3
            stop
        end if
            write(fi,"(I6,1X)",advance="no") iter
            do prep=1,num_n
                write(fi,"(A5,1X)",advance="no") init_rzm_to_write(idx(prep,1))
            end do
            write(fi,"(A1,1X)",advance="no") "|"
            do prep=1,num_p
                write(fi,"(A5,1X)",advance="no") init_rzm_to_write(idx(prep,2))
            end do
        write(fi,*)
        close(fi)
    end subroutine write_prep_sort
    
    subroutine write_wavefunction(num)
        use global_variables,only:wf_plus,wf_minus,r_vec,z_vec
        use constants_and_parameters,only:Nr,Nz,num_n
        implicit none
        integer(int32),intent(in) :: num
        integer(int32) :: i,j,k,is=1,fi
        character(:), allocatable :: filename
        character(10) :: char4
        write (char4, "(i0)") num
        
        allocate(filename, source="wavefunc")
        open(newunit=fi,file=foldername//"orbital_wavefunction/"//filename//trim(char4)//".txt"&
        ,action="write",status="new",iostat=is)
        if(is/=0) then 
            print*, "Error: cannot open file.",filename//char4
            stop
        end if

        write(fi,*) "r z wfn1u wfn1d wfn2u wfn2d wfp1u wfp1d wfp2u wfp2d"
        
        do i=1,Nr
            do j=1,Nz
                ! ここで先に i*dr と j*dz を出力
                write(fi, '(F8.3, F8.3)', advance='no') r_vec(i), z_vec(j)
                ! ここで5次元配列の指定された要素を出力
                do k=1,2
                    write(fi, '(4F24.17)', advance='no')&
                     wf_plus(i,j,k), wf_minus(i,j,k)
                end do
                do k=1,2
                    write(fi, '(4F24.17)', advance='no')&
                     wf_plus(i,j,num_n+k), wf_minus(i,j,num_n+k)
                end do
                ! 新しい行に進むための改行
                write(fi, *)
            end do
        end do
        close(fi)
        deallocate(filename)
    end subroutine write_wavefunction

    subroutine write_d_wavefunction(num)
        use global_variables,only:r_vec,z_vec,dr_wf_plus,dr_wf_minus,dz_wf_plus,dz_wf_minus
        use constants_and_parameters,only:Nr,Nz,num_n,num_p
        implicit none
        integer(int32),intent(in) :: num
        integer(int32) :: i,j,k,is=1,fi
        character(:), allocatable :: filename
        character(10) :: char4
        write (char4, "(i0)") num
        
        allocate(filename, source="wavefunc_d")
        open(newunit=fi,file=foldername//"orbital_wavefunction/"//filename//trim(char4)//".txt"&
            ,action="write",status="new",iostat=is)
        if(is/=0) then 
            print*, "Error: cannot open file.",filename//char4
            stop
        end if

        write(fi,*) "r z dr_wfn1u dz_wfn1u dr_wfn1d dz_wfn1d &
                         dr_wfn2u dz_wfn2u dr_wfn2d dz_wfn2d &
                         dr_wfp1u dz_wfp1u dr_wfp1d dz_wfp1d &
                         dr_wfp2u dz_wfp2u dr_wfp2d dz_wfp2d "
        
        do i=1,Nr
            do j=1,Nz
                ! ここで先に i*dr と j*dz を出力
                write(fi, '(F8.3, F8.3)', advance='no') r_vec(i), z_vec(j)
                ! ここで5次元配列の指定された要素を出力
                do k=1,2
                    write(fi, '(8F24.17)', advance='no')&
                     dr_wf_plus(i,j,k), dz_wf_plus(i,j,k), dr_wf_minus(i,j,k), dz_wf_minus(i,j,k)
                end do
                do k=1,2
                    write(fi, '(8F24.17)', advance='no')&
                     dr_wf_plus(i,j,num_n+k), dz_wf_plus(i,j,num_n+k), dr_wf_minus(i,j,num_n+k), dz_wf_minus(i,j,num_n+k)
                end do
                ! 新しい行に進むための改行
                write(fi, *)
            end do
        end do
        close(fi)
        deallocate(filename)
    end subroutine write_d_wavefunction


    subroutine write_density_and_derivative(num)
        use constants_and_parameters,only:Nr,Nz
        use global_variables,only:rho_n,rho_p,lap_rho_n,lap_rho_p, dr_rho_n,dr_rho_p,dz_rho_n,dz_rho_p &
                                 ,ddr_rho_n,ddr_rho_p,ddz_rho_n,ddz_rho_p &
                                 ,tau_n=>kinetic_dens_n,tau_p=>kinetic_dens_p &
                                 ,thomas_fermi_n,thomas_fermi_p &
                                 ,j_n=>spin_orbit_dens_n,j_p=>spin_orbit_dens_p &
                                 ,divj_n=>div_spin_orbit_dens_n,divj_p=>div_spin_orbit_dens_p, r_vec, z_vec

        implicit none
        character(:), allocatable :: filename
        integer(int32),intent(in) :: num
        integer(int32) :: is=1,fi
        integer(int32) :: i,j
        character(10) :: char4
        write (char4, "(i0)") num
        allocate(filename, source="density_and_derivative")
        open(newunit=fi,file=foldername//"densities/"//filename//trim(char4)//".txt"&
        ,action="write",status="new",iostat=is)
        if(is/=0) then
             print*, "Error: cannot open file.",filename//char4
             stop
        end if

        write(fi,*) "r z &
        rho_n dr_rho_n dz_rho_n ddr_rho_n ddz_rho_n lap_rho_n &
        rho_p dr_rho_p dz_rho_p ddr_rho_p ddz_rho_p lap_rho_p &
        tau_n tau_p jr_n jz_n divj_n jr_p jz_p divj_p thor_n thor_p"
        
        do i=1,Nr
            do j=1,Nz 
                !write(fi,"(F8.3, F8.3, G20.10, G20.10, G20.10, G20.10, &
                !                       G20.10, G20.10, &
                !                       G20.10, G20.10, G20.10, G20.10, G20.10, G20.10)") &
                write(fi,*) &
                r_vec(i), z_vec(j)&
                          , rho_n(i,j), dr_rho_n(i,j), dz_rho_n(i,j), ddr_rho_n(i,j), ddz_rho_n(i,j), lap_rho_n(i,j)&
                          , rho_p(i,j), dr_rho_p(i,j), dz_rho_p(i,j), ddr_rho_p(i,j), ddz_rho_p(i,j), lap_rho_p(i,j) & 
                          , tau_n(i,j), tau_p(i,j) &
                          , j_n(i,j,1), j_n(i,j,3), divj_n(i,j), j_p(i,j,1), j_p(i,j,3), divj_p(i,j) &
                          , thomas_fermi_n(i,j), thomas_fermi_p(i,j)
            end do
        end do
        close(fi)
        deallocate(filename)
    end subroutine write_density_and_derivative


    
    subroutine write_pot(num)
        use constants_and_parameters,only:Nr,Nz,bkin=>HBAR2_over_2m
        use global_variables,only:U_n,U_p,B_n,B_p,W_n,W_p,Cou=>direct_Coulomb_pot,r_vec,z_vec
        implicit none
        
        integer(int32),intent(in) :: num
        integer(int32) :: i,j,is=1,fi
        character(:), allocatable :: filename
        character(10) :: char4
        write (char4, "(i0)") num
        
        allocate(filename, source="pot")
        open(newunit=fi,file=foldername//"potentials/"//filename//trim(char4)//".txt"&
        ,action="write",status="new",iostat=is)
        if(is/=0) then 
            print*, "Error: cannot open file.",filename//char4
            stop
        end if
        
        write(fi,*) "r z Coulomb U_n U_p B_n B_p Wr_n Wz_n, Wr_p Wz_p"
        
        do i=1,Nr
            do j=1,Nz
               ! write(fi, '(F8.3, F8.3, G20.10, G20.10, G20.10, G20.10, G20.10, G20.10, G20.10, G20.10, G20.10)') &
                write(fi, *) &
                r_vec(i), z_vec(j), Cou(i,j),U_n(i,j), U_p(i,j), B_n(i,j)+bkin, B_p(i,j)+bkin&
                , W_n(i,j,1), W_n(i,j,3), W_p(i,j,1), W_p(i,j,3)
            end do
        end do
        close(fi)
        deallocate(filename)
    end subroutine write_pot

    subroutine write_U_term(b0,b0p,b1,b1p,b2,b2p,b3,b3p,b3p2,b4,b4p,dir,ex,np_str)
        use constants_and_parameters,only:Nr,Nz,max_iter
        use global_variables,only:r_vec,z_vec
        implicit none
        
        real(real64),intent(in) :: b0(:,:),b0p(:,:),b1(:,:),b1p(:,:),b2(:,:),b2p(:,:)&
                                    ,b3(:,:),b3p(:,:),b3p2(:,:),b4(:,:),b4p(:,:),dir(:,:),ex(:,:)
        character(len=1),intent(in) :: np_str
        integer(int32),save :: num = 0
        integer(int32) :: i,j,is=1,fi
        character(:), allocatable :: filename
        character(10) :: char4
        write (char4, "(i0)") num
        if(num/=0.or.num/=max_iter)then
            return
        end if
        
        allocate(filename, source="U_term")
        open(newunit=fi,file=foldername//"potentials/"//filename//np_str//char4//".txt",action="write",status="new",iostat=is)
        if(is/=0) then 
            print*, "Error: cannot open file.",filename//np_str//char4
            stop
        end if

        write(fi,*) "r z b0 b0p b1 b1p b2 b2p b3 b3p b3p2 b4 b4p dir ex"
        do j=1,Nz
            do i=1,Nr
                write(fi,*)&
                r_vec(i),z_vec(j), b0(i,j), b0p(i,j), b1(i,j), b1p(i,j), b2(i,j), b2p(i,j)&
                                             , b3(i,j), b3p(i,j), b3p2(i,j), b4(i,j), b4p(i,j), dir(i,j), ex(i,j)
            end do
        end do
        close(fi)
        deallocate(filename)
        if(np_str=="p")then
            num = num + 1
        end if
    end subroutine write_U_term

    subroutine write_iter_Energy(iter)
        use global_variables,only:Skyrme_Energy, Kinetic_Energy, Coulomb_Energy, CM_Energy, Total_Energy, &
                                   b_0_term,b_3_term,b_2_term,b_1_term,b_4_term,c_1_term
        implicit none
        integer(int32) :: fi,is
        integer(int32),intent(in) :: iter
        logical,save :: make_index = .true.
        character(:), allocatable :: filename


        allocate(filename, source="iter_Energy.txt")
        open(newunit=fi,file=foldername//filename,action="write",position="append",status="old",iostat=is)
        if(is/=0) print*, "Error: cannot open file.",filename
        if(make_index) then
            write(fi,*) "Total Kinetic Coulomb CM Skyrme b_0 b_1 b_2 b_3 b_4 c_1"
            make_index = .false.
        end if
        !write(fi,'(f10.3,f10.3,f10.3,f10.3,f10.3,f10.3,f10.3,f10.3,f10.3,f10.3,f10.3)')  &
        write(fi,*) &
        Total_Energy, Kinetic_Energy, Coulomb_Energy, CM_Energy, Skyrme_Energy, &
                    b_0_term,b_1_term,b_2_term,b_3_term,b_4_term,c_1_term

        close(fi)
    end subroutine write_iter_Energy
    
    
    function overlap_to_file(idx_1,idx_2) result(ovlp)
        use math_integrate,only: volume_integrate
        use global_variables,only:wf_plus,wf_minus
        implicit none
        integer(int32),intent(in) :: idx_1
        integer(int32),intent(in) :: idx_2
        real(real64)              :: ovlp
        real(real64)              :: term_plus,term_minus
        
        term_plus = volume_integrate(f=wf_plus(:,:,idx_1)*wf_plus(:,:,idx_2))
        term_minus= volume_integrate(f=wf_minus(:,:,idx_1)*wf_minus(:,:,idx_2))
        
        ovlp = term_plus + term_minus
    end function overlap_to_file
end module write_to_file