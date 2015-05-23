program coursework
    
    use my_module    
5   format(5(2x, f15.10))
    
    ! x : (start_x, end_x)
    step = (end_x - start_x) / (Nx);
    do i = 1, Nx
        x(i) = start_x + i * step
    end do
    
    ! y : (start_y, end_y)
    step = (end_y - start_y) / (Ny);
    do i = 1, Ny
        y(i) = start_y + i * step
    end do
    
    ! count U, V, T
    do i = 1, Nx    
        do j = 1, Ny
            r = sqrt((x(i) - start_x)**2 + (y(j) - start_y)**2)
            U_xy(i, j) = ( complex_i * gamma * (x(i) - start_x) ) / (4d0 * kappa * (lamda + 2 * mu) * (kappa1 - complex_i * k1) * r)
            V_xy(i, j) = ( complex_i * gamma * (y(i) - start_y) ) / (4d0 * kappa * (lamda + 2 * mu) * (kappa1 - complex_i * k1) * r)
            T_xy(i, j) = -(complex_i) / (4d0 * kappa)
            
            call S17DLF(1, 0d0, r * sqrt(kappa1),         2, 'u', Hankel_1, NZ, IFAIL)
            call S17DLF(1, 0d0, r * sqrt(k1 * complex_i), 2, 'u', Hankel_2, NZ, IFAIL)
            
            U_xy(i, j) = U_xy(i, j) * (sqrt(kappa1) * Hankel_1(2) - sqrt(complex_i * k1) * Hankel_2(2))
            V_xy(i, j) = V_xy(i, j) * (sqrt(kappa1) * Hankel_1(2) - sqrt(complex_i * k1) * Hankel_2(2))
            
            call S17DLF(1, 0d0, r * sqrt(k1 * complex_i), 1, 'u', Hankel_2, NZ, IFAIL)
            
            T_xy(i, j) = T_xy(i, j) * Hankel_2(1)
        end do
    end do
    
    ! write U_xy
    open(6, file = 'output_u.dat')
    do i = 1, Nx
        do j = 1, Ny
            write(6, 5) U_xy(i, j)
            print 5, U_xy(i, j)         ! debug
        end do
    end do
    
    ! write V_xy
    open(6, file = 'output_v.dat')
    do i = 1, Nx
        do j = 1, Ny
            write(6, 5) V_xy(i, j)
            print 5, V_xy(i, j)         ! debug
        end do
    end do
    
    ! write T_xy
    open(6, file = 'output_t.dat')
    do i = 1, Nx
        do j = 1, Ny
            write(6, 5) T_xy(i, j)
            print 5, T_xy(i, j)         ! debug
        end do
    end do
    
    
    pause 1
end program