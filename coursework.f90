program coursework
    
    use my_module    
    
1   format(5(2x, f20.15))
    open(2, file = 'x.dat')
    open(3, file = 'y.dat')
    open(4, file = 'output_u.dat')
    open(5, file = 'output_v.dat')
    open(6, file = 'output_t.dat')
    
    ! x : (start_x, end_x)
    step = (end_x - start_x) / Nx
    x(1) = start_x
    do i = 2, Nx
        x(i) = x(i-1) + step
    end do
    
    ! y : (start_y, end_y)
    step = (end_y - start_y) / Ny
    y(1) = 10d0!start_y
    do i = 2, Ny
        y(i) = 0d0!y(i-1) + step
    end do
    
    j=1
    ! count U, V, T
    do i = 1, Nx    
!        do j = 1, Ny
            r = sqrt((x(i) - x0)**2 + (y(j) - y0)**2)
            U_xy(i, j) = ( complex_i * gamma * (x(i) - x0) ) / (4d0 * kappa * (lamda + 2 * mu) * (kappa1 - complex_i * k1) * r)
            V_xy(i, j) = ( complex_i * gamma * (y(j) - y0) ) / (4d0 * kappa * (lamda + 2 * mu) * (kappa1 - complex_i * k1) * r)
            T_xy(i, j) = -complex_i / (4d0 * kappa)
            
            call S17DLF(1, 0d0, r * sqrt(kappa1),         2, 'u', Hankel_1, NZ, IFAIL)
            call S17DLF(1, 0d0, r * sqrt(k1 * complex_i), 2, 'u', Hankel_2, NZ, IFAIL)
            
            U_xy(i, j) = U_xy(i, j) * (sqrt(kappa1) * Hankel_1(2) - sqrt(complex_i * k1) * Hankel_2(2))
            V_xy(i, j) = V_xy(i, j) * (sqrt(kappa1) * Hankel_1(2) - sqrt(complex_i * k1) * Hankel_2(2))
            
            call S17DLF(1, 0d0, r * sqrt(k1 * complex_i), 1, 'u', Hankel_2, NZ, IFAIL)
            
            
            T_xy(i, j) = T_xy(i, j) * Hankel_2(1)
 !       end do
    end do
    
    ! write x, y, U_xy, V_xy, T_xy
    do i = 1, Nx
        write(2, 1) x(i)
        write(3, 1) y(i)
        j=1
  !      do j = 1, Ny
            write(4, 1) abs(U_xy(i, j))
            write(5, 1) abs(V_xy(i, j))
            write(6, 1) abs(T_xy(i, j))
  !      end do
    end do    
    
    pause 1
end program