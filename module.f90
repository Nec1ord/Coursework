MODULE my_module   
    
    ! variables
    integer     NZ, IFAIL             ! for S17DLF - Hankel fns
    integer*2   i
    integer*2   j
    real*8      step
    real*8      r
    
    ! const
    integer*2,  parameter:: Nx          = 500
    integer*2,  parameter:: Ny          = 500
    real*8,     parameter:: start_x     = -100d0
    real*8,     parameter:: start_y     = -100d0
    real*8,     parameter:: end_x       = 100d0
    real*8,     parameter:: end_y       = 100d0
    real*8,     parameter:: x0          = 0d0
    real*8,     parameter:: y0          = 0d0
    real*8,     parameter:: omega       = 2.5d0
    real*8,     parameter:: p           = 1.21d0
    real*8,     parameter:: v           = 11.21d0
    real*8,     parameter:: E           = 53.21d0
    real*8,     parameter:: lambda      = (v*E) / ((1+v)*(1-2*v))!10d0
    real*8,     parameter:: mu          = E / (2*(1+v))!86.2d0
    real*8,     parameter:: gamma       = (3 * lambda + 2 * mu) * 2d0
    complex*16, parameter:: complex_i   = (0d0, 1d0)
    complex*16, parameter:: kappa       = 30d0
    complex*16, parameter:: k           = p * omega**2 
    complex*16, parameter:: kappa1      = k / (lambda + 2 * mu)
    complex*16, parameter:: k1          = omega / kappa
    
    ! arrays
    real*8,     dimension(Nx):: x
    real*8,     dimension(Ny):: y
    complex*16, dimension(2):: Hankel_1
    complex*16, dimension(2):: Hankel_2
    complex*16, dimension(Nx, Ny):: U_xy
    complex*16, dimension(Nx, Ny):: V_xy
    complex*16, dimension(Nx, Ny):: T_xy
    
END MODULE my_module