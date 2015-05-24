MODULE my_module   
    
    ! variables
    integer     NZ, IFAIL             ! for S17DLF - Hankel fns
    integer*2   i
    integer*2   j
    real*8      step
    real*8      r
    
    ! const
    integer*2,  parameter:: Nx      = 200
    integer*2,  parameter:: Ny      = 200
    real*8,     parameter:: start_x = -50d0
    real*8,     parameter:: start_y = -50d0
    real*8,     parameter:: end_x   = 50d0
    real*8,     parameter:: end_y   = 50d0
    real*8,     parameter:: x0      = 20.11d0
    real*8,     parameter:: y0      = 20.43d0
    real*8,     parameter:: omega   = 1d0
    real*8,     parameter:: lambda  = 1d0
    real*8,     parameter:: mu      = 1d0
    real*8,     parameter:: gamma   = 1d0
    complex*16, parameter:: complex_i   = (0d0, 1d0)
    complex*16, parameter:: kappa       = 1d0
    complex*16, parameter:: k           = 1d0
    
    complex*16, parameter:: kappa1      = (k**2) / (lambda + 2 * mu)
    complex*16, parameter:: k1          = (omega) / (kappa)
    
    ! arrays
    real*8,     dimension(Nx):: x
    real*8,     dimension(Ny):: y
    complex*16, dimension(2):: Hankel_1
    complex*16, dimension(2):: Hankel_2
    complex*16, dimension(Nx, Ny):: U_xy
    complex*16, dimension(Nx, Ny):: V_xy
    complex*16, dimension(Nx, Ny):: T_xy
    
END MODULE my_module