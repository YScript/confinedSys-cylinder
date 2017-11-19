MODULE global_para

IMPLICIT NONE

integer :: N_monomerT, N_chain, N_monomer, N_monomerA, N_monomerB, N_solvent, Ntot
integer ::  start_step, stop_step, stepcr, N_snapshot

double precision, parameter :: PI=3.14159265358979323846264338327950d0  

double precision :: Lx, Ly, Lz, Lxi, Lyi, Lzi, rho
double precision, dimension(:), allocatable :: Rxx, Ryy, Rzz

!ix,iy,iz are the numbers of penetrating through x`y`z direction separately
character, dimension(:), allocatable :: type_monomer 









END MODULE global_para
