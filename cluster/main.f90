program main

use global_para 
use cluster_mod
use xyz

IMPLICIT NONE

include 'parameter.h'

Lxi=1.0d0/Lx
Lyi=1.0d0/Ly
Lzi=1.0d0/Lz


write(*,*)N_monomerT, N_chain*(N_monomerA+N_monomerB)+N_solvent

ALLOCATE(Rxx(Ntot))
ALLOCATE(Ryy(Ntot))
ALLOCATE(Rzz(Ntot))
ALLOCATE(type_monomer(Ntot))


call cluster_domain            !调用cluster子程序
!call xyz_out(100)


end program main
