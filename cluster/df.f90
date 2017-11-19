MODULE df
USE global_para
USE xyz

IMPLICIT NONE

CONTAINS

subroutine rdf(rcut, N_bin)
!  Radial distribution function calculation.
!  Reads in the configuration snapshots, and computes g(r) of the snapshots. Need xyz_in
!  g(r) is the ratio between the average number density rho(r) at a distance r from any given atom and the density at
!  a distance r from an atom in an ideal gas at the same overall density.
!  Any deviation of g(r) from unity reflects correlations between the particles due to the intermolecule interactions.
integer :: istep
double precision :: rcut, rcut2, dr,r, vbin
integer :: N_bin, bin
double precision :: r2, dx, dy, dz, hL
integer :: i, j, Ntemp
integer, dimension(:), allocatable :: H

ALLOCATE(H(N_bin))
! compute the size of bins.
dr = rcut/dble(N_bin)
hL = 0.50d0*Lx

if (rcut.gt.hL)then
  write(100,*)'the maximum radius you have requested is greater than one-half the box length!'
  write(100,*)'Reducing rcut to L/2'
  call flush(100)
  rcut = hL
endif

rcut2 = rcut*rcut  
N_snapshot = 0
H(:) = 0
do 100 istep = start_step, stop_step, stepcr
    CALL xyz_in(istep)
    N_snapshot = N_snapshot + 1
!An N^2 algorithm for computing interparticle separations 
!and updating the radial distribution function histogram. 
!Implements parts of Algorithm 7 in F&S.
   do i=1, N_monomerT-1
      do j=i+1,N_monomerT
!        if(type_monomer(i).eq.'C'.and.type_monomer(j).eq.'C')then
         dx = Rxx(i)-Rxx(j)
         dy = Ryy(i)-Ryy(j)
         dz = Rzz(i)-Rzz(j)
         dx = dx - Lx*dnint(dx/Lx)
         dy = dy - Ly*dnint(dy/Ly)
         dz = dz - Lz*dnint(dz/Lz)
         r2 = dx*dx + dy*dy + dz*dz
         if (r2<rcut2) then
            bin=nint(dsqrt(r2)/dr+0.5)
            H(bin)=H(bin)+2
!         endif
        endif
       enddo
    enddo
100 continue

!Normalize and output g(r) to the file rdf.data */
open(unit = 91, file = 'rdf.data', status= 'replace')
do i=1, N_bin
   r=dr*(dble(i)-0.50d0)
   vbin=(4.0d0/3.0d0)*PI*dble((i)**3-(i-1)**3)*dr*dr*dr
   write(91,*) r, dble(H(i))/vbin/dble(N_monomerT)/rho/dble(N_snapshot)
enddo
close (91)

DEALLOCATE(H)
return
end subroutine rdf
!**************************************************************************
subroutine bdf(rcut, N_bin)
! calculate the distribution function of bond length
integer :: istep, i, j, Ngrow, Ngrowminus1
double precision :: rcut, rcut2, dr, r2, r
integer :: N_bin, bin
double precision ::  disp_x, disp_y, disp_z, hL
!double precision, dimension(:), allocatable :: Rxx_dpl, Ryy_dpl, Rzz_dpl
integer, dimension(:), allocatable :: H

!ALLOCATE (Rxx_dpl(Npart))
!ALLOCATE (Ryy_dpl(Npart))
!ALLOCATE (Rzz_dpl(Npart))
ALLOCATE (H(N_bin))

! compute the size of bins.
dr = rcut/dble(N_bin)
hL = 0.50d0*Lx

if (rcut>hL)then
  write(100,*)'the maximum radius you have requested is greater than one-half the box length!'
  write(100,*)'Reducing rcut to L/2'
  call flush(100)
  rcut = hL
endif

rcut2 = rcut*rcut
N_snapshot = 0
H(:) = 0
do 100 istep = start_step, stop_step, stepcr
    CALL xyz_wraped_in(istep)
    N_snapshot = N_snapshot + 1
     do j = 1, N_chain
! set the postion of first monomer is zero.
!       Rxx_dpl((j-1)*N_monomer+1) = 0.d0
!       Ryy_dpl((j-1)*N_monomer+1) = 0.d0
!       Rzz_dpl((j-1)*N_monomer+1) = 0.d0

       do i = 1, N_monomer-1
          Ngrow = (j-1)*N_monomer + (i+1)
          Ngrowminus1 = Ngrow - 1
! Determine the positions of the monomers in a chain with respect to the first monomer
          disp_x = Rxx(Ngrow) - Rxx(Ngrowminus1)
          disp_y = Ryy(Ngrow) - Ryy(Ngrowminus1)
          disp_z = Rzz(Ngrow) - Rzz(Ngrowminus1)
          disp_x = disp_x - Lx*dnint( disp_x/Lx )
          disp_y = disp_y - Ly*dnint( disp_y/Ly )
          disp_z = disp_z - Lz*dnint( disp_z/Lz )
          r2 = disp_x*disp_x + disp_y*disp_y + disp_z*disp_z
          if (r2<rcut2)then
            bin=nint(dsqrt(r2)/dr+0.50d0)
            H(bin) = H(bin)+ 1
          endif
        enddo 
      enddo
100 continue

!Normalize and output p(b) to the file pb.data */
open(unit = 90, file = 'pb.data', status= 'replace')
do i=1, N_bin
   r=dr*(dble(i)-0.50d0)
   write(90,*) r, dble(H(i))/dble(N_monomer-1)/dble(N_chain)/dble(N_snapshot)
enddo
close (90)

!DEALLOCATE (Rxx_dpl)
!DEALLOCATE (Ryy_dpl)
!DEALLOCATE (Rzz_dpl)
DEALLOCATE (H)

end subroutine bdf
!****************************************************************
subroutine bdf_n(rcut,N_bin)

integer :: i, istep
integer :: bin, N_bin
!bin是区间的序号eg.bin=1即第1号区间，N_bin是区间的个数
integer :: id1, id2
double precision :: rcut, rcut2, dr, r2, r
!rcut是区间总长度，dr是每个区间的长度，
double precision ::  disp_x, disp_y, disp_z, hL
integer, dimension(:), allocatable :: H  

ALLOCATE (H(N_bin))    !H是一维数组，代表键长落在第bin个区间内的键的个数

! compute the size of bins.
dr = rcut/dble(N_bin)  !dr是每个区间的长度
hL = 0.50d0*Lx

if (rcut>hL)then
  write(100,*)'the maximum radius you have requested is greater than one-half the box length!'
  write(100,*)'Reducing rcut to L/2' !模拟盒子每个方向的边长都要大于2Rc，这样保证盒子中每个粒子只与 盒子内的另外其他粒子 或 这些粒子的最近邻映象发生作用 
  call flush(100)
  rcut = hL
endif

rcut2 = rcut*rcut
N_snapshot = 0                              
H(:) = 0

call bondinfo_in

do 100 istep = start_step, stop_step, stepcr
!在parameter.h中已经定义过，体系在运行结束之后开始抽样，start_step,开始步，stop_step,结束步，stepcr步间隔
! 要对全部的抽样步进行统计，
    CALL xyz_wraped_in(istep) !将原子的坐标输出并且折到盒子里
    N_snapshot = N_snapshot + 1   !要对从start_step到stop_step，每隔stepcr的步数作统计
    do i=1, N_bond
       id1 = bondinfo(i)%backward
       id2 = bondinfo(i)%forward
       disp_x = Rxx(id2) - Rxx(id1)
       disp_y = Ryy(id2) - Ryy(id1)
       disp_z = Rzz(id2) - Rzz(id1)
       disp_x = disp_x - Lx*dnint( disp_x/Lx )
       disp_y = disp_y - Ly*dnint( disp_y/Ly )
       disp_z = disp_z - Lz*dnint( disp_z/Lz )
       r2 = disp_x*disp_x + disp_y*disp_y + disp_z*disp_z
       if (r2<rcut2)then
          bin=nint(dsqrt(r2)/dr+0.50d0)
          H(bin) = H(bin) + 1
       endif
    enddo

100 continue

!Normalize and output p(b) to the file pb.data */
open(unit = 90, file = 'pbn.data', status= 'replace')
do i=1, N_bin
   r=dr*(dble(i)-0.50d0)
   write(90,*) r, dble(H(i))/dble(N_bond)/dble(N_snapshot)
enddo
close (90)

DEALLOCATE (H)

end subroutine bdf_n

!*******************************************************************************
subroutine bdf_AB(rcut,N_bin)

integer :: i, istep
integer :: bin, N_bin
integer :: id1, id2
double precision :: rcut, rcut2, dr, r, r2
integer, dimension(:),allocatable :: HAA, HBB, HAB
character :: CFGFILE*10
double precision ::  disp_x, disp_y, disp_z, hL

ALLOCATE(HAA(N_bin)) 
ALLOCATE(HBB(N_bin))
ALLOCATE(HAB(N_bin))

dr=rcut/dble(N_bin)
hL=0.50d0*Lx

if (rcut>hL)then
  write(100,*)'the maximum radius you have requested is greater than one-half the box length!'
  write(100,*)'Reducing rcut to L/2'
  call flush(100)
  rcut = hL
endif

rcut2 = rcut*rcut
N_snapshot = 0
HAA(:) = 0
HBB(:) = 0
HAB(:) = 0

CALL bondinfo_in

do 100 istep = start_step, stop_step, stepcr
   CALL xyz_wraped_in(istep)
   N_snapshot = N_snapshot + 1
   do i = 1, N_bond
      id1 = bondinfo(i)%backward
      id2 = bondinfo(i)%forward
      disp_x = Rxx(id2) - Rxx(id1)
      disp_y = Ryy(id2) - Ryy(id1)
      disp_z = Rzz(id2) - Rzz(id1)
      disp_x = disp_x - Lx*DNINT(disp_x/Lx)
      disp_y = disp_y - Ly*DNINT(disp_y/Ly)
      disp_z = disp_z - Lz*DNINT(disp_z/Lz)
      r2 = disp_x*disp_x + disp_y*disp_y + disp_z*disp_z
      if (r2<rcut2) then
         if (type_monomer(id1) .EQ. 'C' .and. type_monomer(id2) .EQ. 'C') then
         bin = nint(dsqrt(r2)/dr+0.50d0)
         HAA(bin) = HAA(bin) + 1
         elseif (type_monomer(id1) .EQ. 'O' .and. type_monomer(id2) .EQ. 'O') then
         bin = nint(dsqrt(r2)/dr+0.50d0)
         HBB(bin) = HBB(bin) + 1
         elseif (type_monomer(id1) .NE. type_monomer(id2)) then
         bin = nint(dsqrt(r2)/dr+0.50d0)
         HAB(bin) = HAB(bin) + 1
         endif
      endif
   enddo
100 continue

open(70,file='PAA(b).data',status='replace')
open(80,file='PBB(b).data',status='replace')
open(90,file='PAB(b).data',status='replace')

do i = 1, N_bin
   r=dr*(dble(i)-0.50d0)
   write(70,*) r, dble(HAA(i))/dble(N_bond)/dble(N_snapshot)
   write(80,*) r, dble(HBB(i))/dble(N_bond)/dble(N_snapshot)
   write(90,*) r, dble(HAB(i))/dble(N_bond)/dble(N_snapshot)
enddo
close(70)
close(80)
close(90)

DEALLOCATE(HAA)
DEALLOCATE(HBB)
DEALLOCATE(HAB)

end subroutine bdf_AB
!**********************************************************
SUBROUTINE alb_homopolymer(rcut,N_bin) 
! calculate the distribution function of bond length
integer :: istep, i, j, Ngrow, Ngrowminus1
double precision :: rcut, rcut2, dr, r2, r
integer :: N_bin, bin
double precision ::  disp_x, disp_y, disp_z, hL
!double precision, dimension(:), allocatable :: Rxx_dpl, Ryy_dpl, Rzz_dpl
integer, dimension(:), allocatable :: H
double precision :: rp, sum_of_rp, b_ave
integer ::  n_temp
character :: xx 

!ALLOCATE (Rxx_dpl(Npart))
!ALLOCATE (Ryy_dpl(Npart))
!ALLOCATE (Rzz_dpl(Npart))
ALLOCATE (H(N_bin))

! compute the size of bins.
dr = rcut/dble(N_bin)
hL = 0.50d0*Lx

if (rcut>hL)then
  write(100,*)'the maximum radius you have requested is greater than one-half the box length!'
  write(100,*)'Reducing rcut to L/2'
  call flush(100)
  rcut = hL
endif

rcut2 = rcut*rcut
N_snapshot = 0
H(:) = 0
do 100 istep = start_step, stop_step, stepcr
    CALL xyz_wraped_in(istep)
    N_snapshot = N_snapshot + 1
     b_ave = 0.0d0
     n_temp = 0
     do j = 1, N_chain
! set the postion of first monomer is zero.
!       Rxx_dpl((j-1)*N_monomer+1) = 0.d0
!       Ryy_dpl((j-1)*N_monomer+1) = 0.d0
!       Rzz_dpl((j-1)*N_monomer+1) = 0.d0

       do i = 1, N_monomer-1
          Ngrow = (j-1)*N_monomer + (i+1)
          Ngrowminus1 = Ngrow - 1
! Determine the positions of the monomers in a chain with respect to the first
! monomer
          disp_x = Rxx(Ngrow) - Rxx(Ngrowminus1)
          disp_y = Ryy(Ngrow) - Ryy(Ngrowminus1)
          disp_z = Rzz(Ngrow) - Rzz(Ngrowminus1)
          disp_x = disp_x - Lx*dnint( disp_x/Lx )
          disp_y = disp_y - Ly*dnint( disp_y/Ly )
          disp_z = disp_z - Lz*dnint( disp_z/Lz )
          r2 = disp_x*disp_x + disp_y*disp_y + disp_z*disp_z
          if (r2<rcut2)then
            b_ave = b_ave + dsqrt(r2)
            n_temp = n_temp + 1
            bin=nint(dsqrt(r2)/dr+0.5)
            H(bin) = H(bin)+ 1
          endif
        enddo
       enddo

!write(*,*)n_temp
!read(*,*)xx

b_ave = b_ave/dble(n_temp)
open(unit = 10, file='bave.data',access = 'append')
write(10,*)  b_ave

100 continue

close(10)
stop

!Normalize and output p(b) to the file pb.data */
!open(unit = 90, file = 'pb.data', status= 'replace')
!open(unit = 208, file = 'sum_of_rp.data', status= 'new')
!sum_of_rp=0.0d0
!do i=1, N_bin
   !r=dr*(dble(i)-0.50d0)
   !rp=r*(dble(H(i))/dble(N_monomer-1)/dble(N_chain)/dble(N_snapshot))
   !sum_of_rp=sum_of_rp+rp
  ! write(90,*) r, dble(H(i))/dble(N_monomer-1)/dble(N_chain)/dble(N_snapshot)
  ! write(208,*) sum_of_rp
!enddo

!close (90)
!close (208)


!DEALLOCATE (Rxx_dpl)
!DEALLOCATE (Ryy_dpl)
!DEALLOCATE (Rzz_dpl)
DEALLOCATE (H)
END SUBROUTINE alb_homopolymer
!*********************************************************
subroutine ave_bond

integer :: i, istep
integer :: id1, id2
integer :: N_temp
double precision :: rcut, rcut2, dr, r2, r
double precision ::  disp_x, disp_y, disp_z, hL
double precision :: sum_of_r=0,A



hL = 0.50d0*Lx

if (rcut>hL)then
  write(100,*)'the maximum radius you have requested is greater than one-half the box length!'
  write(100,*)'Reducing rcut to L/2'
  call flush(100)
  rcut = hL
endif
N_temp = 0
N_snapshot=0
rcut2=rcut*rcut
call bondinfo_in

do 100 istep = start_step, stop_step, stepcr
    CALL xyz_wraped_in(istep)
    N_snapshot = N_snapshot + 1
    open(unit = 90, file = 'ave_bond.data', status = 'replace')
    do i=1, N_bond
       id1 = bondinfo(i)%backward
       id2 = bondinfo(i)%forward
       disp_x = Rxx(id2) - Rxx(id1)
       disp_y = Ryy(id2) - Ryy(id1)
       disp_z = Rzz(id2) - Rzz(id1)
       disp_x = disp_x - Lx*dnint( disp_x/Lx )
       disp_y = disp_y - Ly*dnint( disp_y/Ly )
       disp_z = disp_z - Lz*dnint( disp_z/Lz )

       r2 = disp_x*disp_x + disp_y*disp_y + disp_z*disp_z
       if (r2<rcut2)then
          
       r = dsqrt(r2)
       sum_of_r=sum_of_r+r
       N_temp = N_temp + 1
        endif
    enddo

A = sum_of_r/dble(N_temp)
write(90,*) A
100 continue
close (90)



end subroutine ave_bond
!############################################################
subroutine aver_bondlen(rcut)
integer :: i, istep
integer :: id1, id2
integer :: N_temp
double precision :: rcut, rcut2, dr, r, r2
double precision :: tot_r
double precision ::  disp_x, disp_y, disp_z, hL


hL=0.50d0*Lx

if (rcut>hL)then
  write(100,*)'the maximum radius you have requested is greater than one-half the box length!'
  write(100,*)'Reducing rcut to L/2'
  call flush(100)
  rcut = hL
endif

rcut2 = rcut*rcut
N_temp = 0
tot_r = 0.0d0

CALL bondinfo_in

do 100 istep = start_step, stop_step, stepcr
   call xyz_wraped_in(istep)
   N_snapshot = N_snapshot + 1
   do i = 1, N_bond
      id1 = bondinfo(i)%backward
      id2 = bondinfo(i)%forward
      disp_x = Rxx(id2) - Rxx(id1)
      disp_y = Ryy(id2) - Ryy(id1)
      disp_z = Rzz(id2) - Rzz(id1)
      disp_x = disp_x - Lx*DNINT(disp_x/Lx)
      disp_y = disp_y - Ly*DNINT(disp_y/Ly)
      disp_z = disp_z - Lz*DNINT(disp_z/Lz)
      r2 = disp_x*disp_x + disp_y*disp_y + disp_z*disp_z
      if (r2<rcut2)then
         tot_r = tot_r + dsqrt(r2)
         N_temp = N_temp + 1
      endif
    enddo
100 continue

tot_r = tot_r/dble(N_temp)

open(80, file = 'aver_bond.data', status = 'replace')
write(80,*)  tot_r

close(80)
stop

end subroutine aver_bondlen

!**********************************************************
END MODULE df
