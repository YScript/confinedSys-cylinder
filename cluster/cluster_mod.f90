MODULE cluster_mod
USE global_para

IMPLICIT NONE

TYPE :: cluster
     INTEGER :: last             !胶束中链子或者单体的数目.
     INTEGER, DIMENSION(:), ALLOCATABLE :: molecules   !胶束中链子或者单体的编号.
END TYPE cluster

TYPE :: nlist
     INTEGER :: last            !单体近邻的数目
     INTEGER :: neighbors(200)  !单体近邻的编号
END TYPE nlist

TYPE(cluster),DIMENSION (:), ALLOCATABLE  :: cluster_monomer !胶束中的单体
TYPE(cluster),DIMENSION (:), ALLOCATABLE  :: cluster_chain   !胶束中的链子
integer, DIMENSION(:), ALLOCATABLE :: isVisited
integer, dimension(:), allocatable :: clusterNumber, clusterNumber_chain  !chusternumber是胶束的编号，chusternumber_chain是链子对应的胶束的编号
integer :: N_cluster
double precision, parameter :: d_samecluster=1.50d0

CONTAINS


!************************************************************************************************************************************
SUBROUTINE cluster_domain
integer :: istep, c, i,j
integer :: i_u, i_d, i_t
double precision :: instan_phi_unimer, ave_phi_unimer, instan_Mn, ave_Mn, instan_Mw, ave_Mw
ALLOCATE (cluster_monomer(N_chain))
ALLOCATE (cluster_chain(N_chain))
ALLOCATE (isVisited(N_monomerT))
ALLOCATE (clusterNumber(N_monomerT)) !单体对应的胶束的编号.
ALLOCATE (clusterNumber_chain(N_chain)) !第i条链对应的胶束的编号.

do i = 1, N_chain
  ALLOCATE(cluster_monomer(i)%molecules(N_monomerT))
  ALLOCATE(cluster_chain(i)%molecules(N_chain))
enddo

ave_phi_unimer = 0.0d0
ave_Mn = 0.0d0
ave_Mw = 0.0d0
N_snapshot = 0

open(10, file = 'instan_info.dat')               !,status='replace', access='append')
open(20, file = 'ave_Mn_Mw.dat')

do 1000 istep = start_step, stop_step, stepcr   
        N_snapshot = N_snapshot  + 1      !计算抽样的次数

   call step_cluster(istep, N_cluster, cluster_monomer, clusterNumber, cluster_chain, clusterNumber_chain) 

!------------------------------------------------------------------------------------------------------------------
do c =1 ,N_cluster
write(*,*)c, cluster_monomer(c)%last
enddo

stop

     i_u = 0
   do c = 1, N_cluster
     if(cluster_chain(c)%last.EQ.1) then  !!!若第c个胶束中链子的数目为1
	  i_u = i_u + 1   !!!!计算自由链的数目
	 endif
   enddo
write(*,*)N_cluster,i_u
   instan_phi_unimer = dble(i_u)*N_monomer*PI/Lx/Ly/Lz/6.0d0   !!!!自由链的体积分数
   ave_phi_unimer = ave_phi_unimer + instan_phi_unimer

   instan_Mn = 0.0d0  !瞬时数均聚集数
   instan_Mw = 0.0d0  !瞬时重均聚集数
   do c = 1, N_cluster
      instan_Mn = instan_Mn + dble(cluster_chain(c)%last)  !!! cluster_chain(c)%last第c个胶束中链子的数目
      instan_Mw = instan_Mw + dble(cluster_chain(c)%last) * dble(cluster_chain(c)%last)
   enddo
   instan_Mn = instan_Mn/dble(N_cluster) 
   instan_Mw = instan_Mw/dble(N_chain)   
   ave_Mn = ave_Mn + instan_Mn
   ave_Mw = ave_Mw + instan_Mw

   write(10,*)istep, sngl(instan_phi_unimer), sngl(instan_Mn), sngl(instan_Mw)
   write(*,*)istep, N_cluster, sngl(instan_phi_unimer), sngl(instan_Mn), sngl(instan_Mw)

1000 enddo !snapshot loop.
close(10)

write(20,*)'average=', sngl(ave_phi_unimer/dble(N_snapshot)), sngl(ave_Mn/dble(N_snapshot)), sngl(ave_Mw/dble(N_snapshot))
write(*,*)'average=', sngl(ave_phi_unimer/dble(N_snapshot)), sngl(ave_Mn/dble(N_snapshot)), sngl(ave_Mw/dble(N_snapshot))
close(20)


RETURN
END SUBROUTINE cluster_domain




!************************************************************************************************************************************
SUBROUTINE step_cluster(nstep, N_cc, cluster_monomer, clusterNumber, cluster_chain, clusterNumber_chain)
USE xyz
USE QUEUE_UTILITY

integer:: nstep, N_cc
TYPE(cluster) :: cluster_monomer(*), cluster_chain(*)
integer :: clusterNumber(*), clusterNumber_chain(*)
TYPE(nlist), DIMENSION (:), ALLOCATABLE :: neighborlist
integer :: i, startIndex
integer :: j, c, m, id_chain
integer, dimension(:),allocatable :: chain_selected
ALLOCATE (neighborlist(N_monomerT))

call xyz_in(nstep)
call neighbor_list(N_monomerT, neighborlist)

!以下程序是去找胶束中的单体：可以知道胶束的个数N_cc，每个胶束中单体的个数cluster_monomer(i)%last，胶束中每个单体的编号cluster_monomer(i)%molecules(j)
do i = 1, N_monomerT
   isVisited(i) = 0
   clusterNumber(i) = 0
enddo
do i = 1, N_chain
   cluster_monomer(i)%last = 0 !!!
   cluster_monomer(i)%molecules(:)=0
enddo
N_cc = 0   !!! N_cc用来统计胶束的数目
do i = 1, N_monomerT
   if((isVisited(i).EQ.0).and.(type_monomer(i).EQ.'C'))then 
     startIndex = i
     N_cc = N_cc + 1
     if(N_cc.GT.N_chain)then
            write(6,*)'step=',nstep, 'N_cluster=', N_cc
            write(6,*)'the number of clusters is too large!'
     stop
     endif
     call cluster_by_bfs(startIndex,N_cc,cluster_monomer, clusterNumber,neighborlist)
   endif
enddo
DEALLOCATE (neighborlist)


!以下程序是去寻找胶束中的链子
do i = 1, N_chain
  cluster_chain(i)%molecules(:) = 0
  clusterNumber_chain(i) = 0
enddo
ALLOCATE(chain_selected(N_chain))
chain_selected(:) = 0
do 2000 c=1, N_cc
   call que_clear

  do i = 1, cluster_monomer(c)%last
     j = cluster_monomer(c)%molecules(i)   !get the id of ith monomer in c-th cluster
     if(type_monomer(j).NE.'C')then
          write(*,*)'type of monomer error!'
          stop
     endif  
     m = MOD(j,N_monomer)                
     if(m.EQ.0) m = N_monomer
     id_chain = (j-m)/N_monomer + 1       !get the id of chain in c-th cluster
         if(chain_selected(id_chain).eq.0)then
           call que_push(id_chain)
           chain_selected(id_chain) = 1
         endif
  enddo
  cluster_chain(c)%last = que_size()
  do i = 1, cluster_chain(c)%last
     cluster_chain(c)%molecules(i)= que_front()
     clusterNumber_chain(cluster_chain(c)%molecules(i)) = c
     call que_pop
  enddo
2000 enddo 
DEALLOCATE(chain_selected)
RETURN
END SUBROUTINE step_cluster




!**************************************************************************************************
SUBROUTINE cluster_by_bfs(startI, N_c, cluster_monomer, clusterNumber,neighborlist)
USE QUEUE_UTILITY

TYPE(cluster) :: cluster_monomer(*)
integer :: clusterNumber(*)
TYPE(nlist) ::  neighborlist(*)
integer :: startI, N_c, walker, i, neighbor
integer  pushed(N_monomerT)

call que_clear

pushed(:) = 0
call que_push(startI)
pushed(startI) = 1    !把startI放到队列中
do while(que_empty()==0) !!!!若队列不是空的
   walker = que_front()   !!!把队列中的第一个值赋给Walker
   call que_pop
  if((isVisited(walker).EQ.0).and.(type_monomer(walker).EQ.'C'))then
   clusterNumber(walker) = N_c    !get the ID of cluster for the walker monomer
   isVisited(walker) = 1
   cluster_monomer(N_c)%last = cluster_monomer(N_c)%last + 1
   cluster_monomer(N_c)%molecules(cluster_monomer(N_c)%last) = walker
   do i = 1, neighborlist(walker)%last
     neighbor = neighborlist(walker)%neighbors(i)
     if((isVisited(neighbor).EQ.0).and.(type_monomer(neighbor).EQ.'C').and.(pushed(neighbor).eq.0))then
       call que_push(neighbor)
       pushed(neighbor)=1
     endif
   enddo
  endif 
enddo

RETURN
END SUBROUTINE cluster_by_bfs

!**********************************************************************************************
subroutine neighbor_list(iNtemp,neighborlist)
TYPE(nlist) ::  neighborlist(*)
integer:: iNtemp
integer:: i, ii,jj
double precision :: rxi, ryi, rzi, rxj, ryj, rzj, dr2, drx, dry, drz

do i = 1, iNtemp
neighborlist(i)%last = 0
neighborlist(i)%neighbors(:) = 0
enddo

do ii = 1, iNtemp-1
   rxi = Rxx(ii)
   ryi = Ryy(ii)
   rzi = Rzz(ii)
do jj= ii+1,  iNtemp
   rxj = Rxx(jj)
   ryj = Ryy(jj)
   rzj = Rzz(jj)
           
    drx = rxi-rxj
    dry = ryi-ryj
    drz = rzi-rzj
! apply the minimum image convention
    drx = drx - DNINT(drx*Lxi)*Lx
    dry = dry - DNINT(dry*Lyi)*Ly
    drz = drz - DNINT(drz*Lzi)*Lz
    dr2 = drx*drx + dry*dry + drz*drz     
    if(dsqrt(dr2).LE.d_samecluster)then
      neighborlist(ii)%last = neighborlist(ii)%last + 1
      neighborlist(ii)%neighbors(neighborlist(ii)%last) = jj
      neighborlist(jj)%last = neighborlist(jj)%last + 1
      neighborlist(jj)%neighbors(neighborlist(jj)%last) = ii              
    endif

enddo !jj
enddo  !ii

return
end subroutine neighbor_list


END MODULE cluster_mod
