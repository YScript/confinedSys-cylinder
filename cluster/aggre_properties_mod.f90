MODULE aggre_propertise_mod
USE global_para
USE cluster_mod

IMPLICIT NONE

CONTAINS


SUBROUTINE aggregate_Rg
integer :: i,j,jj,k,c, ngrow,nn,n,mi
double precision :: Rgx2, Rgy2, Rgz2, Rg2
double precision :: Rcx, Rcy, Rcz, dx, dy, dz, sumx, sumy, sumz
call cluster_domain
 ! ALLOCATE (cluster_chain(N_chain))
 ! ALLOCATE (clusterNumber(N_chain)) !the No. of cluster for the ith chain.
 ! do i = 1, N_chain
  !   ALLOCATE(cluster_chain%molecules(N_chain))
  !enddo

!  call step_cluster(1600000, N_cluster, cluster_chain, clusterNumber)

!write(*,*)N_cluster_1
!write(*,*) cluster_chain_1(:)%last
!stop

c = 1 ! the id of the largest aggregate

! calculate the center of mass for the largest aggregate
Rcx = 0.0d0
Rcy = 0.0d0
Rcz = 0.0d0
mi=0
do i = 1, cluster_chain(c)%last
   j =  cluster_chain(c)%molecules(i)
   do k = 1, N_monomer
       n=(j-1)*N_monomer+k 
	   nn=(cluster_chain(c)%molecules(1)-1)*N_monomer+1
          Dx=Rxx(n)-Rxx(nn)
          Dy=Ryy(n)-Ryy(nn)
	      Dz=Rzz(n)-Rzz(nn)
	    if(Dx.gt.lx/2)then
          Rxx(n)=Rxx(n)-lx
	    else if(Dx.lt.-lx/2)then
          Rxx(n)=Rxx(n)+lx
	    endif
	    if(Dy.gt.ly/2)then
          Ryy(n)=Ryy(n)-ly
	    else if(Dy.lt.-ly/2)then
          Ryy(n)=Ryy(n)+ly
	    endif
       if(Dz.gt.lz/2)then
         Rzz(n)=Rzz(n)-lz
	   else if(Dz.lt.-lz/2)then
         Rzz(n)=Rzz(n)+lz
	   endif 
   Rcx =Rcx+Rxx(n)
   Rcy =Rcy+Ryy(n)
   Rcz =Rcz+Rzz(n)
   mi=mi+1
   enddo
enddo
Rcx =Rcx/mi
Rcy =Rcy/mi
Rcz =Rcz/mi
write(*,*)rcx,rcy,rcz
do i = 1, cluster_chain(c)%last
     j =  cluster_chain(c)%molecules(i)
   do k = 1, N_monomer
      jj = (j-1)*N_monomer + k  !ÊòØÂê¶Â∫îËØ•Âå∫ÂàÜÁ±ªÂûãÔºüÔºüÔºü
      dx = Rxx(jj)-Rcx
      dy = Ryy(jj)-Rcy
      dz = Rzz(jj)-Rcz
    !  dx = dx - Lx*dnint( dx*Lxi )
     ! dy = dy - Ly*dnint( dy*Lyi )
     ! dz = dz - Lz*dnint( dz*Lzi )
      Rgx2 =  dx**2
      Rgy2 =  dy**2
      Rgz2 =  dz**2
      Rg2 = Rg2 + Rgx2 + Rgy2 + Rgz2 
   enddo
enddo
Rg2 = Rg2/dble(cluster_chain(c)%last*N_monomer)
write(*,*)c,sngl(dsqrt(Rg2))
RETURN
END SUBROUTINE aggregate_Rg
!*************************************************************

SUBROUTINE aggregate_eccentricity
USE math_mod
INTEGER :: i,j,k,jj,info1,ngrow,c,n,nn,mi
integer, parameter :: rank = 3
DOUBLE PRECISION :: work(3*rank)
DOUBLE PRECISION :: TLNA(2*rank*rank),TLNB(2*rank*rank)
DOUBLE PRECISION :: TLNA1(2*rank*rank),TLNB1(2*rank*rank)
DOUBLE PRECISION :: AH(rank**2+rank)
double precision :: H(rank,rank), eA(rank)
double precision :: Rcx, Rcy, Rcz, dx, dy, dz
double precision :: sum1, sum2, sum3, sum4, sum5, sum6
call cluster_domain
open(50, file = 'aggregate_eccentricity.dat')
!ALLOCATE (cluster_chain(N_chain))
!ALLOCATE (clusterNumber(N_chain)) !the No. of cluster for the ith chain.
!  do i = 1, N_chain
 !    ALLOCATE(cluster_chain(i)%molecules(N_chain))
 ! enddo

!call step_cluster(1600000, N_cluster, cluster_chain, clusterNumber)


!write(*,*)N_cluster_1
!write(*,*) cluster_chain_1(:)%last
!stop

c=19 ! the id of the largest aggregate
! calculate the center of mass for the largest aggregate
Rcx = 0.0d0
Rcy = 0.0d0
Rcz = 0.0d0
do i = 1, cluster_chain(c)%last !!!!!!!!!!!µ⁄c∏ˆΩ∫ ¯÷–µƒ¡¥◊” ˝
   j =  cluster_chain(c)%molecules(i)   !!!!!!µ⁄c∏ˆΩ∫ ¯÷–µ⁄iÃı¡¥◊”µƒ±‡∫≈
   do k = 1, N_monomer
       n=(j-1)*N_monomer+k 
	   nn=(cluster_chain(c)%molecules(1)-1)*N_monomer+1
          Dx=Rxx(n)-Rxx(nn)
          Dy=Ryy(n)-Ryy(nn)
	      Dz=Rzz(n)-Rzz(nn)
	    if(Dx.gt.lx/2)then
          Rxx(n)=Rxx(n)-lx
	    else if(Dx.lt.-lx/2)then
          Rxx(n)=Rxx(n)+lx
	    endif
	    if(Dy.gt.ly/2)then
          Ryy(n)=Ryy(n)-ly
	    else if(Dy.lt.-ly/2)then
          Ryy(n)=Ryy(n)+ly
	    endif
       if(Dz.gt.lz/2)then
         Rzz(n)=Rzz(n)-lz
	   else if(Dz.lt.-lz/2)then
         Rzz(n)=Rzz(n)+lz
	   endif 
   Rcx =Rcx+Rxx(n)
   Rcy =Rcy+Ryy(n)
   Rcz =Rcz+Rzz(n)
   mi=mi+1
   enddo
enddo
Rcx =Rcx/mi
Rcy =Rcy/mi
Rcz =Rcz/mi
!write(*,*)rcx,rcy,rcz
! construct the real-symmetric H matrix
H(:,:) = 0.0d0
sum1 = 0.0d0
sum2 = 0.0d0
sum3 = 0.0d0
sum4 = 0.0d0
sum5 = 0.0d0
sum6 = 0.0d0
do i = 1, cluster_chain(c)%last
   j =  cluster_chain(c)%molecules(i)
   do k = 1, N_monomer
      jj = (j-1)*N_monomer + k  !ÊòØÂê¶Â∫îËØ•Âå∫ÂàÜÁ±ªÂûãÔºüÔºüÔºü
      dx = Rxx(jj)-Rcx
      dy = Ryy(jj)-Rcy
      dz = Rzz(jj)-Rcz
      dx = dx - Lx*dnint( dx*Lxi )
      dy = dy - Ly*dnint( dy*Lyi )
      dz = dz - Lz*dnint( dz*Lzi )  
      sum1 = sum1 + dy**2 + dz**2
      sum2 = sum2 + dx**2 + dz**2
      sum3 = sum3 + dx**2 + dy**2
      sum4 = sum4 - dx*dy
      sum5 = sum5 - dx*dz
      sum6 = sum6 - dy*dz
   enddo
enddo
H(1,1) = sum1
H(2,2) = sum2
H(3,3) = sum3
H(1,2) = sum4
H(1,3) = sum5
H(2,3) = sum6
!=============================================
! Find the eigenvalues and eigenfunctions of
! H using Eigch.f
!=============================================
        k=1
        DO j=1,rank
        DO i=1,j
        AH(2*k-1)=H(i,j)
        AH(2*k)=0.0d0
        k=k+1
        End Do
        End Do
    eA(:) = 0.0d0 
    call EIGCH(AH,rank,1,eA,TLNA,rank,work,info1)
        do j=1,2*rank*rank,2
        TLNA1((j+1)/2)=TLNA(j)
        enddo

        do 111 j=1,rank
        do i=1,rank
        k=rank*(j-1)+i
        H(i,j)=TLNA1(k)
        end do
111     continue
write(*,*)eA(:)
write(50,*)eA(:)
write(*,*)H(:,1)
write(*,*)H(:,2)
write(*,*)H(:,3)
RETURN
close(50)
END SUBROUTINE aggregate_eccentricity     
END MODULE aggre_propertise_mod

