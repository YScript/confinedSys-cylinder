MODULE cluster_analise
USE global_para
USE cluster_mod
IMPLICIT NONE
CONTAINS

!!!!!!!!!!!!!!******************************************************************以下是计算胶束的平均直径
SUBROUTINE cluster_diameter
integer :: i,mi,j,n,Dx,Dy,Dz,nn,DA,K,m,mm,mmo,tt
DOUBLE PRECISION :: rg22, rr,rg22A,xi,yi,zi,rrz,rrx,rry
DOUBLE PRECISION, ALLOCATABLE :: rg2(:),rg2A(:)
open(40, file = 'cluster_diameter.dat')
call cluster_domain
write(*,*) N_cluster
ALLOCATE (rg2(N_cluster),rg2A(N_cluster))
do i=1,N_cluster
  rg2(i)=0
  rg2A(i)=0
end do
rg22=0
rg22A=0
do i=1,N_cluster

   xi=0
   yi=0
   zi=0  !!!!!!(xi,yi,zi表示第N个胶束的质心坐标)
   mi=0
   do 11 j=1,cluster_chain(i)%last
     tt=cluster_chain(i)%molecules(j) !!!n表示第i个胶束上第j个单体的编号
	   do k=1,N_monomer
          n=(tt-1)*N_monomer+k 
		  nn=(cluster_chain(i)%molecules(1)-1)*N_monomer+1
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
     xi=xi+Rxx(n)
	 yi=yi+Ryy(n)
	 zi=zi+Rzz(n)
	 mi=mi+1
	 end do
 11 enddo
     xi=xi/mi
	 yi=yi/mi
	 zi=zi/mi  !!!!!!!此时已经求出质心坐标（xi,yi,zi）   
!	 write(*,*)i,xi,yi,zi   
 do 33 j=1,cluster_chain(i)%last
       m=cluster_chain(i)%molecules(j)  !!!n表示第i个胶束上第j条链子的编号
	   do k=1,N_monomer
	      mm=(m-1)*N_monomer+k         !!!nn表示第i个胶束上第j条链子上第k个单体的编号
          mmo=(cluster_chain(i)%molecules(1)-1)*N_monomer+1 
          Dx=Rxx(mm)-Rxx(mmo)
          Dy=Ryy(mm)-Ryy(mmo)
	      Dz=Rzz(mm)-Rzz(mmo)
	   if(Dx.gt.lx/2)then
         Rxx(mm)=Rxx(mm)-lx
	   else if(Dx.lt.-lx/2)then
         Rxx(mm)=Rxx(mm)+lx
	   endif
	   if(Dy.gt.ly/2)then
         Ryy(mm)=Ryy(mm)-ly
	   else if(Dy.lt.-ly/2)then
         Ryy(mm)=Ryy(mm)+ly
	   endif
       if(Dz.gt.lz/2)then
         Rzz(mm)=Rzz(mm)-lz
	   else if(Dz.lt.-lz/2)then
         Rzz(mm)=Rzz(mm)+lz
	   endif
        rrx=Rxx(mm)-xi
   	    rry=Ryy(mm)-yi
	    rrz=Rzz(mm)-zi
	    rr=rrx*rrx+rry*rry+rrz*rrz
	    rg2(i)=rg2(i)+rr
	  enddo	 
 33 enddo
    rg2(i)=rg2(i)/((cluster_chain(i)%last)*N_monomer)
  !  rg2A(i)=rg2A(i)/cluster_monomer(i)%last
    rg22=rg22+rg2(i)
   ! rg22A=rg22A+rg2A(i)
   write(*,*) i, rg2(i),sqrt(rg2(i))
enddo
rg22=rg22/N_cluster
!rg22A=rg22A/N_cluster
write(40,*) rg22,sqrt(rg22)	 
close(40)
END SUBROUTINE cluster_diameter





!!!!!!!!!!!!!!******************************************************************以下是计算胶束的径向密度分布
SUBROUTINE radical_density_distribution
Integer :: i,mi,j,m,k,mm,mmo,Dx,Dy,Dz,t,tt,u,num,cc,tt0
Double precision :: pa,pb,rr2,rr,rx,ry,rz,xi,yi,zi
INTEGER, ALLOCATABLE :: numa(:),numb(:) 
ALLOCATE (numa(15),numb(15))
open(50, file = 'radical_density_distributiion.dat')
call cluster_domain

do i=1,N_cluster
   if(i.eq.26) then    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!按照挑出来的胶束来决定t的具体值
   xi=0
   yi=0
   zi=0  !!!!!!(xi,yi,zi表示第N个胶束的质心坐标)
   mi=0
    do 33 j=1,cluster_chain(i)%last
          m=cluster_chain(i)%molecules(j)  !!!m表示第i个胶束上第j条链子的编号
	  do 44 k=1,N_monomer
	        mm=(m-1)*N_monomer+k         !!!mm表示第i个胶束上第j条链子上第k个单体的编号
            mmo=(cluster_chain(i)%molecules(1)-1)*N_monomer+1 
            Dx=Rxx(mm)-Rxx(mmo)
            Dy=Ryy(mm)-Ryy(mmo)
	        Dz=Rzz(mm)-Rzz(mmo)
	      if(Dx.gt.lx/2)then
            Rxx(mm)=Rxx(mm)-lx
	      else if(Dx.lt.-lx/2)then
            Rxx(mm)=Rxx(mm)+lx
	      endif
	      if(Dy.gt.ly/2)then
            Ryy(mm)=Ryy(mm)-ly
	      else if(Dy.lt.-ly/2)then
            Ryy(mm)=Ryy(mm)+ly
	      endif
          if(Dz.gt.lz/2)then
            Rzz(mm)=Rzz(mm)-lz
	      else if(Dz.lt.-lz/2)then
            Rzz(mm)=Rzz(mm)+lz
	      endif   
         xi=xi+Rxx(mm)
	     yi=yi+Ryy(mm)
	     zi=zi+Rzz(mm)
	     mi=mi+1
 44    enddo
 33  enddo
     xi=xi/mi
	 yi=yi/mi
	 zi=zi/mi  !!!!!!!此时已经求出质心坐标（xi,yi,zi）
   do cc=1,15
      numa(cc)=0
      numb(cc)=0
   enddo
   do 55 j=1,cluster_chain(i)%last
         t=cluster_chain(i)%molecules(j)  !!!m表示第i个胶束上第j条链子的编号		 
	  do 66 k=1,N_monomer
	        tt=(t-1)*N_monomer+k 			    
            rx=Rxx(tt)-xi
			ry=Ryy(tt)-yi
			rz=Rzz(tt)-zi
			rr=rx*rx+ry*ry+rz*rz		
			rr=sqrt(rr)		
			u=1
			cc=floor(rr/u)+1   !!!!!!!!!!!floor(x)只保留取整的函数
		  if (cc.le.15) then
		       if(type_monomer(tt).EQ.'C')then 
		         numa(cc)=numa(cc)+1
			   else if (type_monomer(tt).EQ.'O')then
			     numb(cc)=numb(cc)+1
			   end if
		  end if
66     enddo
55 enddo
else
endif
enddo
do cc=1,15
   num=numa(cc)+numb(cc)
   pa=1.0*numa(cc)/num
   pb=1.0*numb(cc)/num
   write(*,*) cc, numa(cc),numb(cc)
   write(50,*)pa,pb
enddo
close(50)
END SUBROUTINE radical_density_distribution






!!!!!!!!!!!!!!***********************************************************以下是统计链子构象的改变
SUBROUTINE conformational_change
INTEGER :: Dx,Dy,Dz,i,j,k,mm,mmc,m,rrx,rry,rrz,mmo,Dxo,Dyo,Dzo,N_bin
DOUBLE PRECISION :: rr,p,dr,bin,r
DOUBLE PRECISION, ALLOCATABLE :: num(:)
ALLOCATE (num(50))
call cluster_domain
open(60, file = 'conformational_change.dat')
   do bin=1,50
      num(bin)=0
   end do
 do 11 i=1,N_cluster
   do 22 j=1,cluster_chain(i)%last
         m=cluster_chain(i)%molecules(j) !!!!m表示第i个胶束上第j个链子的编号
      do 33 k=1,N_monomer
         mm=(m-1)*N_monomer+k         !!!mm表示第i个胶束上第j条链子上第k个单体的编号
      !   mmo=(cluster_chain(i)%molecules(1)-1)*N_monomer+1   
		 mmc=(m-1)*N_monomer+11    !!!!!!!!!!! mmc要随着具体的链子长度而改变     
            rrx=Rxx(mm)-Rxx(mmc)
            rry=Ryy(mm)-Ryy(mmc)
	        rrz=Rzz(mm)-Rzz(mmc)			
	      if(rrx.gt.lx/2) then
			rrx=-lx+rrx
			else if(rrx.lt.-lx/2) then
			rrx=lx+rrx
			end if
		    if(rry.gt.ly/2) then
		    rry=-ly+rry
		   	else if(rry.lt.-ly/2) then
			rry=ly+rry
	    	end if
		    if(rrz.gt.lz/2) then      
			rrz=-lz+rrz
			else if(rrz.lt.-lz/2) then
			rrz=lz+rrz
		    end if
		!	write(*,*)rrx,rry,rrz 		 
            rrx=ABS(rrx)                    !!!!!!!!!!!!!ABS(x)表示取绝对值的函数
            rry=ABS(rry)
		    rrz=ABS(rrz)	     
          rr=rrx*rrx+rry*rry+rrz*rrz
		  rr=sqrt(rr)
!		write(*,*)rr
		  dr=0.4 
		  N_bin=10/dr                            !!!!!!!!!!!!!!! dr表示间隔
		  bin=nint(rr/dr+0.5)                 !!!!!!!!!!!!!!!bin 用来统计间隔的个数，也就是说处于第几个间隔内
		  if(bin.le.50.0)then               !!!!!!!!!!!!!!间隔要取到0.1，还未解决？？？？？
		  num(bin)=num(bin)+1
		  end if
      33 end do
  22 end do
11 end do
do i=1,N_bin
r=dr*(dble(i)-0.5)
p=num(i)/(N_chain*N_monomer)
 write(*,*) r, p
write(60,*)r,p
end do
close (60)
END SUBROUTINE conformational_change


!!!!!!!!!!!!!******************************************************************以下是计算链子的伸展程度
SUBROUTINE stretch_hydrophobic_chains
INTEGER :: n,nn,k,m,mm,j,i
Double precision :: avrg2A,avrg2B,avrg2,rry,rrx,rrz
DOUBLE PRECISION, ALLOCATABLE :: dsrg2A(:), dsrg2B(:),dsrg2AB(:), dsrg2(:)
DOUBLE PRECISION, ALLOCATABLE :: rg2A(:), rg2B(:),rg2(:)
CALL cluster_domain
ALLOCATE (dsrg2A(N_cluster), dsrg2B(N_cluster),dsrg2AB(N_cluster),dsrg2(N_cluster))
ALLOCATE (rg2A(N_cluster), rg2B(N_cluster),rg2(N_cluster))
open(66, file = 'stretch_hydrophobic_chains.dat')

avrg2A=0
avrg2B=0
avrg2=0
do 22 i=1,N_cluster    
	  rg2A(i)=0
      rg2B(i)=0
      rg2(i)=0
	  dsrg2(i)=0
	  dsrg2A(i)=0
      dsrg2B(i)=0
      dsrg2AB(i)=0 
  do 33 j=1,cluster_chain(i)%last
        m=cluster_chain(i)%molecules(j)  !!!m表示第i个胶束上第j条链子的编号	    
	do 44 k=1,N_monomer
	      mm=(m-1)*N_monomer+k   
      do 55 n=1,N_monomer
            nn=(m-1)*N_monomer+n 
         if (k.le.N_monomerA.and.n.le.N_monomerA) then
            rrx=Rxx(mm)-Rxx(nn)
            rry=Ryy(mm)-Ryy(nn)
            rrz=Rzz(mm)-Rzz(nn)
            if(rrx.gt.lx/2) then
			rrx=-lx+rrx
			else if(rrx.lt.-lx/2) then
			rrx=lx+rrx
			end if
		    if(rry.gt.ly/2) then
		    rry=-ly+rry
		   	else if(rry.lt.-ly/2) then
			rry=ly+rry
	    	end if
		    if(rrz.gt.lz/2) then      
			rrz=-lz+rrz
			else if(rrz.lt.-lz/2) then
			rrz=lz+rrz
		    end if
	 	    dsrg2A(i)=dsrg2A(i)+rrx*rrx+rry*rry+rrz*rrz
		
         else if (k.gt.N_monomerA.and. n.gt.N_monomerA) then
            rrx=Rxx(mm)-Rxx(nn)
            rry=Ryy(mm)-Ryy(nn)
            rrz=Rzz(mm)-Rzz(nn)
            if(rrx.gt.lx/2) then
			rrx=-lx+rrx
			else if(rrx.lt.-lx/2) then
			rrx=lx+rrx
			end if
		    if(rry.gt.ly/2) then
		    rry=-ly+rry
		   	else if(rry.lt.-ly/2) then
			rry=ly+rry
	    	end if
		    if(rrz.gt.lz/2) then      
			rrz=-lz+rrz
			else if(rrz.lt.-lz/2) then
			rrz=lz+rrz
		    end if
            dsrg2B(i)=dsrg2B(i)+rrx*rrx+rry*rry+rrz*rrz
         else
            rrx=Rxx(mm)-Rxx(nn)
            rry=Ryy(mm)-Ryy(nn)
            rrz=Rzz(mm)-Rzz(nn)
            if(rrx.gt.lx/2) then
			rrx=-lx+rrx
			else if(rrx.lt.-lx/2) then
			rrx=lx+rrx
			end if
		    if(rry.gt.ly/2) then
		    rry=-ly+rry
		   	else if(rry.lt.-ly/2) then
			rry=ly+rry
	    	end if
		    if(rrz.gt.lz/2) then      
			rrz=-lz+rrz
			else if(rrz.lt.-lz/2) then
			rrz=lz+rrz
		    end if
            dsrg2AB(i)=dsrg2AB(i)+rrx*rrx+rry*rry+rrz*rrz
          end if
     55 end do
   44 end do     
 33 end do
 write(*,*)dsrg2AB(i), dsrg2B(i), dsrg2A(i)
    rg2A(i)=dsrg2A(i)/(2.d0*N_monomerA*N_monomerA*cluster_chain(i)%last)  !!!!!!!!!!!此处不能除以链子的总个数N_chain,而是cluster_chain(i)%last，应为是一个个的胶束计算的，并且，胶束中的链子数也不同
	rg2B(i)=dsrg2B(i)/(2.d0*N_monomerB*N_monomerB*cluster_chain(i)%last)	
	rg2(i)=(dsrg2A(i)+dsrg2B(i)+dsrg2AB(i))/(2.d0*N_monomer*N_monomer*cluster_chain(i)%last)
	write(*,*)i, rg2A(i),rg2B(i),rg2(i),N_chain
    avrg2A=avrg2A+rg2A(i)
	avrg2B=avRG2B+rg2B(i)	
	avrg2=avrg2+rg2(i)  
22 end do
   avrg2A=avrg2A/N_cluster
   avrg2B=avrg2B/N_cluster
   avrg2=avrg2/N_cluster
write(66,*)avrg2A,avrg2B,avrg2
close(66)
END SUBROUTINE stretch_hydrophobic_chains





!!!!!!!!!!!!!******************************************************************以下是计算单个链子的伸展程度
SUBROUTINE single_chains_stretching
INTEGER :: n,nn,k,m,mm,j,i,N_bin,bin
Double precision :: avrg2A,avrg2B,avrg2,rry,rrx,rrz,dr,dsrg2,dsrg2A,dsrg2B,dsrg2AB,sum,r,p,rg2,rg2A,rg2B
DOUBLE PRECISION, ALLOCATABLE :: num(:)
CALL cluster_domain
ALLOCATE (num(50))
open(66, file = 'single_chains_stretching.dat')
do bin=1,50
   num(bin)=0
enddo
do 22 i=1,N_cluster    	 
  do 33 j=1,cluster_chain(i)%last
        m=cluster_chain(i)%molecules(j)  !!!m表示第i个胶束上第j条链子的编号	
		rg2A=0
        rg2B=0
        rg2=0
		dsrg2=0
	    dsrg2A=0
        dsrg2B=0
        dsrg2AB=0    
	do 44 k=1,N_monomer
	      mm=(m-1)*N_monomer+k   
      do 55 n=1,N_monomer
            nn=(m-1)*N_monomer+n 
         if (k.le.N_monomerA.and.n.le.N_monomerA) then
            rrx=Rxx(mm)-Rxx(nn)
            rry=Ryy(mm)-Ryy(nn)
            rrz=Rzz(mm)-Rzz(nn)
            if(rrx.gt.lx/2) then
			rrx=-lx+rrx
			else if(rrx.lt.-lx/2) then
			rrx=lx+rrx
			end if
		    if(rry.gt.ly/2) then
		    rry=-ly+rry
		   	else if(rry.lt.-ly/2) then
			rry=ly+rry
	    	end if
		    if(rrz.gt.lz/2) then      
			rrz=-lz+rrz
			else if(rrz.lt.-lz/2) then
			rrz=lz+rrz
		    end if
	 	    dsrg2A=dsrg2A+rrx*rrx+rry*rry+rrz*rrz
		
         else if (k.gt.N_monomerA.and. n.gt.N_monomerA) then
            rrx=Rxx(mm)-Rxx(nn)
            rry=Ryy(mm)-Ryy(nn)
            rrz=Rzz(mm)-Rzz(nn)
            if(rrx.gt.lx/2) then
			rrx=-lx+rrx
			else if(rrx.lt.-lx/2) then
			rrx=lx+rrx
			end if
		    if(rry.gt.ly/2) then
		    rry=-ly+rry
		   	else if(rry.lt.-ly/2) then
			rry=ly+rry
	    	end if
		    if(rrz.gt.lz/2) then      
			rrz=-lz+rrz
			else if(rrz.lt.-lz/2) then
			rrz=lz+rrz
		    end if
            dsrg2B=dsrg2B+rrx*rrx+rry*rry+rrz*rrz
         else
            rrx=Rxx(mm)-Rxx(nn)
            rry=Ryy(mm)-Ryy(nn)
            rrz=Rzz(mm)-Rzz(nn)
            if(rrx.gt.lx/2) then
			rrx=-lx+rrx
			else if(rrx.lt.-lx/2) then
			rrx=lx+rrx
			end if
		    if(rry.gt.ly/2) then
		    rry=-ly+rry
		   	else if(rry.lt.-ly/2) then
			rry=ly+rry
	    	end if
		    if(rrz.gt.lz/2) then      
			rrz=-lz+rrz
			else if(rrz.lt.-lz/2) then
			rrz=lz+rrz
		    end if
            dsrg2AB=dsrg2AB+rrx*rrx+rry*rry+rrz*rrz
          end if
     55 end do
   44 end do
     rg2A=dsrg2A/(2.d0*N_monomerA*N_monomerA)
	 rg2B=dsrg2B/(2.d0*N_monomerB*N_monomerB) 
	 rg2=(dsrg2A+dsrg2B+dsrg2AB)/(2.d0*N_monomer*N_monomer)  
	 dr=1.0
	 N_bin=ceiling(18/dr) !!! N_bin此时必须是一个最大的值
	 bin=nint(rg2/dr+0.5)   !!!!此时是计算每条链子的伸展程度 
	 if (bin.le.50)then
	 num(bin)=num(bin)+1
	 end if
 33 end do 
22 end do
      sum=0
   do i=1,N_bin
      sum=sum+num(i)
   end do
	do i=1,N_bin
       r=dr*(dble(i)-0.5) 
       p=num(i)/sum
       write(66,*)r,p
	end do
close(66)
END SUBROUTINE single_chains_stretching

END MODULE cluster_analise