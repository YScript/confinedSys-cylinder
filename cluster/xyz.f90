MODULE xyz
USE global_para

IMPLICIT NONE


CONTAINS

SUBROUTINE xyz_in(istep)

integer :: istep
character  CFGFILE1*128, CFGFILE2*128
INTEGER ::  nchain(N_chain,N_monomer)
INTEGER :: kk,gg,i,kkk, i1, i2 ,i3,k,j,n
character  st*10
INTEGER, ALLOCATABLE :: atom(:,:)
integer :: buf_type
character, dimension(:), allocatable :: type_monomer_temp

ALLOCATE(type_monomer_temp(Ntot))

write(st,'(i10)')istep
ALLOCATE  (atom(Ntot,3))
CFGFILE1 = trim(adjustl(st))//'.snap'
CFGFILE2 = trim(adjustl(st))//'.chain'

 open(unit=1,FILE=CFGFILE1)
 open(unit=2,FILE=CFGFILE2)

    Rxx(1)=0
    Ryy(1)=0
    Rzz(1)=0
 kkk=0
do i1=1,Lx
  do i2=1,Ly
    do i3=1,Lz
	  kkk=kkk+1
	  Rxx(kkk)=Rxx(1)+(i1-1)
	  Ryy(kkk)=Ryy(1)+(i2-1)
      Rzz(kkk)=Rzz(1)+(i3-1)
	enddo
  enddo
enddo



do kkk=1,Ntot
  read(1,*) atom(kkk,1), atom(kkk,2), atom(kkk,3), buf_type
  if(buf_type==1)then
   type_monomer_temp(kkk) = 'C'
  elseif(buf_type==2)then
   type_monomer_temp(kkk) = 'O'
  elseif(buf_type==3)then
   type_monomer_temp(kkk) = 'W'
  else
    write(*,*)'monomer type error!'
  endif
end do

do kk = 1, N_chain
		read(2,*) (Nchain(kk,gg),gg=1,N_monomer)
end do

 n=0    
do i = 1, N_chain
  do j=1,N_monomer
     n=n+1          !!! n表示新格点的编号
     k=nchain(i,j)  !!! k表示旧格点的编号
	 Rxx(n)=atom(k,1)
     Ryy(n)=atom(k,2)
	 Rzz(n)=atom(k,3)
     type_monomer(n)=type_monomer_temp(k) !!!对格点的属性重新赋值

  enddo 
enddo
!do n=1,24
!write(*,*) n,Rxx(n),Ryy(n),Rzz(n), type_monomer(n)
!end do

END SUBROUTINE xyz_in

SUBROUTINE xyz_out(istep)
integer :: i,istep
character  CFGFILE*128
character  st*10


call xyz_in(istep)

write(st,'(i10)')istep
CFGFILE = trim(adjustl(st))//'.xyz'

open(11,file=CFGFILE)
write(11,*)N_monomerT
write(11,*)
do i = 1, N_monomerT
 write(11,*)type_monomer(i), sngl(Rxx(i)), sngl(Ryy(i)), sngl(Rzz(i))
enddo 
close(11)

END SUBROUTINE xyz_out   





END MODULE xyz



