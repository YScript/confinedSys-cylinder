            !   60%
   	Dimension atom(216000,3),as(216000,3),nna(216000,18),nn(216000),icha(216000),ichaP(216000)
	Dimension eab(4,4),nchain(3985,12),ichapp(13),nempty(216000),mc(216000),mw(216000),mne(216000)
	DOUBLE PRECISION E,EP,ep0,epp,de,randNum
		  
	
	INTEGER :: idum, ISEED,zl
 	COMMON /CSEED/ ISEED
	REAL :: ranf
	integer,parameter::typeA = 1,typeB=2,typeS = 3,typeW = 4
	integer::len_A_segment,len_B_segment
	double precision::concentration,r0,r1,dwall
	real(kind=4):: time_tic,time_toc,time_toc1,time_toc2,time_toc3

	! 	concentration=0.34480 float-point exception
! 	concentration = 0.1500
 	concentration = 0.500

	OPEN (UNIT=2,FILE='atom.txt')
	OPEN (UNIT=3,FILE='chain.txt')  
	OPEN (UNIT=4,FILE='energy.txt')
	OPEN (UNIT=7,FILE='v61.txt')
	OPEN (UNIT=8,FILE='v62.txt')
	OPEN (UNIT=9,FILE='v63.txt')
	OPEN (UNIT=10,FILE='v64.txt') 
	OPEN (UNIT=11,FILE='v65.txt')
		 	
	len_A_segment = 10
	len_B_segment = 2
	nnd1=10
	nnd=12

	zl=60			
	r0=12.0
	dwall=2.0
	r1=r0+dwall
	llongth=2*r1+1
		
	nnstep=2 !5000
		
	idum = -191
	ISEED= -137

	do i = 1, 4, 1
		do j = 1, 4, 1
			eab(i,j) = 0.0
		end do
	end do
	eab(typeA,typeB) = 1.0
	eab(typeA,typeS) = -1.0
	eab(typeA,typeW) = -1.0
	eab(typeB,typeS) = 1.0
	eab(typeB,typeW) = 1.0
	do i = 1, 4, 1
		do j = 1, 4, 1
			eab(j,i) = eab(i,j)
		end do
	end do
	as(1,1)=-r1
	as(1,2)=-r1
	as(1,3)=0.
	kkk=0
   	zll=zl/2
	do 1 i1=1,llongth   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		do 2 i2=1,llongth
		   	do 387 i3=1,zl	 
			   	kkk=kkk+1
			   	as(kkk,1)=as(1,1)+(i1-1)
			   	as(kkk,2)=as(1,2)+(i2-1)
			   	as(kkk,3)=as(1,3)+(i3-1)
		   387 continue
		2 continue
	1 continue
	write(*,*) 'kkk=',kkk

	k0=0
	do 3 i=1,kkk
		rr=sqrt(as(i,1)*as(i,1)+as(i,2)*as(i,2))         
		if(rr.le.r0)then
			k0=k0+1	
			atom(k0,1)=as(i,1)
			atom(k0,2)=as(i,2)
			atom(k0,3)=as(i,3)
		else
		endif
	3 continue
	ntot=k0
	k1=k0
	do 4 i=1,kkk
		rr=sqrt(as(i,1)*as(i,1)+as(i,2)*as(i,2))		  
	   	if(rr.gt.r0.and.rr.le.r1)then
			k1=k1+1 
        	atom(k1,1)=as(i,1)
			atom(k1,2)=as(i,2)
			atom(k1,3)=as(i,3)
		else
		endif
	4 continue
	write(*,*) 'k0,k1=',k0,k1

	call cpu_time(time_tic)
	do 5 i=1,ntot
		nn(i)=0
5 	continue
  	do 6 i=1,ntot
    	j=i
		do while(nn(i).lt.18)
	        if(j.lt.ntot)then
				j=j+1		
			    rrx=atom(i,1)-atom(j,1)
			    rry=atom(i,2)-atom(j,2)
				rrz=atom(i,3)-atom(j,3)
				if(rrz.gt.zll)then
				 	rrz=-zl+rrz
				else if(rrz.lt.-zll)then
					rrz=zl+rrz
				else
				endif
		     	rr=rrx*rrx+rry*rry+rrz*rrz
		   		if(rr.lt.2.5)then
		           	nn(i)=nn(i)+1
		           	nn(j)=nn(j)+1
		           	nna(i,nn(i))=j
		           	nna(j,nn(j))=i
		        else
		        endif
		    else
                if(j.lt.k1)then
			        j=j+1		
		           	rrx=atom(i,1)-atom(j,1)
		           	rry=atom(i,2)-atom(j,2)
		           	rrz=atom(i,3)-atom(j,3)
					if(rrz.gt.zll)then
						rrz=-zl+rrz
				   	else if(rrz.lt.-zll)then
						rrz=zl+rrz
	              	else
				  	endif
		     		rr=rrx*rrx+rry*rry+rrz*rrz
		            if(rr.lt.2.5)then
		                nn(i)=nn(i)+1
		                nna(i,nn(i))=j
		  		    else
                    endif
               	else
		        endif			  		  
		    endif
    	enddo
6 	continue
	call cpu_time(time_toc)
	print*,'time on the neighbours:',time_toc - time_tic

write(*,*) 'after0'
 !  pause
    do 7 i1=1,ntot
       icha(i1)=3
    7 continue
    do 8 i2=ntot+1,k1
        icha(i2)=4
    8 continue
    !write(*,*) 'after1'
	k=int(zl/nnd)
	ntotc=int(ntot/zl)*k
	nntotc=int(ntot*concentration/nnd)

  	nndx=nnd-nnd1
  	nndxp=int(nnd/2)
    nx=0
	i=1
	do while(i.lt.ntot)
		!write(*,*) k
		!	pause	
    	j=i
		do 9 ii=1,k		    
			if(icha(j).eq.3)then
				nx=nx+1
				nn0=1
				nchain(nx,nn0)=j
				icha(j)=1
				do while(nn0.ne.nnd1)
					nx2=j+1
				    if(icha(nx2).eq.3)then
			            icha(nx2)=1
						nn0=nn0+1
						nchain(nx,nn0)=nx2
						j=nx2
					else
					endif
				enddo
	    		do while(nn0.ne.nnd)
			    	nx2=j+1
					if(icha(nx2).eq.3)then
						icha(nx2)=2
						nn0=nn0+1
						nchain(nx,nn0)=nx2
						j=nx2
			  		else
				  		write(*,*) 'error'
				 		! pause
				  	endif
				enddo
			else
			endif
			j=j+1
		9 continue
		i=i+zl
	enddo
	!print*,'chains:',nchain(1,1),nchain(1,2),nchain(1,3),nchain(1,4)
!------------------------------------------------
 	ntrnt=ntotc-nntotc
   	write (*,*)ntotc,ntrnt,nntotc
	j=1
do 88 i=1,ntrnt
	nn1=3
	do while(nn1.eq.3)
		j=int(1.0+ntotc*ranf(idum))
		nx=nchain(j,1)
		nn1=icha(nx)
		if(nn1.ne.3)then
			icha(nx)=3
		else
		endif
	enddo
	nn0=1
	do while(nn0.ne.nnd)
		if(nn1.ne.3)then
			nn0=nn0+1
			nx=nchain(j,nn0)
			icha(nx)=3
		else
		endif
	enddo
88 continue

	k=0
   	do 888 i=1,ntotc
	   	nx=nchain(i,1)
	   	nn1=icha(nx)
	   	if(nn1.ne.3)then
	   		k=k+1
	   	else
   		endif
       	do 999 j=1,nnd
		   	if(nn1.ne.3)then
			   	nchain(k,j)=nchain(i,j)
			   	mc(nchain(i,j))=k
			   	mw(nchain(i,j))=j
		   	else
		   	endif
	999 continue
888 continue

   	ntotc=nntotc
   	percent=float(ntotc*nnd)/ntot
 	write(*,*)ntotc,percent 
 	write(*,*)'after enddo 9'	
	
	i=0
    DO  900 kk=1,ntot
     	if(icha(kk).eq.3)then
     		i=i+1
	 		nempty(i)=kk
			mne(kk)=i
 		!    write(3,*)'i,nempty(i)=',i,kk,icha(kk)
     	else
  	 	endif
900	continue
 	ntote=i 
 
 	write(*,*) ntote

	e=0.
	do 118 ii=1,ntot
		do 118 JJP=1,18
			nnii=icha(nna(ii,JJP))
		 	if(ICHA(Ii).lt.nnii)then	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!			energy calculate
  				e=e+eab(icha(ii),nnii)
				!	elseif(icha(i1).eq.nnii)then
				!ee1=eab(icha(ii),nnii)/2 
				!ee1=e+ee1
			else
		  	endif
			e=e !+ee1
118 continue 
   	write(*,*)'e0=',e
	do 18 i1=1,k1
		ichaP(i1)=ICHA(I1)
18	continue 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 nnum=2
    ep0=10000000.
	epp=e
	cons1=30.
	
	call cpu_time(time_tic)
	do 50 jj=1,60
		if(ep0-epp.gt.500.)then
	 		cons1=cons1*.95
	 	else
	  		cons1=cons1*.92
	  	endif
		EP=0.
		do 5002 ii02=1,nnstep
		 	nrrend=0
		 	nmontes=0
			do 5003 ii03=1,ntot  	 
	 			cons=1./cons1
    			ncc=1
	  
	    		nnx=INT(1.0+ntotc*ranf(idum))
		
				do ii=1,nndx
					nndxx=ii+nnd1
					ichap(nchain(nnx,ii))=2
					ichap(nchain(nnx,nndxx))=1
		        enddo
				de=0.
				do 13 ii=1,nndx
					nndxx=ii+nnd1
					do 14 kk=1,18
						mii=nna(nchain(nnx,ii),kk)
						mii1=nna(nchain(nnx,nndxx),kk)
						de=de+eab(2,ichap(mii))-eab(1,icha(mii))		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!			energy calculate
				   		de=de+eab(1,ichap(mii1))-eab(2,icha(mii1))
				14   continue		
		  13   	continue

  				if(de.le.0)then
	     			e=e+de	
					do 122 ii=1,nndx
							nndxx=ii+nnd1
							icha(nchain(nnx,ii))=2
							icha(nchain(nnx,nndxx))=1
					122  continue		
					do 111 ii=1,nndx
							nndxx=ii+nnd1
							nchax=nchain(nnx,ii)
							nchain(nnx,ii)=nchain(nnx,nnd+1-ii)
							nchain(nnx,nnd+1-ii)=nchax
							mw(nchain(nnx,nnd+1-ii))=nnd+1-ii
							mw(nchain(nnx,ii))=ii
					111 continue		
					do 133 ii=nndx+1,nndxp
							nchax=nchain(nnx,ii)
							nchain(nnx,ii)=nchain(nnx,nnd-ii+1)
							nchain(nnx,nnd-ii+1)=nchax
							mw(nchain(nnx,nnd+1-ii))=nnd+1-ii
				     		mw(nchain(nnx,ii))=ii
					133 continue		
		  			nmontes=nmontes+1
				else
	  				pp=exp(-de*cons)
				 	if(pp.ge.ranf(idum))then
				   		e=e+de	
						do 16 ii=1,nndx
							nndxx=ii+nnd1
							icha(nchain(nnx,ii))=2
							icha(nchain(nnx,nndxx))=1
					16  continue		
						do 17 ii=1,nndx
							nndxx=ii+nnd1
							nchax=nchain(nnx,ii)
							nchain(nnx,ii)=nchain(nnx,nnd+1-ii)
							nchain(nnx,nnd+1-ii)=nchax
							mw(nchain(nnx,nnd+1-ii))=nnd+1-ii
							mw(nchain(nnx,ii))=ii
					17  continue		
				  		do 134 ii=nndx+1,nndxp
							nchax=nchain(nnx,ii)
							nchain(nnx,ii)=nchain(nnx,nnd-ii+1)
							nchain(nnx,nnd-ii+1)=nchax
							mw(nchain(nnx,nnd+1-ii))=nnd+1-ii
							mw(nchain(nnx,ii))=ii
				  	134 continue		
	  				nmontes=nmontes+1
   					else
						do ii=1,nndx
							nndxx=ii+nnd1
							ichap(nchain(nnx,ii))=1
							ichap(nchain(nnx,nndxx))=2
				        enddo
    				endif
				endif
		 		EP=EP+e
				nrrend=nrrend+1
        
 	!! Chain movement  only need this part------------------------------------------------------------------------------------------------------------------------------
				randNum=ranf(idum)
				if (randNum.ge.0.0.and.randNum.lt.0.50)then
	 				nnx=INT(1.0+ntotc*ranf(idum))
			 		nxx=INT(1.0+nnd*ranf(idum))
				  	nx=nchain(nnx,nxx)
		        	nx1=INT(1.0+18*ranf(idum))
				  	nx2=nna(nx,nx1)

 
  					if((icha(nx2).eq.3) .or. (icha(nx2).eq.7))then
						nny=mne(nx2)
  						ichg=1
!write(*,*)'nmontes=',nmontes,nxx,nx,nx1,nx2
					   	if (nxx.eq.1)then
						   	nstar=nxx+1
						    nend=nstar
						else if(nxx.eq.nnd)then
						   	nstar=nxx-1
						    nend=nstar
					   	else
						   	nstar=nxx-1
						    nend=nxx+1
					   	endif
 
						do ii=nstar,nend,nnum
							nni=nchain(nnx,ii)
							rrx=atom(nni,1)-atom(nx2,1)
							rry=atom(nni,2)-atom(nx2,2)
							rrz=atom(nni,3)-atom(nx2,3)
							rr=rrx*rrx+rry*rry+rrz*rrz
	
							if(rr.gt.2.5)then
							   	ichg=ichg+1
						      	nii=ii
							  	njj1=nni
							else
							endif
						enddo
	  					njj=nx

				 		if(ichg.eq.2)then
							if(nii.lt.nxx)then
								mm1=-1
							else
								mm1=1
							endif

							do while(ichg.eq.2)
								nii=nii+Mm1
								 	  !if 3
								if(nii.ge.1.and.nii.le.nnd)then 
									nni=nchain(nnx,nii)
									rrx=atom(nni,1)-atom(njj,1)
									rry=atom(nni,2)-atom(njj,2)
									rrz=atom(nni,3)-atom(njj,3)
							 
									rr=rrx*rrx+rry*rry+rrz*rrz
									if(rr.gt.2.5)then
								      	njj=njj1
									  	njj1=nni

							     	else
									   	ichg=1
								    	njj=njj1
							   		endif
							    else
									ichg=1
									njj=njj1
								endif
								ncc=ncc+1
							enddo
						else
						endif ! end if ichg .eq.2
						if(ichg.eq.1)then
							de=0.
						  	if(ncc.eq.1)then
						    	ichanx=icha(nx)
								ichanx2=icha(nx2)
								do ii=1,18
									iia=nna(nx,ii)
									iib=nna(nx2,ii)
									if(iia.ne.nx2)then   
					 					ichaa=icha(iia)
										de=de+eab(ICHANX2,ichaa)-eab(ichanx,ichaa)	!energy calculate
									else
									endif

									if(iib.ne.nx)then   
							     	ichab=icha(iib)
							 
									de=de+eab(ichanx,ichab)-eab(ICHANX2,ichab) !	energy calculate
									else
									endif
								enddo
  								if(de.le.0)then
									e=e+de	
							   		icha(nx2)=ichanx
							        icha(nx)=ICHANX2
							 		nchain(nnx,nxx)=nx2
							 		mc(nx2)=nnx
							 		mw(nx2)=nxx
							 		nempty(nny)=nx
							  		mne(nx)=nny
							   		ichap(nx2)=ichanx
							   		ichap(nx)=ICHANX2
							  		nmontes=nmontes+ncc
						  		else
		  							pp=exp(-de*cons)
	   
									if(pp.ge.ranf(idum))then
										e=e+de	
										icha(nx2)=ichanx
										icha(nx)=ICHANX2
										ichap(nx2)=ichanx
									   	ichap(nx)=ICHANX2
									 	nchain(nnx,nxx)=nx2
									 	mc(nx2)=nnx
									 	mw(nx2)=nxx
									  	nempty(nny)=nx
									  	mne(nx)=nny
									 	nmontes=nmontes+1

								    else 
								    endif !enf if pp.ge.ranf(idum)
								endif
							else !	else ->ncc.eq.1
								NCC1=NCC-1
								NCC2=NCC+1
								nx1=nx
								ichnx1=ichap(nx1)
								ichapp(1)=nx2

								DO  JCC=2,NCC2
									ichapp(jcc)=nx1
									if(jcc.ne.ncc2)then
										mmp=nxx+(jcc-1)*mm1
										nx1=NCHAIN(NNX,mmP)
									else
									endif
									enddo
									DO JCC=1,NCC
										jcc1=jcc+1
										njcc=ichapp(jcc)
										njcc1=ichapp(jcc1)
										ichap(njcc)=icha(njcc1)
									enddo
									ichap(ichapp(ncc2))=ICHA(NX2)
									DO  JCC=1,NCC
										jcc1=jcc+1
										njcc=ichapp(jcc)
										njcc1=ichapp(jcc1)
										de=de-eab(ichAP(NJCC),ichAP(NJCC1))+eab(ichA(NJCC),ichA(NJCC1))				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!			energy calculate
									ENDDO

									DO JCC=1,NCC2
										nccp=ichapp(jcc)
										ichnccp=icha(nccp)
										ichncc=ichap(nccp)
										do  ii=1,18
										  	iia=nna(nccp,ii)	!; write(100,*)nccp,ii,iia
										   	ichaap=icha(iia)   
										    ichaa=ichap(iia)  ! ; write(100,*)ichncc,ichaa,ichnccp,ichaap
											de=de+eab(ichncc,ichaa)-eab(ichnccp,ichaap)								!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!			energy calculate
										enddo
									enddo

									DO JCC=1,NCC
										nccp=ichapp(jcc)
										do  ii=1,18
										  	iia=nna(nccp,ii)
											DO IIP=JCC+2,NCC2
										  		NCCPP=ichapp(IIP)
										   		IF(IIA.EQ.NCCPP)THEN	
										   			de=de-eab(ichAP(NCCP),ichAP(NCCPP))+Eab(ichA(NCCP),ichA(NCCPP))				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!			energy calculate
										   		ELSE
									   			ENDIF
									   		ENDDO
										ENDDO
									ENDDO

									if(de.le.0)then
									   e=e+de	
										DO  JCC=1,NCC
											nccp=ichapp(jcc)
											icha(nccp)=ichap(nccp)
											mmp=nxx+(jcc-1)*mm1
											nchain(nnx,mmp)=nccp
											mc(nccp)=nnx
											mw(nccp)=mmp
										enddo
										nccp=ichapp(ncc2)
										icha(nccp)=ichap(nccp)
										nempty(nny)=nccp
										mne(nccp)=nny
		  								nmontes=nmontes+1
								  	else
									  	pp=exp(-de*cons)
								 		if(pp.gt.ranf(idum))then
								   			e=e+de	
											DO  JCC=1,NCC
												nccp=ichapp(jcc)
												icha(nccp)=ichap(nccp)
												mmp=nxx+(jcc-1)*mm1
												nchain(nnx,mmp)=nccp
												mc(nccp)=nnx
												mw(nccp)=mmp
											enddo
											nccp=ichapp(ncc2)
											icha(nccp)=ichap(nccp)
											nempty(nny)=nccp
											mne(nccp)=nny
										   	nmontes=nmontes+1
	    								else
											DO  JCC=1,NCC2
												njcc=ichapp(jcc)
												ichap(njcc)=icha(njcc)
											enddo
	    								endif
									endif
								continue

							endif
	
						else
						endif

					else
					endif	
		
!-------------------------------------			
!---------------------------------------------------------------------------------------------------------------------------------------------------------------------
				else if (randNum.ge.0.50.and.randNum.lt.1.0)then 
					! solvent movement
		   			ncc=1
		         	nny=INT(1.0+ntote*ranf(idum))
				  	nx2=nempty(nny)
				  	nx1=int(1.0+18*ranf(idum))
				  	nx=nna(nx2,nx1)
 
 					if(icha(nx).lt.3)then
			  		  	nnx=mc(nx)
					  	nxx=mw(nx)
  						ichg=1

					   	if (nxx-1.eq.0)then
						   	nstar=nxx+1
						    nend=nstar
						else if(nxx.eq.nnd)then
						   	nstar=nxx-1
						    nend=nstar
					   	else
						   	nstar=nxx-1
						    nend=nxx+1
					   	endif
					 
	do 10 ii=nstar,nend,nnum

	nni=nchain(nnx,ii)
!		write(*,*)'ii,nni=',ii,nni
		rrx=atom(nni,1)-atom(nx2,1)
		rry=atom(nni,2)-atom(nx2,2)
		rrz=atom(nni,3)-atom(nx2,3)
 
	rr=rrx*rrx+rry*rry+rrz*rrz
	
	if(rr.gt.2.5)then
	   ichg=ichg+1
      nii=ii
	  njj1=nni
	  else
	endif

10 continue
	  njj=nx

! write(2,*)'nx2=',nx2

!	mm1=0
	 !do2
  if(ichg.eq.2)then
	 if(nii.lt.nxx)then
	mm1=-1
	else
	mm1=1
	endif

	 do while(ichg.eq.2)
	nii=nii+Mm1
	 	  !if 3
	if(nii.ge.1.and.nii.le.nnd)then 
		nni=nchain(nnx,nii)
!		write(*,*)'ii,nni=',ii,nni
		rrx=atom(nni,1)-atom(njj,1)
		rry=atom(nni,2)-atom(njj,2)
		rrz=atom(nni,3)-atom(njj,3)
 
	rr=rrx*rrx+rry*rry+rrz*rrz
		if(rr.gt.2.5)then
      njj=njj1
	  njj1=nni

     	   else
	   ichg=1
    	njj=njj1
      		endif
    	else
	ichg=1
	njj=njj1
        
	endif
!endif 3
ncc=ncc+1
	enddo
else
endif
  !if2,enddo2
!write(*,*)'jj=',jj,'nrrend=',nrrend,'nmontes=',nmontes,ncc,de
	if(ichg.eq.1)then
      !if 1
	de=0.
  if(ncc.eq.1)then
    	ichanx=icha(nx)
ichanx2=icha(nx2)
		do 15 ii=1,18
    iia=nna(nx,ii)
	iib=nna(nx2,ii)
 

if(iia.ne.nx2)then   
   ichaa=icha(iia)
    
de=de+eab(ICHANX2,ichaa)-eab(ichanx,ichaa)			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!			energy calculate
else
endif

if(iib.ne.nx)then   
     ichab=icha(iib)
 
de=de+eab(ichanx,ichab)-eab(ICHANX2,ichab)				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!			energy calculate
else
endif
15  continue
!if 0
  if(de.le.0)then
  	
		 e=e+de	

	   icha(nx2)=ichanx
	        icha(nx)=ICHANX2
	 nchain(nnx,nxx)=nx2
	 mc(nx2)=nnx
	 mw(nx2)=nxx
	 nempty(nny)=nx
	 
	  mne(nx)=nny
	 	 

	   ichap(nx2)=ichanx
	   ichap(nx)=ICHANX2

	  nmontes=nmontes+1
	  else
	  pp=exp(-de*cons)
   
	 if(pp.ge.ranf(idum))then
	  e=e+de	

	 icha(nx2)=ichanx
	        icha(nx)=ICHANX2
ichap(nx2)=ichanx
	   ichap(nx)=ICHANX2
	 nchain(nnx,nxx)=nx2
	 mc(nx2)=nnx
	 mw(nx2)=nxx
	 nempty(nny)=nx
mne(nx)=nny

	  nmontes=nmontes+1

    else
    endif
endif
!endif0
else 
!write(*,*)'jj=',jj,'nrrend=',nrrend,'nmontes=',nmontes,ncc,de
NCC1=NCC-1
NCC2=NCC+1
nx1=nx
ichnx1=ichap(nx1)
ichapp(1)=nx2

 DO 28 JCC=2,NCC2
ichapp(jcc)=nx1
if(jcc.ne.ncc2)then
mmp=nxx+(jcc-1)*mm1
 nx1=NCHAIN(NNX,mmP)

else
endif
28 continue

 DO 38 JCC=1,NCC
 jcc1=jcc+1
njcc=ichapp(jcc)
njcc1=ichapp(jcc1)
ichap(njcc)=icha(njcc1)
38 continue
ichap(ichapp(ncc2))=ICHA(NX2)
DO  JCC=1,NCC
 jcc1=jcc+1
njcc=ichapp(jcc)
njcc1=ichapp(jcc1)
de=de-eab(ichAP(NJCC),ichAP(NJCC1))+eab(ichA(NJCC),ichA(NJCC1))				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!			energy calculate
  ENDDO
DO 48 JCC=1,NCC2
nccp=ichapp(jcc)
ichnccp=icha(nccp)
ichncc=ichap(nccp)
do 115 ii=1,18
  iia=nna(nccp,ii)
   ichaap=icha(iia)   
    ichaa=ichap(iia)   
de=de+eab(ichncc,ichaa)-eab(ichnccp,ichaap)				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!			energy calculate
 115  continue

48  continue

DO JCC=1,NCC
nccp=ichapp(jcc)
do  ii=1,18
  iia=nna(nccp,ii)
DO IIP=JCC+2,NCC2
  NCCPP=ichapp(IIP)
   IF(IIA.EQ.NCCPP)THEN	
   de=de-eab(ichAP(NCCP),ichAP(NCCPP))+Eab(ichA(NCCP),ichA(NCCPP))			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!			energy calculate
   ELSE
   ENDIF
   ENDDO
ENDDO
ENDDO

 

!if 00
  if(de.le.0)then
 
   e=e+de	

DO 58 JCC=1,NCC
nccp=ichapp(jcc)
icha(nccp)=ichap(nccp)
mmp=nxx+(jcc-1)*mm1
nchain(nnx,mmp)=nccp
mc(nccp)=nnx
mw(nccp)=mmp
58 continue
nccp=ichapp(ncc2)
icha(nccp)=ichap(nccp)
nempty(nny)=nccp
	 
  mne(nccp)=nny
	  nmontes=nmontes+1

	  else
	  pp=exp(-de*cons)
 

	 if(pp.gt.ranf(idum))then
	    e=e+de	

		DO 68 JCC=1,NCC
nccp=ichapp(jcc)
icha(nccp)=ichap(nccp)
mmp=nxx+(jcc-1)*mm1
nchain(nnx,mmp)=nccp
mc(nccp)=nnx
mw(nccp)=mmp
68 continue
nccp=ichapp(ncc2)
icha(nccp)=ichap(nccp)
nempty(nny)=nccp
 mne(nccp)=nny
	   nmontes=nmontes+1


    else
 DO 78 JCC=1,NCC2
njcc=ichapp(jcc)
ichap(njcc)=icha(njcc)
78 continue
!write(*,*)'jj=',jj,'nrrend=',nrrend,'nmontes=',nmontes,ncc,de
    endif
 
endif
!endif00

endif
!endif1
else
endif
!endif2
!write(*,*)'jj=',jj,'nrrend=',nrrend,'nmontes=',nmontes,ncc,de

!!!!!!!!!!!!!!!!

else
endif
!end
 EP=EP+E

	 nrrend=nrrend+1
else
endif
!End solvent movement
!write(*,*)'jj=',jj,'nrrend=',nrrend,'nmontes=',nmontes

!-----------------------------------	

	


5003 continue
5002 continue 
!endo1
EP=EP/ntot
EP=EP/nnstep
write(*,*)'jj=',jj,'nrrend=',nrrend,'nmontes=',nmontes,cons1,cons,EP
write(4,*) jj,nmontes,cons,EP
  	
  	ep0 = epp
  	epp=ep


   if(jj.eq.30)then
   	 do 801 kkk=1,ntotc    	
	   write(7,*)(nchain(kkk,j),j=1,nnd)
	  801 continue
        else if(jj.eq.40)then
	 do 802 kkk=1,ntotc    	
	  write(8,*)(nchain(kkk,j),j=1,nnd)
	  802 continue	
	           else if(jj.eq.45)then
			   do 803 kkk=1,ntotc    	
	           write(9,*)(nchain(kkk,j),j=1,nnd)
	            803 continue
				   else if(jj.eq.50)then
				    do 804 kkk=1,ntotc    	
	           write(10,*)(nchain(kkk,j),j=1,nnd)
	            804 continue
				         else if(jj.eq.55)then
						 do 805 kkk=1,ntotc    	
	           write(11,*)(nchain(kkk,j),j=1,nnd)
	            805 continue
          endif
			  
 
 50 continue
   ! write(*,*) ntotc
!	pause
       do 9091 i=1,ntotc
		 do 9092 j=1,nnd-1
		 rrx=atom(nchain(i,j),1)-atom(nchain(i,j+1),1)
		 rry=atom(nchain(i,j),2)-atom(nchain(i,j+1),2)
		 rrz=atom(nchain(i,j),3)-atom(nchain(i,j+1),3)
			 if(rrz.gt.zll)then
			 rrz=-zl+rrz
			    else if(rrz.lt.-zll)then
				rrz=zl+rrz
            else
			  endif
		     rr=rrx*rrx+rry*rry+rrz*rrz
		 if(rr.gt.2.5)then
		 write(*,*) 'error'
	!	 pause
		 else
		 endif
		 9092 continue
       9091 continue

2222     	do 11 kkk=1,k1    	
	     write(2,*)atom(kkk,1),atom(kkk,2),atom(kkk,3),icha(kkk),ichap(kkk)
        11	   continue
		 do 12 kkk=1,ntotc    	
	   write(3,*)(nchain(kkk,j),j=1,nnd)	 
           12   continue
		
write(*,*)'e=',e
   
    e=0.
   	do 128 ii=1,ntot
	do 128 JJP=1,18
	nnii=icha(nna(ii,JJP))
 if(ICHA(Ii).lt.nnii)then
 
  e=e+eab(icha(ii),nnii)				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!			energy calculate
  else
  endif

  128   continue
 
  write(*,*)'ep=',e 
10000 print*,'++'
1111	   stop 
 end


		
		FUNCTION ranf(idum)
IMPLICIT NONE
INTEGER, PARAMETER ::  K4B=selected_int_kind(9)
INTEGER(K4B), INTENT(INOUT) :: idum
REAL :: ranf
INTEGER(K4B), PARAMETER :: IA=16807, IM=2147483647, IQ=127773, IR=2836
REAL, SAVE :: am
INTEGER(K4B), SAVE :: ix=-1, iy=-1, k
if (idum <= 0 .or. iy < 0) then					! Initialize 
	am=nearest(1.0,-1.0)/IM
	iy=ior(ieor(888889999, abs(idum)),1)
	ix=ieor(777755555, abs(idum))
	idum=abs(idum)+1							! Set idum positive
end if
ix=ieor(ix, ishft(ix, 13))						! Marsaglia shift sequence with period 2^^32-1
ix=ieor(ix, ishft(ix, -17))
ix=ieor(ix, ishft(ix, 5))
k=iy/IQ											! Park-Miller sequence by Schrage's method, period 2^^31-1
iy=IA*(iy-k*IQ)-IR*k
if (iy < 0 ) iy=iy+IM
ranf=am*ior(iand(IM, ieor(ix,iy)),1)			! Combine the two generators with masking to ensure nonzero value

END FUNCTION ranf
		
		
		
		
		
	
