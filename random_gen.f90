module random
	use  global_parameter
	implicit none

	contains

	double precision function ranf( )


        if (iseed <= 0 .or. iy < 0) then                                        ! Initialize
           am=nearest(1.0,-1.0)/IM
           iy=ior(ieor(888889999, abs(iseed)),1)
           ix=ieor(777755555, abs(iseed))
           iseed=abs(iseed)+1                                                   ! Set idum positive
	    end if
	    ix=ieor(ix, ishft(ix, 13))                                                ! Marsaglia shift sequence with period 2^^32-1
	    ix=ieor(ix, ishft(ix, -17))
	    ix=ieor(ix, ishft(ix, 5))
	    kx=iy/IQ                                                                                   ! Park-Miller sequence by Schrage's method, period 2^^31-1
	    iy=IA*(iy-kx*IQ)-IR*kx

	    if (iy < 0 ) iy=iy+IM
	    ranf=am*ior(iand(IM, ieor(ix,iy)),1)                          ! Combine the two generators with masking to ensure nonzero value
    
	end function 
end module random

!
! readme:

!this module with a function used to generator random number,you do not need to gasp ,but 
! you should know how to get a return value in function; 
!caution the format of function arguments and dereference; 

!********************************************************************************************************
! From Chap. B7 in Numerical Recipes F90.
! "Minimal" random number generator of Park and Miller combined with a Marsaglia shift sequence.
! Returns a uniform random deviate between 0.0 and 1.0 (exclusive of the end point values).
! Call with idum a negative integer to initialize; thereafter, do not alter idum except to reinitialize.
! The period of this generotor is about 3.1*10^18.
! The random number sequence is controlled by idum, ixran, iyran and kran.
!*********************************************************************************************************
