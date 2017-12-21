program main

	real:: a
	call subprog_1()
	write(*,*)a

end program main


subroutine subprog_1()
	implicit none
	integer::a = 1
	call subprog_1_1()

end subroutine subprog_1


subroutine subprog_1_1()
	implicit none
	print*,'hello world'
! 	print*,a

end subroutine subprog_1_1