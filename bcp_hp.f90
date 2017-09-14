program main
	use global_parameter
	use system_init
	use motion
	implicit none

	open(unit = 2, file = 'atom.txt	')
	open(unit = 3, file = 'chains.txt')
	open(unit = 4, file = 'annealing.txt')

	open(unit = 7, file ='v61.txt')
	open(unit = 8, file ='v62.txt')
	open(unit = 9, file ='v63.txt')
	open(unit = 10, file ='v64.txt')
	open(unit = 11, file ='v65.txt')

	call interaction()
	call box_building()

end program main
