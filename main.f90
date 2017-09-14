program main
	use global_parameter
	use random_gen
	use system_init
	use motion
	implicit none
!the storge device serial count print*, from number 3
	double precision r0, len_NP ,concentration
	double precision eas,ebs,eaw,ebw
! 	double precision init_config_energy,e
	integer hz,num_NP,len_segmentA,len_segmentB,len_BCP,box_len
	integer	number_atom,number
	_wall,n_atom_free
	open(unit=3, file=' annealing.txt', status="new", action="write")
	print*,'open file annealing.txt'
	
	iseed = -137
! 	allocate(box(lxyz:3))
	eas = 1.0
	ebs = -3.0
	eaw = 3.0
	ebw = -1.0
	call interaction_init(eas,ebs,eaw,ebw)

	r0 = g_confined_r0
	hz = g_cylinder_height
	call box_cylinder_confined_init(r0,hz)
	
	len_NP = 2.5
	num_NP = 18
	box_len = g_box_len
	call near_point_init(len_NP,num_NP,box_len,hz)

	concentration = g_concentration
	len_segmentA = g_length_Ablock
	len_segmentB = g_length_Bblock
	len_BCP = g_length_BCP
	call polymer_chains_init(len_segmentA,len_segmentB,concentration)

	call Marking_empty_particle()

	number_atom = g_num_atom
	number_wall = g_num_wall
	call calculate_config_E(number_atom,num_NP)
	init_config_energy = config_energy
	print*,'init_config_energy=',init_config_energy

	call tmp_character(number_atom,number_wall)
	!****************************************************************!
	!			end the initialization partion						 !
	!****************************************************************!



end program main
