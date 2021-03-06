!	Author: Fengyuan
!	Institution: Physics institution of NanKai Unversity; Room -106;
!	E-mail:	fengy104@mail.nankai.edu.cn

module globalparameters
	implicit none
	integer,parameter:: long_int = SELECTED_INT_KIND(8)
	integer,parameter:: dble_real = SELECTED_REAL_KIND(8)
	integer(kind=8)::i,j,k
	real(kind=8),parameter::g_concentration = 0.600
	real(kind=4),parameter::g_cylinder_r0 = 6.0,g_wall_width = 3.0
	integer,parameter::g_len_ASegment = 10,g_len_BSegment = 2
	real(kind=8),parameter::g_eas = -0.5,g_eaw = -3.0,g_ebs = 0.5,g_ebw= 3.0
	real(kind=8),parameter::g_init_temp = 30
	integer(kind=4),parameter::g_annealing_steps = 60
	integer(kind=8),parameter::g_standard_MCS = 2 !5000

	integer(kind=4),parameter::g_lx = 60, g_ly = 60, g_lz = 60
	integer(kind=4),parameter::g_lxyz = g_lx*g_ly*g_lz
	integer(kind=4)::g_num_box,g_num_atom,g_total_parti,g_num_chains,g_num_solvents
	real(kind=4),parameter::g_len_near_points = 2.5
	integer(kind=4),parameter::g_num_neighbours = 18
	integer,parameter::g_typeA = 1,g_typeB = 2,g_typeS=3,g_typeW= 4
	integer(kind=8)::g_num_accepted_steps,g_num_trial_steps
	integer(kind=8)::iseed
	real(kind=8)::g_config_e,g_sum_var_e,g_sum_confe_PerTemper

	real(kind=8),dimension(:,:)::eab(4,4)
	integer(kind=4),allocatable,dimension(:,:)::box,atom,id_of_neigh,nchain
	integer(kind=4),dimension(:,:)::neighbours(g_num_neighbours,3)
	integer(kind=4),allocatable,dimension(:)::postn2chnOrder,postn2MrOrder
	integer(kind=4),allocatable,dimension(:)::num_neigh,icha,ichatmp,solvents2Position,position2Solvents
end module globalparameters

program main
	use globalparameters
	implicit none
	real(kind=4)::eas,eaw,ebs,ebw
	real(kind=8)::concentration,config_e,init_config_e
	real(kind=8)::time_tic,time_toc0,time_toc1,time_toc2,time_toc3
	integer(kind=4)::lx,ly,lz,len_ASegment,len_BSegment
	
	real(kind=8)::configE_calculating
	real(kind=8)::ranf	

	call cpu_time(time_tic)
	iseed = -77
	eas = g_eas
	eaw = g_eaw
	ebs = g_ebs
	ebw = g_ebw

	call interaction_init(eas,eaw,ebs,ebw)

	call latticePoints_init(g_lx,g_ly,g_lz)

	call Parti_coor_init()
	print*,'config_e after lattice initialization:',configE_calculating()

	concentration = g_concentration
	len_ASegment = g_len_ASegment
	len_BSegment = g_len_BSegment
	call AB_diblock_copolymer_init(len_ASegment,len_BSegment)
	print*,'config_e after arrange chains:',configE_calculating()

	call export2cc1_Figure(0)

	call chains_remove_pro(len_ASegment,len_BSegment,concentration)
	config_e = configE_calculating()
	init_config_e = configE_calculating()
	print*,'config_e after remove chains/-the initial config_e:',config_e,init_config_e
	
	call count_num_solvents()
	call tempType_init()
	call cpu_time(time_toc0)
	print*,'The total time on the initialization process:',time_toc0 - time_tic
	!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
	!	the following partion is the simulated Annealing partion    !
	print*,'_____________________________________________________________________________________________'
	call simulatedAnnealing_plan(init_config_e,g_init_temp,g_annealing_steps)
	
	call cpu_time(time_toc1)
	print*,'The total time on the program main:',time_toc1- time_tic
	stop "end and exit"
end program main
!_____________________________________________________________________________________________
!_____________________________________________________________________________________________
!_____________________________________________________________________________________________
!_____________________________________________________________________________________________

subroutine simulatedAnnealing_plan(init_config_e,init_temp,annealing_steps)
	!	the plan of simulated Annealing
	use globalparameters
	implicit none
	real(kind=8),intent(in)::init_temp,init_config_e
	integer(kind=4),intent(in)::annealing_steps
	integer(kind=4):: step
	real(kind=8)::annealing_temp,annealing_factor
	integer(kind=8)::mcs,order,standard_MCS,loopNum,num_trial_steps,num_accepted_steps
	real(kind=8)::configE_average_before,configE_average_current,ep
	real(kind=8)::current_configE,changed_configE
	real(kind=8)::configE_calculating

	annealing_temp = init_temp
	standard_MCS = g_standard_MCS
	loopNum = g_num_atom
	configE_average_before = 10000000.0
	!___________________________________________
	current_configE = init_config_e
	changed_configE = 0.0
	g_config_e = init_config_e
	configE_average_current = init_config_e
	print*,'test floating-point exceptions'
	
	!-----------------------------------------main_annealing_plan------------------!
	main_annealing_plan:do step = 1, annealing_steps, 1		
		if ( (configE_average_before - configE_average_current) .gt. 500 ) then
			annealing_temp = annealing_temp * 0.95
		else
			annealing_temp = annealing_temp * 0.92
		endif

		annealing_factor = 1.0/ annealing_temp
		g_num_trial_steps = 0
		g_num_accepted_steps = 0
		ep = 0.0
		g_sum_var_e = 0.0
		g_sum_confe_PerTemper = 0.00

		standard_MCS_loop:do mcs = 1, standard_MCS, 1
			g_num_atom_loop:do order = 1, loopNum, 1

				current_configE = g_config_e
				call chains_reversing(annealing_factor,current_configE)
				g_num_trial_steps = g_num_trial_steps + 1

				current_configE = g_config_e
				call chains_snake(annealing_factor,current_configE,g_concentration)
				g_num_trial_steps = g_num_trial_steps + 1

			end do g_num_atom_loop
		end do standard_MCS_loop

		call figure(step)

		configE_average_before = configE_average_current
 		!configE_average_current = g_sum_confe_PerTemper/(g_num_accepted_steps)
		configE_average_current = g_sum_confe_PerTemper/(g_num_trial_steps)

 		write(*,*)step,g_num_accepted_steps,g_num_trial_steps,g_sum_var_e,g_config_e
 		
 		
	end do	main_annealing_plan
	print*,'last g_config_e',configE_calculating()
end subroutine simulatedAnnealing_plan

subroutine chains_reversing(annealing_factor,init_config_e)
	!	chains reversig motion;
	!	delete the dummy variables of deta_config_e, num_trial_steps,num_accepted_steps
	!	set them as a globalparameters;
	use globalparameters
	implicit none
	real(kind=8),intent(in)::init_config_e,annealing_factor
	integer(kind=4)::moveChains_chosen
	integer(kind=4)::half_copolymer,len_ASegment,len_BSegment,len_copolymer
	integer(kind=4)::atom_X_position,atom_Y_position,neighs_of_A,neighs_of_B
	real(kind=8)::possibility,config_e,deta_config_e
	real(kind=8)::ranf

	config_e = init_config_e
	len_ASegment = g_len_ASegment
	len_BSegment = g_len_BSegment
	len_copolymer = len_ASegment + len_BSegment
	half_copolymer = int(len_copolymer/2)
	moveChains_chosen = int(1.0+g_num_chains*ranf())

	do i = 1, len_BSegment, 1
		atom_X_position = i
		atom_Y_position = i + len_ASegment
		ichatmp(nchain(moveChains_chosen,atom_X_position)) = g_typeB
		ichatmp(nchain(moveChains_chosen,atom_Y_position)) = g_typeA
	end do
	deta_config_e = 0.000
	outer:do i = 1, len_BSegment, 1
		atom_X_position = i
		atom_Y_position = i + len_ASegment
		
		inner:do j = 1, g_num_neighbours, 1
			neighs_of_A = id_of_neigh(nchain(moveChains_chosen,atom_X_position),j)
			neighs_of_B = id_of_neigh(nchain(moveChains_chosen,atom_Y_position),j)
			deta_config_e = deta_config_e+ eab(g_typeB,ichatmp(neighs_of_A)) -eab(g_typeA,icha(neighs_of_A))
			deta_config_e = deta_config_e+ eab(g_typeA,ichatmp(neighs_of_B)) -eab(g_typeB,icha(neighs_of_B))
			!	record the neighbours of atom A as 'neighs_of_A',so as to atom B;
			!	record the changed config_e of remove an A by B as a assume;
			!	so do it by remove the B atom;
		end do inner
	end do outer
	!	the following codes: sum the configuration Energy, would be the NOTE:
	!	Note: The following floating-point exceptions are signalling: IEEE_DENORMAL
	!	because the float precision of config_e & possibility is real(kind=8),not .eq. 0.0 and ranf()
	if ( deta_config_e .le. 0.0)then
		g_num_accepted_steps = g_num_accepted_steps + 1
		config_e = config_e + deta_config_e
		call chaRev_pro(moveChains_chosen,len_ASegment,len_BSegment) ! call the process of the exchange of each pair of atoms	
	else if ( deta_config_e .gt. 0.0)then
		possibility = exp(-deta_config_e *annealing_factor )
		if ( possibility .ge. ranf() ) then
			g_num_accepted_steps = g_num_accepted_steps + 1
			config_e = config_e + deta_config_e
			call chaRev_pro(moveChains_chosen,len_ASegment,len_BSegment)
		else
			do i = 1, len_BSegment, 1
				atom_X_position = i
				atom_Y_position = i + len_ASegment
				ichatmp(nchain(moveChains_chosen,atom_X_position)) = g_typeA
				ichatmp(nchain(moveChains_chosen,atom_Y_position)) = g_typeB
			end do
		endif	
	else 	!the else belong the outer if statement
	endif	!the endif belong the outer if statement
	g_sum_confe_PerTemper = g_sum_confe_PerTemper + config_e
	g_config_e = config_e 		! return current configuration Energy;
	!	return tge g_sum_var_e as the deta_config_e;
	g_sum_var_e = g_sum_var_e + deta_config_e 
end subroutine chains_reversing

subroutine chaRev_pro(moveChains_chosen,len_ASegment,len_BSegment)
	! to achieve the process of chains reversing;
	use globalparameters
	implicit none
	integer(kind=4),intent(inout)::moveChains_chosen,len_ASegment,len_BSegment
	integer(kind=4)::atom_X_position,atom_Y_position,id_of_positionX,id_of_positionY
	integer(kind=4)::len_copolymer,half_copolymer
	len_copolymer = len_ASegment + len_BSegment
	half_copolymer = int(len_copolymer/2)

	do i = 1, len_BSegment, 1
		atom_X_position = i 
		atom_Y_position = i + len_ASegment
		icha(nchain(moveChains_chosen,atom_X_position)) = g_typeB
		icha(nchain(moveChains_chosen,atom_Y_position)) = g_typeA

	end do
	!	exchange the atoms for length of the B segments
	do i = 1, len_BSegment, 1
		atom_X_position = i
		atom_Y_position = len_copolymer + 1-i
		id_of_positionX = nchain(moveChains_chosen,atom_X_position)
		id_of_positionY = nchain(moveChains_chosen,atom_Y_position)
		nchain(moveChains_chosen,atom_X_position) = id_of_positionY
		nchain(moveChains_chosen,atom_Y_position) = id_of_positionX
		postn2MrOrder(id_of_positionY) = atom_X_position
		postn2MrOrder(id_of_positionX) = atom_Y_position
	end do
	!	exchange the atoms remain
	do i = 1+len_BSegment, half_copolymer, 1
		id_of_positionX = nchain(moveChains_chosen,i)
		id_of_positionY = nchain(moveChains_chosen,len_copolymer+1-i)
		nchain(moveChains_chosen,i) = id_of_positionY
		nchain(moveChains_chosen,len_copolymer+1-i) = id_of_positionX
		postn2MrOrder(id_of_positionX) = len_copolymer+1-i
		postn2MrOrder(id_of_positionY) = i
	end do
end subroutine chaRev_pro



subroutine chains_snake(annealing_factor,init_config_e,concentration)
	!	chains_snake moving;
	use globalparameters
	implicit none
	real(kind=8),intent(in)::annealing_factor,init_config_e,concentration
	! 	call movingAtom_chosen(concentration)
	integer(kind=4)::len_copolymer
	integer(kind=4)::moveChains_chosen,Monomer_order,neigh_order
	integer(kind=4)::selectedSolvents_chosen
	integer(kind=4)::id_moveAtom,id_neigh_of_moveAtom,id_of_solvent
	real(kind=8)::ranf
	len_copolymer = g_len_ASegment + g_len_BSegment
	neigh_order = int(1.0+g_num_neighbours*ranf())
	! unit the variables of id_moveAtom and id_neigh_of_moveAtom;
	if ( concentration .ge.0.500 ) then
		moveChains_chosen = int(1.0+g_num_chains*ranf())
		Monomer_order = int(1.0+len_copolymer*ranf())
		id_moveAtom = nchain(moveChains_chosen,Monomer_order)
		id_neigh_of_moveAtom = id_of_neigh(id_moveAtom,neigh_order)
	else
		selectedSolvents_chosen = int(1.0+g_num_solvents*ranf())
		id_of_solvent = solvents2Position(selectedSolvents_chosen)
		id_neigh_of_moveAtom = selectedSolvents_chosen
		id_moveAtom = id_of_neigh(id_of_solvent,neigh_order)
	endif
end subroutine chains_snake

!_____________________________________________________________________________________________
!_____________________________________________________________________________________________

subroutine tempType_init()
	!	the initialization of temporary array of particles type 
	use globalparameters
	implicit none
	integer(kind=4)::err
	allocate(ichatmp(g_total_parti), stat=err)
	if (err /= 0) print *, "ichatmp: Allocation request denied"
	
	do i = 1, g_total_parti, 1
		ichatmp(i) = icha(i)
	end do
end subroutine tempType_init

subroutine count_num_solvents()
	!	its the counter of number of solvents;
	use globalparameters
	implicit none
	integer(kind=4)::counter,err
	counter = 0
	allocate(solvents2Position(g_num_atom), stat=err)
	if (err /= 0) print *, "solvents2Position: Allocation request denied"
	
	allocate(position2Solvents(g_num_atom), stat=err)
	if (err /= 0) print *, "position2Solvents: Allocation request denied"
	
	do i = 1, g_num_atom, 1
		if ( icha(i) .eq. g_typeS )then
			counter = counter + 1
			solvents2Position(counter) = i
			position2Solvents(i) = counter
		else
		endif
	end do
	g_num_solvents = counter
	print*,'number of solvents:',g_num_solvents
end subroutine count_num_solvents

subroutine chains_remove_pro(len_ASegment,len_BSegment,concentration)
	!	remove the chains after initial configuration;
	use globalparameters
	implicit none
	integer(kind=4),intent(in)::len_ASegment,len_BSegment
	real(kind=8),intent(in)::concentration

	integer(kind=4):: type_of_currAtom
	integer(kind=4)::len_copolymer
	integer:: num_remove_chains,num_chains,num_init_chains
	integer(kind=4)::current_chain,current_atom,Monomer_order,chain_new_order
	integer(kind=4)::err,counter,firstMonomer_of_chain,type_of_firstMonomer
	real(kind=8)::ranf

	len_copolymer = len_ASegment + len_BSegment
	num_init_chains = g_num_chains
	num_chains = int(g_num_atom*concentration/len_copolymer)
	num_remove_chains =num_init_chains - num_chains

	!print*,'num_init_chains:',num_init_chains
	!print*,'num_chains:',num_chains
	!print*,'num_remove_chains:',num_remove_chains
	
	remove_chains:do i = 1, num_remove_chains, 1

		type_of_currAtom = g_typeS
		Monomer_order = 1	!---2
		do while ( type_of_currAtom .eq.g_typeS )
			current_chain = int(1+num_init_chains*ranf())
			current_atom = nchain(current_chain,Monomer_order)
			type_of_currAtom = icha(current_atom)
			if ( icha(current_atom).ne.g_typeS ) then
				icha(current_atom) = g_typeS
			else
			endif
		end do ! dowhile_loop type_of_currAtom .eq. g_typeS
		do while ( Monomer_order .ne. len_copolymer )
			if ( type_of_currAtom .ne. g_typeS) then 
				Monomer_order = Monomer_order + 1
				current_atom = nchain(current_chain,Monomer_order)
				icha(current_atom) = g_typeS
			else
			endif
		end do ! dowhile_loop Monomer_order .lt. len_copolymer		
	end do remove_chains
	allocate(postn2chnOrder(g_num_atom), stat=err)
		if (err /= 0) print *, "postn2chnOrder: Allocation request denied"
	allocate(postn2MrOrder(g_num_atom), stat=err)
		if (err /= 0) print *, "postn2MrOrder: Allocation request denied"
	
	!	re_identify the order of chains and the monomer of a chain after remove the chains;
	counter = 0
	outer:do i = 1, g_num_chains, 1
		firstMonomer_of_chain = nchain(i,1)
		type_of_firstMonomer = icha(firstMonomer_of_chain)
		if ( type_of_firstMonomer .ne. g_typeS ) then
			counter = counter + 1
		else
		endif
		inner:do j = 1, len_copolymer, 1
			if ( type_of_firstMonomer .ne. g_typeS) then
				nchain(counter,j) = nchain(i,j)
				postn2chnOrder(nchain(i,j)) = counter
				postn2MrOrder(nchain(i,j)) = j
			else
			endif
		end do inner
	end do outer
	print*,'end subroutine chains_remove_pro(), the g_num_chains:',g_num_chains
end subroutine chains_remove_pro                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        

subroutine AB_diblock_copolymer_init(len_ASegment,len_BSegment)
	!	subroutine to initialise the AB BCP chains;
	use globalparameters
	implicit none
	integer,intent(in)::len_ASegment,len_BSegment
	integer(kind=4)::num_init_chains,height_of_cylinder,perChain_one_line
	integer(kind=4)::len_copolymer,id_of_chains,id_of_atom
	integer(kind=4)::current_atom,next_atom,current_chain,next_chain

	integer::err

	len_copolymer = len_ASegment + len_BSegment
	!print*,'len_ASegment:',len_ASegment,'len_BSegment',len_BSegment
	!print*,'len_copolymer:',len_copolymer

	height_of_cylinder = g_lz
	
	num_init_chains = int(height_of_cylinder/len_copolymer)*g_num_atom/height_of_cylinder
	perChain_one_line = int(height_of_cylinder/len_copolymer)

	!print*,'num_init_chains:',num_init_chains
	!print*,'perChain_one_line:',perChain_one_line

	allocate(nchain(num_init_chains,len_copolymer), stat=err)
	if (err /= 0) print *, "nchain: Allocation request denied"
	
	i = 1
	current_chain = 0
	do while ( i .lt. g_num_atom )
		j = i 
		do k = 1, perChain_one_line, 1
			if ( icha(j) .eq. g_typeS ) then
				current_chain = current_chain + 1
				current_atom = 1
				nchain(current_chain,current_atom) = j
				icha(j) = g_typeA
				do while ( current_atom .ne. len_ASegment )
					next_atom = j + 1
					if ( icha(next_atom).eq.g_typeS ) then
						icha(next_atom) = g_typeA
						current_atom = current_atom + 1
						nchain(current_chain,current_atom) = next_atom
						j = next_atom
					else
					endif
				end do !	len_ASegment_loop
				do while ( current_atom .ne. len_copolymer )
					next_atom = j + 1
					if ( icha(next_atom).eq.g_typeS )then
						icha(next_atom) = g_typeB
						current_atom = current_atom + 1
						nchain(current_chain,current_atom) = next_atom
						j = next_atom
					else
						print*,'error in building chains array'
					endif 
				end do !	len_ASegment_loop
			else
			endif
			j = j + 1
		end do
		i = i + height_of_cylinder
	end do  !	dowhile_i_loop
	if ( current_chain .ne. num_init_chains ) print*,'error in building chains array and order'
	g_num_chains = num_init_chains
	print*,'end subroutine AB_diblock_copolymer_init()'
end subroutine AB_diblock_copolymer_init

subroutine interaction_init(eas,eaw,ebs,ebw)
	!	the subroutine is used to build the interaction array between each particles
	use globalparameters
	implicit none
	real(kind=4), intent(in) :: eas,eaw,ebs,ebw
	do i = 1, 4, 1
		do j = 1, 4, 1
			eab(i,j) = 0.0
		end do
	end do
	eab(g_typeA,g_typeB) = 1.0
	eab(g_typeA,g_typeS) = eas
	eab(g_typeA,g_typeW) = eaw
	eab(g_typeB,g_typeS) = ebs
	eab(g_typeB,g_typeW) = ebw
	do i = 1, 4, 1
		do j = 1, 4, 1
			eab(j,i) = eab(i,j)
		end do
	end do
	print*,'end subroutine interaction_init()'
end subroutine interaction_init

subroutine latticePoints_init(lx,ly,lz)
	!	the subroutine of latticePoints_init used to build a array of box
	use globalparameters
	implicit none
	integer, intent(in) :: lx,ly,lz
	integer::height_of_cylinder,counter
	integer::err
	height_of_cylinder = lz
	! 	print*,'g_lxyz=',g_lxyz
	! 	print*,'lx,ly,lz:',lx,ly,lz
	allocate(box(g_lxyz,3), stat=err)
	if (err /= 0) print *, "box: Allocation request denied"
	counter = 0
	do i = 1, lx, 1
		do j = 1, ly, 1
			do k = 0, height_of_cylinder-1, 1
				counter = counter + 1
				box(counter,1) = i 
				box(counter,2) = j
				box(counter,3) = k
			end do
		end do
	end do
	g_num_box = counter
	!print*,'g_num_box',g_num_box,counter
	print*,'end subroutine latticePoints_init()'
end subroutine latticePoints_init

subroutine Parti_coor_init()
	! 	build the particles coordinates array and the order to every atom
	use globalparameters
	implicit none
	real(kind=4)::r0,dwall,len_near_points
	integer(kind=4)::height_of_cylinder
	integer(kind=4)::err,num_neighbours,num_atom,num_cylinder

	allocate(atom(g_lxyz,3), stat=err)
	if (err /= 0) print *, "atom: Allocation request denied"

	r0 = g_cylinder_r0
	dwall = g_wall_width
	height_of_cylinder = g_lz

	call confined_system_init(r0,dwall,height_of_cylinder)
	num_atom = g_num_atom
	num_cylinder = g_total_parti

	num_neighbours = g_num_neighbours
	len_near_points = g_len_near_points
	!print*,num_atom,num_cylinder,num_neighbours,len_near_points
	!call neighbourPoints_init(num_neighbours,len_near_points)
	call neighPts_init_distCalcu(num_atom,num_cylinder,num_neighbours,len_near_points)
	allocate(icha(num_cylinder), stat=err)
	if (err /= 0) print *, "icha: Allocation request denied"

	do i = 1, num_atom, 1
		icha(i) = g_typeS
	end do
	do i = num_atom+1, num_cylinder, 1
		icha(i) = g_typeW
	end do
	print*,'end subroutine Parti_coor_init()'
end subroutine Parti_coor_init

subroutine confined_system_init(r0,dwall,height_of_cylinder)
	!	build the confined system enviroment 
	!	the center of the cylinder is (half_lx,half_ly)
	use globalparameters
	implicit none
	real(kind=4),intent(in)::r0,dwall ! r0 is the confined radius
	integer(kind=4),intent(in)::height_of_cylinder
	integer(kind = 4)::counter,cylinder_radius
	real(kind=4)::x,y,z,rr
	integer(kind=4)::centerX,centerY
	centerX = g_lx /2
	centerY = g_ly /2

	counter = 0
	cylinder_radius = int(r0+dwall)

	!print*,'r0=',r0
	do i = 1, g_lxyz, 1
		x = box(i,1) - centerX
		y = box(i,2) - centerY
		rr = sqrt(dble(x*x + y*y))
		if ( rr .le. r0 ) then
			counter = counter + 1
			atom(counter,1) = box(i,1)
			atom(counter,2) = box(i,2)
			atom(counter,3) = box(i,3)
		else
		end if	
	end do
	g_num_atom = counter

	do i = 1, g_lxyz, 1
		x = box(i,1) -centerX
		y = box(i,2) -centerY
		rr = sqrt(dble(x*x + y*y))
		if ( rr .gt. r0 .and. rr .le. (r0+dwall))then
			counter =counter +1
			atom(counter,1) = box(i,1)
			atom(counter,2) = box(i,2)
			atom(counter,3) = box(i,3)
		else
		end if 
	end do
	print*,'subroutine-confined_system_init:'
	print*,'r0,dwall=',r0,dwall
	g_total_parti = counter
	print*,'g_num_atom=',g_num_atom
	print*,'g_total_parti=',g_total_parti
	print*,'end subroutine confined_system_init()'
end subroutine confined_system_init

subroutine neighbourPoints_init(num_neighbours,len_near_points)
	!	the subroutine make to build the neighbourPoints array
	use globalparameters
	implicit none
	real(kind=4), intent(in) :: len_near_points
	integer(kind=4),intent(in):: num_neighbours
	integer(kind=4)::x,y,z,counter
	real(kind=4)::rr
	do i = -1, 1, 1
		do j = -1, 1, 1
			do k = -1, 1, 1
				x = i
				y = j
				z = k
				rr = (x*x+ y*y+ z*z)
				if ( rr .le. len_near_points .and. rr .gt.0) then
					counter = counter  + 1
					neighbours(counter,1) = x
					neighbours(counter,2) = y
					neighbours(counter,3) = z

					! 	Program received signal SIGSEGV: Segmentation fault - invalid memory reference.
					!	Backtrace for this error
				else
				end if		
			end do
		end do
	end do
	!print*,'len_near_points:',len_near_points
	!print*,num_neighbours,counter
	if ( counter .ne. num_neighbours ) then
		print*,'error in build the array of neighbours(num_neighbours,3)'
	else
	end if
	do i = 1, g_num_atom, 1
		
	end do
	print*,'end subroutine neighbourPoints_init()'
end subroutine neighbourPoints_init

subroutine neighPts_init_distCalcu(num_atom,num_cylinder,num_neighbours,len_near_points)
	!	this coordinates of neighbours is calculated by the divided of the distance of double atoms;
	use globalparameters
	implicit none
	integer(kind=4),intent(in)::num_neighbours,num_atom,num_cylinder
	real(kind=4),intent(in):: len_near_points
	integer(kind=4)::x,y,z,counter,err,height_of_cylinder
	real(kind =4)::rr2,time_tic,time_toc

	allocate(num_neigh(num_atom), stat=err)
	if (err /= 0) print *, "num_neigh: Allocation request denied"

	do i = 1, num_atom, 1
			num_neigh(i) = 0
	end do

	allocate(id_of_neigh(num_atom,num_neighbours), stat=err)
	if (err /= 0) print *, "id_of_neigh: Allocation request denied"

	call cpu_time(time_tic)
	height_of_cylinder = g_lz
  	call find_id_of_neigh_met2(num_atom,num_cylinder)
	call cpu_time(time_toc)

	print*,'cpu_time spand on neighbourPoints:',time_toc - time_tic
end subroutine neighPts_init_distCalcu

subroutine find_id_of_neigh_met1(num_atom,num_cylinder)
	!	the first method to find the id_of_neigh(:,:) array;
	use globalparameters
	implicit none
	integer(kind=4)::x,y,z,height_of_cylinder,num_neighbours
	integer(kind=4),intent(in)::num_atom,num_cylinder
	integer(kind=4)::rr2
	integer(kind=4)::counter
	height_of_cylinder = g_lz
	num_neighbours = g_num_neighbours
	!print*,'num_atom:',num_atom,'num_cylinder:',num_cylinder
	!print*,'num_neighbours',num_neighbours
	counter = 0
	do i = 1, num_atom, 1
		num_neigh(i) = 0
		counter = counter + 1
 		do j = 1, num_cylinder, 1
			x = atom(counter,1) - atom(j,1)
			y = atom(counter,2) - atom(j,2)	
			z = atom(counter,3) - atom(j,3)	
			call PBC_criterion(z,height_of_cylinder)
			rr2 = x*x + y*y+ z*z
			if ( rr2 .lt. g_len_near_points .and. rr2 .gt.0.0) then
				num_neigh(counter) = num_neigh(counter) + 1
				id_of_neigh(counter,num_neigh(counter)) = j
			else
			endif
 		end do
  		if ( num_neigh(counter) .ne. num_neighbours) print*,'error in building neigh array'
 	end do
 	!	testing
 	do i = 1, num_atom, 1
 		if ( num_neigh(i) .ne. g_num_neighbours) then
 			!print*,'hello world'
 			goto 1001
 		else
 		endif
 	end do
 1001	print*,'error in find_id_of_neigh_met1'
end subroutine find_id_of_neigh_met1

subroutine find_id_of_neigh_met2(num_atom,num_cylinder)
	!	the Optimized algorithm to find the id_of_neigh();
	use globalparameters
	implicit none
	integer(kind=4),intent(in)::num_atom,num_cylinder
	integer(kind=4)::height_of_cylinder,num_neighbours,x,y,z
	real(kind=4)::rr2
	height_of_cylinder = g_lz
	num_neighbours = g_num_neighbours
	do i = 1, num_atom, 1
		j = i
		do while ( num_neigh(i) .lt.num_neighbours )
			j = j + 1
			x = atom(i,1) -atom(j,1)
			y = atom(i,2) -atom(j,2)
			z = atom(i,3) -atom(j,3)
			call PBC_criterion(z,height_of_cylinder)
			rr2 = x*x + y*y + z*z
			if (rr2 .lt. g_len_near_points .and. rr2 .gt.0.0) then 
				if ( j .le. num_atom ) then
					num_neigh(j) = num_neigh(j) + 1
					id_of_neigh(j,num_neigh(j)) = i
				else
				endif
				num_neigh(i) = num_neigh(i) + 1
				id_of_neigh(i,num_neigh(i)) = j
			else
			endif
			
		end do
	end do
	!	testing
	do i = 1, num_atom, 1
 		if ( num_neigh(i) .ne. g_num_neighbours) then
 			!print*,'hello world'
 			goto 1002
 		else
 		endif
 	end do
 1002 print*,'error in find_id_of_neigh_met2'	
end subroutine find_id_of_neigh_met2

subroutine figure(step)
	use globalparameters
	implicit none
	integer(kind=4),intent(in)::step
end subroutine figure

function configE_calculating()
	! the subroutine is used to calculate the configuration Energy
	use globalparameters
	implicit none
	real(kind=8) :: configE_calculating
	integer(kind=4)::id
	real(kind=8)::config_e
	config_e = 0
	do i = 1, g_num_atom, 1
		do j = 1, g_num_neighbours, 1
			id = id_of_neigh(i,j)
			if ( icha(i) .lt. icha(id) ) then
				config_e = config_e + eab(icha(i),icha(id))
			else
			endif
		end do
	end do
	configE_calculating = config_e
end function configE_calculating

subroutine PBC_criterion(var,period)
	! subroutine of Period Boundary Condition
	implicit none
	integer(kind=4),intent(inout)::var,period
	integer(kind=4)::half_period,var_init
	var_init = var
	half_period = period/2
	if ( var .gt.half_period ) then
		var = var - period
	else if ( var .le. -half_period ) then
		var = var + period
	else
	endif
end subroutine PBC_criterion

function ranf()
	use globalparameters
	implicit none
	real(kind=8) :: ranf
	integer,parameter:: K4B = selected_int_kind(9)
	integer(K4B),parameter:: IA=16807, IM=2147483647, IQ=127773, IR=2836
	real(kind=4),save::am
	integer(K4B),save::ix = -1,iy = -1,kx
	if (iseed <= 0 .or. iy < 0) then                                        ! Initialize
           am=nearest(1.0,-1.0)/IM
           iy=ior(ieor(888889999, abs(iseed)),1)
           ix=ieor(777755555, abs(iseed))
           iseed=abs(iseed)+1                                                   ! Set idum positive
	    end if
	    ix=ieor(ix, ishft(ix, 13))                                                ! Marsaglia shift sequence with period 2^^32-1
	    ix=ieor(ix, ishft(ix, -17))
	    ix=ieor(ix, ishft(ix, 5))
	    kx=iy/IQ                            ! Park-Miller sequence by Schrage's method, period 2^^31-1
	    iy=IA*(iy-kx*IQ)-IR*kx

	    if (iy < 0 ) iy=iy+IM
	    ! Combine the two generators with masking to ensure nonzero value
	    ranf=am*ior(iand(IM, ieor(ix,iy)),1)     
end function ranf

subroutine export2cc1_Figure(step)
	use globalparameters
	implicit none
	integer, intent(in) :: step
	integer(kind=4)::num_atom
	character(len=256) :: esti_morph,temp_Char
	if ( step .lt.10 ) then
		write(temp_Char,'(I1)')step
	else if ( step .ge. 10 .and. step .lt. 100 ) then
		write(temp_Char, '(I2)')step
	else
	endif
		
	write(esti_morph,*)'morph_',trim(temp_Char),'.cc1'

	open(unit = 2,file = esti_morph)
	num_atom = (g_len_ASegment + g_len_BSegment)*g_num_chains
	j = 0
	write(2,*)num_atom
	do i = 1, g_num_atom, 1
		if ( icha(i).eq.g_typeA ) then
			j = j + 1
			write(2,*)'Na',j,atom(i,1),atom(i,2),atom(i,3),190
		else if ( icha(i).eq.g_typeB ) then
			j = j + 1
			write(2,*)'Br',j,atom(i,1),atom(i,2),atom(i,3),110
		else
		endif
	end do
	close(3)
end subroutine export2cc1_Figure


!_____________________________________________________________________________________________

!	ReadMe:

! 	this procedure is used to solve the problem of self-assemly of copolymer confined_system;
!	To simplify the procedure, Making all of the array as dynamic and global array;
!	To make the procedure is easy to transplant, identify the global and local variables;
!	Try to change it ! 
!	And you can contact with me;

!	Attention:

!		It would be error of the numerical value in testing config_e after remove the chains 
!	when using different seed of the generation function of random number ;
!		the Energy of the different initial chains_configuration with the confined pore wall are different;		
!		So you should keep the seed, the concentration, the confined_system radius, the interaction array,
!	and even the length of each segment of chains should all be the same when testing;

!	the following situation would be the NOTE:
	!	Note: The following floating-point exceptions are signalling: IEEE_DENORMAL
	!	because the dummy variables of function is not initialization for a value;