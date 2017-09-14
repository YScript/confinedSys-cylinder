module system_init
!	made by fengyuan 2017/6/8
!
!	set 5 kind of particles for interaction,  A-1 & B-2 & C-3 & Solvent->4  &  Wall->5 ;

	use global_parameter
	implicit none
contains

	subroutine interaction_init(eas,ebs,eaw,ebw)	
		!   the interaction array in confined system is different from others
		integer i,j,k
		double precision,intent(in)::eas,ebs,eaw,ebw
		do i = 1, 5, 1
			do j = 1, 5, 1
				eab(i,j) = 0.0
			end do		 	
		end do
		eab(1,2) = 1.0
		!		eab(1,3) = 1.0			
		eab(1,4) = eas
		eab(1,5) = eaw
		!		eab(2,3) = 1.0
		eab(2,4) = ebs
		eab(2,5) = ebw
				
		do i = 1, 5, 1
			do j = 1, 5, 1
				eab(j,i) = eab(i,j)
			end do
		end do
		print*,'subroutine produce interaction_init()'
	end subroutine interaction_init
! end	subroutine interaction_init-------------------------------------

	subroutine box_cylinder_confined_init(r0,hz)
		
		integer, intent(in) :: hz
		double precision,intent(in)::r0
		!	integer,intent(out)::g_num_box,g_num_atom,g_num_wall
		integer i,j,k,l,m,n,ncount,id,x,y,rx0,ry0
		double precision dwall , total_diameter , rr

		ncount = 0
		dwall = 3.0 	! confined cylinder wall width;
		!	total_diameter = 2*(r0+ dwall) + 1
		do i = 0, g_box_len, 1
			do j = 0, g_box_len, 1
				do k = 0, (hz-1), 1
					ncount = ncount +1
					box(ncount,1) = i
					box(ncount,2) = j
					box(ncount,3) = k
				end do					
			end do
		end do
		! creat a box ,its cordinates be from x(0:diameter);y(0:diameter);z(0:hz-1);
		!	>> (diameter *diameter * hz);
		g_num_box = ncount
		allocate(atom(g_num_box)) !!! array atom record the  identifier of cylinder
		ncount = 0
		rx0 = g_box_len/2
		ry0 = rx0
		!	>>>>> Marking the atoms in the cylinder;
		do i = 1, g_num_box, 1
			x = box(i,1) - rx0
			y = box(i,2) - ry0
			rr = sqrt((x*x) + (y*y) )  !!! is there need to convert type?
			if ( rr .le. r0 ) then
				ncount = ncount + 1
				atom(ncount) = i 	 !  the array record the id of atom,not cordinates;
			else
			endif
		end do
		g_num_atom = ncount
		!	>>>>> Marking the atoms in the cylinder wall;
		do i = 1, g_num_box, 1
			x = box(i,1) - rx0
			y = box(i,2) - ry0
			rr = sqrt((x*x) + (y*y) ) 
			if ( rr .le. (r0+dwall) ) then
				ncount = ncount + 1
				atom(ncount) = i 
			else
			endif
		end do
		g_num_wall = ncount - g_num_atom
		print*,'number_atom',g_num_atom,'number_box',g_num_box
		print*,'finish the box_cylinder_confined_init section'
	end subroutine box_cylinder_confined_init
! 	end	subroutine box_cylinder_confined_init-------------------------------------
	
	subroutine near_point_init(len_NP,num_NP,box_len,hz)
		!	subroutine of build nearest neighbouring lattice points coordinates;
		double precision, intent(in) :: len_NP
		integer,intent(in)::num_NP,box_len,hz
		integer i,j,k,l,m,n,x,y,z,ncount
		double precision rr
		ncount = 0
		allocate(nearP_coor(num_NP,3))
		do i = -1, 1, 1
			do j = -1, 1, 1
				do k = -1, 1, 1
					x = i
					y = j
					z = k
					rr = (x*x)+(y*y) + (z*z)
					if ( (rr .gt.0).and.(rr .le. len_NP) )then
						ncount = ncount +1
						nearP_coor(ncount,1) = x
						nearP_coor(ncount,2) = y 
						nearP_coor(ncount,3) = z
					else
					endif
				end do
			end do
		end do
		ncount_NP = ncount
		
		if ( num_NP .ne. ncount_NP ) then
			print*,'error in building nearP_coor'
			pause
		endif

		do i = 1, n_box, 1
			num_nearpoint(i) = 0
		end do

		allocate(id_NP(n_box,ncount_NP))
		do i = 1, n_box, 1
			do j = 1, ncount_NP, 1
				x = box(i,1) + nearP_coor(j,1)
				y = box(i,2) + nearP_coor(j,2)
				z = box(i,3) + nearP_coor(j,3)

				if ( x .gt.box_len )then 
					x = -box_len + x
				else if ( x .lt.0 ) then
					x = box_len + x
				else
				endif
				if ( y .gt. box_len ) then
					y = -box_len + y
				else if ( y .lt. 0 ) then
					y = box_len + x
				else
				endif
				if ( z .gt.hz ) then
					z = -hz + z
				else if ( z .lt. 0 ) then
					z = hz + z
				else
				endif

				id_NP(i,j) = x*box_len*hz + y*hz + z + 1		! record the id of j nearst points for i atom;
				num_nearpoint(i) = num_nearpoint(i) + 1
			end do
		end do
		do i = 1, g_num_atom, 1
			if ( num_nearpoint(i).ne.num_NP)then
				print*,'error in building nearest neighbouring points'
			endif
		end do
		print*,'finish subroutine -- near_point_init'
	end subroutine near_point_init
!	end subroutine near_point_init---------------------------------------------------------------

	subroutine polymer_chains_init(len_segmentA,len_segmentB,concentration)
		!	subroutine for building polymer_chains solution with concentration in cylinder;
		!	this subroutine depends on the arguments and parameters::g_num_atom,g_cylinder_height;
		integer, intent(in) ::len_segmentA,len_segmentB 
		double precision,intent(in)::concentration
		integer,intent(out)::num_BCP,len_BCP,
		integer num_build_polymer_chains,num_chain_per_line,num_real_polymer_chains,num_remove_polymer_chains
		integer i,j,k,l
		integer len_arr
		integer target_atom_id,next_target_atom_id
		integer target_chain_id
		integer type_of_monomer_A = 1,type_of_monomer_B = 2,type_of_monomer_C = 3
		integer type_of_monomer_S = 4,type_of_monomer_W = 5
		integer target_units_of_chain,type_of_monomer

		len_BCP = len_segmentA + len_segmentB
		allocate(icha(g_num_atom+g_num_wall))
		do i = 1, g_num_atom, 1
			icha(i) = type_of_monomer_S		! set all the atoms type is solvents; 
		end do		
		do i = g_num_atom, n_atom+n_wall, 1
			icha(i) = type_of_monomer_W		! set the wall type;
		end do

		num_chain_per_line = int( g_cylinder_height /len_BCP)
		num_real_polymer_chains = int(concentration *g_num_atom/len_BCP)
		num_build_polymer_chains = g_num_atom*num_chain_per_line/g_cylinder_height
		g_num_chains = num_real_polymer_chains

		num_remove_polymer_chains = num_build_polymer_chains - num_real_polymer_chains	
		print*,'num_real_polymer_chains',num_real_polymer_chains
		print*,'num_build_polymer_chains',num_build_polymer_chains
		print*,'num_remove_polymer_chains',num_remove_polymer_chains

		i = 0
		target_atom_id = 0
		target_chain_id = 0
		next_target_atom_id = 0
		do while(i .lt. g_num_atom)
			!	change one solvent atom to start building chains array;
			!	the number of atoms is larger than that of chains no doubt;
			target_atom_id = i
			!	record this atom
			do l = 1, num_chain_per_line, 1
				!	undergo all the polymer_chains;
				if ( icha(target_atom_id).eq.type_of_monomer_S ) then
					target_chain_id = target_chain_id + 1
					len_arr = 1
					!	len_arr is the length of chains arrangement;	
					nchain(target_chain_id,len_arr) = target_atom_id
					icha(target_atom_id) = type_of_monomer_A
					do while ( len_arr .ne. len_segmentA )
						next_target_atom_id = target_atom_id + 1
						
						if ( icha(next_target_atom_id) .eq. type_of_monomer_S ) then
							len_arr = len_arr + 1
							nchain(target_chain_id,len_arr) = next_target_atom_id
							icha(next_target_atom_id) = type_of_monomer_A
							target_atom_id = next_target_atom_id
						else
						endif
					end do
					!	end building A segments
					do while ( len_arr .ne. len_BCP )
						next_target_atom_id = target_atom_id + 1

						if ( icha(next_target_atom_id) .eq.type_of_monomer_S ) then
							len_arr = len_arr + 1
							nchain(target_chain_id,len_arr) = next_target_atom_id
							icha(next_target_atom_id) = type_of_monomer_B
							target_atom_id = next_target_atom_id
						else
						endif
					end do
					!	end building B segments
				else
				endif
				target_atom_id = target_atom_id + 1
			end do
			i = i + g_cylinder_height
			!building num_build_polymer_chains until enough;
		end
		print*,'end build chains'

		target_chain_id = 0
		do i = 1, num_remove_polymer_chains, 1
			!	all the chains that needs to be remove;
			target_units_of_chain = 1
			type_of_monomer = type_of_monomer_S 
			! ensure the produce would be excuted for each loop;
			do while ( type_of_monomer .eq. type_of_monomer_S )
				target_chain_id = int(1.0+g_num_chains*ranf())
				target_atom_id = nchain(target_chain_id,target_units_of_chain)
				type_of_monomer = icha(target_atom_id)
				if ( type_of_monomer .ne. type_of_monomer_S ) then
					icha(target_atom_id) = type_of_monomer_S
				else
				endif
			end do			
			do while ( target_units_of_chain .ne. len_BCP )
				if ( type_of_monomer .ne. type_of_monomer_S ) then
					target_units_of_chain = target_units_of_chain + 1
					target_atom_id = nchain(target_chain_id,target_units_of_chain)
					icha(target_atom_id) = type_of_monomer_S
				else
				endif
			end do
		end do
		print*,'end remove chains'

		target_chain_id = 0
		allocate(id_chains(g_num_atom))
		allocate(id_chain_units(g_num_atom))
		do i = 1, g_num_chains, 1
			target_atom_id = nchain(i,1)
			type_of_monomer = icha(target_atom_id)
			if ( type_of_monomer .ne. type_of_monomer_S ) then
				target_chain_id = target_chain_id + 1
			else
			endif
			do j = 1, len_BCP, 1
				if ( type_of_monomer .ne. type_of_monomer_S ) then
					nchain(target_chain_id,j) = nchain(i,j)
					id_chains(nchain(i,j)) = target_chain_id
					id_chain_units(nchain(i,j)) = j
				else
				endif
			end do
		end do
		print*,'concentration',concentration
		print*,"end subroutine polymer_chains_init"
	end subroutine polymer_chains_init
!	end subroutine polymer_chains_init---------------------------------------------
	
	subroutine Marking_empty_particle()
		!	the subroutine used to Mark the empty particle
		integer i,j,ntotempty
		integer type_of_monomer_S = 4
		i = 0
		allocate(id_empty(g_num_atom))
		do j = 1, g_num_atom, 1
			if ( icha(j).eq. type_of_monomer_S) then
				i = i + 1
				id_empty(i)=j
			!	mne(j)=i
			else
			endif
		end do
		ntotempty = i
		print*,'number of empty',ntotempty
		print*,'end subroutine Marking_empty_particle'
	end subroutine Marking_empty_particle
!	end subroutine Marking_empty_particle

	subroutine calculate_config_E(ntot,num_NP)
		!	calculate the configation Energy with ntot atoms in cylinder and with num_NP neighbour points;
		! 	return the double precision e
		integer, intent(in) :: ntot,num_NP
		integer i,j
		integer type_of_nearest_Point
		double precision, intent(out) :: e 
		e = 0.0
		do i = 1, ntot, 1
			do j = 1, num_NP, 1
				type_of_nearest_Point = icha(id_NP(i,j))
				if ( icha(i) .le. type_of_nearest_Point ) then
					e=e+eab(icha(i),type_of_nearest_Point)
				else
				endif
			end do
		end do
		config_energy =  e
		print*,'config-E=',e
		print*,"end subroutine calculate_config_E"
	end subroutine calculate_config_E
!	end subroutine calculate_config_E

	subroutine tmp_character(number_atom,num_wall)
		integer, intent(inout) :: number_atom,num_wall
		integer total_number = number_atom + num_wall
		allocate(ichap(total_number))
		do i = 1, total_number, 1
			ichap(i) = icha(i)
		end do
		print*,"end subroutine tmp_character"
	end subroutine tmp_character
!	end subroutine tmp_character

end module system_init
