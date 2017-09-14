module motion
	use global_parameter
	use random_gen
	use system_init

	implicit none
contains
	subroutine annealing_plan(annealing_init_temp,annealing_steps)
		double precision, intent(inout) :: annealing_init_temp
		integer, intent(in):: annealing_steps
		double precision:: ann_coef1,ann_coef2,epp,ep0,ep,e2p
		integer::nmontes,num_esti,esti_num,MCS_per_T_counter,MCS_per_T
		integer:: chain_re_num.chain_re_acc,chain_mv_num,chain_mv_acc
		integer:: solvent1_num,solvent2_num,solvent1_acc,solvent2_acc

		ann_coef1 = 0.985
		ann_coef2 = 0.955
		
		do i = 1, annealing_steps, 1
			if (ep0-epp .gt. 4000. ) then
				annealing_init_temp = annealing_init_temp * ann_coef1
			else
				annealing_init_temp = annealing_init_temp * ann_coef2
			endif

			ep0 = epp  !! old status energy
			ep = 0.		!! new status energy
			e2p = 0.	!! energy^2
			nmontes =  0 !! monomer move counter
			MCS_per_T_counter = 0		!!	Monte Carlo step
			esti_num = 0
			num_esti = 0

			chain_re_num = 0		 !! chain reverse number
			chain_re_acc = 0		 !! accepted chain reverse counter
			chain_mv_num = 0		!! chain move number
			chain_mv_acc = 0		!! accepted chain move counter
			Probabi_OverTurn = 1 !! float(1)/float(jj) !! reverse possbility
			inv_T = 1./annealing_init_temp
			do while ( MCS_per_T_counter .lt. MCS_per_T ) !! 	loop
				chain_selected = int(1.0 + num_chain * ranf())
				
			end do

			print*,annealing_steps,config_e,tmp_e,de,annealing_init_temp
		end do

		print*,'end subroutine annealing_plan'
	end subroutine annealing_plan



	subroutine chain_snake()
	print*,'end subroutine chain_snake'		
	end subroutine chain_snake

	subroutine solvent_exchang()
	print*,'end subroutine solvent_exchang'
	end subroutine solvent_exchang




	subroutine chains_reserve()
		! 
		print*,'end subroutine chains_reserve'
	end subroutine chains_reserve

	subroutine chains_snake_move()
		! 
		print*,'end subroutine chains_move'
	end subroutine chains_snake_move

	subroutine solvent_chosen()
		!, intent() :: 
		print*,'end subroutine solvent_chosen'
	end subroutine solvent_chosen

	subroutine monomer_chosen()
		!, intent(inout) :: 
		print*,'end subroutine monomer_chosen'
	end subroutine monomer_chosen

	subroutine exchange_units()
		
		print*,'end subroutine exchange_units'
	end subroutine exchange_units

end module motion
