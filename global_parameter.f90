module global_parameter

	implicit none
! the golbal parameter used universally !
	integer, parameter:: g_length_Ablock = 10,g_length_Bblock = 2
	integer,parameter:: g_length_BCP = g_length_Ablock + g_length_Bblock
	integer,parameter::g_box_len = 60, g_cylinder_height = 60

	double precision,parameter:: g_confined_r0 = 10.0
	integer,parameter:: iseed = -191

	double precision,parameter:: g_concentration = 0.45
	integer g_num_box,g_num_atom,g_num_wall
!	end global_parameter

!	system initialization variable !
	integer,allocatable,dimension(:,:)::box,nearP_coor
	integer,allocatable,dimension(:)::atom,num_nearpoint,icha,icha_tmp
	double precision,dimension(5,5)::eab
	integer,allocatable,dimension(:,:)::id_of_NP
	integer n_box,n_atom
	integer,allocatable,dimension(:)::id_empty
	integer,allocatable,dimension(:)::id_chains,id_chain_units
	double precision config_energy
! 	end system initialization		

!	random number generator !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	
	INTEGER, PARAMETER :: K4B=selected_int_kind(9)
    INTEGER(K4B) iseed
    INTEGER(K4B), PARAMETER :: IA=16807, IM=2147483647, IQ=127773, IR=2836
    REAL, SAVE :: am
    INTEGER(K4B), SAVE :: ix=-1, iy=-1, kx
!	end random number generator
end module global_parameter



! Readme for parameters :
! g_length_BCP,g_length_Ablock -the length of A segment and BCP;
! g_box_len,g_cylinder_height		-the wdith of box and the confined cylinder height;
! g_confined_r0		-the confined cylinder radius;
! g_concentration 	-the BCP concentration in the solution;
! iseed	-the random number generation function seed;


! Readme for arrays:
! box  -the coordinates records of box;
! nearP_coor 	-the near points coordinates;
! atom -recording the moving atom identifier;
! num_nearpoint - number of near points for a site;
! icha icha_tmp 	- the type and temporary type of an atom;
! eab 	-the interaction between two atoms;
! id_NP	- the identifier of an near point;



! Readme for variables:
! g_num_box	-the total number of box sites;
! g_num_atom	-the total number of the confined system,include the wall;
! g_num_wall 	-the total number of the atoms in the cylinder wall;