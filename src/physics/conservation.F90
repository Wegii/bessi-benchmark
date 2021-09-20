module conservation
	! Check conservation of different values

	use bessi_defs

	implicit none

	! Mass conservation checks
    real(kind=8), dimension(ndays) :: check_init_mass_snow  !kg/m2
    real(kind=8), dimension(ndays) :: check_init_mass_water !kg/m2
    real(kind=8), dimension(ndays) :: check_accum_snow      !kg/m2
    real(kind=8), dimension(ndays) :: check_accum_water     !kg/m2
    real(kind=8), dimension(ndays) :: check_melted_snow     !kg/m2
    real(kind=8), dimension(ndays) :: check_runoff_water    !kg/m2
    real(kind=8), dimension(ndays) :: check_refreeze_water  !kg/m2
    real(kind=8), dimension(ndays) :: check_end_mass_snow   !kg/m2
    real(kind=8), dimension(ndays) :: check_end_mass_water  !kg/m2
    real(kind=8), dimension(ndays) :: check_ice2ice         !kg/m2

    ! Energy conservation checks
    real(kind=8), dimension(ndays) :: check_init_energy_snow    !J/m2
    real(kind=8), dimension(ndays) :: check_init_energy_water   !J/m2
    real(kind=8), dimension(ndays) :: check_end_energy_snow     !J/m2
    real(kind=8), dimension(ndays) :: check_end_energy_water    !J/m2

    real(kind=8), dimension(ndays) :: check_surface_e_flux_obs  !J/m2
    real(kind=8), dimension(ndays) :: check_surface_e_flux_diag !J/m2
    real(kind=8), dimension(ndays) :: check_e_accum_tot         !J/m2

    real(kind=8), dimension(ndays) :: check_e_melt_tot      !J/m2
    real(kind=8), dimension(ndays) :: check_e_melt_heat     !J/m2
    real(kind=8), dimension(ndays) :: check_e_melt_qq       !J/m2
    real(kind=8), dimension(ndays) :: check_e_melt_runoff   !J/m2

    real(kind=8), dimension(ndays) :: check_e_freeze_tot    !J/m2
    real(kind=8), dimension(ndays) :: check_e_freeze_heat   !J/m2

    real(kind=8), dimension(ndays) :: check_e_perc_runoff   !J/m2
    real(kind=8), dimension(ndays) :: check_e_ice2ice       !J/m2

    character(1) :: tab = char(09)

	public :: check_init
	public :: check_start_loop
	public :: check_accumulation
	public :: get_check_energy
	public :: check_flux
	public :: check_melt
	public :: check_runoff
	public :: check_refreeze
	public :: check_end_loop
	public :: check_smb
	public :: check_end
contains
	subroutine check_init()
		! Initialize variables with zero

		implicit none

        check_init_mass_snow = 0.       
        check_init_mass_water = 0.      
        check_accum_snow  = 0.          
        check_accum_water = 0.          
        check_melted_snow = 0.          
        check_runoff_water = 0.         
        check_end_mass_snow = 0.        
        check_end_mass_water = 0.       
        check_refreeze_water = 0.       
        check_ice2ice = 0.


        check_init_energy_snow  = 0.    
        check_end_energy_snow  = 0.     
        check_init_energy_water = 0.    
        check_end_energy_water = 0.     

        check_surface_e_flux_obs = 0.   
        check_surface_e_flux_diag = 0.  
        check_e_accum_tot = 0.          

        check_e_freeze_tot = 0.         
        check_e_freeze_heat = 0.        

        check_e_melt_tot = 0.           
        check_e_melt_heat = 0.          
        check_e_melt_qq = 0.            
        check_e_melt_runoff = 0.        

        check_e_perc_runoff = 0.        
        check_e_ice2ice = 0.    
	end subroutine check_init

	subroutine check_start_loop(ix, iy, time)
		! Accumulate variables

		implicit none
        
		integer, intent(in) :: ix
		integer, intent(in) :: iy
		integer, intent(in) :: time

        check_init_mass_snow(time) = check_init_mass_snow(time) + sum(snowman(ix, iy, :))     
        check_init_mass_water(time) = check_init_mass_water(time) + sum(lwmass(ix, iy, :))
        check_init_energy_snow(time) = check_init_energy_snow(time) + sum(snowman(ix, iy, :)*snow_temp(ix, iy, :))*c_i
        check_init_energy_water(time) = check_init_energy_water(time) + sum(lwmass(ix, iy, :))*(kelvin*c_i + L_lh)
	end subroutine check_start_loop


	subroutine check_accumulation(ix, iy, time, accum, rainman)
		! Accumulate accumulation variables

		implicit none

		integer, intent(in) :: ix
		integer, intent(in) :: iy
		integer, intent(in) :: time
		real(kind=8), intent(in) :: rainman
		real(kind=8), intent(in) :: accum

		check_accum_snow(time) = check_accum_snow(time) + accum*dt_firn
        check_accum_water(time) = check_accum_water(time) + rainman*dt_firn
        check_e_accum_tot(time) = check_e_accum_tot(time) + dt_firn*accum*c_i*snow_temp(ix, iy, 1) + &
        							dt_firn*rainman*(c_i*kelvin + L_lh)

	end subroutine check_accumulation

	pure function get_check_energy(ix, iy) result(rv)
		! Calculate energy

		implicit none

		integer, intent(in) :: ix
		integer, intent(in) :: iy
		real(8) :: rv

		rv = sum(snowman(ix, iy, :)*snow_temp(ix, iy, :))*c_i + sum(lwmass(ix, iy, :))*(kelvin*c_i + L_lh)
	end function get_check_energy

	subroutine check_flux(ix, iy, time, dummy_energy, dummy_heat)
		! Accumulate flux variables

		implicit none

		integer, intent(in) :: ix
		integer, intent(in) :: iy
		integer, intent(in) :: time
	    real(kind=8), intent(in) :: dummy_energy
		real(kind=8), intent(in) :: dummy_heat
		
		check_surface_e_flux_obs(time) = check_surface_e_flux_obs(time) - dummy_energy + get_check_energy(ix, iy)
		check_surface_e_flux_diag(time) = check_surface_e_flux_diag(time) + dummy_heat
	end subroutine check_flux


	subroutine check_melt(ix, iy, time, dummy_melt, dummy_energy, dummy_heat, dummy_runoff, dummy_e_qq)
		! Accumulate melt variables

		implicit none
		
		integer, intent(in) :: ix
		integer, intent(in) :: iy
		integer, intent(in) :: time
		real(kind=8), intent(in) :: dummy_melt
	    real(kind=8), intent(in) :: dummy_energy
		real(kind=8), intent(in) :: dummy_heat
		real(kind=8), intent(in) :: dummy_runoff
		real(kind=8), intent(in) :: dummy_e_qq

        check_melted_snow(time) = check_melted_snow(time) + dummy_melt
        check_runoff_water(time) = check_runoff_water(time) + dummy_runoff
        check_e_melt_tot(time) = check_e_melt_tot(time) - dummy_energy + get_check_energy(ix, iy)                       
        check_e_melt_heat(time) = check_e_melt_heat(time) + dummy_heat
        check_e_melt_qq(time) = check_e_melt_qq(time) + dummy_e_qq
        check_e_melt_runoff(time) = check_e_melt_runoff(time) + dummy_runoff*(c_i*kelvin + L_lh)
	end subroutine check_melt

	subroutine check_runoff(time, dummy_runoff)
		! Accumulate runoff variables

		implicit none

		integer, intent(in) :: time
		real(kind=8), intent(in) :: dummy_runoff

		check_runoff_water(time) = check_runoff_water(time) + dummy_runoff
		check_e_perc_runoff(time) = check_e_perc_runoff(time) + dummy_runoff*(c_i*kelvin + L_lh)		
	end subroutine check_runoff

	subroutine check_refreeze(ix, iy, time, dummy_heat, dummy_refreeze, dummy_energy)
		! Accumulate refreeze variables

		implicit none

		integer, intent(in) :: ix
		integer, intent(in) :: iy
		integer, intent(in) :: time
		real(kind=8), intent(in) :: dummy_heat
		real(kind=8), intent(in) :: dummy_refreeze
		real(kind=8), intent(in) :: dummy_energy

		check_e_freeze_heat(time) = check_e_freeze_heat(time) + dummy_heat
		check_refreeze_water(time) = check_refreeze_water(time) + dummy_refreeze
		check_e_freeze_tot(time) = check_e_freeze_tot(time) - dummy_energy + get_check_energy(ix, iy)
	end subroutine check_refreeze

	subroutine check_end_loop(ix, iy, time)
		! Accumulate variables

		implicit none

		integer, intent(in) :: ix
		integer, intent(in) :: iy
		integer, intent(in) :: time

		check_end_mass_snow(time) = check_end_mass_snow(time) + sum(snowman(ix, iy, :))
		check_end_mass_water(time) = check_end_mass_water(time) + sum(lwmass(ix, iy, :))

		check_end_energy_snow(time) = check_end_energy_snow(time) + sum(snowman(ix, iy, :)*snow_temp(ix, iy, :))*c_i 
		check_end_energy_water(time) = check_end_energy_water(time) + sum(lwmass(ix, iy, :))*(kelvin*c_i + L_lh)
	end subroutine check_end_loop

	subroutine check_smb(dummy_ice2ice, dummy_energy, dummy_runoff)
		! Accumulate smb variables

		implicit none

		real(kind=8), intent(in) :: dummy_ice2ice
		real(kind=8), intent(in) :: dummy_energy
		real(kind=8), intent(in) :: dummy_runoff

		check_ice2ice(ndays) = check_ice2ice(ndays) + dummy_ice2ice
		check_e_ice2ice(ndays) = check_e_ice2ice(ndays) + dummy_energy
		check_runoff_water(ndays) = check_runoff_water(ndays) + dummy_runoff
		check_e_perc_runoff(ndays) = check_e_perc_runoff(ndays) + dummy_runoff*(c_i*kelvin + L_lh)
	end subroutine check_smb

	subroutine check_end(ix, iy)
		! Accumulate variables

		implicit none

		integer, intent(in) :: ix
		integer, intent(in) :: iy

        check_end_mass_snow(ndays) = check_end_mass_snow(ndays) + sum(snowman(ix, iy, :))
        check_end_mass_water(ndays) = check_end_mass_water(ndays) + sum(lwmass(ix, iy, :))

        check_end_energy_snow(ndays) = check_end_energy_snow(ndays) + sum(snowman(ix, iy, :)*snow_temp(ix, iy, :))*c_i
        check_end_energy_water(ndays) = check_end_energy_water(ndays) + sum(lwmass(ix, iy, :))*(kelvin*c_i + L_lh)
	end subroutine check_end
end
