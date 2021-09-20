module precipitation_accumulate
	! Routine for calculating the amount of new snow

	use bessi_defs

	implicit none

	private
	public :: accumulate
contains
	subroutine accumulate(ix, iy, time, rainman, accum, nday_snowfall, H_lh, K_lh, air_temp_ice, precip_ice)
		! Determines the amount [kg/m2] of snow that falls on a gridcell
		! and adjusts the inner energy of the gridcell i.e. T. in case of rain the
		! fraction that can freeze shall be determined and inner energy will be
		! adjusted too.
		! In case of snowfall the albedo is increased.

		implicit none

		integer, 		intent(in) 		:: ix
		integer, 		intent(in) 		:: iy 
		integer, 		intent(in) 		:: time
		real(kind=8), 	intent(inout) 	:: rainman
		real(kind=8), 	intent(inout) 	:: accum
		integer, 		intent(inout)	:: nday_snowfall
		real(kind=8), 	intent(inout) 	:: H_lh
		real(kind=8), 	intent(inout) 	:: K_lh
		real(kind=8), 	intent(in) 		:: air_temp_ice(nx, ny, ndays)
    	real(kind=8), 	intent(in) 		:: precip_ice(nx, ny, ndays) 

		real(kind=8) :: massum

		H_lh = 0.
        K_lh = 0.
        accum = 0.
        rainman = 0.
        massum = 0.

		if ((air_temp_ice(ix, iy, time) - kelvin <= snow_fall_temperature) .and. &
	        (precip_ice(ix, iy, time) > precip_cutoff/3600./24./1000.)) then
	        ! Snowing

	        ! Precipitation that falls at temperatures smaller than snow_fall_temperature C

	        ! Seeding temperature for first snow
	        if(snowman(ix, iy, 1) == 0) then 
	            snow_temp(ix, iy, 1) = air_temp_ice(ix, iy, time)
	        end if

	        ! Accumulated snow during one timestep in kg/m2/s
	        accum = precip_ice(ix, iy, time)*rho_w
	        
	        nday_snowfall = time
	        if(albedo_module > 0) then
	        	! Based on oerlemans depth scale inverse, but not entirely similar
	            
	            albedo_dynamic(ix, iy) = min(albedo_snow_new, albedo_dynamic(ix, iy) + (albedo_snow_new - albedo_snow_wet)* &
	                                    	(1 - exp(-accum*dt_firn/3.)))
	                        
	            ! Alternative also strange
	            ! albedo_dynamic(ix, iy, 1) = min(albedo_snow_new, albedo_dynamic(ix, iy, 1) + &
	            !								(albedo_snow_new - albedo_dynamic(ix, iy, 1))* &
	            !                            	(1 - exp(-accum*dt_firn))
	        end if
	        
	        rainman = 0.
	        massum = snowman(ix, iy, 1) + accum*dt_firn
	        rho_snow(ix, iy, 1) = massum / (snowman(ix, iy, 1)/rho_snow(ix, iy, 1) + &
	                                accum*dt_firn/rho_s)
	        
	        snowman(ix, iy, 1) = massum

	        ! Latent heat flux for Energybalance
	        H_lh = precip_ice(ix, iy, time)*rho_w*c_i
	        K_lh = precip_ice(ix, iy, time)*rho_w*c_i*air_temp_ice(ix, iy, time)

	    else if(snowman(ix,iy,1) .gt. 0.) then
	        ! Raining

	        accum = 0.
	        rainman = precip_ice(ix, iy, time)*rho_w ! kg/m2/s
	        lwmass(ix, iy, 1) = lwmass(ix, iy, 1) + rainman*dt_firn

	        ! Latent heat flux for Energybalance
	        H_lh = 0.
	        K_lh = precip_ice(ix, iy, time)*rho_w*c_w*(air_temp_ice(ix, iy, time) - kelvin)
	    end if
	end subroutine accumulate
end module precipitation_accumulate