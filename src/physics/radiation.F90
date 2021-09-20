module radiation
	! Routine for calculating the incoming radiation

	use bessi_defs

	implicit none

	private
	public :: radiate_point
contains
	subroutine radiate_point(ix, iy, nday_snowfall, time, K_sw, lwc)
		! Calculated the incoming solar radiation based on the current
		! surface temperature and wetness (liquid water content) and
		! calculation the snow albedo for dry and wet snow

		! TODO: What does happen to the albedo if the snow became wet, but it got colder again, but no precipitation
		! TODO: Remove duplicate lines (e.g. K_sw = ...) -> There is a lot

		implicit none

		integer, 		intent(in) 		:: ix
	   	integer, 		intent(in) 		:: iy
		integer, 		intent(in)		:: nday_snowfall
		integer, 		intent(in) 		:: time
	    real(kind=8), 	intent(inout) 	:: K_sw !short-wave radiation absorbed by the snow cover (first cell)
	    real(kind=8), 	intent(in) 		:: lwc
	    
	    real(16), parameter :: PI_8 = 4 * atan (1.0_8) 
	    real(8) :: tempa

	    tempa = albedo_dynamic(ix, iy)

        if (albedo_module == 1) then
        	! Classic calculation

            if (snow_temp(ix, iy, 1) < kelvin) then
                albedo_dynamic(ix, iy) = albedo_snow_new
                K_sw = P_sun(ix, iy, time)*(1 - albedo_dynamic(ix, iy))
            else
                albedo_dynamic(ix, iy) = albedo_snow_wet
                K_sw = P_sun(ix, iy, time)*(1 - albedo_dynamic(ix, iy))
    
            end if  
        elseif (albedo_module == 2) then
        	! Harmonic, for now simple with a sinus and maxium

            if (snow_temp(ix, iy, 1) < kelvin) then
                albedo_dynamic(ix, iy) = min(albedo_dynamic(ix, iy), 0.9 - sin(time/360.*PI_8)*0.25) !avg of 0.74albedo 
                K_sw = P_sun(ix, iy, time)*(1 - albedo_dynamic(ix, iy))
            else
                albedo_dynamic(ix, iy) = albedo_snow_wet
                K_sw = P_sun(ix, iy, time)*(1 - albedo_dynamic(ix, iy))
            end if   
        elseif (albedo_module == 3) then 
        	! Temporal decay temperature INdependent Oerlemans and Kap 1998 neglecting the depth scale, which on longer
        	! timescales will be ingnored

            if (snow_temp(ix, iy, 1) < kelvin) then 
                albedo_dynamic(ix, iy) = min(tempa,albedo_snow_wet + (albedo_snow_new - albedo_snow_wet)* &
                							exp((nday_snowfall - time)/3.))
                K_sw = P_sun(ix, iy, time)*(1 - albedo_dynamic(ix, iy))
            else
                albedo_dynamic(ix, iy) = min(tempa,albedo_snow_wet + (albedo_snow_new - albedo_snow_wet)* &
                							exp((nday_snowfall - time)/5.))
                K_sw = P_sun(ix, iy, time)*(1 - albedo_dynamic(ix, iy))
            end if   
        elseif (albedo_module == 5) then
        	! Temporal decay Bougamont (Oerlemans and Kap 1998) neglecting the depth scale, LWC based reduced albedo

        	! Uses a linear increase in reduction rate depending on the lwc of the snowlayer, exponential possible for
            if (snow_temp(ix, iy, 1) < kelvin - 10) then 
                albedo_dynamic(ix, iy) = min(tempa,albedo_snow_wet + (albedo_snow_new - albedo_snow_wet)* &
                							exp((nday_snowfall-time)/100.))
                K_sw = P_sun(ix, iy, time)*(1 - albedo_dynamic(ix, iy))
            elseif (snow_temp(ix, iy, 1) < kelvin) then
                albedo_dynamic(ix, iy) = min(tempa,albedo_snow_wet + (albedo_snow_new - albedo_snow_wet)* &
                							exp((nday_snowfall - time)/(30. + 7*abs(snow_temp(ix, iy, 1) - kelvin)))) 
                K_sw = P_sun(ix, iy, time)*(1 - albedo_dynamic(ix, iy))
            else
                albedo_dynamic(ix, iy) = min(tempa, albedo_snow_wet + (albedo_snow_new - albedo_snow_wet)* &
                							exp((nday_snowfall - time)/(15. - 14*(lwc/max_lwc)))) 
                ! Option use even lower value for wet snow with max lwc, problem for lwc/max>1
                K_sw = P_sun(ix, iy, time)*(1 - albedo_dynamic(ix, iy))
            end if     
        elseif (albedo_module == 4) then
        	! Calculate dry snow albedo reduction based on Aoki 2003 suggest weaker effect due to fewer impurities 
        	! in Antarctica and remote sea ice

            albedo_dynamic(ix, iy) = min(tempa, tempa - ((snow_temp(ix, iy, 1) - kelvin)*1.35e-3 + 0.0278))
            if (albedo_dynamic(ix, iy) < albedo_snow_wet) then
                albedo_dynamic(ix, iy) = albedo_snow_wet
            end if

            if (lwc > 0) then
	            tempa = albedo_dynamic(ix, iy)
	            albedo_dynamic(ix, iy) = max(albedo_snow_wet, min(tempa, tempa - (tempa - albedo_snow_wet)*(lwc/max_lwc)))
	            K_sw = P_sun(ix, iy, time)*(1 - albedo_dynamic(ix, iy))
            end if
            K_sw = P_sun(ix, iy, time)*(1 - albedo_dynamic(ix, iy))
        elseif(albedo_module==0) then
            K_sw = P_sun(ix, iy, time)*(1 - albedo_dynamic(ix, iy))
        end if
	end subroutine radiate_point
end module radiation