module melting
	! Melting routines for ice and snow, and calculation routines
	! for the surface mass balance
	
	use bessi_defs
	use regridding

	implicit none

	private
	public :: melt_snow
	public :: melt_ice
	public :: smb_point
contains
	subroutine melt_snow(ix, iy, time, T_air, K_sw, K_lh, H_lh, Q_heat, melted_snow, runoff_water, ice_melt, &
							vaporflux, D_lf, dummy_melt_ice, dummy_regrid, used_q, l_heat)
		implicit none

	    integer, 		intent(in) 		:: ix
	    integer, 		intent(in) 		:: iy
	    integer, 		intent(in) 		:: time
	    real(kind=8), 	intent(in) 		:: H_lh
	    real(kind=8), 	intent(in) 		:: T_air
	    real(kind=8), 	intent(in) 		:: K_sw
	    real(kind=8), 	intent(in) 		:: K_lh
	    real(kind=8), 	intent(in) 		:: D_lf
	    real(kind=8), 	intent(inout) 	:: Q_heat
	    real(kind=8), 	intent(inout) 	:: melted_snow
	    real(kind=8), 	intent(inout) 	:: runoff_water
	    real(kind=8), 	intent(inout) 	:: ice_melt
	    real(kind=8), 	intent(inout) 	:: used_q
	    real(kind=8), 	intent(inout) 	:: l_heat
	    real(kind=8), 	intent(inout) 	:: vaporflux
	    real(kind=8), 	intent(inout) 	:: dummy_melt_ice
	    integer, 		intent(inout) 	:: dummy_regrid

	    ! Local variables
	    real(kind=8) :: QQ
	    real(kind=8) :: Qp_lw
	    real(kind=8) :: Qp_sh
	    real(kind=8) :: Qp_lh
	    real(kind=8) :: Q_v
	    real(kind=8) :: dT
	    real(kind=8) :: dm
	    real(kind=8) :: DewpT
	    real(kind=8) :: lwrd_l

	    melted_snow = 0.
	    runoff_water = 0.
	    QQ = 0.
	    l_heat = 0.
	    used_q = 0.

	    ! Melt the snow depending on energy upptaken by the snowcover if it
	    ! reaches 273K
	    if(longwave_from_air_temp) then
	    	lwrd_l = sigma*eps_air*(T_air)**4.
	    else 
	    	lwrd_l = lwrd(ix, iy, time)
	    end if
	    
	    Qp_lw =lwrd_l - sigma*eps_snow*(kelvin)**4.
	    Qp_sh = D_sf*(T_air-kelvin)
	    Qp_lh = K_lh - H_lh*snow_temp(ix, iy, 1)
	    Q_v = vaporflux
	    
	    QQ = max((K_sw + Qp_lw + Qp_sh + Qp_lh + Q_v)*dt_firn - Q_heat, real(0.))
	    used_q = QQ
	    
	    Q_heat = 0.
	    do while (QQ .gt. 0.)
	    	! Maximum possible heating
	        dT = QQ/c_i/snowman(ix, iy, 1)

	        if (kelvin - snow_temp(ix, iy, 1) .gt. dT) then
	        	! QQ too small to melt something
	            
	            snow_temp(ix, iy, 1) = snow_temp(ix, iy, 1) + dT
	            QQ = 0.

	        else if(snowman(ix, iy, 1) .gt. c_i*snowman(ix, iy, 1)*(dT - (kelvin - snow_temp(ix, iy, 1)))/L_lh) then 
	        	! Snow is heated to 0C and partially melted

	            dm = QQ/L_lh - c_i*snowman(ix, iy, 1)*(kelvin - snow_temp(ix, iy, 1))/L_lh
	            l_heat = l_heat + dm*L_lh
	            snow_temp(ix, iy, 1) = kelvin
	            snowman(ix, iy, 1) = snowman(ix, iy, 1) - dm
	            lwmass(ix, iy, 1)  = lwmass(ix, iy, 1) + dm
	            melted_snow = melted_snow + dm
	            QQ = 0.
	        
	        else if(snowman(ix, iy, 2) .gt. 0) then
	        	! The entire grid cell melts 

	            QQ = QQ - c_i*snowman(ix, iy, 1)*(kelvin - snow_temp(ix, iy, 1)) - snowman(ix, iy, 1)*L_lh
	            lwmass(ix, iy, 1)  = lwmass(ix, iy, 1) + snowman(ix, iy, 1)
	            melted_snow = melted_snow + snowman(ix, iy, 1)
	            l_heat = l_heat + snowman(ix, iy, 1)*L_lh
	            snowman(ix, iy, 1) = 0.
	            snow_temp(ix, iy, 1) = snow_temp(ix, iy, 2)
	            rho_snow(ix, iy, 1) = rho_snow(ix, iy, 2)
	            call regrid_point(1,1, dummy_regrid)

	        else
	        	! Entire grid cell melts and no second layer to melt. i.e. reset grid cell

	            QQ = QQ - c_i*snowman(ix, iy, 1)*(kelvin - snow_temp(ix, iy, 1)) - snowman(ix, iy, 1)*L_lh
	            lwmass(ix, iy, 1) = lwmass(ix, iy, 1) + snowman(ix, iy, 1)
	            melted_snow = melted_snow + snowman(ix, iy, 1)
	            
	            ! TODO: if balance is not ok, think about this line again
	            runoff_water = runoff_water + lwmass(ix, iy, 1)    
	            l_heat = l_heat + snowman(ix, iy, 1)*L_lh
	            snowman(ix, iy, 1) = 0.
	            snow_temp(ix, iy, 1) = snow_temp(ix, iy, 2)
	            rho_snow(ix, iy, 1) = rho_snow(ix, iy, 2)
	            
	            dummy_melt_ice= max(-QQ/L_lh,0.)

	            ! This underestimates the melt since the ice is darker than the snow
	            ice_melt= ice_melt + min(-QQ/L_lh/rho_ice/seconds_per_year, 0.)
	            used_q =  used_q - QQ

	            snowman(ix, iy, :) = 0.
	            lwmass(ix, iy, :) = 0.
	            ! Temperature to 273K
	            snow_temp(ix, iy, :) = 0. 
	            rho_snow(ix, iy, :) = rho_s
	            QQ = 0.
	        end if
	    end do
	end subroutine melt_snow

	subroutine melt_ice(ix, iy, time, ice_melt, T_air, precipitation, vaporflux, D_lf, dummy_melt_ice, dummy_rain_ice, &
						p_air)
		! Melt ice in grid point without snow. It is assumed that blank ice has a surface temperature 0C.
		! Water can melt the ice but doesn't freeze, rain and meltwater run off

	    implicit none
	    
	    integer, intent(in) :: ix
	    integer, intent(in) :: iy
	    integer, intent(in) :: time
	    ! Temperature at the level of the ice elevation
	    real(kind=8), intent(in) :: T_air
	    real(kind=8), intent(in) :: precipitation
	    real(kind=8), intent(in) :: p_air
	    real(kind=8), intent(inout) :: ice_melt
	    real(kind=8), intent(inout) :: D_lf
	    real(kind=8), intent(inout) :: dummy_melt_ice
	    real(kind=8), intent(inout) :: dummy_rain_ice
	   	real(kind=8), intent(out) :: vaporflux
	    
	    ! Local variables
		real(kind=8) :: ea
	    real(kind=8) :: es
	    real(kind=8) :: lwrd_l
	    real(kind=8) :: dQ_lh
	    real(kind=8) :: dQ_sh
	    real(kind=8) :: dQ_sw
	    real(kind=8) :: dQ_lw
	    real(kind=8) :: dQ_v
	    real(kind=8) :: dQ_tot

		dQ_lh = 0
		dQ_sh = 0
		dQ_sw = 0
		dQ_lw = 0
		dQ_v = 0
		dQ_tot = 0
	    

	    ea = 610.8*exp(17.27*DewpT(ix, iy, time)/(DewpT(ix, iy, time) + 237.3))
	    es = 611.2*exp(22.46*(kelvin - 273)/(272.62 + (kelvin - 273))) 
	    
	    vaporflux = D_lf/p_air*(ea - es)

	    if (longwave_from_air_temp) then
	        lwrd_l = sigma*eps_air*(T_air)**4.
	    else 
	        lwrd_l = lwrd(ix, iy, time)
	    end if
	    
	    dQ_lw = lwrd_l - sigma*eps_snow*(kelvin)**4.
	    
	    dQ_lh = (T_air - kelvin)*precipitation*c_w*rho_w
	    dQ_sh = D_sf*(T_air - kelvin)
	    dQ_sw = P_sun(ix, iy, time)*(1. - albedo_ice)
	    dQ_v = vaporflux
	    vaporflux = dQ_v
	    dQ_tot = dQ_lw + dQ_sh + dQ_lh + dQ_sw + dQ_v
	    ! Use ice density of ice model to get the right height
	    ice_melt = ice_melt - max(dQ_tot,real(0))/L_lh*dt_firn/rho_ice/seconds_per_year + dQ_v/(L_v + L_lh)* &
	    			dt_firn/rho_ice/seconds_per_year

	    dummy_melt_ice = max(dQ_tot,real(0))/L_lh*dt_firn
	    dummy_rain_ice = precipitation*rho_w

        albedo_dynamic(ix, iy) = albedo_ice
	end subroutine melt_ice

	subroutine smb_point(ix, iy, dummy_runoff, dmass, dummy_ice2ice, dummy_energy, smb_ice)
		! Calculate smb for ice model for a given grid point

		implicit none

		integer, intent(in) :: ix
	    integer, intent(in) :: iy
	    real(kind=8), intent(in) :: dmass
	    real(kind=8), intent(inout) :: dummy_runoff
	    real(kind=8), intent(inout) :: dummy_ice2ice
	    real(kind=8), intent(inout) :: dummy_energy
	    real(kind=8), intent(inout) :: smb_ice(nx, ny)

		smb_ice(ix, iy) = max(dmass, 0.)/rho_ice/seconds_per_year

	    call continuous_snow_removal(ix, iy, dummy_runoff, dmass, dummy_ice2ice, dummy_energy)
	end subroutine smb_point

	subroutine continuous_snow_removal(ix, iy, runoff, dmass, ice2ice, energy_ice2ice)
		! Remove snow at the lower end of the snowmodel and convert it a smb for the ice model

	    implicit none

	    integer, intent(in) :: ix
	    integer, intent(in) :: iy
	    real(kind=8), intent(in) :: dmass ! TODO
	    real(kind=8), intent(inout) :: runoff
	    real(kind=8), intent(inout) :: ice2ice
	    real(kind=8), intent(inout) :: energy_ice2ice

	    ! local variables
	    integer :: nn
	    real(kind=8) :: d_m
	    real(kind=8) :: d_lw
	    real(kind=8) :: perc

	    nn = n_snowlayer
	    perc = 0.
	    d_m = dmass
	    ice2ice = 0.
	    runoff = 0.
	    energy_ice2ice = 0.

	    do while (d_m .gt. 0.)
	        if (d_m .gt. snowman(ix, iy, nn) ) then
	            d_m = d_m - snowman(ix, iy, nn)
	            ice2ice = ice2ice + snowman(ix, iy, nn)
	            runoff = runoff + lwmass(ix, iy, nn)
	            energy_ice2ice = energy_ice2ice + snowman(ix, iy, nn)*c_i*snow_temp(ix, iy, nn)

	            rho_snow(ix, iy, nn) = rho_s
	            snow_temp(ix, iy, nn) = 0.
	            lwmass(ix, iy, nn) = 0.
	            snowman(ix, iy, nn) = 0.
	        else
	            d_lw = d_m * lwmass(ix, iy, nn) / snowman(ix, iy, nn)
	            snowman(ix, iy, nn) = snowman(ix, iy, nn) - d_m
	            ice2ice = ice2ice + d_m
	            energy_ice2ice = energy_ice2ice + d_m*c_i*snow_temp(ix, iy, nn)
	            runoff = runoff + d_lw

	            lwmass(ix, iy, nn) = lwmass(ix, iy, nn) - d_lw
	            d_m = 0.
	        end if

	        nn = nn - 1
	    end do  
	end subroutine continuous_snow_removal

	subroutine go_depleet_snowman(snowman, rho_snow, snow_temp, lwmass, smb_ice, passon_mass)
	    ! Remove the lowest gridboxes with densities larger thant rho_e and add the snow to the smb_ice
	    ! Currently not in use

	    implicit none
	    real(kind=8), intent(inout) :: snowman(nx,ny,n_snowlayer)
	    real(kind=8), intent(inout) :: lwmass(nx,ny,n_snowlayer)
	    real(kind=8), intent(inout) :: snow_temp(nx,ny,n_snowlayer)
	    real(kind=8), intent(inout) :: rho_snow(nx,ny,n_snowlayer)
	    real(kind=8), intent(inout) :: smb_ice(nx,ny)
	    real(kind=8), intent(inout) :: passon_mass

	    ! local variables
	    integer :: mm
	    integer :: ix
	    integer :: iy
	    logical :: move
	    logical :: hit_snow

	    passon_mass = 0
		move = .true.
		hit_snow = .false.

	    do ix = 1, nx, 1
	        do iy = 1, ny, 1
	            move = .true.
	            hit_snow = .false.
	            mm = n_snowlayer + 1

	            ! Empty lowest grid box
	            smb_ice(ix,iy) = smb_ice(ix,iy) + snowman(ix,iy,n_snowlayer)/rho_ice/seconds_per_year

	            ! Reset snowman cell
	            snowman(ix, iy, n_snowlayer) = 0.
	            snow_temp(ix, iy, n_snowlayer) = 0.
	            rho_snow(ix, iy, n_snowlayer) = rho_s
	            lwmass(ix, iy, n_snowlayer) = 0.

	            do while ((move) .and. (mm .gt. 1) .and. (snowman(ix, iy, 1) .gt. 0.))
	                mm = mm - 1

	                if (snowman(ix, iy, mm) .ne. 0.) hit_snow = .true.


	                if(hit_snow) then
	                    if (rho_snow(ix,iy,mm) .ge. rho_pass_on) then
	                        ! Add ice to icesheet
	                        smb_ice(ix,iy) = smb_ice(ix,iy) + snowman(ix,iy,mm)/rho_ice/seconds_per_year
	                                    
	                        ! Reset snowman cell
	                        snowman(ix, iy, mm) = 0.
	                        snow_temp(ix, iy, mm) = 0.
	                        rho_snow(ix, iy, mm) = rho_s
	                        lwmass(ix, iy, mm) = 0.
	                    else
	                        move = .false.
	                    end if
	                end if
	            end do 
	        end do 
	    end do 
	end subroutine go_depleet_snowman
end module melting