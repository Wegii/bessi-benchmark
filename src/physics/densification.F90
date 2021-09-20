module densification
	! Densify layers depending on the densification_model and
	! layer content

	use bessi_defs

	implicit none

	private
	public :: densify
contains
	subroutine densify(ix, iy, At, dz)
		! Densify layers

	    implicit none

	    integer, 		intent(in) 		:: ix
	   	integer, 		intent(in) 		:: iy
	    real(kind=8), 	intent(in) 		:: At
    	real(kind=8), 	intent(inout) 	:: dz(n_snowlayer)

	    real(kind=8) :: ddens
	    real(kind=8) :: f
	    real(kind=8) :: columnsnow
	    real(kind=8) :: P_ice
	    real(kind=8) :: P_bubble
	    real(kind=8) :: dp
	    integer :: mm
	    integer :: im

	    if (densification_model .and. (snowman(ix, iy,  3) .gt. 0.)) then
		    do mm = 1, n_snowlayer, 1
		        if(rho_snow(ix, iy, mm) .lt. 550.) then
		            ! densification  H-L
		            if(snowman(ix, iy, mm) > 0.) then

		            	! is of order 7e-7
		                ddens = 0.011*exp(-10160.0/8.13/snow_temp(ix, iy, mm))*(rho_i - rho_snow(ix, iy, mm))*max(At, real(0))

		                rho_snow(ix, iy, mm) = max(rho_snow(ix, iy, mm), rho_snow(ix, iy, mm) + ddens*dt_firn)
		                rho_snow(ix, iy, mm) = min(rho_snow(ix, iy, mm), rho_i)
		            end if
		        else if (hl == 1) then ! H-L for rho>550
		            if(mm > 1) then
		                if(snowman(ix, iy, mm)>0.) then

		                	! is of order 2e-7
		                    ddens = 0.575*exp(-21400.0/8.13/snow_temp(ix, iy, mm))*(rho_i - rho_snow(ix, iy, mm))* &
		                    		(1000./3600./24./365.)**(0.5)*(max(real(0), At))**0.5

		                    rho_snow(ix, iy, mm) = max(rho_snow(ix, iy, mm), rho_snow(ix, iy, mm) + ddens*dt_firn)
		                    rho_snow(ix, iy, mm) = min(rho_snow(ix, iy, mm), rho_i)
		                end if
		            end if

		        else if((rho_snow(ix, iy, mm) .lt. 800.) .and. (rho_snow(ix, iy, mm) .ge. 550.) ) then
		            ! Densification time for rho 550-800 P-B schwander et all
		            f = 10.**(-29.166*(rho_snow(ix, iy, mm)/rho_i)**3. + 84.422*(rho_snow(ix, iy, mm)/rho_i)**2.-87.425* &
		            	(rho_snow(ix, iy, mm)/rho_i) + 30.673)
		            
		            columnsnow = 0
		            if (mm > 1) then
		                do im = 2, mm, 1
		                    columnsnow = columnsnow + snowman(ix, iy, im - 1)
		                end do
		            end if

		            P_ice = 9.81*(columnsnow+snowman(ix, iy, mm)/2.)/1e6
		            if(rho_snow(ix, iy, mm)>rho_e) then
		                P_bubble = P_atm*((1./rho_e-1./rho_i)/(1./rho_snow(ix, iy, mm) - 1./rho_i) - 1.)/1e6
		            else
		                P_bubble = 0
		            end if
		            
		            dp = P_ice - P_bubble
		            ddens = 25400.*exp(-60000./8.13/snow_temp(ix, iy, mm))*(rho_snow(ix, iy, mm))*f *(dp)**3.
		            
		            rho_snow(ix, iy, mm) = max(rho_snow(ix, iy, mm), rho_snow(ix, iy, mm) + ddens*dt_firn)
		            rho_snow(ix, iy, mm) = min(rho_snow(ix, iy, mm), rho_i)

		        else if(rho_snow(ix, iy, mm) .gt. 800.) then
		            ! Densification time for rho > 800 P-B schwander et all
		            f = 3./16.*(1. - (rho_snow(ix, iy, mm)/rho_i))/(1. - (1. - (rho_snow(ix, iy, mm)/rho_i))**(1./3.))**3.
		            
		            columnsnow = 0.
		            if (mm > 1) then
		                do im = 2, mm, 1
		                    columnsnow = columnsnow + snowman(ix, iy, im - 1)
		                end do
		            end if
		            
		            P_ice = 9.81*(columnsnow + snowman(ix, iy, mm)/2.)/1e6
		            if (rho_snow(ix, iy, mm) > rho_e) then
		                P_bubble = P_atm*((1./rho_e-1./rho_i)/(1./rho_snow(ix, iy, mm) - 1./rho_i) - 1.)/1e6
		            else
		                P_bubble = 0.
		            end if

		            dp = P_ice - P_bubble
		            ddens = 25400.*exp(-60000./8.13/snow_temp(ix, iy, mm))*(rho_snow(ix, iy, mm))*f *(dp)**3.

		            rho_snow(ix, iy, mm) = max(rho_snow(ix, iy, mm), rho_snow(ix, iy, mm) + ddens*dt_firn)
		            rho_snow(ix, iy, mm) = min(rho_snow(ix, iy, mm), rho_i)
		        end if
		    end do
		end if

		dz(:) = snowman(ix, iy, :)/rho_snow(ix, iy, :)
	end subroutine densify
end module densification