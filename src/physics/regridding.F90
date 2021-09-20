module regridding
	! Routines for regridding, fusing and splitting of boxes

	use bessi_defs

	implicit none

	private
	public :: regrid
	public :: regrid_point
contains
	subroutine regrid(ix, iy, dummy_regrid)
		! Fuse and split boxes if necessary

		implicit none

		integer, intent(in) :: ix
	    integer, intent(in) :: iy
	    integer, intent(inout) :: dummy_regrid

	    if ((snowman(ix, iy, 1) .gt. upper_massbound) .or. (snowman(ix, iy, 1) .lt. lower_massbound)) then
          	! Fuse and split boxes
            call regrid_point(ix, iy, dummy_regrid)

            do while((snowman(ix, iy, 1) .gt. upper_massbound))
                call regrid_point(ix, iy, dummy_regrid)
            end do
        end if
	end subroutine regrid

	subroutine regrid_point(ix, iy, dummy_regrid)
		! Adjust boxsizes and move boxes arround, 
		! depending on the soll_mass.

	    implicit none
	    
	    integer, intent(in) :: ix
	    integer, intent(in) :: iy
	    integer, intent(inout) :: dummy_regrid

	    real(kind=8) :: masssum
	    integer :: ii
	    integer :: kk

	    dummy_regrid = dummy_regrid + 1
	    
	    if (snowman(ix, iy, 1) .gt. upper_massbound) then
	    	! Split top box (accumulating)
	    	! Fuse the two lowest boxes
	        masssum = snowman(ix, iy, n_snowlayer) + snowman(ix, iy, n_snowlayer-1)

	        if (masssum == 0) then
	            rho_snow(ix, iy, n_snowlayer) = rho_s
	            snow_temp(ix, iy, n_snowlayer) = 0.
	        else if(snowman(ix, iy, n_snowlayer) == 0) then
	            rho_snow(ix, iy, n_snowlayer) = rho_snow(ix, iy, n_snowlayer - 1)
	            snow_temp(ix, iy, n_snowlayer) = snow_temp(ix, iy, n_snowlayer - 1)
	        else
	            rho_snow(ix, iy, n_snowlayer) = masssum / (snowman(ix, iy, n_snowlayer)/rho_snow(ix, iy, n_snowlayer) + &
	            						snowman(ix, iy, n_snowlayer - 1)/rho_snow(ix, iy, n_snowlayer - 1))

	            snow_temp(ix, iy, n_snowlayer) = (snowman(ix, iy, n_snowlayer)*snow_temp(ix, iy, n_snowlayer) + &
	            							snowman(ix, iy, n_snowlayer - 1)*snow_temp(ix, iy, n_snowlayer - 1))/masssum
	        end if

	        snowman(ix, iy, n_snowlayer) = masssum
	        lwmass(ix, iy, n_snowlayer) = lwmass(ix, iy, n_snowlayer) + lwmass(ix, iy, n_snowlayer - 1)

	        ! Push all layers down, except the toplayer
	        kk = n_snowlayer - 1;
	        do while (kk .gt. 2)
	            snowman(ix, iy, kk) = snowman(ix, iy, kk - 1)
	            rho_snow(ix, iy, kk) = rho_snow(ix, iy, kk - 1)
	            snow_temp(ix, iy, kk) = snow_temp(ix, iy, kk - 1)
	            lwmass(ix, iy, kk) = lwmass(ix, iy, kk - 1)

	            kk = kk - 1
	        end do

	        ! Split top box
	        ! Move proper amount of content down to second box
	        snowman(ix, iy, 2) = soll_mass
	        rho_snow(ix, iy, 2) = rho_snow(ix, iy, 1)
	        snow_temp(ix, iy, 2) = snow_temp(ix, iy, 1)
	        lwmass(ix, iy, 2) = lwmass(ix, iy, 1)*soll_mass/snowman(ix, iy, 1)

	        ! Adjust content in topbox (except density and temperature)
	        snowman(ix, iy, 1) = snowman(ix, iy, 1) - soll_mass
	        lwmass(ix, iy, 1) = lwmass(ix, iy, 1) - lwmass(ix, iy, 2)
	    else if ((snowman(ix, iy, 1) .lt. lower_massbound) .AND. (snowman(ix, iy, 2) .gt. 0.)) then
	    	! Fusing top boxes (melting)
	        ! Splitting lowest box in a meaningfull way and fuse both uppermost boxes if there are two

	        masssum = snowman(ix, iy, 1) + snowman(ix, iy, 2)

	        rho_snow(ix, iy, 1) = masssum / (snowman(ix, iy, 1)/rho_snow(ix, iy, 1) + snowman(ix, iy, 2)/rho_snow(ix, iy, 2))
	        snow_temp(ix, iy, 1) = (snowman(ix, iy, 1)*snow_temp(ix, iy, 1) + snowman(ix, iy, 2)*snow_temp(ix, iy, 2))/masssum
	        snowman(ix, iy, 1) = masssum
	        lwmass(ix, iy, 1) = lwmass(ix, iy, 1) + lwmass(ix, iy, 2)

	        ! Drag all layers up, except the toplayer
	        do ii = 2, n_snowlayer-1, 1
	            snowman(ix, iy, ii) = snowman(ix, iy, ii + 1)
	            snow_temp(ix, iy, ii) = snow_temp(ix, iy, ii + 1)
	            lwmass(ix, iy, ii) = lwmass(ix, iy, ii + 1)
	            rho_snow(ix, iy, ii) = rho_snow(ix, iy, ii + 1)
	        end do

	        ! Split second lowest box in 2. The upper should have the
	        ! soll_mass, the lowest not less than the soll_mass
	        if(snowman(ix, iy, n_snowlayer - 1) .ge. 2*soll_mass) then
	            ! Move down to lowest box (density and temperature do not change)
	            snowman(ix, iy, n_snowlayer) = snowman(ix, iy, n_snowlayer - 1) - soll_mass
	            lwmass(ix, iy, n_snowlayer)= (1. - soll_mass/snowman(ix, iy, n_snowlayer - 1))*lwmass(ix, iy, n_snowlayer - 1)

	            ! Adjust content in second lowest box (except density and temperature)
	            snowman(ix, iy, n_snowlayer - 1) = soll_mass
	            lwmass(ix, iy, n_snowlayer - 1)= lwmass(ix, iy, n_snowlayer - 1) - lwmass(ix, iy, n_snowlayer)

	        else
	            ! Reset lowest box
	            rho_snow(ix, iy, n_snowlayer) = rho_s
	            snowman(ix, iy, n_snowlayer) = 0.
	            snow_temp(ix, iy, n_snowlayer) = 0.
	            lwmass(ix, iy, n_snowlayer) = 0.
	        end if
	    end if
	end subroutine regrid_point
end module regridding