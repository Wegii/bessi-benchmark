module percolation
	! Routines for percolation and refreezing of water
	
	use bessi_defs

	implicit none

	private
	public :: percolate
	public :: refreeze
contains
	subroutine percolate(ix, iy, runoff)
		! Calculate water percolating through the firn layers.
	    ! If the whatercontent of one gridcell has grown over a certain age
	    ! of the free bubble space, push everything above downwards.

		implicit none

	    integer, intent(in) :: ix
	    integer, intent(in) :: iy
	    real(kind=8), intent(inout) :: runoff

	    ! Local variables
	    real(kind=8) :: lwc
	    real(kind=8) :: percolating
	    integer :: ii

	   	runoff = 0.
	    ii = 1

	    ! TODO: maybe its better to use a while loop, one could abbort earlier	    
	    ! do ii = 1, n_snowlayer, 1

	   	! Calculate water content [fraction of free volume]
	    do while (ii <= n_snowlayer)
	        if (snowman(ix, iy, ii) > 0.) then
	            if (rho_snow(ix, iy, ii) > rho_i - 10.) then 
	                ! In case of very dense snow percolate all water downwards. 
	                ! we do so in order to avoid division by zero

	                percolating = lwmass(ix, iy, ii)
	                lwmass(ix, iy, ii) = 0.
	                !runoff = percolating
	                
	                if (ii < n_snowlayer) then
	                    if (snowman(ix, iy, ii + 1) > 0.) then
	                        lwmass(ix, iy, ii + 1) = lwmass(ix, iy, ii + 1) + percolating
	                    else
	                        runoff = runoff + percolating
	                    end if
	                else
	                    runoff = runoff + percolating
	                end if

	            else
	                lwc = lwmass(ix, iy, ii)/snowman(ix, iy, ii)/rho_w/(1./rho_snow(ix, iy, ii) - 1./rho_i)
	                
	                if (lwc > max_lwc) then
	                    ! Perform percolation
	                    
	                    percolating = (lwc - max_lwc)*rho_w*snowman(ix, iy, ii)*(1./rho_snow(ix, iy, ii) - 1./rho_i)
	                    lwmass(ix, iy, ii) = lwmass(ix, iy, ii) - percolating
	                    !runoff = percolating
	                    
	                    if (ii < n_snowlayer) then
	                        if (snowman(ix, iy, ii + 1) > 0) then
	                            lwmass(ix, iy, ii + 1) = lwmass(ix, iy, ii + 1) + percolating
	                        else
	                            runoff = runoff + percolating
	                        end if
	                    else
	                        runoff = runoff + percolating
	                    end if
	                end if
	            end if
	        else
	            ii = n_snowlayer
	        end if

	        ii = ii + 1
	    end do
	end subroutine percolate

	subroutine refreeze(ix, iy, dummy_refreeze, heat_fusion)
		! Refreeze water in firn layers

		implicit none

		integer, intent(in) :: ix
	    integer, intent(in) :: iy
	    real(kind=8), intent(inout) :: dummy_refreeze
	    real(kind=8), intent(inout) :: heat_fusion

	    ! Local variables
	    integer :: ii
	    real(kind=8) :: icecube

	    dummy_refreeze = 0.
	    heat_fusion = 0.

	    do ii = 1, n_snowlayer, 1
	        if ((snowman(ix, iy, ii) .gt. 0.) .and. (lwmass(ix, iy, ii) .gt. 0.)) then
	            if ((kelvin - snow_temp(ix, iy, ii))*c_i*snowman(ix, iy, ii) .lt. lwmass(ix, iy, ii)*L_lh) then
	                ! Water freezes partly

	                icecube = (kelvin - snow_temp(ix, iy, ii))*c_i*snowman(ix, iy, ii)/L_lh
	                heat_fusion = heat_fusion + icecube*L_lh

	                snow_temp(ix, iy, ii) = kelvin
	                rho_snow(ix, iy, ii) = rho_snow(ix, iy, ii)*(icecube + snowman(ix, iy, ii)) / snowman(ix, iy, ii)
	                snowman(ix, iy, ii) = snowman(ix, iy, ii) + icecube
	                lwmass(ix, iy, ii) = lwmass(ix, iy, ii) - icecube
	                dummy_refreeze = dummy_refreeze + icecube
	            else
	                ! All water freezes
	                snow_temp(ix, iy, ii) = (lwmass(ix, iy, ii)*L_lh/c_i + lwmass(ix, iy, ii)*kelvin + &
	                						snow_temp(ix, iy, ii)*snowman(ix, iy, ii) )/(lwmass(ix, iy, ii) + snowman(ix, iy, ii))
	                    
	                rho_snow(ix, iy, ii) = rho_snow(ix, iy, ii)*(lwmass(ix, iy, ii) + snowman(ix, iy, ii)) / snowman(ix, iy, ii)
	                snowman(ix, iy, ii) = lwmass(ix, iy, ii) + snowman(ix, iy, ii)
	                dummy_refreeze = dummy_refreeze + lwmass(ix, iy, ii)
	                heat_fusion = heat_fusion + lwmass(ix, iy, ii)*L_lh

	                lwmass(ix, iy, ii) = 0.
	            end if
	        end if
	    end do
	end subroutine refreeze
end module percolation
