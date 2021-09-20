module bessi_data
	! Data and writing module for annual, monthly, and daily variables.
	! Module stores all variables to be written to disk on a time basis,
	! and contains routines for writing the variables to disk 

	use bessi_defs
	use io
	use conservation

	implicit none

	! Annual data and variable IDs
    real(kind=4) :: m_snowman_annual(nx, ny, n_snowlayer)
    real(kind=4) :: m_lwmass_annual(nx, ny, n_snowlayer)
    real(kind=4) :: m_rho_snow_annual(nx, ny, n_snowlayer)
    real(kind=4) :: m_snow_temp_annual(nx, ny, n_snowlayer)
    real(kind=4) :: m_albedo_dynamic_annual(nx, ny)
    real(kind=4) :: m_snow_temp_ave_surface_annual(nx, ny)
    real(kind=4) :: m_vaporflux_annual(nx, ny)
    real(kind=4) :: m_accum_annual(nx, ny)
    real(kind=4) :: m_rain_annual(nx, ny)
    real(kind=4) :: m_refreezing_annual(nx, ny)
    real(kind=4) :: m_melt_annual(nx, ny)
    real(kind=4) :: m_snow_annual(nx, ny)
    real(kind=4) :: m_runoff_annual(nx, ny)
    real(kind=4) :: m_melt_ice_annual(nx, ny)
    integer      :: snow_mask_annual(nx, ny)
    integer      :: m_regridding_annual(nx, ny)
    real(kind=4) :: m_real_mass_balance_annual(nx, ny)
    integer :: filehandle_netcdf_annual
    character(128) :: netcdf_output_filename_annual
    ! Variable for the netCDF snow mass, parameter
    integer :: snowman_annual_varid
    integer :: rho_snow_annual_varid
    integer :: snow_temp_annual_varid
    integer :: lwmass_annual_varid
    integer :: accum_annual_varid
    ! Variable for the netCDF snow mass, parameter
    integer :: rain_annual_varid
    ! Variable for the netCDF density of snow, parameter
    integer :: refreezing_annual_varid
    ! Variable for the netCDF density of snow, parameter
    integer :: melt_annual_varid
    ! Variable for the netCDF snow, parameter
    integer :: snow_annual_varid
    ! Variable for the netCDF runoff of snow, parameter
    integer :: runoff_annual_varid
    integer :: latent_heat_annual_varid
    integer :: albedo_dynamic_annual_varid
    ! Variable snow_mask
    integer :: snow_mask_annual_varid
    integer :: regridding_annual_varid
    integer :: melt_ice_annual_varid
    integer :: snow_temp_ave_surface_annual_varid
    integer :: real_mass_balance_annual_varid
    
    ! Monthly data and variable IDs
    real(kind=4) :: m_snowman_monthly(12, n_snowlayer)
    real(kind=4) :: m_lwmass_monthly(12, n_snowlayer)
    real(kind=4) :: m_rho_snow_monthly(12, n_snowlayer)
    real(kind=4) :: m_snow_temp_monthly(12, n_snowlayer)
    real(kind=4) :: m_albedo_dynamic_monthly(12)
    real(kind=4) :: m_snow_temp_ave_surface_monthly(12)
    real(kind=4) :: m_vaporflux_monthly(12)
    real(kind=4) :: m_accum_monthly(12)
    real(kind=4) :: m_rain_monthly(12)
    real(kind=4) :: m_refreezing_monthly(12)
    real(kind=4) :: m_melt_monthly(12)
    real(kind=4) :: m_snow_monthly(12)
    real(kind=4) :: m_runoff_monthly(12)
    real(kind=4) :: m_melt_ice_monthly(12)
    integer     :: m_snow_mask_monthly(12)
    integer     :: m_regridding_monthly(12)
    real(kind=4) :: m_real_mass_balance_monthly(12)
    integer :: filehandle_netcdf_monthly
    character(128) :: netcdf_output_filename_monthly
    integer :: snowman_monthly_varid
    integer :: rho_snow_monthly_varid
    integer :: snow_temp_monthly_varid
    integer :: lwmass_monthly_varid
    integer :: accum_monthly_varid
    ! Variable for the netCDF snow mass, parameter
    integer :: rain_monthly_varid
    ! Variable for the netCDF density of snow, parameter
    integer :: refreezing_monthly_varid
    ! Variable for the netCDF density of snow, parameter
    integer :: melt_monthly_varid
    ! Variable for the netCDF snow, parameter
    integer :: snow_monthly_varid
    ! Variable for the netCDF runoff of snow, parameter
    integer :: runoff_monthly_varid
    integer :: latent_heat_monthly_varid
    integer :: albedo_dynamic_monthly_varid
    ! Variable snow_mask
    integer :: snow_mask_monthly_varid
    integer :: regridding_monthly_varid
    integer :: melt_ice_monthly_varid
    integer :: snow_temp_ave_surface_monthly_varid
    integer :: real_mass_balance_monthly_varid
    integer, dimension(12) ::  A = (/ 1, 32, 60, 91, 121, 152, 182, 213, 244, 274,305, 335 /)
    integer, dimension(12) ::  B = (/ 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365 /)

    ! Daily data and variable IDs
    real(kind=4) :: m_snowman_daily(ndays, n_snowlayer)
    real(kind=4) :: m_lwmass_daily(ndays, n_snowlayer)
    real(kind=4) :: m_rho_snow_daily(ndays, n_snowlayer)
    real(kind=4) :: m_snow_temp_daily(ndays, n_snowlayer)
    real(kind=4) :: m_albedo_dynamic_daily(ndays)
    real(kind=4) :: m_snow_temp_ave_surface_daily(ndays)
    real(kind=4) :: m_vaporflux_daily(ndays)
    real(kind=4) :: m_accum_daily(ndays)
    real(kind=4) :: m_rain_daily(ndays)
    real(kind=4) :: m_refreezing_daily(ndays)
    real(kind=4) :: m_melt_daily(ndays)
    real(kind=4) :: m_snow_daily(ndays)
    real(kind=4) :: m_runoff_daily(ndays)
    real(kind=4) :: m_melt_ice_daily(ndays)
    integer      :: m_snow_mask_daily(ndays)
    integer      :: m_regridding_daily(ndays)
    real(kind=4) :: m_real_mass_balance_daily(ndays)
    integer :: filehandle_netcdf_daily
    character(128) :: netcdf_output_filename_daily
    ! Variable for the netCDF snow mass, parameter    
    integer :: snowman_daily_varid
    integer :: rho_snow_daily_varid
    integer :: snow_temp_daily_varid
    integer :: lwmass_daily_varid
    integer :: accum_daily_varid
    ! Variable for the netCDF snow mass, parameter
    integer :: rain_daily_varid
    ! Variable for the netCDF density of snow, parameter
    integer :: refreezing_daily_varid
    ! Variable for the netCDF density of snow, parameter
    integer :: melt_daily_varid
    ! Variable for the netCDF snow, parameter
    integer :: snow_daily_varid
    ! Variable for the netCDF runoff of snow, parameter
    integer :: runoff_daily_varid
    integer :: latent_heat_daily_varid
    integer :: albedo_dynamic_daily_varid
    ! Variable snow_mask
    integer :: snow_mask_daily_varid
    integer :: regridding_daily_varid
    integer :: melt_ice_daily_varid
    integer :: snow_temp_ave_surface_daily_varid
    integer :: real_mass_balance_daily_varid
    
	private
	public :: bessi_write_init
	public :: bessi_store_values
	public :: bessi_store_smb
	public :: bessi_write_values_annual
	public :: bessi_write_values_monthly
	public :: bessi_write_values_daily
	public :: write_check_mass
	public :: write_check_energy
contains

	subroutine bessi_write_init(year)
		! Initialize annual, monthly, and daily ncdf file and
		! allocate respective arrays

		implicit none

        integer, intent(in) :: year

	    ! Allocate annual output variables and netcdf output file
	    if (annual_data .and. (((year .ge. daily_data_start) .and. (year .le. daily_data_end)) &
			.or.(mod(year, annual_data_frequency) == 0))  ) then

	        print*,'*** Inititalizing files for annual data', year
	        write (netcdf_output_filename_annual, '( "ANNUAL_", I7.7, ".nc" )') year

	        call init_netcdf_snow_file(TRIM(adjustl(output_directory)) // netcdf_output_filename_annual, &
	            filehandle_netcdf_annual, ny, nx, n_snowlayer, 1, int(dy), int(dx), 1, snowman_annual_varid, &
	            lwmass_annual_varid, rho_snow_annual_varid, snow_temp_annual_varid, albedo_dynamic_annual_varid, &
	            latent_heat_annual_varid, accum_annual_varid, rain_annual_varid, melt_annual_varid, &
	            refreezing_annual_varid, snow_annual_varid, runoff_annual_varid, melt_ice_annual_varid, &
	            snow_mask_annual_varid, snow_temp_ave_surface_annual_varid, regridding_annual_varid, &
	            real_mass_balance_annual_varid, new_input_file)
	        
			m_snowman_annual = 0.
			m_lwmass_annual = 0.
			m_rho_snow_annual = 0.
			m_snow_temp_annual = 0.
			m_albedo_dynamic_annual = 0.
			m_snow_temp_ave_surface_annual = 0.
			m_vaporflux_annual = 0.
			m_accum_annual = 0.
			m_rain_annual = 0.
			m_refreezing_annual = 0.
			m_melt_annual = 0.
			m_snow_annual = 0.
			m_runoff_annual = 0.
			m_melt_ice_annual = 0.
			snow_mask_annual = 0
			m_regridding_annual = 0
			m_real_mass_balance_annual = 0.
	    end if

	    ! Allocate monthly output variables and netcdf output file
    	if (monthly_data .and. (((year .ge. daily_data_start) .and. (year .le. daily_data_end)) &
			.or. (mod(year, monthly_data_frequency) == 0))) then
	        
	        print*,'*** Inititalizing files for monthly data',year
	        write (netcdf_output_filename_monthly, '( "MONTHLY_", I7.7, ".nc" )') year
	   
	        call init_netcdf_snow_file_time_first(TRIM(adjustl(output_directory)) // netcdf_output_filename_monthly, &
	            filehandle_netcdf_monthly, ny, nx, n_snowlayer, 12, int(dy), int(dx), 1, snowman_monthly_varid, &
	            lwmass_monthly_varid, rho_snow_monthly_varid, snow_temp_monthly_varid, albedo_dynamic_monthly_varid, &
	            latent_heat_monthly_varid, accum_monthly_varid, rain_monthly_varid, melt_monthly_varid, &
	            refreezing_monthly_varid, snow_monthly_varid, runoff_monthly_varid, melt_ice_monthly_varid, &
	            snow_mask_monthly_varid, snow_temp_ave_surface_monthly_varid, regridding_monthly_varid, &
	            real_mass_balance_monthly_varid, new_input_file)
	        
	        m_snowman_monthly = 0.
			m_lwmass_monthly = 0.
			m_rho_snow_monthly = 0.
			m_snow_temp_monthly = 0.
			m_albedo_dynamic_monthly = 0.
			m_snow_temp_ave_surface_monthly = 0.
			m_vaporflux_monthly = 0.
			m_accum_monthly = 0.
			m_rain_monthly = 0.
			m_refreezing_monthly = 0.
			m_melt_monthly = 0.
			m_snow_monthly = 0.
			m_runoff_monthly = 0.
			m_melt_ice_monthly = 0.
			m_snow_mask_monthly = 0
			m_regridding_monthly = 0
			m_real_mass_balance_monthly = 0.
	    end if

	    ! Allocate daily output variables and netcdf output file
		if ((monthly_data .or. daily_data) .and. (((year .ge. daily_data_start) .and. (year .le. daily_data_end)) &
			.or. (mod(year, monthly_data_frequency) == 0 .or. mod(year, daily_data_frequency) == 0))) then

	        print*,'*** Inititalizing files for daily data', year
	        write (netcdf_output_filename_daily, '( "DAILY_", I7.7, ".nc" )') year
	    
	        if (daily_data .and. (((year .ge. daily_data_start) .and. (year .le. daily_data_end)) &
	            .or. (mod(year,daily_data_frequency) == 0))) then
	            call init_netcdf_snow_file_time_first(TRIM(adjustl(output_directory)) // netcdf_output_filename_daily, &
	                filehandle_netcdf_daily, ny, nx, n_snowlayer, 365, int(dy), int(dx), 1, snowman_daily_varid, &
	                lwmass_daily_varid, rho_snow_daily_varid, snow_temp_daily_varid, albedo_dynamic_daily_varid, &
	                latent_heat_daily_varid, accum_daily_varid, rain_daily_varid, melt_daily_varid, &
	                refreezing_daily_varid, snow_daily_varid, runoff_daily_varid, melt_ice_daily_varid, &
	                snow_mask_daily_varid, snow_temp_ave_surface_daily_varid, regridding_daily_varid, &
	                real_mass_balance_daily_varid, new_input_file)
	        end if
	        
	        m_snowman_daily = 0.
			m_lwmass_daily = 0. 
			m_rho_snow_daily = 0.
			m_snow_temp_daily = 0.
			m_albedo_dynamic_daily = 0.
			m_snow_temp_ave_surface_daily = 0.
			m_vaporflux_daily = 0.
			m_accum_daily = 0.
			m_rain_daily = 0.
			m_refreezing_daily = 0.
			m_melt_daily = 0.
			m_snow_daily = 0.
			m_runoff_daily = 0.
			m_melt_ice_daily = 0.
			m_snow_mask_daily = 0
			m_regridding_daily = 0
			m_real_mass_balance_daily = 0.
	    end if
	end subroutine bessi_write_init

	subroutine bessi_store_values(year, ix, iy, time, accum, vaporflux, mass0, dummy_runoff, dummy_refreeze, &
									dummy_melt, dummy_melt_ice, dummy_regrid, dummy_rain_ice, smb_ice, rainman, &
									albedo_runtime)
		! Store annual and daily data.
		! All values are accumulated and stored in their respective variable

		implicit none

		integer, 		intent(in) 	:: year
		integer, 		intent(in) 	:: ix
		integer, 		intent(in) 	:: iy 
		integer, 		intent(in) 	:: time
		real(kind=8), 	intent(in) 	:: accum
		real(kind=8), 	intent(in) 	:: vaporflux
		real(kind=8), 	intent(in) 	:: mass0
		real(kind=8), 	intent(in) 	:: dummy_runoff
		real(kind=8), 	intent(in) 	:: dummy_refreeze
		real(kind=8), 	intent(in) 	:: dummy_melt
		real(kind=8), 	intent(in) 	:: dummy_melt_ice
		real(kind=8), 	intent(in) 	:: dummy_rain_ice
		integer, 		intent(in) 	:: dummy_regrid
	    real(kind=8),	intent(in) 	:: rainman
		real(kind=8), 	intent(in) 	:: smb_ice(nx, ny)
		real(kind=8),	intent(inout)  :: albedo_runtime(nx, ny) 

		! Accumulate annual values
		if (annual_data .and. (((year .ge. daily_data_start) .and. (year .le. daily_data_end)) &
        	.or. (mod(year, annual_data_frequency) == 0))) then

	        m_snow_annual(ix, iy) = m_snow_annual(ix, iy) + real(accum*dt_firn)
	        m_rain_annual(ix, iy) = m_rain_annual(ix, iy) + real(rainman*dt_firn)   !rain on ice?
	        m_vaporflux_annual(ix, iy) = m_vaporflux_annual(ix,iy) + real(vaporflux*dt_firn/(L_v + L_lh))
	        m_runoff_annual(ix, iy) = m_runoff_annual(ix, iy) + real(dummy_runoff)
	        m_refreezing_annual(ix, iy) = m_refreezing_annual(ix, iy) + real(dummy_refreeze)
	        m_melt_annual(ix, iy) = m_melt_annual(ix, iy) + real(dummy_melt) !melt on ice
	        m_melt_ice_annual(ix, iy) = m_melt_ice_annual(ix, iy) + real(dummy_melt_ice)
	        albedo_runtime(ix, iy) = albedo_runtime(ix, iy) + albedo_dynamic(ix, iy)
	        m_snow_temp_ave_surface_annual(ix, iy) = m_snow_temp_ave_surface_annual(ix, iy) + snow_temp(ix, iy, 1)
	        m_regridding_annual(ix, iy) = m_regridding_annual(ix, iy) + dummy_regrid

	        if (include_values_over_ice) then
	            m_runoff_annual(ix, iy) = m_runoff_annual(ix, iy) + real(dummy_melt_ice)
	            m_rain_annual(ix, iy) = m_rain_annual(ix, iy) + real(dummy_rain_ice*dt_firn)
	        end if

	        if (snowman(ix, iy, 1) .gt. 0.) snow_mask_annual(ix,iy) = snow_mask_annual(ix, iy) + 1
	    end if

	    ! Store daily values
	    if ((monthly_data .or. daily_data) .and. (((year .ge. daily_data_start) .and. (year .le. daily_data_end)) &
            .or. (mod(year, monthly_data_frequency) == 0 .or. mod(year, daily_data_frequency) == 0))) then

            m_snowman_daily(time, 1:15) = real(snowman(ix, iy, :))
            m_lwmass_daily(time, 1:15) = real(lwmass(ix, iy, :))
            m_rho_snow_daily(time, 1:15) = real(rho_snow(ix, iy, :))
            m_snow_temp_daily(time, 1:15) = real(snow_temp(ix, iy, :))
            
            if (m_snow_temp_daily(time, 1) == 0) m_snow_temp_daily(time, 1:15) = 273.16

            m_albedo_dynamic_daily(time) = real(albedo_dynamic(ix, iy))
            m_vaporflux_daily(time) = real(vaporflux)
            m_accum_daily(time) = smb_ice(ix, iy)*rho_ice*seconds_per_year
            m_rain_daily(time) = real(rainman*dt_firn)
            m_refreezing_daily(time) = real(dummy_refreeze)
            m_melt_daily(time) = real(dummy_melt)
            m_snow_daily(time) = real(accum*dt_firn)
            m_runoff_daily(time) = real(dummy_runoff)
            m_melt_ice_daily(time) = real(dummy_melt_ice)

            if (time .eq. 1) then
                m_real_mass_balance_daily(time) = real(sum(snowman(ix, iy, :))) - mass0 - real(dummy_melt_ice)
            else
                m_real_mass_balance_daily(time) = real(sum(snowman(ix, iy, :))) - sum(m_snowman_daily(time - 1, :)) - &
                									real(dummy_melt_ice)
            end if

            m_regridding_daily(time) = dummy_regrid
            m_snow_temp_ave_surface_daily(time) = sum(m_real_mass_balance_daily(1:time))
            
            ! Extra for wet snow else useless as snowman, general not very useful
            if (snowman(ix, iy, 1) .gt. 0.) then 
                m_snow_mask_daily(time) = 1
            end if

            if (include_values_over_ice) then
                m_runoff_daily(time) = m_runoff_daily(time) + real(dummy_melt_ice)
                m_rain_daily(time) = m_rain_daily(time) + real(dummy_rain_ice*dt_firn)
            end if
        end if

	end subroutine bessi_store_values

	subroutine bessi_store_smb(year, ix, iy, mass0)
		! Calculate and store surface mass balance

		implicit none

		integer, 		intent(in) 	:: year
		integer, 		intent(in) 	:: ix
		integer, 		intent(in) 	:: iy 
		real(kind=8), 	intent(in) 	:: mass0

		if (annual_data .and. (((year .ge. daily_data_start) .and. (year .le. daily_data_end)) &
			.or.(mod(year, annual_data_frequency) == 0))) then
                        
        	m_real_mass_balance_annual(ix,iy) = real(sum(snowman(ix, iy, :)) - mass0 - m_melt_ice_annual(ix, iy))
        end if
	end subroutine bessi_store_smb


	subroutine bessi_write_values_annual(year, albedo_runtime, smb_ice)
		! Write annual data to disk

		implicit none

		integer, 		intent(IN) :: year
		real(kind=8), 	intent(IN) :: smb_ice(nx, ny)
		real(kind=8),	intent(IN) :: albedo_runtime(nx, ny) 


		if (annual_data .and. (((year .ge. daily_data_start) .and. (year .le. daily_data_end)) &
	        .or. (mod(year,annual_data_frequency) == 0))) then

	        !conversion of annual data which are continued to be used on runtime
	        m_snowman_annual(:, :, :) = real(snowman(:, :, :))
	        m_lwmass_annual = real(lwmass)
	        m_snow_temp_annual = real(snow_temp)
	        m_rho_snow_annual = real(rho_snow)
	        
	        !averageing of annual surface data
	        m_albedo_dynamic_annual = real(albedo_runtime/ndays)
	        m_snow_temp_ave_surface_annual = m_snow_temp_ave_surface_annual/ndays
	        
	        m_accum_annual = m_accum_annual + smb_ice*rho_ice * seconds_per_year
	        

	        print *,'writing annual output netcdf file'        
	        call writeNCDFSNOW3DValues(filehandle_netcdf_annual, 1, snowman_annual_varid, m_snowman_annual(:,:,:), ny, nx,&
	                 n_snowlayer)
	        call writeNCDFSNOW3DValues(filehandle_netcdf_annual, 1, lwmass_annual_varid, m_lwmass_annual(:,:,:), ny, nx,&
	                 n_snowlayer)
	        call writeNCDFSNOW3DValues(filehandle_netcdf_annual, 1, snow_temp_annual_varid, m_snow_temp_annual(:,:,:), ny, nx,&
	                 n_snowlayer)        
	        call writeNCDFSNOW3DValues(filehandle_netcdf_annual, 1, rho_snow_annual_varid, m_rho_snow_annual(:,:,:), ny, nx,&
	                 n_snowlayer)
	                        
	        !Write 2D values per grid point           
	        call writeNCDFSNOW2DValues(filehandle_netcdf_annual, 1, &
	            accum_annual_varid, real(m_accum_annual(:,:)), ny, nx,1)
	        !sum of the latent heat flux  annually
	        call writeNCDFSNOW2DValues(filehandle_netcdf_annual, 1, snow_annual_varid, real(m_snow_annual(:,:)), ny, nx,1)
	        !sum of accumulation annually
	        call writeNCDFSNOW2DValues(filehandle_netcdf_annual, 1, rain_annual_varid, real(m_rain_annual(:,:)), ny, nx,1)
	        !sum of liquid precip  annually
	        call writeNCDFSNOW2DValues(filehandle_netcdf_annual, 1, latent_heat_annual_varid, &
	            real(m_vaporflux_annual(:,:)), ny, nx,1)
	        !sum of melt annually
	        call writeNCDFSNOW2DValues(filehandle_netcdf_annual, 1, albedo_dynamic_annual_varid, &
	            real(m_albedo_dynamic_annual(:,:)), ny, nx,1)
	        !runoff
	        call writeNCDFSNOW2DValues(filehandle_netcdf_annual, 1, runoff_annual_varid, real(m_runoff_annual(:,:)), ny, nx,1)
	        !sum of the refreeze  annually
	        call writeNCDFSNOW2DValues(filehandle_netcdf_annual, 1, melt_annual_varid, real(m_melt_annual(:,:)), ny, nx,1)
	        !sum of liquid solid annually
	        call writeNCDFSNOW2DValues(filehandle_netcdf_annual, 1, refreezing_annual_varid, &
	            real(m_refreezing_annual(:,:)), ny, nx,1)  
	        call writeNCDFSNOW2DValues(filehandle_netcdf_annual, 1, melt_ice_annual_varid, &
	            real(m_melt_ice_annual(:,:)), ny, nx,1)
	        call writeNCDFSNOW2DValues(filehandle_netcdf_annual, 1, snow_mask_annual_varid, real(snow_mask_annual(:,:)), ny, nx,1)
	        call writeNCDFSNOW2DValues(filehandle_netcdf_annual, 1, regridding_annual_varid, &
	            real(m_regridding_annual(:,:)), ny, nx,1)
	        call writeNCDFSNOW2DValues(filehandle_netcdf_annual, 1, &
	            snow_temp_ave_surface_annual_varid, real(m_snow_temp_ave_surface_annual(:,:)), ny, nx,1) 
	        call writeNCDFSNOW2DValues(filehandle_netcdf_annual, 1, &
	            real_mass_balance_annual_varid, &
	            real(m_real_mass_balance_annual(:,:)), ny, nx,1)
	        call closeNCDFFile(filehandle_netcdf_annual)
	    end if
	end subroutine bessi_write_values_annual

	subroutine bessi_write_values_monthly(year, ix, iy)
		! Write monthly data to disk

		implicit none

		integer, intent(IN) :: year
		integer, intent(IN) :: ix
		integer, intent(IN) :: iy 

		integer :: aa
		integer :: bb
		integer :: l

		! TODO: Is it even necessary to store the monthly values if the are immediately written to disk?
        if (monthly_data .and. (((year .ge. daily_data_start) .and. (year .le. daily_data_end)) &
            .or. (mod(year, monthly_data_frequency) == 0))) then
               
            	! Duplicate of loop below?
                l = 1
                aa = A(1)
                bb = B(1)
                m_snowman_monthly(l, 1:n_snowlayer) = sum(m_snowman_daily(aa:bb, 1:n_snowlayer), 1)
                ! TODO: Mistake?
                m_lwmass_monthly(l, 1:n_snowlayer) = sum(m_snowman_daily(aa:bb, 1:n_snowlayer), 1)
                m_snow_temp_monthly(l, 1:n_snowlayer) = sum(m_snow_temp_daily(aa:bb, 1:n_snowlayer), 1)
                m_rho_snow_monthly(l, 1:n_snowlayer) = sum(m_rho_snow_daily(aa:bb, 1:n_snowlayer), 1)
                m_albedo_dynamic_monthly(l) = sum(m_albedo_dynamic_daily(aa:bb))/(bb - aa)
                m_vaporflux_monthly(l) = sum(m_vaporflux_daily(aa:bb))*dt_firn/(L_v + L_lh)
                m_accum_monthly(l) = sum(m_accum_daily(aa:bb))
                m_refreezing_monthly(l) = sum(m_refreezing_daily(aa:bb))
                m_melt_monthly(l) = sum(m_melt_daily(aa:bb))
                m_snow_monthly(l) = sum(m_snow_daily(aa:bb))
                m_rain_monthly(l) = sum(m_rain_daily(aa:bb))
                m_runoff_monthly(l) = sum(m_runoff_daily(aa:bb))
                m_snow_mask_monthly(l) = sum(m_snow_mask_daily(aa:bb))
                m_melt_ice_monthly(l) = sum(m_melt_ice_daily(aa:bb))
                m_snow_temp_ave_surface_monthly(l) = sum(m_snow_temp_ave_surface_daily(aa:bb))/(bb - aa)
                m_regridding_monthly(l) = sum(m_regridding_daily(aa:bb))
                m_real_mass_balance_monthly(l) = sum(m_real_mass_balance_daily(aa:bb))
                
                do l = 2, 12, 1
                    aa = A(l)
                    bb = B(l)
                    m_snowman_monthly(l, 1:n_snowlayer) = sum(m_snowman_daily(aa:bb, 1:n_snowlayer), 1)
                    m_lwmass_monthly(l, 1:n_snowlayer) = sum(m_lwmass_daily(aa:bb, 1:n_snowlayer), 1)
                    m_snow_temp_monthly(l, 1:n_snowlayer) = sum(m_snow_temp_daily(aa:bb, 1:n_snowlayer), 1)
                    m_rho_snow_monthly(l, 1:n_snowlayer) = sum(m_rho_snow_daily(aa:bb, 1:n_snowlayer), 1)
                    m_albedo_dynamic_monthly(l) = sum(m_albedo_dynamic_daily(aa:bb))/(bb - aa)
                    m_vaporflux_monthly(l) = sum(m_vaporflux_daily(aa:bb))*dt_firn/(L_v + L_lh)
                    m_accum_monthly(l) = sum(m_accum_daily(aa:bb))
                    m_refreezing_monthly(l) = sum(m_refreezing_daily(aa:bb))
                    m_melt_monthly(l) = sum(m_melt_daily(aa:bb))
                    m_snow_monthly(l) = sum(m_snow_daily(aa:bb))
                    m_rain_monthly(l) = sum(m_rain_daily(aa:bb))
                    m_runoff_monthly(l) = sum(m_runoff_daily(aa:bb))
                    m_snow_mask_monthly(l) = sum(m_snow_mask_daily(aa:bb))
                    m_melt_ice_monthly(l) = sum(m_melt_ice_daily(aa:bb))
                    m_snow_temp_ave_surface_monthly(l) = sum(m_snow_temp_ave_surface_daily(aa:bb)/(bb - aa))
                    m_regridding_monthly(l) = sum(m_regridding_daily(aa:bb))
                    m_real_mass_balance_monthly(l) = sum(m_real_mass_balance_daily(aa:bb))
                end do

                call writeNCDFSNOW3D_pointwise_t_zxy(filehandle_netcdf_monthly, 1, snowman_monthly_varid, m_snowman_monthly, &
                										iy, ix, n_snowlayer, 12)
                call writeNCDFSNOW3D_pointwise_t_zxy(filehandle_netcdf_monthly, 1, lwmass_monthly_varid, m_lwmass_monthly, &
                										iy, ix, n_snowlayer, 12)
                call writeNCDFSNOW3D_pointwise_t_zxy(filehandle_netcdf_monthly, 1, rho_snow_monthly_varid, m_rho_snow_monthly, &
                										iy, ix, n_snowlayer, 12)
                call writeNCDFSNOW3D_pointwise_t_zxy(filehandle_netcdf_monthly, 1, snow_temp_monthly_varid, m_snow_temp_monthly, &
                										iy, ix, n_snowlayer, 12)
                call writeNCDFSNOW2D_pointwise_t_xy(filehandle_netcdf_monthly, 1, albedo_dynamic_monthly_varid, &
                									m_albedo_dynamic_monthly, iy, ix, 12)
                call writeNCDFSNOW2D_pointwise_t_xy(filehandle_netcdf_monthly, 1, latent_heat_monthly_varid, &
                									m_vaporflux_monthly, iy, ix, 12)
                call writeNCDFSNOW2D_pointwise_t_xy(filehandle_netcdf_monthly, 1, accum_monthly_varid, &
                									m_accum_monthly, iy, ix, 12)
                call writeNCDFSNOW2D_pointwise_t_xy(filehandle_netcdf_monthly, 1, rain_monthly_varid, &
                									m_rain_monthly, iy, ix, 12)
                call writeNCDFSNOW2D_pointwise_t_xy(filehandle_netcdf_monthly, 1, refreezing_monthly_varid, &
                									m_refreezing_monthly, iy, ix, 12)
                call writeNCDFSNOW2D_pointwise_t_xy(filehandle_netcdf_monthly, 1, melt_monthly_varid, &
                									m_melt_monthly, iy, ix, 12)
                call writeNCDFSNOW2D_pointwise_t_xy(filehandle_netcdf_monthly, 1, snow_monthly_varid, &
                									m_snow_monthly, iy, ix, 12)
                call writeNCDFSNOW2D_pointwise_t_xy(filehandle_netcdf_monthly, 1, runoff_monthly_varid, &
                									m_runoff_monthly, iy, ix, 12)
                ! TODO: Cast m_snow_mask_monthly to integer instead of real?
                call writeNCDFSNOW2D_pointwise_t_xy(filehandle_netcdf_monthly, 1, snow_mask_monthly_varid, &
                									real(m_snow_mask_monthly), iy, ix, 12)
                call writeNCDFSNOW2D_pointwise_t_xy(filehandle_netcdf_monthly, 1, melt_ice_monthly_varid, &
                									m_melt_ice_monthly, iy, ix, 12)
                call writeNCDFSNOW2D_pointwise_t_xy(filehandle_netcdf_monthly, 1, real_mass_balance_monthly_varid, &
                									m_real_mass_balance_monthly, iy, ix, 12)
                call writeNCDFSNOW2D_pointwise_t_xy(filehandle_netcdf_monthly, 1, snow_temp_ave_surface_monthly_varid, &
                									m_snow_temp_ave_surface_monthly, iy, ix,12)
                call writeNCDFSNOW2D_pointwise_t_xy(filehandle_netcdf_monthly, 1, regridding_monthly_varid, &
                									real(m_regridding_monthly), iy, ix,12)
                
                m_snowman_monthly = 0.
                m_lwmass_monthly = 0.
                m_rho_snow_monthly = 0.
                m_snow_temp_monthly = 0.
                m_albedo_dynamic_monthly = 0.
                m_snow_temp_ave_surface_monthly = 0.
                m_vaporflux_monthly = 0.
                m_accum_monthly = 0.
                m_rain_monthly = 0.
                m_refreezing_monthly = 0.
                m_melt_monthly = 0.
                m_snow_monthly = 0.
                m_runoff_monthly = 0.
                m_snow_mask_monthly = 0
                m_melt_ice_monthly = 0.
                m_real_mass_balance_monthly = 0.
                !m_regridding_monthly = 0 TODO: This too, correct?
        end if
	end subroutine bessi_write_values_monthly

	subroutine bessi_write_values_daily(year, ix, iy)
		! Write daily data to disk

		implicit none

		integer, intent(IN) :: year
		integer, intent(IN) :: ix
		integer, intent(IN) :: iy 

		! Write daily data
        if (daily_data .and. (((year .ge. daily_data_start) .and. (year .le. daily_data_end)) &
            .or. (mod(year, daily_data_frequency) == 0))) then

                call writeNCDFSNOW3D_pointwise_t_zxy(filehandle_netcdf_daily, 1, snowman_daily_varid, m_snowman_daily, &
                										iy, ix, n_snowlayer, ndays)
                call writeNCDFSNOW3D_pointwise_t_zxy(filehandle_netcdf_daily, 1, lwmass_daily_varid, m_lwmass_daily, &
                										iy, ix, n_snowlayer, ndays)
                call writeNCDFSNOW3D_pointwise_t_zxy(filehandle_netcdf_daily, 1, rho_snow_daily_varid, m_rho_snow_daily, &
                										iy, ix, n_snowlayer, ndays)
                call writeNCDFSNOW3D_pointwise_t_zxy(filehandle_netcdf_daily, 1, snow_temp_daily_varid, m_snow_temp_daily, &
                										iy, ix, n_snowlayer, ndays)
                call writeNCDFSNOW2D_pointwise_t_xy(filehandle_netcdf_daily, 1, albedo_dynamic_daily_varid, &
                									m_albedo_dynamic_daily, iy, ix, ndays)
                call writeNCDFSNOW2D_pointwise_t_xy(filehandle_netcdf_daily, 1, latent_heat_daily_varid, &
                									m_vaporflux_daily, iy, ix, ndays)
                call writeNCDFSNOW2D_pointwise_t_xy(filehandle_netcdf_daily, 1, accum_daily_varid, &
                									m_accum_daily, iy, ix, ndays)
                call writeNCDFSNOW2D_pointwise_t_xy(filehandle_netcdf_daily, 1, rain_daily_varid, &
                									m_rain_daily, iy, ix, ndays)
                call writeNCDFSNOW2D_pointwise_t_xy(filehandle_netcdf_daily, 1, refreezing_daily_varid, &
                									m_refreezing_daily, iy, ix, ndays)
                call writeNCDFSNOW2D_pointwise_t_xy(filehandle_netcdf_daily, 1, melt_daily_varid, &
                									m_melt_daily, iy, ix, ndays)
                call writeNCDFSNOW2D_pointwise_t_xy(filehandle_netcdf_daily, 1, snow_daily_varid, &
                									m_snow_daily, iy, ix, ndays)
                call writeNCDFSNOW2D_pointwise_t_xy(filehandle_netcdf_daily, 1, runoff_daily_varid, &
                									m_runoff_daily, iy, ix, ndays)
                ! TODO: Cast m_snow_mask_monthly to integer instead of real?
                call writeNCDFSNOW2D_pointwise_t_xy(filehandle_netcdf_daily, 1, snow_mask_daily_varid, &
                									real(m_snow_mask_daily), iy, ix, ndays)
                call writeNCDFSNOW2D_pointwise_t_xy(filehandle_netcdf_daily, 1, melt_ice_daily_varid, &
                									m_melt_ice_daily, iy, ix, ndays)
                call writeNCDFSNOW2D_pointwise_t_xy(filehandle_netcdf_daily, 1, real_mass_balance_daily_varid, &
                									m_real_mass_balance_daily, iy, ix, ndays)
                call writeNCDFSNOW2D_pointwise_t_xy(filehandle_netcdf_daily, 1, snow_temp_ave_surface_daily_varid, &
                									m_snow_temp_ave_surface_daily, iy, ix, ndays)
        end if    
        
        if ((monthly_data .or. daily_data) .and. (((year .ge. daily_data_start) .and. (year .le. daily_data_end)) &
                .or. (mod(year, monthly_data_frequency) == 0 .or. mod(year, daily_data_frequency) == 0))) then                 
                
                m_snowman_daily = 0.
                m_lwmass_daily = 0.
                m_rho_snow_daily = 0.
                m_snow_temp_daily = 0.
                m_albedo_dynamic_daily = 0.
                m_snow_temp_ave_surface_daily = 0. !delete -> TODO: Why delete?
                m_vaporflux_daily = 0.
                m_accum_daily = 0.
                m_rain_daily = 0.
                m_refreezing_daily = 0.
                m_melt_daily = 0.
                m_snow_daily = 0.
                m_runoff_daily = 0.
                m_snow_mask_daily = 0
                m_melt_ice_daily = 0.
                m_real_mass_balance_daily = 0.
                m_regridding_daily = 0
        end if
	end subroutine bessi_write_values_daily

	subroutine alloc_val_3D(var, value, d1, d2, d3)
		! Allocate a 3D array and initialize it with a
		! value

		implicit none

		real(kind=4), dimension(:, :, :), allocatable, intent(inout) :: var
		real, intent(in) :: value
		integer, intent(in) :: d1
		integer, intent(in) :: d2	
		integer, intent(in) :: d3	

		! TODO: Consider aligned allocation
		if (.not. allocated(var)) allocate(var(d1, d2, d3))
		var = value
	end subroutine

	subroutine alloc_val_2D_r(var, value, d1, d2)
		! Allocate a 2D array and initialize it with a
		! value

		implicit none

		real(kind=4), dimension(:, :), allocatable, intent(inout) :: var
		real, intent(in) :: value
		integer, intent(in) :: d1
		integer, intent(in) :: d2		

		if (.not. allocated(var)) allocate(var(d1, d2))
		var = value
	end subroutine

	subroutine alloc_val_2D_i(var, value, d1, d2)
		! Allocate a 2D array and initialize it with a
		! value

		implicit none

		integer, dimension(:, :), allocatable, intent(inout) :: var
		integer, intent(in) :: value
		integer, intent(in) :: d1
		integer, intent(in) :: d2		

		if (.not. allocated(var)) allocate(var(d1, d2))
		var = value
	end subroutine

	subroutine alloc_val_1D_r(var, value, d1)
		! Allocate a 1D array and initialize it with a
		! value

		implicit none

		real(kind=4), dimension(:), allocatable, intent(inout) :: var
		real, intent(in) :: value
		integer, intent(in) :: d1		

		if (.not. allocated(var)) allocate(var(d1))
		var = value
	end subroutine

	subroutine alloc_val_1D_i(var, value, d1)
		! Allocate a 1D array and initialize it with a
		! value

		implicit none

		integer, dimension(:), allocatable, intent(inout) :: var
		integer, intent(in) :: value
		integer, intent(in) :: d1	

		if (.not. allocated(var)) allocate(var(d1))
		var = value
	end subroutine

	subroutine write_check_mass(year)
		! Create a txt file with all mass fluxes

		implicit none

		integer, intent(in) :: year
		integer :: ii

		if(mass_checking) then
            if (year == 1) then
                open(42, file=trim(adjustl(output_directory)) // 'diagnosed_massflux_Snowman3D.txt', &
                         status="new", action="write")
                write(42,*) 'This file contains the total mass of water and snow as well as massfluxes of snow and water '
                write(42,*) 'inside the model Snoman3D. Total masses are in kg and the fluxes are in kg per timestep.'
                write(42,*) 'initial mass snow', tab, 'end mass snow', tab,'initial mass water', tab, &
                         	'end mass water', tab, 'accumulation snow', tab, 'accumulation water', tab, 'snow melt', &
                            tab, 'water runoff', tab, 'refreeze', tab, 'snow to ice model'
            else
                open(42, file=trim(adjustl(output_directory)) // 'diagnosed_massflux_Snowman3D.txt',&
                         status="old", position="append", action="write")
            end if

            check_init_mass_snow = check_init_mass_snow * dx * dy
            check_end_mass_snow = check_end_mass_snow * dx * dy

            check_init_mass_water = check_init_mass_water * dx * dy
            check_end_mass_water = check_end_mass_water * dx * dy

            check_accum_snow  = check_accum_snow * dx * dy
            check_accum_water = check_accum_water * dx * dy

            check_melted_snow = check_melted_snow * dx * dy
            check_runoff_water = check_runoff_water * dx * dy
            check_refreeze_water = check_refreeze_water * dx * dy
            check_ice2ice = check_ice2ice * dx * dy

            ! write fluxes to txt file
            do ii = 1, ndays, 1
                write(42,*) check_init_mass_snow(ii), tab, check_end_mass_snow(ii), tab, check_init_mass_water(ii), tab, &
                     check_end_mass_water(ii), tab, check_accum_snow(ii), tab, check_accum_water(ii), tab, &
                    check_melted_snow(ii), tab, check_runoff_water(ii), tab, check_refreeze_water(ii), tab, check_ice2ice(ii)
            end do

            close(42)
        end if
	end subroutine write_check_mass

	subroutine write_check_energy(year)
		! Creating a txt file with all energy fluxes
        
		implicit none

		integer, intent(in) :: year
		integer :: ii

		if(energy_checking) then
            if (year == 1) then
                open(420, file=trim(adjustl(output_directory)) // 'diagnosed_energyflux_Snowman3D.txt',&
                         status="new", action="write")
                write(420,*) 'This file contains the total energy stored in water and snow as well as energyfluxes into and out '
                write(420,*) 'of the model Snoman3D. Total energies are in J and the fluxes are in J per timestep.'
                write(420,*) 'initial energy snow', tab, 'end energy snow', tab, 'initial energy water', tab, 'end energy water', &
                				tab, 'obs. surface energy flux', tab, 'diag. surface energy flux', tab, 'accum tot energy flux', &
                				tab, 'energy runoff by percolation', tab, 'energy passed to ice', tab, 'freezing energy diag', &
                				tab, 'freezing energy heat', tab, 'melting energi diag', tab, 'melting energy heat', tab, &
                        		'Energy entering snow for melt', tab, 'energy runoff by melt-up'
            else
                open(420, file=trim(adjustl(output_directory)) // 'diagnosed_energyflux_Snowman3D.txt',&
                         status="old", position="append", action="write")
            end if

            check_init_energy_snow = check_init_energy_snow * dx * dy 
            check_init_energy_water = check_init_energy_water * dx * dy 
            check_end_energy_snow = check_end_energy_snow * dx * dy 
            check_end_energy_water = check_end_energy_water * dx * dy 
            check_surface_e_flux_obs = check_surface_e_flux_obs * dx * dy 
            check_surface_e_flux_diag = check_surface_e_flux_diag * dx * dy 
            check_e_accum_tot = check_e_accum_tot * dx * dy 

            check_e_freeze_tot  = check_e_freeze_tot * dx * dy 
            check_e_freeze_heat = check_e_freeze_heat * dx * dy 

            check_e_perc_runoff = check_e_perc_runoff * dx * dy 
            check_e_ice2ice = check_e_ice2ice * dx * dy

            check_e_melt_tot = check_e_melt_tot * dx * dy 
            check_e_melt_heat = check_e_melt_heat * dx * dy 
            check_e_melt_qq = check_e_melt_qq * dx * dy 
            check_e_melt_runoff = check_e_melt_runoff * dx * dy 

            ! write fluxes to txt file
            do ii=1,ndays,1
                write(420,*) check_init_energy_snow(ii), tab, check_end_energy_snow(ii), tab, check_init_energy_water(ii), tab, &
                     			check_end_energy_water(ii), tab,check_surface_e_flux_obs(ii), tab, check_surface_e_flux_diag(ii), &
                     			tab, check_e_accum_tot(ii), tab, check_e_perc_runoff(ii), tab, check_e_ice2ice(ii), tab, &
                     			check_e_freeze_tot(ii), tab, check_e_freeze_heat(ii), tab, check_e_melt_tot(ii), tab, &
                     			check_e_melt_heat(ii), tab, check_e_melt_qq(ii), tab, check_e_melt_runoff(ii)
            end do

            close(420)
        end if
	end subroutine write_check_energy
end module bessi_data
