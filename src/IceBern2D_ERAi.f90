program IceModel

    !=========================
    ! Include Own Modules
    !=========================
    use bessi_defs       ! Module with all variables
    !use variables_snow
    use io              ! Own module to read the values
    !use smb_pdd         ! Module for positive degree day surface mass balance (Basil Neff)
    use smb_emb         ! Module for energy mass balance  (Michael Imhof)
    !use smb_troll       ! Module troll model (Michael Imhof)

    ! TODO: PERFORMANCE MEASUREMENT
    CHARACTER(100) :: num1char
    CHARACTER(100) :: num2char

    !use OMP_LIB
    INTEGER :: c1,c2,cr,cm,c_start,c_end, reorder_year
    REAL :: rate
        print *,"Spam"
    ! output directory
    if (store_input_netcdf .or. write_netcdf .or. annual_data .or. daily_data .or. monthly_data) then
        call init_output_directory(output_directory, TRIM(adjustl(experiment_description)))
    endif


    !=========================
    ! initialze Variable
    !=========================
    if(debug > 0 ) then
        print *, "Initialize Variables from external files"
        print *, "----------------------------------------"
    end if
if (erai_reorder_climate) then
    print*,"reading order vector"
    OPEN(UNIT=11, FILE=year_vector_file)
    do kk=1,maxyears_spam 
	 read(11,*) erai_vector(kk)
        if( iostat < 0 )then
        print*,'Warning: Year order file shorter than simulation period'
            print*,kk
            exit
        else if( iostat > 0 )then
            print*,'Error: error reading file'
            stop
        end if
    end do
!     do kk=1,234,1
!         erai_vector(kk)=1979
!     end do
!     do kk=10,50,1
!         erai_vector(kk)=2011
!     end do
!     erai_vector=(/ /)
    print*,erai_vector
end if

    ! for parallelisazion
    !result = system(TRIM(adjustl(set_core_command)))
    !result = system('echo $OMP_NUM_THREADS')

    CALL system_clock(count_rate=cr)
    CALL system_clock(count_max=cm)
    rate = REAL(cr)
    WRITE(*,*) "system_clock rate ",rate


    ! If the sea level gets not adjusted, set the offset to 0
    if (.not.adjust_sea_level) then
        sea_level_offset = 0.
    end if

    ! Read the relaxed Bedrock from a NetCDF file
    !--------------------------------------------
    if(debug > 1 ) then
        print *, "*** Read '", TRIM(adjustl(netcdf_input_bedrock_variable)) ,"' from file: ", &
                    TRIM(adjustl(netcdf_input_bedrock))
    end if
    ! Bedrock_initial is how the bedrock would be without iceload.
    Bedrock_initial = read_variable(netcdf_input_bedrock, nx, ny, TRIM(adjustl(netcdf_input_bedrock_variable)))

    ! set the initial value for ice free hemisphere
    elevation_netcdf = Bedrock_initial
    bedrock_netcdf = Bedrock_initial


    if(store_input_netcdf) then
        if(debug > 0) then
            print *, "*** Copy bedrock input (", TRIM(adjustl(netcdf_input_bedrock)) , &
                     ") file to the output directory: ", TRIM(adjustl(output_directory))
        end if
        call copy_to_output_directory(netcdf_input_bedrock, output_directory)
    end if


    ! 1 = water, 0 = normal case, 3 = unstable grid points
    ! assignment_mask = read_watermask(Bedrock_initial, nx, ny, (sea_level + sea_level_offset), Bedrock_netcdf, ice_thickness)




    ! load a present day climate icesheet
	if(initial_gis) then
		! set bedrock to state with ice		
		bedrock_netcdf = read_variable(netcdf_input_bedrock, nx, ny, netcdf_input_pd_bedrock_variable) ! TODO remove _special
		
		! set surface to state with ice
		elevation_netcdf = read_variable(netcdf_input_bedrock, nx, ny, netcdf_input_icesurf_variable) ! TODO remove _special
		print*,'*** PD ice sheets loaded'
	end if

    ! load a LGM climate icesheet
	if(initial_lgm_ice) then
		!load LGM topography
		bedrock_netcdf = read_variable(netcdf_input_lgm_ice, nx, ny, initial_bedrock_variable_name )
		elevation_netcdf = read_variable(netcdf_input_lgm_ice, nx, ny, initial_height_variable_name )
		print*,'*** LGM ice sheets loaded'
	end if


	! avoid negative ice thicknesses
	Ice_thickness = elevation_netcdf - bedrock_netcdf
	where(Ice_thickness(:,:).lt. 0.)
		elevation_netcdf(:,:) = bedrock_netcdf(:,:)
	end where	

	!elevation_netcdf = elevation
	!bedrock_netcdf = Bedrock

	Ice_thickness = elevation_netcdf - bedrock_netcdf


	! Restore Sea level
	!sea_level = get_sea_level(Ice_thickness, nx, ny, dx, dy, ocean_area)

	! Assign the assignment_mask
	if(debug > 0) then
		print *, "*** Get assignment_mask from the initial bedrock."
	end if
	assignment_mask = read_watermask(Bedrock_initial, nx, ny, (sea_level + sea_level_offset),&
	 bedrock_netcdf, Ice_thickness) ! 1 = water, 0 = normal case, 3 = unstable grid points
	if(debug > 0) then
		print *, "*** Created assignment_mask from the initial bedrock."
	end if

	! Set the height of the bedrock on the water to the sea level, Otherwise there will be an unstable integration,
	! cause the flux at the coast will get very high

	elevation = elevation_netcdf
	bedrock = bedrock_netcdf

	where(assignment_mask(:,:) .eq. 1)
		Bedrock(:,:) = sea_level + sea_level_offset
		elevation(:,:) = sea_level + sea_level_offset
		! integrated mass balance calculated before
		ice_thickness(:,:) = 0d0
		elevation_netcdf(:,:) = bedrock_netcdf(:,:)
	end where

	! Elevation of ice surface above sea level (over the water: Bedrock and ice_thickness are 0 -> elevation = 0)
	!elevation = Bedrock + ice_thickness

    	! Distance of each grid box from bottom domain boundary
    	y = y * dy
    	! Distance of each grid box from left domain boundary
    	X = x * dx



    ! Read the initial state out of an NetCDF File
    ! For this the watermask from the bedrock is used!
    ! THIS HAS NOT BEEN TESTED!!!!	
    if (read_initial_state) then
        ! Store the input files
        if(store_input_netcdf) then
            call copy_to_output_directory(initial_netcdf_file, output_directory)
        end if

        ! Bedrock_initial is how the bedrock would be without iceload.
        ! already done before, so not needed twice
        ! Bedrock_initial = read_variable(netcdf_input_bedrock, nx, ny, TRIM(adjustl(netcdf_input_bedrock_variable)))

        if(debug > 0 ) then
            print *, '*** Read bedrock state from file: ', TRIM(adjustl(initial_netcdf_file))
        end if
        Bedrock = read_variable(initial_netcdf_file, nx, ny, initial_bedrock_variable_name)

        elevation_netcdf = read_variable(initial_netcdf_file, nx, ny, initial_height_variable_name)
        elevation = elevation_netcdf

        ice_thickness = elevation - Bedrock

        ! Restore Sea level
        sea_level = get_sea_level(ice_thickness, nx, ny, dx, dy, ocean_area)

        ! Assign the assignment_mask
        assignment_mask = read_watermask(Bedrock_initial, nx, ny, (sea_level + sea_level_offset), bedrock_netcdf, ice_thickness) ! 1 = water, 0 = normal case, 3 = unstable grid points

        ! Set the height of the bedrock on the water to the sea level, Otherwise there will be an unstable integration,
        ! cause the flux at the coast will get very high
        do ix=1,nx,1
            do iy=1,ny,1
                if (assignment_mask(ix,iy) == 1) then
                    Bedrock(ix,iy) = sea_level + sea_level_offset
                    ice_thickness(ix,iy) = 0d0
                endif
            end do
        end do
    end if  ! ENDIF: Read initial state




    ! ELRA - Kugelfunktion
    if (active_elra) then
        open(5, iostat=ios, file=TRIM(adjustl(elra_kei_file)), status='old')
        if (ios /= 0) stop ' Error when opening the kei file!'
        do elra_ii=1,1059,1
            read(5,*) (kei(elra_ii,jj), jj=1,2)
        end do
        close(5, status='keep')
        print *, "*** Read kei file from ", TRIM(adjustl(elra_kei_file)), " successful."
    end if



	! READING THE CLIMATE INPUT DATA
	!-------------------------------


	! read climate iterim base data. 
	! potential_temperature and initial_climate_precipitation contain the base climate of eraiterim. 
	! temerature and precipitation can be adaptet to topography and time.

	if((eraiterim_climate).or.(erai_backandforth_climate).or.(erai_reorder_climate))then	
		! Load climate reference elevation (includes ice)
	    	initial_climate_elevation = read_variable(netcdf_input_eraiterim_initial_climate_elevation, nx, ny, &
					     TRIM(adjustl(netcdf_input_eraiterim_initial_climate_elevation_variable)))
					     print*,'here'
	
	end if

	if(eraiterim_climate)then
		inp_temp =   read_climate_long(netcdf_input_eraiterim_climate, nx, ny, TRIM(adjustl(netcdf_input_eraiterim_temp_variable)) )
		inp_precip = read_climate_long(netcdf_input_eraiterim_climate, nx, ny, TRIM(adjustl(netcdf_input_eraiterim_precip_variable)) )
        inp_dewpT = read_climate_long(netcdf_input_eraiterim_climate, nx, ny, TRIM(adjustl(netcdf_input_eraiterim_dewpT_variable)) )
        inp_wind = read_climate_long(netcdf_input_eraiterim_climate, nx, ny, TRIM(adjustl(netcdf_input_eraiterim_wind_variable)) )
        inp_lwrd = read_climate_long(netcdf_input_eraiterim_climate, nx, ny, TRIM(adjustl(netcdf_input_eraiterim_lwrd_variable)) )
    
		! Remove errorous negative precipitation
		where(inp_precip(:,:,:) .lt. 0.)
			inp_precip(:,:,:) = 0.
		end where

		if(ndays==365) then
			do id=1,ndays,1
				temperature(:,:,id)= inp_temp(:,:,id)+kelvin ! TODO TODO remove test cooling
				precipitation(:,:,id)= inp_precip(:,:,id)/3600./24./1000. ! from mm/day to m/sec
                DewpT(:,:,id)= inp_dewpT(:,:,id)
                wind(:,:,id)= inp_wind(:,:,id)
                if (longwave_downscaling) then
                        lwrd(:,:,id)= inp_lwrd(:,:,id)/temperature(:,:,id)**4*(temperature(:,:,id) +&
                        ((initial_climate_elevation-elevation) * temperature_lapse_rate))**4
                        else 
                        lwrd(:,:,id)= inp_lwrd(:,:,id)
                    end if 
                    if (lwrd_unit==2) then
                        lwrd(:,:,id)= lwrd(:,:,id)/3600./24.
                end if
				! Calculate potential temperature (at sea level)
				potential_temperature(:,:,id) = temperature(:,:,id) + (initial_climate_elevation * temperature_lapse_rate)

			end do

			!print*,'test temp',temperature(placex,placey,3)
			initial_climate_precipitation = precipitation

			! load short wave radiation
			inp_temp =   read_climate_long(netcdf_input_eraiterim_climate, nx, ny, TRIM(adjustl(netcdf_input_eraiterim_swradboa_variable)) )

			! reduce to 96 Timesteps
			do id=1,ndays,1
				P_sun0(:,:,id)= inp_temp(:,:,id)
				P_sun(:,:,id) = P_sun0(:,:,id)
			end do
		end if

		if(short_wave_damping) then
			call swrad_damping(P_sun,P_sun0, elevation, initial_climate_elevation)
		end if
		print*,'*** ERA-I climate loaded'

	end if
	



            ! Calculate Elevation Desertification (Budd and Smith (1979)), logic from Vizcaino et al. (2009)
		precipitation(:,:,:) = initial_climate_precipitation(:,:,:)



	! calculate elevation feedback on short wave radiation
	P_sun = P_sun0
	if(short_wave_damping) then
		call swrad_damping(P_sun, P_sun0, elevation, initial_climate_elevation )
	end if

	! add short wave radiation deviation
    	do id=1,ndays,1
		P_sun(:,:,id) = P_sun(:,:,id) + deviation_P_sun(:,:,id)

		where(P_sun(:,:,id).lt. 0.)
			P_sun(:,:,id) = 0.
		end where
	end do


    if(restart)then
        snowman=read_snow_data(restart_file, nx, ny,n_snowlayer, TRIM(adjustl('snowmass')))
        lwmass(1:nx,1:ny,1:n_snowlayer)=read_snow_data(restart_file, nx, ny, n_snowlayer, TRIM(adjustl('lwmass')) )
        snow_temp(1:nx,1:ny,1:n_snowlayer)=read_snow_data(restart_file, nx, ny, n_snowlayer, TRIM(adjustl('snowtemp')) )
        rho_snow(1:nx,1:ny,1:n_snowlayer)=read_snow_data(restart_file, nx, ny, n_snowlayer, TRIM(adjustl('snowdensity')) )
    end if

	! save the climate if wanted
	if((save_climate_forcing).or.(save_initial_climate)) then
		write (netcdf_output_climate_filename, '( "Climate_", I7.7, ".nc" )') 0
		call init_netcdf_climate_file(TRIM(adjustl(output_directory)) // netcdf_output_climate_filename,&
 				filehandle_netcdf3, ny, nx, int(dy), int(dx), precip_varid, temp_varid, swrad_varid,&
				 dev_precip_varid, dev_temp_varid, dev_swrad_varid )
		do id=1,ndays,1
			call writeNCDFGridValues(filehandle_netcdf3, id, precip_varid, real(precipitation(:,:,id)), ny, nx)
			call writeNCDFGridValues(filehandle_netcdf3, id, temp_varid, real(temperature(:,:,id)), ny, nx)
			call writeNCDFGridValues(filehandle_netcdf3, id, dev_precip_varid, real(deviation_precip(:,:,id)) , ny, nx)
			call writeNCDFGridValues(filehandle_netcdf3, id, dev_temp_varid, real(deviation_temp(:,:,id)), ny, nx)
			call writeNCDFGridValues(filehandle_netcdf3, id, swrad_varid, real(P_sun(:,:,id)), ny, nx)
			call writeNCDFGridValues(filehandle_netcdf3, id, dev_swrad_varid, real(deviation_P_sun(:,:,id)), ny, nx)
		end do
		call closeNCDFFile(filehandle_netcdf3)
	end if


    ! initialize netcdf output files
    !--------------------------------

    ! store the data
    if(store_input_netcdf) then
        call copy_to_output_directory(netcdf_input_eraiterim_climate, output_directory)
        !call copy_to_output_directory(netcdf_input_precip_calib, output_directory)
        !call copy_to_output_directory(netcdf_input_swradboa, output_directory)
        !call copy_to_output_directory(netcdf_input_albedo, output_directory)
    end if
    
    ! Integrated Mass Balance Output
    if (store_integrated_mass_balance) then
        open (unit=imb_filehandle,file=TRIM(adjustl(output_directory)) // imb_output_filename,action="write",status="replace")
        ! CSV Header
        ! year, accumulation, ablation, calving, isostaticmelt, totalmass, mass_change, net_change
        write (imb_filehandle,*) "year;accumulation;ablation;calving;isostaticmelt;totalmass;mass_change,net_change;sealevel"
    end if




    !=========================
    ! Lets go, do the loop
    !=========================
    if(debug > 0) then
        print *, "Loop over time steps"
        print *, "--------------------"
    endif
	call cpu_time(clock_start)
	CALL SYSTEM_CLOCK(c_start)
    do 
        ! Time series of diagnostics
        myyear(it) = it * int(real(dt)/(3600.*24.*365.)) ! the year that will be calculatet now
        ! Check if the loop conditions are at the end
        if ((myyear(it) > maxyears).or.((erai_backandforth_climate).and.(erai_year == erai_end_year+1)) ) then
            print *, 'The end is near (last year calculated): ', myyear(it)-1
            exit ! Jumps out of the loop. Does not exit the application (call exit(1)), otherwise the script is not terminated correctly
        end if



        !=====================================
        ! UPDATING CLIMATE AND MAP
        !=====================================
        ! this includes also elevation feedbacks and loading of new data

        ! Adjust Sea Level every 50 years
        if (store_integrated_mass_balance) then
            imb_isostaticmelt = 0d0
        end if


        ! load new climate data and update climate and swradboa
	!---------------------------------------------------------------

		if(erai_backandforth_climate) then
	        
			! creat string for the next year
			CALL SYSTEM_CLOCK(c1)
			
			write (spec_format, '("(A", I0,",I",I0,".", I0,",A", I0,")")') &
			len_trim(netcdf_input_name_leading_string), netcdf_input_digit_specification, &
			netcdf_input_digit_specification , len_trim(netcdf_input_name_end_string) 
			if(debug > 1) then
	            print*,'input_file_format',spec_format
	        end if
			
			
			write (new_input_file,spec_format) &
			netcdf_input_name_leading_string,erai_year,netcdf_input_name_end_string
			
	        print*,"*** Calculating mass balance with backandforth climate"
			!write (new_input_file , '( "ERAinterim_",I4.4,".interp.cdf" )') erai_year
	! 		write (new_input_file , '( "LGM_global",I4.4,".interp.cdf" )') erai_year

			new_input_file = TRIM(adjustl(netcdf_input_eraiterim_directory)) // TRIM(adjustl(new_input_file))

			! set number of next input file
			! 1979 1980 ... 2015 2016 2016 2015 ... 1980 1979 1979 1980 ...
			if (erai_climate_backwards)then
				erai_year = erai_year-1
			else
				erai_year = erai_year+1
			end if

			if((erai_year == erai_year_turn+1).and.(erai_current_iteration .le. erai_iterations_max ))then
				erai_climate_backwards = .true.
				erai_year = erai_year_turn
			end if

			if(erai_year == erai_year_begin-1)then
				erai_current_iteration = erai_current_iteration+1
				erai_climate_backwards = .false.
				erai_year = erai_year_begin
			end if

			if (( (      erai_climate_backwards).and.(erai_year .ne. erai_year_turn-1 ) ) .or. &
			    ( (.not. erai_climate_backwards).and.(erai_year .ne. erai_year_begin+1) ) .or. &
			    (myyear(it) == 1)	) then

				! load new year
	            call read_climate_once(new_input_file, nx, ny, TRIM(adjustl(netcdf_input_eraiterim_temp_variable)), &
	             TRIM(adjustl(netcdf_input_eraiterim_precip_variable)), TRIM(adjustl(netcdf_input_eraiterim_dewpT_variable)), &
	            TRIM(adjustl(netcdf_input_eraiterim_wind_variable)), TRIM(adjustl(netcdf_input_eraiterim_lwrd_variable)), &
	            TRIM(adjustl(netcdf_input_eraiterim_swradboa_variable)),&
	            inp_temp, inp_precip, inp_dewpT, inp_wind, inp_lwrd, inp_swrd, rate)	

			!	! Remove errorous negative precipitation
			!	where(inp_precip(:,:,:) .lt. 0.)
			!		inp_precip(:,:,:) = 0.
			!	end where
			
				if(ndays==365) then
					do id=1,ndays,1
						temperature(:,:,id)= inp_temp(:,:,id)+kelvin ! TODO TODO remove test cooling
						precipitation(:,:,id)= inp_precip(:,:,id)/3600./24./1000. ! from mm/day to m/sec
	                    DewpT(:,:,id)= inp_dewpT(:,:,id)
	                    wind(:,:,id)= inp_wind(:,:,id)
	                    if (longwave_downscaling) then
	                        lwrd(:,:,id)= inp_lwrd(:,:,id)/temperature(:,:,id)**4*(temperature(:,:,id) +&
	                        ((initial_climate_elevation-elevation) * temperature_lapse_rate))**4
	                        else 
	                        lwrd(:,:,id)= inp_lwrd(:,:,id)
	                    end if 
	                    if (lwrd_unit==2) then
	                        lwrd(:,:,id)= lwrd(:,:,id)/3600./24.
	                    end if
	                    P_sun0(:,:,id)= inp_swrd(:,:,id)
						P_sun(:,:,id) = P_sun0(:,:,id)
						! Calculate potential temperature (at sea level)
						potential_temperature(:,:,id) = temperature(:,:,id) + (initial_climate_elevation * temperature_lapse_rate)

					end do

					initial_climate_precipitation = precipitation

				
				end if

				if(short_wave_damping) then
					call swrad_damping(P_sun,P_sun0, elevation, initial_climate_elevation)
				end if
				print*,'*** ERA-I ',trim(adjustl(new_input_file)),' climate loaded'

			end if

	        CALL SYSTEM_CLOCK(c2)
	        WRITE(*,*) "system_clock : ",(c2 - c1)/rate
		end if
	
		if(erai_reorder_climate) then
	        
	        
	        reorder_year=erai_vector(erai_year)
	        print*,"*** Reoder climate forcing using year order vector:", reorder_year
			! creat string for the next year
			CALL SYSTEM_CLOCK(c1)
			
			write (spec_format, '("(A", I0,",I",I0,".", I0,",A", I0,")")') &
			len_trim(netcdf_input_name_leading_string), netcdf_input_digit_specification, &
			netcdf_input_digit_specification , len_trim(netcdf_input_name_end_string) 
			if(debug > 1) then
	            print*,'input_file_format',spec_format
	        end if
			
			write (new_input_file,spec_format) &
			netcdf_input_name_leading_string,reorder_year,netcdf_input_name_end_string
			!write (new_input_file , '( "ERAinterim_",I4.4,".interp.cdf" )') reorder_year
	! 		write (new_input_file , '( "LGM_global",I4.4,".interp.cdf" )') reorder_year

			new_input_file = TRIM(adjustl(netcdf_input_eraiterim_directory)) // TRIM(adjustl(new_input_file))

			! set number of next input file
			! 1979 1980 ... 2015 2016 2016 2015 ... 1980 1979 1979 1980 ...
			erai_year=erai_year+1


			if (( (      erai_climate_backwards).and.(erai_year .ne. erai_year_turn-1 ) ) .or. &
			    ( (.not. erai_climate_backwards).and.(erai_year .ne. erai_year_begin+1) ) .or. &
			    (myyear(it) == 1)	) then

				! load new year
	            call read_climate_once(new_input_file, nx, ny, TRIM(adjustl(netcdf_input_eraiterim_temp_variable)), &
	             TRIM(adjustl(netcdf_input_eraiterim_precip_variable)), TRIM(adjustl(netcdf_input_eraiterim_dewpT_variable)), &
	            TRIM(adjustl(netcdf_input_eraiterim_wind_variable)), TRIM(adjustl(netcdf_input_eraiterim_lwrd_variable)), &
	            TRIM(adjustl(netcdf_input_eraiterim_swradboa_variable)),&
	            inp_temp, inp_precip, inp_dewpT, inp_wind, inp_lwrd, inp_swrd, rate)
	    
	            ! Corrections
				if(ndays==365) then
					do id=1,ndays,1
						temperature(:,:,id)= inp_temp(:,:,id)+kelvin ! TODO TODO remove test cooling
						precipitation(:,:,id)= inp_precip(:,:,id)/3600./24./1000. ! from mm/day to m/sec
	                    DewpT(:,:,id)= inp_dewpT(:,:,id)
	                    wind(:,:,id)= inp_wind(:,:,id)
	                    if (longwave_downscaling) then
	                        lwrd(:,:,id)= inp_lwrd(:,:,id)/temperature(:,:,id)**4*(temperature(:,:,id) +&
	                        ((initial_climate_elevation-elevation) * temperature_lapse_rate))**4
	                        else 
	                        lwrd(:,:,id)= inp_lwrd(:,:,id)
	                    end if 
	                    if  (lwrd_unit==2) then
	                        lwrd(:,:,id)= lwrd(:,:,id)/3600./24.
	                    end if 
	                    P_sun0(:,:,id)= inp_swrd(:,:,id)
						P_sun(:,:,id) = P_sun0(:,:,id)
						! Calculate potential temperature (at sea level)
						potential_temperature(:,:,id) = temperature(:,:,id) + (initial_climate_elevation * temperature_lapse_rate)
					end do

					!print*,'test temp',temperature(placex,placey,3)
					initial_climate_precipitation = precipitation
				end if

				if(short_wave_damping) then
					call swrad_damping(P_sun,P_sun0, elevation, initial_climate_elevation)
				end if
				print*,'*** ERA-I ',trim(adjustl(new_input_file)),' climate loaded'

			end if

	        CALL SYSTEM_CLOCK(c2)
	        WRITE(*,*) "system_clock : ",(c2 - c1)/rate
		end if
    
        ! Calculate the elevation feedback every 20 years.
	if ( (mod(myyear(it), calc_rate_climate) == 0) .or. (erai_backandforth_climate) .or. &
			(erai_reorder_climate))then  		!TODO set back to 50 years and remove that heating
		
!		if(adjust_sea_level) then
!			sea_level = get_sea_level(ice_thickness, nx, ny, dx, dy, ocean_area)
!			! TODO: Do not take the initial bedrock, use the current one but forget the ice above it (and it should still increase, even if it is below the water
!			! But this could get to complicated, cause we do not want any sea in the middle of america.
!			assignment_mask = read_watermask(Bedrock_initial, nx, ny, (sea_level + sea_level_offset), &
!							! (use bedrock_netcdf, cause the bedrock is not equal to the sea level)
!							 Bedrock_netcdf, ice_thickness)
!		        
!			! integrated mass balance
!			! Die IF-Abfrage kann nicht innerhalb vom where integriert werden, daher einzeln
!			if (store_integrated_mass_balance) then
!				do ix=1,nx,1
!					do iy=1,ny,1
!						if ((assignment_mask(ix,iy) .eq. 1) .and. (ice_thickness(ix,iy) .gt. 0d0)) then
!							imb_isostaticmelt = imb_isostaticmelt - ice_thickness(ix,iy)
!						end if
!					end do
!				end do
!		        end if
!		        
!			! Set the sea level in the netcdf file for every water point
!			! 0 = ice, 1 = water, 2 = no ice, 3 = unstable integration
!			where(assignment_mask(:,:) .eq. 1)
!				Bedrock(:,:) = sea_level + sea_level_offset
!				elevation(:,:) = sea_level + sea_level_offset
!				! integrated mass balance calculated before
!				ice_thickness(:,:) = 0d0
!			end where
!			sunken_snow=0
!			do ix=1,nx,1
!				do iy=1,ny,1
!			    		! reset snow where there is sea
!			    		if((sum(snowman(ix,iy,:)).gt. 0.).and.(assignment_mask(ix,iy) .eq. 1)) then
!						sunken_snow = sunken_snow + sum(snowman(ix,iy,:)) + sum(lwmass(ix,iy,:))
!						do ii=1,n_snowlayer,1
!			    				snowman(ix,iy,ii)=0.
!							lwmass(ix,iy,ii)=0.
!							rho_snow(ix,iy,ii)=rho_s
!							snow_temp(ix,iy,ii)=0.
!						end do
!				    	end if
!				end do
!			end do
!		end if ! end adjust_sea_level


		! add the lapse rate for every day
		if(active_hysteresis) then
			do id = 1,ndays,1
			
				! It's still getting colder
				if ( myyear(it) .lt. hysteresis_return_point_in_time ) then
					! Check if the temperature is hold constant
					if((hysteresis_stop) .and. (hysteresis_stop_year .lt. myyear(it))) then
						! Do nothing, hysteresis_stop_temperature is reached
					else
					! 10 because it is calculated every 10 years
						potential_temperature(:,:,id) = potential_temperature(:,:,id) - (10 * hysteresis_temperature_factor)
					end if
					! Its getting warmer again
				else
					!if((hysteresis_stop_at_specific_temperature) == 2 .and. &
					!    (hysteresis_stop_temperature <= ((hysteresis_inital_temperature_offset - hysteresis_temperature_delta) &
					!    + ((myyear(it) - hysteresis_return_point_in_time) * hysteresis_temperature_factor)) ) ) then
					if((hysteresis_stop) .and. (hysteresis_stop_year .lt. myyear(it))) then
						! Do nothing, hysteresis_stop_temperature is reached
					else
						! 10 because it is calculated every 10 years
						potential_temperature(:,:,id) = potential_temperature(:,:,id) + (10. * hysteresis_temperature_factor)
					end if
				end if
			end do

		end if


		! calculate temperature lapsrate
        inp_temp(:,:,:)=temperature(:,:,:)
		
		
		if(active_temperature_lapsrate) then
		temperature(:,:,:) = potential_temperature(:,:,:)
			do id = 1,ndays,1
				temperature(:,:,id) = potential_temperature(:,:,id)+ deviation_temp(:,:,id) - &
							(elevation * temperature_lapse_rate)
                if(active_dewpoint_lapserate)then
                    dewpT(:,:,id) = dewpT(:,:,id)+ (initial_climate_elevation-elevation) * dewpoint_lapse_rate
                end if
			end do
		end if
        

		! calculate elevation feedback on short wave radiation
		P_sun = P_sun0
		if(short_wave_damping) then
			call swrad_damping(P_sun, P_sun0, elevation, initial_climate_elevation )
		end if
		! add short wave radiation deviation
            	!do id=1,ndays,1
!			P_sun(:,:,id) = P_sun(:,:,id) + deviation_P_sun(:,:,id)
!			where(P_sun(:,:,id).lt. 0.)
!				P_sun(:,:,id) = 0.
!			end where
!		end do
		


            ! Calculate Elevation Desertification (Budd and Smith (1979)), logic from Vizcaino et al. (2009)
	    precipitation(:,:,:) = initial_climate_precipitation(:,:,:)
            if(active_elevation_desertification) then
                do id=1,ndays,1
                    precipitation(:,:,id) = (initial_climate_precipitation(:,:,id) )
		    ! where the initial elevation is below the threshold
                    where (initial_climate_elevation(:,:) .lt. precipitation_threshold)
                        ! Only where the elevation is above the threshold, otherwise precipitation stays the same
                        where (elevation(:,:) .gt. precipitation_threshold)
                            precipitation(:,:,id) = (initial_climate_precipitation(:,:,id) )&
                                * exp( - precipitation_lapse_rate &
                                * (elevation(:,:) - precipitation_threshold) )
                        end where
                    ! where the initial climate elevation is above the threshold
                    ! (in this case, the precipitation can get amplified in lower areas)
                    elsewhere
                        ! Where the elevation is above the threshold, the precipitation gets lower
                        where (elevation(:,:) .gt. precipitation_threshold)
                            precipitation(:,:,id) = (initial_climate_precipitation(:,:,id)  )&
                                * exp( - precipitation_lapse_rate &
                                * (elevation(:,:) - initial_climate_elevation(:,:)) )
                        ! If the elevation falls below the threshold, amplify the precipitation
                        elsewhere ! turn of precipitation amplification
                            precipitation(:,:,id) = (initial_climate_precipitation(:,:,id) )&
                                * exp( - precipitation_lapse_rate &
                                * (precipitation_threshold - initial_climate_elevation(:,:)) )
                        end where
                    end where

			precipitation(:,:,id) = precipitation(:,:,id) + deviation_precip(:,:,id)

		    where(precipitation(:,:,id).lt. 0.)
			precipitation(:,:,id) = 0.
		    end where
                end do
 	        !print*,'minimum precip',minval(precipitation)
            end if

            ! pdd accumulation model
            ! the accumulation/ablation i.e. surfacebalance is calculatet for one year
            ! the pdd model smb is then used for 10 years
            ! pdd model
            if(smb_model==1)then
                accumulation = temperature_dependent_accumulation(nx, ny, temperature, precipitation, &
                    precipitation_unit, accumulation_daily_temperature_threshold, ndays)
                ablation = temperature_dependent_ablation(nx, ny, temperature, beta, ndays)
               

                surface_mass_balance = (accumulation - ablation)/seconds_per_year

		where( ( elevation_netcdf(:,:).lt. sea_level + sea_level_offset ) )
			surface_mass_balance(:,:) = min( surface_mass_balance(:,:),0.)
		end where


            end if

                ! fun model
            if(smb_model==3)then
                ! insert yearly data
                surface_mass_balance  = get_troll_accumulation(nx,ny,ndays,temperature, precipitation, assignment_mask, P_sun)
                !print *, "can you see me? "
                !print *, 'Year: ', myyear(it)

		where( ( elevation_netcdf(:,:).lt. sea_level + sea_level_offset ) )
			surface_mass_balance(:,:) = min( surface_mass_balance(:,:),0.)
		end where

                accumulation = surface_mass_balance
                surface_mass_balance = surface_mass_balance/seconds_per_year
            end if
		
        end if ! mod 20 year if for elevation feedbacks
          


	! write climate forcing data to netcdf
	!-------------------------------------
        if((mod(myyear(it), monthly_data_freq)==0).and.(save_climate_forcing)) then 
		write (netcdf_output_climate_filename, '( "Climate_", I7.7, ".nc" )') myyear(it)
		call init_netcdf_climate_file(TRIM(adjustl(output_directory)) // netcdf_output_climate_filename, &
 				filehandle_netcdf3, ny, nx, int(dy), int(dx), precip_varid, temp_varid, swrad_varid,&
				 dev_precip_varid, dev_temp_varid, dev_swrad_varid )
		do id=1,ndays,1
			call writeNCDFGridValues(filehandle_netcdf3, id, precip_varid, real(precipitation(:,:,id)), ny, nx)
			call writeNCDFGridValues(filehandle_netcdf3, id, temp_varid, real(temperature(:,:,id)), ny, nx)
			call writeNCDFGridValues(filehandle_netcdf3, id, dev_precip_varid, real(deviation_precip(:,:,id)) , ny, nx)
			call writeNCDFGridValues(filehandle_netcdf3, id, dev_temp_varid, real(deviation_temp(:,:,id)), ny, nx)
			call writeNCDFGridValues(filehandle_netcdf3, id, swrad_varid, real(P_sun(:,:,id)), ny, nx)
! 			call writeNCDFGridValues(filehandle_netcdf3, id, dev_swrad_varid, real(deviation_P_sun(:,:,id)), ny, nx)
		end do
		call closeNCDFFile(filehandle_netcdf3)
	end if

	
        !=====================================
        ! SNOW MODEL OR SMB STARTS HERE
        !=====================================

        ! the emb model accumulation must be calculatet for every year seperatly due to theyr descrete accumulation
        ! the models calculate the cumulative smb for one year and then the icesheet is adjusted
        ! emb model

        if(smb_model==2)then
            !print *,'myyear(it)=', myyear(it)
		call get_accumulation_snowman(ndays, temperature, precipitation, assignment_mask,&
		        	elevation_netcdf, sea_level + sea_level_offset, &
				surface_mass_balance, nc_counter, myyear(it), adaptive_timestep) ! topography and S_BOA
        ! assignment_mask from yelmo or via elevation
        ! snowman = "snowmass"
        ! lwmass = "water"
        ! surface_mass_balance output for yelmo (smb_ice in smb_emb/BESSI)
        ! nc_counter counter variable for netcdf files
        ! albedo_dynamic = "albedo from last year"
        ! elevation_netcdf "climate model elevation"
        ! elevation "effective elevation"

		sunken_snow=0.
		
		! surface_mass_balance in m_ice/second
		accumulation = surface_mass_balance*seconds_per_year ! in m_ice/ year
		    	!ablation = sum(snowman,3) ! snow in kg/m2 

        end if  

                print*,'*** end of snow part of ice model'

        !------------------------smb part ends here-------------------------------------

        ! Print Heartbeat
        if ((mod(myyear(it), netcdf_timesteps) == 0) .and. (debug > 0)  ) then !.and. (last_netcdf_year .lt. myyear(it))
            !call ETIME(execution_time, runtime)
	     call cpu_time(clock_end)
	     CALL SYSTEM_CLOCK(c_end)
	    runtime = runtime + clock_end - clock_start
            Write( heartbeat, '(i8)' ) myyear(it)
            print *, 'Year: ', TRIM(adjustl(heartbeat)), ', Runtime [s]: ', int(runtime), &
                !', End in [s]: ', int(( (runtime/myyear(it)) * maxyears) - runtime), &
                !', s/1000yr: ', (runtime/myyear(it) * 1000), &
                ', current yr/hour: ', (netcdf_timesteps/(clock_end - clock_start)*3600.), &		
                ', ave yr/hour: ', (myyear(it)/runtime * 60*60), &
                ', Sea level: ', (sea_level + sea_level_offset), &
                ', NetCDF TS: ', nc_counter-1 !((myyear(it)/netcdf_timesteps) + 1)
                print*,'system_time',(c_end-c_start)/rate
          CALL SYSTEM_CLOCK(c_start)
	     call cpu_time(clock_start)
        end if

	    call cpu_time(clock_end)
	    CALL SYSTEM_CLOCK(c_end)
	    runtime = runtime + clock_end - clock_start
            print *, "Runtime [mm:ss]: ", int(runtime/60), ':', mod(int(runtime),60)
            print *, 'Year: ', myyear(it)
 
    if(debug > 0 ) then
        print *, '========================'
        print *, "Successfully terminated!"
        print *, "We are at the end"
        print *, "-----------------"
        print *, "Some statistics:"
        !call ETIME(execution_time, runtime)
	call cpu_time(clock_end)
	CALL SYSTEM_CLOCK(c_end)
	runtime = runtime + clock_end - clock_start
        print *, "Runtime [mm:ss]: ", int(runtime/60), ':', mod(int(runtime),60)
        print *, "Seconds per 1000 years: ", (runtime/maxyears * 1000)
        print *, 'yr/hour: ', (maxyears/runtime * 60*60)
        print *, "-----------------"
        if (store_input_netcdf .or. write_netcdf .or. annual_data .or. daily_data .or. monthly_data) then
            print *, 'The output is stored in the directory: ', TRIM(adjustl(output_directory))
        else
            print *, 'No output was stored from this run!'
        endif
    end if
    print*,it
     it = int(it) + 1
end do

end program
