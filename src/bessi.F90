! Created on: May 12, 2015
!   Author: Imhof Michael
!   Developer: Imhof Michael
!   Mail: imhof@vaw.baug.ethz.ch or imhof@climate.unibe.ch
! Updated Version: September, 2020
!   Author: Tobias Zolles
!   Developer: Tobias Zolles
!   Mail: tobias.zolles@uib.no
! Update Version: September, 2021
!   Author: Peter Wegmann
!   Developer: Peter Wegmann
!   Mail: Peter.Wegmann@student.uib.no

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            
                 !!!!!!!!!  !!!!!!!!  !!!!!!!!  !!!!!!!!  !!
                 !!      !! !!        !!        !!        !!
                 !!!!!!!!!  !!!!!!!!  !!!!!!!!  !!!!!!!!  !! 
                 !!      !! !!              !!        !!  !!
                 !!!!!!!!!  !!!!!!!!  !!!!!!!!  !!!!!!!!  !!

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module smb_emb
    ! This module calculates the mass surface balance with the engergyfluxes in a multilayer firn

    use bessi_defs
    use bessi_data
    use io ! -> TODO: rename to bessi_io
    
    use regridding
    use precipitation_accumulate
    use densification
    use radiation
    use energy_flux
    use melting
    use percolation
    use conservation

    use omp_lib

    implicit none

    private
    public :: get_accumulation_snowman
contains
    subroutine get_accumulation_snowman(ndays, air_temp_ice, precip_ice, landmask, seafloor, sealevel, &
            smb_ice, nc_entry, year, spinup)
        ! Calculating the surface mass balance and energy balance of the firn for every grid point
        ! the whole year is calculated

        implicit none

        integer,        intent(in)      :: ndays
        integer,        intent(in)      :: year
        integer,        intent(in)      :: nc_entry
        logical,        intent(in)      :: spinup
        integer,        intent(in)      :: landmask(nx, ny) ! ice = 0, water  = 1, land with no ice = 3
        real(kind=8),   intent(in)      :: air_temp_ice(nx, ny, ndays) ! Temperature at the level of the ice elevation
        real(kind=8),   intent(in)      :: precip_ice(nx, ny, ndays)   ! Precipitation meta data
        real(kind=8),   intent(in)      :: seafloor(nx, ny)   ! TODO: unused
        real(kind=8),   intent(in)      :: sealevel          ! TODO: unused
        real(kind=8),   intent(inout)   :: smb_ice(nx, ny) ! height in ice equivalents for the ice wrapper

        real(kind=8) :: D_lf
        real(kind=8) :: Q_heat
        real(kind=8) :: dummy
        real(kind=8) :: H_lh
        real(kind=8) :: K_lh
        real(kind=8) :: K_sw
        real(kind=8) :: accum
        real(kind=8) :: rainman
        real(kind=8) :: lwc
        real(kind=8) :: vaporflux
        real(kind=8) :: dummy_melt_ice
        real(kind=8) :: dummy_rain_ice
        real(kind=8), dimension(n_snowlayer) :: dz
        real(kind=8), dimension(nx,ny) :: albedo_runtime       
        real(kind=8), dimension(nx,ny) :: p_air(nx,ny)  ! air pressure at the grid points elevation

        integer :: time
        integer :: nday_snowfall
        integer :: dummy_regrid

        real(kind=8) :: dummy_energy  !J/m2
        real(kind=8) :: dummy_heat    !J/m2
        real(kind=8) :: dummy_e_qq    !J/m2
        real(kind=8) :: dummy_melt
        real(kind=8) :: dummy_runoff
        real(kind=8) :: dummy_refreeze
        real(kind=8) :: dummy_ice2ice

        real(kind=8) :: mass0
        real(kind=8) :: dmass
        
        ! Loop variables 
        integer :: ii, ix, iy
        ! Timing variables
        real(kind=8) :: tm, tn

        ! initialize variables
        albedo_runtime = 0.
        mass0 = 0.
        rainman = 0.
        accum = 0.
        K_lh = 0.
        H_lh = 0.
        K_sw = 0.
        dz = 0.
        vaporflux = 0.
        nday_snowfall = 0
        lwc = 0
        dmass = 0
        time = 0

        dummy_energy = 0.          
        dummy_heat = 0.            
        dummy_e_qq = 0.     
        dummy_melt = 0.
        dummy_refreeze = 0.
        dummy_regrid = 0    
        dummy_melt_ice = 0
        dummy_rain_ice = 0
        dummy_runoff = 0.

        ! Calculate pressure coordinates of elevation
        p_air(:, :) = 101325*exp(-9.80665*0.0289644*elevation(:, :)/288.15/8.31447)
        
        ! Possibility to read albedo as input (for example annual averaged)
        if (albedo_module == 0) then
            albedo_dynamic(1:nx, 1:ny) = read_variable(albedo_file_path, nx, ny, albedo_file_variable_name)
        end if
        
        ! Calculate exchange coefficient for turbulent latent heat flux based on ratio and sensible heat flux    
        if (latent_heat_flux_on) then
            if (latent_heat_flux_analog_to_sensible_heat_flux) then
                 D_lf = ratio*D_sf/cp_air*0.622*(L_v + L_lh)
            else
                 D_lf = D_sf/cp_air*0.622*(L_v + L_lh)
            end if
        else
            D_lf = 0
        end if

        ! Intitialize conservation variables
        if (check_conservation) call check_init()
        ! Initialize output files
        call bessi_write_init(year)
        
        if ((year .lt. start_speed_up ) .or. (year .ge. end_speed_up)) fast_calculation = .false.

#ifdef OMPRUN
            tn = omp_get_wtime()
#else
            call cpu_time(tn)
#endif

        !$OMP PARALLEL DO COLLAPSE(2) SCHEDULE(static, 48) NUM_THREADS(48) &
        !$OMP     firstprivate(mass0, time, H_lh, K_lh, accum, rainman, dummy_melt, dummy_runoff, dummy_refreeze, &
        !$OMP         dummy_melt_ice, dummy_rain_ice, dummy_regrid, nday_snowfall, dz, lwc, K_sw, dummy_energy, &
        !$OMP         dummy_ice2ice, dummy_heat, vaporflux, D_lf, Q_heat, china_syndrome, dmass)
        do iy = 1, ny, 1        
            do ix = 1, nx, 1    

                if (landmask(ix,iy) .ne. 1) then
                    if ((fast_calculation(ix, iy) .eqv. .false.) .or. (mod(myyear(it), calc_rate_snow) == 0) .or. &
                        (spinup .eqv. .false.)) then 
                        
                        smb_ice(ix, iy) = 0.
                        mass0 = sum(snowman(ix,iy,:))

                        ! Reset fast calculation everywhere.
                        fast_calculation(ix,iy) = .true.
                        
                        do time = 1, ndays, 1
                            ! Check for energy conservation
                            if (check_conservation) call check_start_loop(ix, iy, time)

                            ! Reset dummy variables
                            dummy_melt = 0.    
                            dummy_runoff = 0.  
                            dummy_refreeze = 0.
                            dummy_melt_ice = 0.
                            dummy_rain_ice = 0.
                            dummy_regrid = 0

                            ! Accumulation
                            call accumulate(ix, iy, time, rainman, accum, nday_snowfall, H_lh, K_lh, air_temp_ice, precip_ice)

                            ! Further calculations only where there is snow
                            if(snowman(ix, iy, 1) > 0.) then 
                                if (check_conservation) call check_accumulation(ix, iy, time, accum, rainman)

                                ! Splitting and fusing of boxes
                                call regrid(ix, iy, dummy_regrid)
                                
                                ! Densification
                                call densify(ix, iy, accum + rainman, dz)
                                
                                ! Calculate liquid water content
                                lwc = lwmass(ix, iy, 1)/snowman(ix, iy, 1)/rho_w/(1./rho_snow(ix, iy, 1) - 1./rho_i)
                                
                                ! Net solar radiation 
                                call radiate_point(ix, iy, nday_snowfall, time, K_sw, lwc)

                                ! Energy Fluxes
                                dummy_energy = get_check_energy(ix, iy)
                                
                                call fluxify(ix, iy, time, air_temp_ice(ix, iy, time), dz, K_sw, H_lh, K_lh, &
                                        Q_heat, dummy_heat, vaporflux, D_lf, p_air(ix, iy), china_syndrome)

                                if (check_conservation) call check_flux(ix, iy, time, dummy_energy, dummy_heat)

                                ! Melt snow
                                if(china_syndrome) then
                                    fast_calculation(ix,iy) = .false.
                                    dummy_energy = get_check_energy(ix, iy)
                                    
                                    call melt_snow(ix,iy,time, air_temp_ice(ix,iy,time), K_sw, K_lh, H_lh, Q_heat, dummy_melt, &
                                                    dummy_runoff, smb_ice(ix,iy), vaporflux, D_lf, dummy_melt_ice, dummy_regrid, &
                                                    dummy_e_qq, dummy_heat)

                                    if (check_conservation) call check_melt(ix, iy, time, dummy_melt, dummy_energy, dummy_heat, &
                                                                                dummy_runoff, dummy_e_qq)
                                end if
                            
                                ! Water Percolation 
                                if(maxval(lwmass(ix, iy, :)) > 0.) then
                                    fast_calculation(ix,iy) = .false.
                                    
                                    call percolate(ix, iy, dummy_runoff)

                                    if (check_conservation) call check_runoff(time, dummy_runoff)
                                    dummy_energy = get_check_energy(ix, iy)
                                    
                                    ! Refreezing
                                    call refreeze(ix, iy, dummy_refreeze, dummy_heat)

                                    if (check_conservation) call check_refreeze(ix, iy, time, dummy_heat, dummy_refreeze, &
                                                                                dummy_energy)
                                end if

                            else
                                ! Melt ice of snow free, but ice containing, grid cells
                                call melt_ice(ix, iy, time, smb_ice(ix, iy), air_temp_ice(ix, iy, time), precip_ice(ix, iy, time), &
                                                vaporflux, D_lf, dummy_melt_ice, dummy_rain_ice, p_air(ix, iy))
                            end if

                            ! Store annual and daily data in variables
                            call bessi_store_values(year, ix, iy, time, accum, vaporflux, mass0, dummy_runoff, dummy_refreeze, &
                                        dummy_melt, dummy_melt_ice, dummy_regrid, dummy_rain_ice, smb_ice, rainman, albedo_runtime)
                            
                            ! TODO: This check can be omitted?!
                            if(time .lt. ndays) then
                                if (check_conservation) call check_end_loop(ix, iy, time)
                            end if
                        end do

                        ! Store and calculate surface mass balance. Calculation of annual mass balance has to be done prior to the
                        ! shifing of the mass balance, and can not be done within the timesteps.
                        call bessi_store_smb(year, ix, iy, mass0)
                        
                        ! Calculate upper_massbound
                        dmass = sum( snowman(ix,iy,:)) - n_snowlayer * soll_mass *1.5
                        ! Calculate mass balance for ice model after each year
                        if ((dmass .gt. 0.)) then
                            ! Maximum mass in one snowcolumn is limited to upper_massbound*n_snowlayer. All surplus snow is moved
                            ! to the ice model
                            call smb_point(ix, iy, dummy_runoff, dmass, dummy_ice2ice, dummy_energy, smb_ice)

                            if (check_conservation) call check_smb(dummy_ice2ice, dummy_energy, dummy_runoff)
                        else
                            ! Slow calculation if snow grid is not yet filled with snow
                            fast_calculation(ix,iy) = .false.
                        end if

                        if (check_conservation) call check_end(ix, iy)

                        ! Fast coupling on land, unless there is ice nearby. 
                        if ((smb_ice(ix,iy).lt. 0.).and.(landmask(ix,iy).eq. 2) )then  
                            fast_calculation(ix,iy) = .true.

                            if( (landmask(max(ix-1,1),iy).eq. 0) .or. (landmask(min(ix+1,nx),iy).eq. 0) .or. &
                                (landmask(ix,max(iy-1,1)).eq. 0) .or. (landmask(ix,min(iy+1,ny)).eq. 0) ) then
                                fast_calculation(ix,iy) = .false.
                            end if
                        end if

                        ! Seasonal snow on ice slow coupling, 
                        if ((smb_ice(ix,iy).lt. 0.).and.(landmask(ix,iy).eq. 0) )then  
                            fast_calculation(ix,iy) = .false.
                        end if

                        ! Collumns with growing mass but no ice yet, slow
                        if ((smb_ice(ix,iy).ge. 0.).and.(landmask(ix,iy).eq. 2) )then  
                            fast_calculation(ix,iy) = .false.
                        end if
                    end if
                end if

                ! Write monthly and daily values to disk
                call bessi_write_values_monthly(year, ix, iy)
                call bessi_write_values_daily(year, ix, iy)
            end do 
        end do
        !$OMP END PARALLEL DO
        
#ifdef OMPRUN
            tm = omp_get_wtime()
#else
            call cpu_time(tm)
#endif
        print *, 'end of xy loop ', tm - tn
        print*,'*** end of year', year

        if (minval(snowman)<0.) print*,'EOY: lightes snow grid cell', minval(snowman),'heaviest snow grid cell', maxval(snowman)
        if (minval(lwmass)<0.) print*,'EOY: negative lw content', minval(lwmass),'wettest snow grid cell', maxval(lwmass)
        
        ! Write and deallocate annual data
        call bessi_write_values_annual(year, albedo_runtime, smb_ice)

        if (check_conservation) then
            ! Creating a txt file with all mass and energy fluxes
            call write_check_mass(year)
            call write_check_energy(year)
        end if
    end subroutine get_accumulation_snowman
end module smb_emb
