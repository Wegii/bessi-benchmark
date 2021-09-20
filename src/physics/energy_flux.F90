module energy_flux
	! Calculate energy fluxes in a multilayer firn

	use bessi_defs

	implicit none

	private
	public :: fluxify
contains
	subroutine fluxify(ix, iy, time, air_temp_ice, dz, K_sw, H_lh, K_lh, Q_heat, dummy_heat, vaporflux, D_lf, p_air, &
						china_syndrome)
		! Calculate energy fluxes

		implicit none

		integer, intent(in) :: ix
		integer, intent(in) :: iy
		integer, intent(in) :: time
		real(kind=8), intent(in) :: air_temp_ice
		real(kind=8), intent(in) :: dz(n_snowlayer)
		real(kind=8), intent(in) :: K_sw
		real(kind=8), intent(in) :: H_lh
		real(kind=8), intent(in) :: K_lh
	    real(kind=8), intent(in) :: p_air
	    real(kind=8), intent(inout) :: D_lf
	    real(kind=8), intent(inout) :: Q_heat
		real(kind=8), intent(inout) :: dummy_heat
		real(kind=8), intent(inout) :: vaporflux
	    logical, intent(inout) :: china_syndrome


	    ! Use complex energy balance in case of latent heat flux on 
	    if (.not. latent_heat_flux_on) then 
	        D_lf = 0.
	    end if

	    call calculate_energy_flux(ix, iy, time, air_temp_ice, dz, K_sw, H_lh, K_lh, &
	    					 		Q_heat, dummy_heat, vaporflux, D_lf, p_air, china_syndrome)

	    if (sublimation_only) then
	        snowman(ix, iy, 1) = max(0., snowman(ix, iy, 1) + vaporflux/(L_v + L_lh)*dt_firn)
	    else
	        ! Sublimation and vaporization are not treated equally
	        if(snow_temp(ix, iy, 1) .lt. kelvin) then   
	            ! Add turbulent latent heat mass flux explicitly
	            snowman(ix, iy, 1) = max(0., snowman(ix, iy, 1) + vaporflux/(L_v + L_lh)*dt_firn)   
	        else 
	            lwmass(ix, iy, 1) = lwmass(ix, iy, 1) + vaporflux/(L_v)*dt_firn
	        end if
	    end if
	end subroutine fluxify

	subroutine calculate_energy_flux(ix, iy, time, T_air, ddz, K_sw, H_lh, K_lh, Q_heat, heating, vaporflux, D_lf, p_air, &
										china_syndrome)
	    ! Routine calculation the energy balance implicitly including latent heat flux
	    ! based on go_energy_flux_new by Michael Imhof
	    !   Author: Tobias Zolles
	    !   Developer: Tobias Zolles
	    !   last interation June 2019
	    
	    ! Tobias Zolles, adapted from Michael Imhof, January 2018

	    ! this subroutine calcualtes the energy fluxes in the snowcap ie the
	    ! new temperatures. an implicit sceme is used for that.
	    ! sublimation was added

	    implicit none
	    
	    integer, intent(in) :: ix
	    integer, intent(in) :: iy
	    integer, 		intent(in) 		:: time
	    real(kind=8), intent(in) :: T_air
	    real(kind=8), intent(in) :: K_sw
	    real(kind=8), intent(in) :: H_lh
	    real(kind=8), intent(in) :: K_lh
	    real(kind=8), intent(in) :: D_lf
	    real(kind=8), intent(in) :: p_air
	    real(kind=8), intent(in) :: ddz(n_snowlayer)
	    real(kind=8), intent(inout) :: Q_heat
	    real(kind=8), intent(inout) :: heating
	    real(kind=8), intent(inout) :: vaporflux
	    logical, intent(inout) :: china_syndrome

	    ! Local variables
	    integer :: nn ! amount of active grid cells in vertical direction
	    integer :: inn
	    integer :: ii
	    integer, dimension(n_snowlayer) ::  istheresnow

	    logical :: no_ground
	    real(kind=8) :: H 
	    real(kind=8) :: dummy
	    real(kind=8) :: BB_mid1_backup
	    real(kind=8) :: lwrd_l
	    
	    real(kind=8) :: ea
	    real(kind=8) :: es
	    real(kind=8) :: RH
	    real(kind=8) :: Q2M
	    
	    ! local allocatable variables
	    real(kind=8) :: tempi(n_snowlayer)
	    real(kind=8) :: new_temp(n_snowlayer)
	    real(kind=8) :: backup(n_snowlayer)
	    real(kind=8) :: K(n_snowlayer)
	    real(kind=8) :: K_snow(n_snowlayer)
	    real(kind=8) :: dz1(n_snowlayer)
	    real(kind=8) :: BB_up(n_snowlayer)
	    real(kind=8) :: BB_mid(n_snowlayer)
	    real(kind=8) :: BB_down(n_snowlayer)

	    nn = 0
		inn = 1
		ii = 0
		H = 0.
		no_ground = .true.

	    ! extract active gridcells and add 2 dummy gridcells at surface and bottom.

	    heating = 0.
	    istheresnow = 0

	    where (snow_temp(ix, iy, :) > 0.)
	        istheresnow(:) = 1
	    end where
	    nn = int(sum(istheresnow))

	    tempi = 0.
		new_temp = 0.
		backup = 0.
		K = 0.
		K_snow = 0.
		dz1 = 0.
		BB_up = 0.
		BB_mid = 0.
		BB_down = 0.

		! in case of melting must occure this will become true
	    china_syndrome = .false. 
	    Q_heat = 0.
	    BB_mid1_backup = 0.

	    ! Initialize grid with variables
	    ! initialize resized dz and temperature

	    tempi(1:nn)= snow_temp(ix, iy, 1:nn)
	    dz1(1:nn) = ddz(1:nn)

	    ! store a backup of the initial temperature
	    backup(1:nn) = tempi(1:nn)

		! theoretical info to the equations
		!!Latent heat flux: H=0.622rho_aL_v*C_h_u(ea-es)P???1,Oerlemans Rolstad   C_h ~2e-3
		!rho_a air densitiy, P air pressure, L_v 2.5e6
		!es= 611.20Pa exp(L_s(2834e6)/R(461.5)*(1/kelvin(273.16)-1/T_s))
		!es=6.112exp(22.46*T/(272.62+T)) T in C, es hPa Guide to Meterological Instruments and Methods of Obersvation (WMO, 2008)
		!RH=ea/ews*100
		!ea=0.6108*exp(17.27*DewpT(ix, iy, time)/(DewpT(ix, iy, time)+237.3))
		!sensible: rho_a cp Ch u =A Ch=A/rho_a/cp


		! TODO: COMPUTATIONALLY EXPENSIVE
	    if (humiditystyle == 'D2M') then
	    	ea = 610.8*exp(17.27*DewpT(ix, iy, time)/(DewpT(ix, iy, time) + 237.3))
	    elseif (humiditystyle == 'RH_water') then
	    	ea = 610.8*exp(17.27*(T_air - kelvin)/(T_air - kelvin + 237.3))*RH/100
	    elseif (humiditystyle =='RH_ice') then
	        if (T_air < kelvin) then
	            ea = 610.8*exp(22.46*(T_air - kelvin)/(T_air - kelvin + 272.62))*RH/100
	        else 
	            ea = 610.8*exp(17.27*(T_air - kelvin)/(T_air - kelvin + 237.3))*RH/100
	        end if
	    elseif (humiditystyle == 'Q2M') then
	        ea=Q2M*461.89/287.058*p_air
	    end if
	    
	    es = 611.2*exp(22.46*(snow_temp(ix, iy, 1) - kelvin)/(272.62 + (snow_temp(ix, iy, 1) - kelvin))) !hpa=6.112
	    
	    ! calculate incoming longwave radiation
	    if (longwave_from_air_temp) then
	    	lwrd_l = sigma*eps_air*(T_air)**4.
	    else 
	    	lwrd_l = lwrd(ix, iy, time)
	    end if
	    

	    ! TODO: COMPUTATIONALLY EXPENSIVE
	    ! set up K vector (Snow Temperature INdependent parts)
	    K(1) = dt_firn/c_i/snowman(ix, iy, 1)*((T_air)*D_sf + lwrd_l+sigma*(eps_snow*3.*(tempi(1))**4.) + &
	    		K_sw + K_lh + D_lf/p_air*(ea - es*1 + es*22.46*272.62*snow_temp(ix, iy, 1)/(272.62 + (snow_temp(ix, iy, 1) - &
	    		kelvin))**2))

	    ! K(nn) = dt_firn*Q_geo/c_i/snowman(ix, iy, 1)    ! This is aobut how the geothermal heatflux would look like. 
	                        ! A corresponding melting routine would be needed.

	    ! set up H in BB_mid (Snow Temperature DEpendent parts)
	    H = dt_firn/c_i/snowman(ix, iy, 1)*(D_sf +sigma*eps_snow*4.*(tempi(1))**3. +H_lh + D_lf/p_air*es*22.46*272.62/(272.62 + &
	    	(snow_temp(ix,iy, 1)-kelvin))**2)
	    
	    ! Set up tridiagonal system with two cases (1 and 2+ boxes)
	    if (nn == 1) then
	        new_temp(1) = (tempi(1) + K(1))/(1 + H)
	        
	        if(new_temp(1) .gt. kelvin)  then
				! If the surface temperature has become larger than 0C, do the
	        	! calculation again without the surface energyfluxes, but set the
	        	! surfacetemperature to 273 at the beginning and at the end. this is
	        	! done this way to avoid unphysical temperature fluxes into the depth of
	        	! the snowcover. the energy used to heat the snow to 0 degrees is subtractet later in the melting
	            china_syndrome = .true.

	            ! Set temperature in topbox to 0C and ignor incomming energy fluxes. keep topbox temperature at 0C and
	            ! calculate energy necessary to rise temperature to 0C
	            Q_heat = (kelvin-backup(1))*c_i*snowman(ix, iy, 1)
	            new_temp(1)=kelvin
	            heating = Q_heat

	        else
	            ! clean the temperature profile. where snow is at melting point over several layers due to water
	            ! it can happen that the temperature of one or more boxes can become 273.00000000000006 K due to
	            ! computational limits.
	        	! where(new_temp(:) .gt. kelvin)
	        	! 	new_temp(:) = kelvin
	        	! end where
	            heating = dt_firn * (new_temp(1)*(-D_sf -sigma*eps_snow*4.*(tempi(1))**3. -H_lh ) + &
	                    ((T_air)*D_sf+sigma*(eps_air*(T_air)**4.+eps_snow*3.*(tempi(1))**4.) + K_sw + K_lh) )

	        end if

	!        vaporflux = D_lf*(6.108*exp(17.27*DewpT(ix, iy, time)/(DewpT(ix, iy, time)+237.3))-&
	!         6.112*exp(22.46*(snow_temp(1)-273)/(272.62+(snow_temp(1)-273))))
	        ! update temperature

	        ! Update temperature
	        snow_temp(ix,iy, 1) = new_temp(1)
	        vaporflux = D_lf/p_air*(6.108*exp(17.27*DewpT(ix, iy, time)/(DewpT(ix, iy, time)+237.3)) - &
	        			6.112*exp(22.46*(snow_temp(ix,iy, 1)-kelvin)/(272.62+(snow_temp(ix,iy, 1)-kelvin))))
	    else
	        ! set up model for the temperature diffusion parameter
	        if(diff_model==1) then
	            K_snow(1:nn) = K_ice*(rho_snow(ix, iy, 1:nn)/1000.)**(1.88) ! yen1981
	        elseif(diff_model==2) then
	            do ii=1,nn,1
	                if (rho_snow(ix, iy, ii).gt. 156.) then
	                    K_snow(ii) = (0.138-1.01e-3*rho_snow(ix, iy, ii)+ 3.233e-6*rho_snow(ix, iy, ii)**2.) ! Sturm 1997
	                else
	                    K_snow(ii) = 0.023+ 0.234e-3*rho_snow(ix, iy, ii)
	                end if
	            end do  
	        else
	            K_snow(1:nn) = (2.1e-2+4.2e-4*rho_snow(ix, iy, 1:nn) + 2.2e-9*rho_snow(ix, iy, 1:nn)**3.) ! Van Dusen (1929)
	        end if

	        ! set up tridiagonal system
	        BB_up(1)   = -2.*dt_firn/rho_snow(ix, iy, 1)/c_i/dz1(1)*( K_snow(1)*dz1(1)+K_snow(1+1)*dz1(1+1)) /(dz1(1)+dz1(1+1))**2.
	        BB_up(nn)  = -999.

	        BB_down(1)  = -999.
	        BB_down(nn) = -2.*dt_firn/rho_snow(ix, iy, nn)/c_i/dz1(nn)*( K_snow(nn)*dz1(nn)+K_snow(nn-1)*dz1(nn-1)) /(dz1(nn)+dz1(nn-1))**2.

	        BB_mid(1)  = 1. - BB_up(1)
	        BB_mid(nn) = 1. - BB_down(nn)

	        if(nn .gt. 2)then
	            do ii=2,nn-1,1
	                BB_down(ii) = -2.*dt_firn/rho_snow(ix, iy, ii)/c_i/dz1(ii)*( K_snow(ii)*dz1(ii)+K_snow(ii-1)*dz1(ii-1)) /(dz1(ii)+dz1(ii-1))**2.
	                BB_up(ii)   = -2.*dt_firn/rho_snow(ix, iy, ii)/c_i/dz1(ii)*( K_snow(ii)*dz1(ii)+K_snow(ii+1)*dz1(ii+1)) /(dz1(ii)+dz1(ii+1))**2.
	                BB_mid(ii)  = 1. - BB_down(ii) - BB_up(ii)
	            end do
	        end if

	        BB_mid1_backup = BB_mid(1)
	        ! Snow Temperature dependent parts (only outgoing sensible heat and longwave rad and latent heatfux)
	        BB_mid(1) = BB_mid(1) + H 

	        ! Tridag solver
	        call tridag(BB_down,BB_mid,BB_up,(tempi+K),new_temp,nn)

	        ! If the surface temperature has become larger than 0C, do the
	        ! calculation again without the surface energyfluxes, but set the
	        ! surfacetemperature to 273 at the beginning and at the end. this is
	        ! done this way to avoid unphysical temperature fluxes into the deep of
	        ! the snowcover. the energy used to heat the snow to 0 degrees is subtractet later in the melting
	        if(new_temp(1) .gt. kelvin)  then

	            china_syndrome = .true.
	            ! Set temperature in topbox to 0C and ignor incomming energy fluxes. keep topbox temperature
	            ! at 0C and calculate energy necessary to rise temperature to 0C
	            Q_heat= (kelvin-backup(1))*c_i*snowman(ix, iy, 1)
	            backup(1)=kelvin
	            BB_mid(1) = BB_mid1_backup

	            ! Tridag solver
	            call tridag(BB_down(1:nn),BB_mid(1:nn),BB_up(1:nn),backup(1:nn),new_temp(1:nn),nn)

	            ! set temperature in topbox to 0C. The energy used for this will be subtracted in the melting routine
	            Q_heat = Q_heat + (kelvin-new_temp(1))*c_i*snowman(ix, iy, 1)
	            new_temp(1)=kelvin

	            ! due to computational limits it is possible that new_temp(x)= 273.00000000000006. this happens especially
	            ! when the entire collumn is at the freezing point. in order to avoid that we set temperatures larger than
	            ! 273.000 K to 273. K
	            where(new_temp(:) .gt. kelvin)
	                new_temp(:) = kelvin
	            end where
	            heating = Q_heat

	            !new_temp(1)=kelvin
	            !if(new_temp(2) .gt. kelvin) then
	            !   print*,'new temp 2', new_temp(2)
	            !end if
	            !QQ=(QQ+(kelvin-tempi(2))*c_i*snowman(ix, iy, 1))/dt_firn ! energy flux due to warming which cant be used for 
	            ! melting subtract energy that has been used for keeping the snow at 0C, note QQ is negative
	        else
	            ! Clean the temperature profile. where snow is at melting point over several layers due to water
	            ! it can happen that the temperature of one or more boxes can become 273.00000000000006 K due to
	            ! computational limits.
	            where(new_temp(:) .gt. kelvin)
	                new_temp(:) = kelvin
	            end where

	            heating = dt_firn * (new_temp(1)*(-D_sf -sigma*eps_snow*4.*(tempi(1))**3. -H_lh ) + &
	                    ((T_air)*D_sf+sigma*(eps_air*(T_air)**4.+eps_snow*3.*(tempi(1))**4.) + K_sw + K_lh) )
	        end if

	        snow_temp(ix,iy, 1:nn) = new_temp(1:nn)
	        vaporflux = D_lf/p_air*(610.8*exp(17.27*DewpT(ix, iy, time)/(DewpT(ix, iy, time)+237.3)) - &
	        			611.2*exp(22.46*(snow_temp(ix,iy, 1)-273)/(272.62+(snow_temp(ix,iy, 1)-273))))
	    end if
	end subroutine calculate_energy_flux

	subroutine tridag(a,b,c,r,u,n)
	    ! from Numerical Recipes in Fortran77 Second Edition 1992
	    ! input
	    implicit none
	    integer, intent(in) :: n

	    real(kind=8), intent(in) :: a(n)
	    real(kind=8), intent(in) :: b(n)
	    real(kind=8), intent(in) :: c(n)
	    real(kind=8), intent(in) :: r(n)

	    ! output
	    real(kind=8), intent(inout) :: u(n)

	    ! locals
	    integer, parameter :: NMAX = 666

	    integer :: j
	    real(kind=8) :: bet
	    real(kind=8) :: gam(NMAX)
	    !Solves for a vector u(1:n) of length n the tridiagonal linear set given by equation (2.4.1).
	    !a(1:n), b(1:n), c(1:n), and r(1:n) are input vectors and are not modified.

	    !Parameter: NMAX is the maximum expected value of n.

	    !One vector of workspace, gam is needed.
	    !if(b(1).eq.0.)pause ???tridag: rewrite equations???
	    !If this happens then you should rewrite your equations as a set of order N ??? 1, with u2
	    !trivially eliminated.

	    u=0.
	    bet=b(1)
	    u(1)=r(1)/bet
	    do j=2,n
	        !Decomposition and forward substitution.
	        gam(j)=c(j-1)/bet
	        bet=b(j)-a(j)*gam(j)

	        u(j)=(r(j)-a(j)*u(j-1))/bet
	    end do
	    do j=n-1,1,-1
	        !Backsubstitution.
	        u(j)=u(j)-gam(j+1)*u(j+1)
	    enddo
	    return
	end subroutine tridag
end module energy_flux
