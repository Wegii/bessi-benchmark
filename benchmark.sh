l=1
k=1
i=1

albedo_new=0.85
albedo_wet=0.45
albedo_ice=0.3
D_sf=6
D_lf=1
e_air=0.8
latent_heat_switch=1
max_lwc=0.1
albedo_module=5
lapse=0.0065
fileshort=${l}_${k}

echo run"$i".a_new"$albedo_new".a_wet"$albedo_wet".a_i"$albedo_ice".D_sf"$D_sf".D_lf"$D_lf".eps"$e_air".latent"$latent_heat_switch".a_mod"$albedo_module".lwc"$max_lwc"
parameters1=\'""eps_"${e_air}"lapse_rate"${lapse}"\'
parameters="${parameters1//./}"
echo $parameters

#convert switch to logic
if [ $latent_heat_switch -eq 1 ]
then
	latent_heat_switch=.true.
else
	latent_heat_switch=.false.
fi

cd /work/pwegmann/BESSI_TEST/bessi_test
outputpath=/work2/pwegmann/output/ERAinterim/grl10_interp/

mkdir bin

cp src/bessi_defs.f90 bin/insert_variables_file.f90
cd bin
sed -i 's/a_new_foo/'${albedo_new}'/' insert_variables_file.f90
sed -i 's/a_wet_foo/'${albedo_wet}'/' insert_variables_file.f90
sed -i 's/a_ice_foo/'${albedo_ice}'/' insert_variables_file.f90
sed -i 's/D_sf_foo/'${D_sf}'/' insert_variables_file.f90
sed -i 's/D_lf_foo/'${D_lf}'/' insert_variables_file.f90
sed -i 's/eps_air_foo/'${e_air}'/' insert_variables_file.f90
sed -i 's/latent_heat_switch_foo/'${latent_heat_switch}'/' insert_variables_file.f90
sed -i 's/albedo_module_foo/'${albedo_module}'/' insert_variables_file.f90
sed -i 's/lwc_foo/'${max_lwc}'/' insert_variables_file.f90
sed -i 's|output_directory_foo|'${outputpath}'|' insert_variables_file.f90
# Manually added
sed -i 's/lapse_rate_foo/'${lapse}'/' insert_variables_file.f90
sed -i 's|erai_vector_foo|'${vector}'|' insert_variables_file.f90
sed -i 's|sim_length_foo|'${l}'|' insert_variables_file.f90

cd ..
cp src/IceBern2D_ERAi.f90 bin/IceBern2D.f90
sed -i 's|sim_length_foo|'${l}'|' bin/IceBern2D.f90

run_str=\'$(printf "%05d" $run)\'
# change path
cp src/ionc.f90 bin/io.f90
# change path
cp src/bessi.F90 bin/bessi.F90
cp src/bessi_data.F90 bin/bessi_data.F90
cp src/physics/regridding.F90 bin/physics/regridding.F90
cp src/physics/precipitation_accumulate.F90 bin/physics/precipitation_accumulate.F90
cp src/physics/densification.F90 bin/physics/densification.F90
cp src/physics/radiation.F90 bin/physics/radiation.F90
cp src/physics/energy_flux.F90 bin/physics/energy_flux.F90
cp src/physics/melting.F90 bin/physics/melting.F90
cp src/physics/percolation.F90 bin/physics/percolation.F90
cp src/physics/conservation.F90 bin/physics/conservation.F90

sed -i 's/run_foo/'${run_str}'/' bin/io.f90
sed -i 's/run_foo/'${run_str}'/' bin/bessi.F90
sed -i 's/parameters_foo/'${parameters}'/' bin/io.f90

cd bin
cp insert_variables_file.f90 variables.f90

export OMP_PLACES=threads


echo "**** Start compling scripts ****"
exename="SMB_simu_${run}_${k}"
rm *.o *.mod

#INPUT='/pwegmann/input'
#Year_start=
#Year_end=

THREADS=(24 48)
#SCHEDULE=(1 2 4 8 12 24 32 48 96)  

for t in "${THREADS[@]}"                                                   
do                                                                                                                                       
    #for s in "${SCHEDULE[@]}"                                              
    #do                                                                     
		gfortran -O3 -mcmodel=medium -c variables.f90 
		gfortran -O3 -mcmodel=medium -c physics/conservation.F90 -ffree-line-length-none #-Wall -Wextra
		gfortran -O3 -mcmodel=medium -c -I/work/zolles/libs/system/include io.f90
		gfortran -O3 -mcmodel=medium -c bessi_data.F90 -ffree-line-length-none #-Wall -Wextra
		gfortran -O3 -mcmodel=medium -c physics/regridding.F90 -ffree-line-length-none #-Wall -Wextra
		gfortran -O3 -mcmodel=medium -c physics/precipitation_accumulate.F90 -ffree-line-length-none #-Wall -Wextra
		gfortran -O3 -mcmodel=medium -c physics/densification.F90 -ffree-line-length-none #-Wall -Wextra
		gfortran -O3 -mcmodel=medium -c physics/radiation.F90 -ffree-line-length-none #-Wall -Wextra
		gfortran -O3 -mcmodel=medium -c physics/energy_flux.F90 -ffree-line-length-none #-Wall -Wextra
		gfortran -O3 -mcmodel=medium -c physics/melting.F90 -ffree-line-length-none #-Wall -Wextra
		gfortran -O3 -mcmodel=medium -c physics/percolation.F90 -ffree-line-length-none #-Wall -Wextra
		gfortran -O3 -mcmodel=medium -c bessi.F90 -fopenmp -DOMPRUN=true -ffree-line-length-none #-Wall -Wextra 
		echo "*** SMB EMB compiled"
		gfortran -O3 -mcmodel=medium -L/work/zolles/libs/system/lib -o ${exename}.x variables.o io.o bessi.o bessi_data.o \
		 		regridding.o precipitation_accumulate.o densification.o radiation.o energy_flux.o melting.o percolation.o \
		 		conservation.o IceBern2D.f90 \
				-lnetcdff -DOMPRUN=true -Wl,-rpath -Wl,/work/zolles/libs/system/lib -fopenmp
		echo "*** Ice Bern compiled"
		
		echo $exename".x"
		#./SMB_simu__1.x ${t} ${s} "/work2/pwegmann/output/bessi_performance/static_t${t}_s${s}"
        ./SMB_simu__1.x ${t} 0 "/work2/pwegmann/output/bessi_performance/new_static_t${t}_s0"        
    #done                                                                                                                              
done     
