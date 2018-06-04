#!/bin/bash

#PBS -N inversionDebug
#PBS -j oe
#PBS -V
#PBS -m n
#PBS -M sbiasse@ntu.edu.sg
#PBS -l nodes=1:ppn=12
#PBS -q q12
#PBS

# Load modules
module load python/2.7.4
module load openmpi/1.4.5-gnu

# Usage on OpenPBS
# qsub -t 9-11 T2inversion_openPBS.sh    Batch, mass between 10^9 and 10^11
# qsub runMass.sh                        Single run

# Sets whether the inversion is ran as a batch or as a single run
# 0: single run
# 1: batch
BATCH=1

# INPUT FILES
# Configuration file for the inversion
configFile=Anei_inversion.conf
# Field input for the inversion - km/m2
inputFile=Anei_input.txt
# Wind profile used for the inversion
windFile=wind.gen
# Calculation grid for the forward solution
gridFile=saku.utm
# Configuration file for the forward solution
configFileF=tephra2.conf

# INPUT RANGES
# Min plume height (m asl)
minHt=15000
# Max plume height (m asl)
maxHt=30000

# In case the inversion is run as a single run, define the mass boundaries below
# Note that these variables specify the exponent of a power 10 (i.e. 10^minMass kg)
# Min mass
minMass=10
# Max mass
maxMass=12

# In case the inversion is run as a batch
# Mass boundaries for each mass increment
deltaMass=1
# Mass increment
incrMass=1
# Height boundaries for each height increment
deltaHt=1000
# Height increment
incrHt=1000


#____________________________________________________________
#
# DO NOT CHANGE ANYTHING BELOW THIS LINE
#____________________________________________________________

function sum() {
perl -e "print $1 + $2"
}

cd $PBS_O_WORKDIR
echo $PBS_O_WORKDIR

chmod 755 genConfig.py

if [ "$BATCH" = "1" ]; then
echo "BATCH"
mass1=$PBS_ARRAYID
mass2=$(($mass1 + deltaMass))

	for ht in `seq $minHt $deltaHt $maxHt`; do
		outDir=res_${mass1}_$(($ht/1000))
		mkdir $outDir
		cp $inputFile $outDir/inversion_input.utm
		cp $configFileF $outDir/tephra2.conf
		cd $outDir

		touch h_q_rmse1.dat
		../genConfig.py ../$configFile tmp.conf $ht `sum $ht $incrHt` "1e$mass1" "1e$mass2"
		date
         	mpirun  ../../tephra2012_inversion tmp.conf ../$inputFile ../$windFile
		date

		rm node_*

		# Write out a configuration file for tephra2
		mpirun -np 1 perl ../../scripts/write_conf.pl parameters.README tephra2.conf 75 75
		# Run the tephra2 forward model
		mpirun -np 1 ../../tephra2-2012 tephra2.conf ../$gridFile wind_levels.out > tephra2.out

		cd ..
	done

else
	echo "SINGLE"
	outDir=res_${minMass}_$(($minHt/1000))
	mkdir $outDir
	cp $inputFile $outDir/inversion_input.utm
	cp $configFileF $outDir/tephra2.conf
	cd $outDir

	touch h_q_rmse1.dat

	../genConfig.py ../$configFile tmp.conf $minHt $maxHt "1e$minMass" "1e$maxMass"
	date
	mpirun -np $NCPUS -machinefile ../../tephra2012_inversion tmp.conf ../$inputFile ../$windFile
	##mpirun -np  $NCPUS -machinefile $PBS_NODEFILE  /home/sbiasse/Inversion/tephra2012_inversion /home/sbiasse/Inversion/Anei/res_10_15/tmp.conf /home/sbiasse/Inversion/Anei/Anei_input.txt
	date

	rm node_*

	# Write out a configuration file for tephra2
	mpirun -np 1 -machinefile $PBS_NODEFILE perl ../../scripts/write_conf.pl parameters.README tephra2.conf 75 75
	# Run the tephra2 forward model
	mpirun -np 1 -machinefile $PBS_NODEFILE ../../tephra2-2012 ../tephra2.conf ../$gridFile wind_levels.out >tephra2.out

	cd ..
fi

