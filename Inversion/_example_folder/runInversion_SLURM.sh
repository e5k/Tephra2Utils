#!/bin/bash

#SBATCH --cpus-per-task=1
#SBATCH -N 1
#SBATCH --output=slurm-%J.out
#SBATCH --job-name=Tephra2Inversion
#SBATCH --mail-user=email@domain.edu
#SBATCH --mail-type=END
#SBATCH --exclusive
module load foss/2016a Python/2.7.11

# Usage on Baobab
# sbatch --array=8-12 -p askja,shared -t 05:00:00 runMass.sh	BATCH
# sbatch -p askja,shared -t 05:00:00 runMass.sh					Single run
# -p: partition name
# -t: maximum run time
# --array: used in case of a batch, represent the minimum and the maximum mass


#____________________________________________________________
#
# DO NOT CHANGE ANYTHING BELOW THIS LINE
#____________________________________________________________

function sum() {
  perl -e "print $1 + $2"
}

source inversionConfig.conf

chmod 755 ../_scripts/genConfig.py
chmod 755 ../_scripts/genConfigForward.py

# Check the seed
if [ "$SEED" = "-1" ]; then
	RANGE=500
	SEED=$RANDOM
	let "SEED %= $RANGE"
fi

if [ "$BATCH" = "1" ]; then
	echo "BATCH"
	mass1=${SLURM_ARRAY_TASK_ID}
	mass2=$(($mass1 + deltaMass))

	for ht in `seq $minHt $deltaHt $maxHt`; do
		outDir=mass${mass1}_ht$(($ht/1000))
		mkdir $outDir
		cp $inputFile $outDir/inversionInput.txt
		cd $outDir

		touch h_q_rmse1.dat

		../../_scripts/genConfig.py $ht `sum $ht $incrHt` "1e$mass1" "1e$mass2" $ventE $ventN $ventA $minDiff $maxDiff $eddy $minMedPhi $maxMedPhi $minSigPhi  $maxSigPhi $minAlpha $maxAlpha $minBeta $maxBeta $minFTT $maxFTT $minWindSpeed $maxWindSpeed $minWindDir $maxWindDir $plumeModel $fixedWind $windLevels $colSteps $partSteps $lithicDensity $pumiceDensity $minPhi $maxPhi $fitTest $SEED
		date
		srun -n $SLURM_TASKS_PER_NODE ../../tephra2012_inversion tmp.conf ../$INPUT ../$WIND
		date
		  
		rm node_*
		  
		# Write out a configuration file for tephra2
		  ../../_scripts/genConfigForward.py
		# Run the tephra2 forward model
		srun -n 1 --exclusive ../../../tephra2-2012 tephra2.conf ../$gridFile wind_levels.out > tephra2.out

		cd ..
	done

elif [ "$BATCH" = "2" ]; then
	echo "SINGLE with variable seed"
	SEED=${SLURM_ARRAY_TASK_ID}

	outDir=mass${minMass}_ht$(($minHt/1000))_$SEED
	mkdir $outDir

	cp $inputFile $outDir/inversionInput.txt
	cd $outDir

	touch h_q_rmse1.dat
	../../_scripts/genConfig.py $minHt $maxHt "1e$minMass" "1e$maxMass" $ventE $ventN $ventA $minDiff $maxDiff $eddy $minMedPhi $maxMedPhi $minSigPhi  $maxSigPhi $minAlpha $maxAlpha $minBeta $maxBeta $minFTT $maxFTT $minWindSpeed $maxWindSpeed $minWindDir $maxWindDir $plumeModel $fixedWind $windLevels $colSteps $partSteps $lithicDensity $pumiceDensity $minPhi $maxPhi $fitTest $SEED
	date
		srun -n $SLURM_TASKS_PER_NODE ../../../tephra2012_inversion tmp.conf ../$inputFile ../$windFile
	date

	rm node_*

	# Write out a configuration file for tephra2
	../../_scripts/genConfigForward.py
	# Run the tephra2 forward model
	srun -n 1 --exclusive ../../../tephra2-2012 tephra2.conf ../$gridFile wind_levels.out > tephra2.out

	cd ..

else
	echo "SINGLE"
	outDir=mass${minMass}_ht$(($minHt/1000))
	mkdir $outDir
	cp $inputFile $outDir/inversionInput.txt
	cd $outDir

	touch h_q_rmse1.dat

	../../_scripts/genConfig.py $minHt $maxHt "1e$minMass" "1e$maxMass" $ventE $ventN $ventA $minDiff $maxDiff $eddy $minMedPhi $maxMedPhi $minSigPhi  $maxSigPhi $minAlpha $maxAlpha $minBeta $maxBeta $minFTT $maxFTT $minWindSpeed $maxWindSpeed $minWindDir $maxWindDir $plumeModel $fixedWind $windLevels $colSteps $partSteps $lithicDensity $pumiceDensity $minPhi $maxPhi $fitTest $SEED
	date
	srun -n 1 --exclusive ../../../tephra2012_inversion tmp.conf ../$inputFile ../$windFile
	date
	  
	rm node_*
	  
	# Write out a configuration file for tephra2
	../../_scripts/genConfigForward.py
	# Run the tephra2 forward model
	srun -n 1 --exclusive ../../../tephra2-2012 tephra2.conf ../$gridFile wind_levels.out > tephra2.out

	cd ..
fi