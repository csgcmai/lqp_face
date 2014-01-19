#!/bin/bash


# Copyright (c) 2013, Sibt ul Hussain <sibt.ul.hussain at gmail dot com>
# All rights reserved.
# Released under BSD License
# -------------------------
# For license terms please see license.lic


# Script for Computing  LBP, LTP & LQP Features...
# feature_type[ lbp or ltp or lqp or lbp+ltp] path_of_file_containing_list_of_images[./face-recog-100x170.txt] dir[directory path for storing feature files.] optional: path_of_a_separate_file_containing_list_of_images_for_codebook_learning> e.g bash computeFeatures.sh lqp ~/dataset/lfw-list.txt ~/experiments/data/ ~/dataset/lfw-list-view1.txt"

count=1
srcdir=${PWD}
exefile=${srcdir}/build/mainFeatures # mainFeatures-large
efname="runExp.sh"
nthreads=1 #number of threads
suffix=""
computeLQP()
{
ftype=14
	addcmd=" ${addcmd} --FeatureType=14 --CodeBook-DMetric=0 --ClusteringRounds=10 --Patch-PruneCount=10  --Patch-PatchStride=1"
#a
#	for dirname in ${ofname[@]}
#	do
		dirname=$ofname
		echo ${dirname}
		for patchsize in ${psize[@]}
		do

			for cbsizevar in ${cbsize[@]}
			do
				cd $mdir
				dname="${dirname}-PatchSize-${patchsize}-CodeType${codetype}-${cbsizevar}-${suffix}"
				echo "Running Dir Name"
				mkdir $dname

				cp  ${exefile} ${srcdir}/${efname} ${dname} # 
				echo "Training $dname " 
				cd $dname
				echo "${cbsizevar} run.sh ${cbsizevar}  ${patchsize} ${ptype[$count]} ${codetype}"
				
				 naddcmd="--Normalization=${norm} --CodeBookSize=${cbsizevar} --PatchType=${ptype}  --Patch-PatchSize=${patchsize}   --Patch-CodingType=${codetype} ${addcmd}"
				
				
				echo  "Running Process Number = $count"
#				tvar=$(($count%$nthreads))
#				echo " tvar="$tvar	

				if [ $(($count%$nthreads)) -eq 0 ]
				then
					echo  "Blocking Call "
 					sh ${efname} $exefile  "${naddcmd}"
 					wait

				else
					 sh ${efname} $exefile "${naddcmd}"&	
				fi
#				sh removefile.sh 	
		 	   	count=$(($count+1))
			done
		done

}

lqp()
{

#Best results on view1 are found using Tolerance=7 and Disk size of 7
#echo "Running Using LTP Features"
#computeLQP #
	for norm in 1
	do
		echo "Running Disc Split Only"
		foldname="Circ-Split-${fname}" # original(HOG
		codetype=4 # code type for circular > 16 = 2
		psize=(7)
#		cbsize=(100 150 200)
		cbsize=(150)
		tcount=0
		for ptype in 2 # horizontal5  diagonal 9 combined
		do
			for tol in 7 5
			do	
				echo "Running Using LQP Features"
				ofname="${foldname[${tcount}]}-norm-${norm}-ptype-${ptype}-tol-${tol}"
				addcmd=" --LTP-Tolerance=${tol} ${taddcmd} "
				computeLQP #
			done
			tcount=$(($tcount+1))
		done        
	      

		# Split-LTP 
	       	echo "Running Horizontal+Vertical+Diag+ADiag Only"
		foldname=Hor+Vert+Diag+ADiag # original(HOG
		codetype=4 # code type for circular > 16 = 2
		psize=7
		cbsize=(150)
		tcount=0
		for ptype in 9 # horizontal5  diagonal 9 combined
		do
			for tol in 5 7
			do	
				echo "Running Using LQP Features"
				ofname="${foldname[${tcount}]}-norm-${norm}-ptype-${ptype}-tol-${tol}"
				addcmd=" --LTP-Tolerance=${tol}"
#				computeLQP #
			done
			tcount=$(($tcount+1))
		done     
     done        
}
lbp()
{
	echo "Feature Destination Dir=${mdir} "
	cd $mdir
	echo  "Running Process Number = $count"
	echo "Running in Dir Name=${dname}"
	dname="${dname}${suffix}"
	mkdir $dname
	cp   ${exefile} ${srcdir}/${efname} ${dname} # 
	echo "Training in $dname " 
	cd $dname
#	tvar=$(($count%$nthreads))
#	echo " tvar="$tvar	
	if [ $(($count%$nthreads)) -eq 0 ]
	then
		echo  "Blocking Call "
		sh ${efname} $exefile "${addcmd}"
		wait

	else
		
		sh ${efname} $exefile "${addcmd}"&	
	fi
   	count=$(($count+1))

}


if [ $# -lt 3 ]
then
	echo "Error Wrong Number of Arguments: <$0  feature_type [lbp or ltp or lqp or lbp+ltp] path_of_file_containing_list_of_images [./face-recog-100x170.txt] dir [directory path for storing feature files.] optional: path_of_a_separate_file_containing_list_of_images_for_codebook_learning> e.g
	bash computeFeatures.sh lqp ~/dataset/lfw-list.txt ~/experiments/data/ ~/dataset/lfw-list-view1.txt"
	exit
#	return
fi

mdir=$3
if [ ! -d $mdir ]
then
	echo "Making Directory ${mdir}"
	mkdir -p ${mdir}
fi	
tfile=$2
if [ $# -eq 3 ]
then
	vfile=${tfile}
else
	vfile=$4
fi

echo "---------------------------------- Computing Features --------------------"
echo "Computing $1 Features .... "
echo "Number of threads = ${nthreads} "
if [ $1 = "lqp" -o $1 = "LQP" ]
then
	
	fname="100x170"
	taddcmd=" --Win-Width=80 --Win-Height=150 --Validation-File=${vfile} --TrainingFile=${tfile} "
	lqp
#call lbp & ltp computation.
elif [ $1 = "lbp" -o $1 = "LBP" ]
then 
	norm=1
	dname="lbp-norm-${norm}"
	addcmd=" --Win-Width=80 --Win-Height=150 --Validation-File=${vfile} --FeatureType=1 "
	lbp
elif [ $1 = "ltp" -o $1 = "LTP" ]
then 
	norm=1
	for tol in 5 7
	do
		dname="ltp-norm-${norm}-tol-${tol}"
		addcmd=" --Win-Width=80 --Win-Height=150 --Validation-File=${vfile} --FeatureType=2  --LTP-Tolerance=${tol} "
		lbp
	done 
elif [ $1 = "lbp+ltp" -o $1 = "LBP+LTP" ]
then 
	norm=1
	for tol in 5 7
	do
		dname="lbp+ltp-norm-${norm}-tol-${tol}"
		addcmd=" --Win-Width=80 --Win-Height=150 --Validation-File=${vfile} --FeatureType=3  --LTP-Tolerance=${tol} "
		lbp
	done 	
fi


