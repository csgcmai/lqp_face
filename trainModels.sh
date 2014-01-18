#!/bin/bash

# Copyright (c) 2013, Sibt ul Hussain <sibt.ul.hussain at gmail dot com>
# All rights reserved.
# Released under BSD License
# -------------------------
# For license terms please see license.lic

#Example bash file for running experiments using python code.

srcdir=${PWD}
outputfile=record_trained_on_view1
outputdir=/home/hussain/research/code-release/new-data-xoffset=1-yoffset=-4/trained_on_view1
date >> ${outputfile}
#python trainAndEvaluateLFW.py  -V view2 -f LQP -t with-pca -p='[1250]'  -d /home/hussain/research/code-release/new-data-xoffset=1-yoffset=-4/tmp/Circ-Split-100x170-norm-1-ptype-2-tol-5-PatchSize-5-CodeType4-150--UOMap >> ${outputfile}

#python trainAndEvaluateLFW.py  -V view2 -f LQP -t with-pca -p='[1250]'  -d /home/hussain/research/code-release/new-data-xoffset=1-yoffset=-4/tmp/Circ-Split-100x170-norm-1-ptype-2-tol-7-PatchSize-5-CodeType4-150--UOMap >> ${outputfile}

#cd /home/hussain/research/code-release/

#bash computeFeatures.sh lqp /home/hussain/datasets/LFW/lfw-list-xoffset=1-yoffset=-4.txt ${outputdir}  /home/hussain/datasets/LFW/view1.txt >> ${outputfile}

#cd ${srcdir}

#python trainAndEvaluateLFW.py  -V view1 -f LQP -t with-pca -p='[100,200,400,500,700,1000,1200,1250,1300, 1400,1600,1800,2000]' -d ${outputdir}/Circ-Split-100x170-norm-1-ptype-2-tol-7-PatchSize-7-CodeType4-150--UOMap >> ${outputfile}

#python trainAndEvaluateLFW.py  -V view1 -f LQP -t with-pca -p='[100,200,400,500,700,1000,1200,1250,1300, 1400,1600,1800,2000]' -d ${outputdir}/Circ-Split-100x170-norm-1-ptype-2-tol-5-PatchSize-7-CodeType4-150--UOMap >> ${outputfile}

python trainAndEvaluateLFW.py  -V view2 -f LQP -t with-pca -p='[1600]' -d ${outputdir}/Circ-Split-100x170-norm-1-ptype-2-tol-5-PatchSize-7-CodeType4-150--UOMap >> ${outputfile}

python trainAndEvaluateLFW.py  -V view2 -f LQP -t with-pca -p='[1600]' -d ${outputdir}/Circ-Split-100x170-norm-1-ptype-2-tol-7-PatchSize-7-CodeType4-150--UOMap >> ${outputfile}

