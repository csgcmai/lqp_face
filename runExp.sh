#!/bin/bash


# Copyright (c) 2013, Sibt ul Hussain <sibt.ul.hussain at gmail dot com>
# All rights reserved.
# Released under BSD License
# -------------------------
# For license terms please see license.lic


exefile=$1
addcmd=$2

pstride=1
outfile=record


echo "Running on `hostname` ">> $outfile
echo "In the Dir $PWD" >> $outfile
echo "On Following Date `date`" >> $outfile
echo " Added CMD = $addcmd "
echo " Added CMD = $addcmd " >> $outfile

cmd="${exefile} --Color-Channel=4 --Add-RGBLBP=1  --Norm-Sep=0 --Cell-Size=10 --LBP-Stride=10  ${addcmd}"

echo $cmd
echo $cmd >> $outfile
$cmd >> $outfile

