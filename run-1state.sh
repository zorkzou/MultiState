#!/bin/bash

# Gaussian 16
module load chemsoft/g16c-avx
export GAUSS_SCRDIR=/tmp/zouwl/GaussianScr

mkdir -p $GAUSS_SCRDIR

export gaussian_ein=$2
export gaussian_eou=$3

# MS@GWEV
export CurrDir=`pwd`
export multistate_dir=$CurrDir/MultiState
export ms_scrdir=$CurrDir/JOB000
mkdir -p $ms_scrdir
export tem_mstate=$multistate_dir/templet-ch4
export inp_state1=$ms_scrdir/state1.gjf
export fch_state1=$ms_scrdir/state1.fch

# generate input files for state 1
$multistate_dir/multistate.exe -gen -nst 1 -gin $gaussian_ein -ctp $tem_mstate \
  -in1 $inp_state1

rm -f $fch_state1

# state 1 calculation
g16 -fchk=$fch_state1 $inp_state1

# write Gaussian's *.EOu file
$multistate_dir/multistate.exe -mix -nst 1 -gin $gaussian_ein -gou $gaussian_eou \
  -fc1 $fch_state1

