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
export ms_scrdir=$CurrDir/JOB001
mkdir -p $ms_scrdir
export tem_mstate=$multistate_dir/templet-o2
export inp_state1=$ms_scrdir/state1.gjf
export inp_state2=$ms_scrdir/state2.gjf
export fch_state1=$ms_scrdir/state1.fch
export fch_state2=$ms_scrdir/state2.fch

# generate input files for states 1 and 2
$multistate_dir/multistate.exe -gen -gin $gaussian_ein -ctp $tem_mstate \
  -in1 $inp_state1 -in2 $inp_state2

rm -f $fch_state1 $fch_state2

# state 1 calculation
g16 -fchk=$fch_state1 $inp_state1

# state 2 calculation
g16 -fchk=$fch_state2 $inp_state2

# write Gaussian's *.EOu file
$multistate_dir/multistate.exe -mix -chi 400 -gin $gaussian_ein -gou $gaussian_eou \
  -fc1 $fch_state1 -fc2 $fch_state2

