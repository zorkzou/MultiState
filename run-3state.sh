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
export ms_scrdir=$CurrDir/JOB002
mkdir -p $ms_scrdir
export tem_mstate=$multistate_dir/templet-fes
export inp_state1=$ms_scrdir/state1.gjf
export inp_state2=$ms_scrdir/state2.gjf
export inp_state3=$ms_scrdir/state3.gjf
export fch_state1=$ms_scrdir/state1.fch
export fch_state2=$ms_scrdir/state2.fch
export fch_state3=$ms_scrdir/state3.fch

# generate input files for states 1, 2, and 3
$multistate_dir/multistate.exe -gen -nst 3 -gin $gaussian_ein -ctp $tem_mstate \
  -in1 $inp_state1 -in2 $inp_state2 -in3 $inp_state3

rm -f $fch_state1 $fch_state2 $fch_state3

# state 1 calculation
g16 -fchk=$fch_state1 $inp_state1

# state 2 calculation
g16 -fchk=$fch_state2 $inp_state2

# state 3 calculation
g16 -fchk=$fch_state3 $inp_state3

# write Gaussian's *.EOu file
$multistate_dir/multistate.exe -mix -nst 3 -chi 400 -gin $gaussian_ein -gou $gaussian_eou \
  -fc1 $fch_state1 -fc2 $fch_state2 -fc3 $fch_state3

