#!/bin/bash

# This is a script for running the "sv-mbirct" program for 
# parallel beam computed tomography. Specifically, it will 
# reconstruct sample data available for download at 
# "http://github.com/HPImaging/mbir-demos.git".
#
# More information on the command line usage can be found
# in the Readme file, and also by running "./mbir_ct -help".

export OMP_NUM_THREADS=20
export OMP_DYNAMIC=true

cd "$(dirname $0)"

### Set executable and data locations
execdir="../bin"

dataDir="../../mbir-demos"
dataName="shepp"
#dataName="xradia"

if [[ "$#" -gt 0 ]]; then
  dataDir="$(dirname $1)"
  dataName="$(basename $1)"
fi

parName="$dataDir/$dataName/par/$dataName"
sinoName="$dataDir/$dataName/sino/$dataName"
wgtName="$dataDir/$dataName/weight/$dataName"
recName="$dataDir/$dataName/recon/$dataName"

matDir="./sysmatrix"
matName="$matDir/$dataName"

# create folders that hold pre-computed items and output if they don't exist
if [[ ! -d "$matDir" ]]; then
  echo "Creating directory $matDir"
  mkdir "$matDir"
fi
if [[ ! -d "$(dirname $recName)" ]]; then
  echo "Creating directory $(dirname $recName)" 
  mkdir "$(dirname $recName)" 
fi

### Compute reconstruction

### Form 1: Reconstruct with a single call (uncomment the next two lines to use)
# $execdir/mbir_ct -i $parName -j $parName -k $parName -s $sinoName \
#    -w $wgtName -r $recName -m $matName -e $matName
# exit 0

### Form 2: Pre-compute system matrix and initial projection and write to file.
###   Then reconstruct. First call only has to be done once for a given set of
###   image/sinogram dimensions--resolution, physical size, offsets, etc.
# $execdir/mbir_ct -i $parName -j $parName -m $matName -f $matName
# $execdir/mbir_ct -i $parName -j $parName -k $parName -s $sinoName \
#    -w $wgtName -r $recName -m $matName -e $matName
# exit 0

### Form 3: The code below checks if the matrix for the input problem 
###   dimensions was previouly computed and saved in the $matDir folder. If no,
###   it computes and saves the matrix; If yes, it skips to the reconstruction.

### PRE-COMPUTE STAGE
# generate the hash value and check the exit status
HASH="$(./genMatrixHash.sh $parName)"
if [[ $? -eq 0 ]]; then
   matName="$matDir/$HASH"
else
   echo "Matrix hash generation failed. Can't read parameter files?"
   [[ -f "$matName.2Dsvmatrix" ]] && /bin/rm "$matName.2Dsvmatrix" 
fi

# check for matrix file, and compute if not present
if [[ ! -f "$matName.2Dsvmatrix" ]]; then
   echo "Generating system matrix file: $matName.2Dsvmatrix"
   echo "Generating projection file: $matName.2Dprojection"
   $execdir/mbir_ct -i $parName -j $parName -m $matName -f $matName
else
   echo "System matrix file found: $matName.2Dsvmatrix"
fi

touch $matName.lastused

### RECONSTRUCTION STAGE

$execdir/mbir_ct -i $parName -j $parName -k $parName -s $sinoName -w $wgtName \
   -r $recName -m $matName -e $matName 
#   2>&1 | tee $(dirname $recName)/out

exit 0


