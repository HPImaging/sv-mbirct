#!/bin/bash

# This is a script for running the "sv-mbirct" program for 
# parallel beam computed tomography. Specifically, it will 
# reconstruct sample data available for download at 
# "http://github.com/HPImaging/mbir-demos.git".
#
# More information on the command line usage can be found
# in the Readme file, and also by running "./mbir_ct -help".

# Set to number of physical cores
export OMP_NUM_THREADS=20
export OMP_DYNAMIC=true

cd "$(dirname $0)"

### Set executable and data locations
execdir="../bin"

dataDir="."
dataName="shepp"
#dataName="xradia"

# sample sub-folder organization
parName="$dataDir/$dataName/par/$dataName"
sinoName="$dataDir/$dataName/sino/$dataName"
#wgtName="$dataDir/$dataName/weight/$dataName"
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

# Run separate utility that generates a hash for identifying a
# matrix file by the problem geometry
HASH="$(./genMatrixHash.sh $parName)"
if [[ $? -eq 0 ]]; then
   matName="$matDir/$HASH"
else
   echo "Matrix hash generation failed. Can't read parameter files?"
   [[ -f "$matName.2Dsvmatrix" ]] && /bin/rm "$matName.2Dsvmatrix" 
fi

### MBIR calls

# Compute system matrix and write to file (only have to run this once for a given geometry)
#   -i specifies image parameter file
#   -j specifies sinogram parameter file
#   -m specifies output system matrix file basename
#   -v verbose level 0=quiet, 1=show progress (default), 2=show even more

if [[ ! -f "$matName.2Dsvmatrix" ]]; then
    $execdir/mbir_ct -i $parName -j $parName -m $matName -v 2
else
   echo "System matrix file found: $matName.2Dsvmatrix"
   touch $matName.2Dsvmatrix  # reset modification time
fi

# Reconstruct using pre-computed matrix (leave off -m to compute matrix internally)
#   -m specifies input output system matrix file basename (optional)
#   -k specifies reconstruction parameter file
#   -s specifies input sinogram file basename
#   -w specifies input sinogram weights (optional)
#   -r specifies output image file basename
$execdir/mbir_ct -m $matName -i $parName -j $parName -k $parName -s $sinoName -r $recName -v 2

# For fun, re-project the reconstruction image
#   -t specifies input image basename
#   -f specifies output projection basename
projName="$dataDir/$dataName/proj/$dataName"
$execdir/mbir_ct -m $matName -i $parName -j $parName -t $recName -f $projName -v 2

# This is reconstruction call that also outputs projection of final image state (usefule for Plug & Play)
#$execdir/mbir_ct -m $matName -i $parName -j $parName -k $parName -s $sinoName -r $recName -f $projName -v 2

# This is reconstruction call that inputs an initial image and initial projection (useful for Plug & Play)
#   -t specifies initial condition for recon
#   -e specifies projection of initial image input
#initImg=$recName
#initProj=$projName
#$execdir/mbir_ct -m $matName -i $parName -j $parName -k $parName -s $sinoName -t $initImg -e $initProj -r $recName -f $projName -v 2

# This computes the backprojection only (no MBIR)
#   -b option flag that calls for back projection
#   no recon parameter file needed
#$execdir/mbir_ct -b -m $matName -i $parName -j $parName -s $sinoName -r $recName -v 2

exit 0


