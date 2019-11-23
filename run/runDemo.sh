#!/bin/bash

export OMP_NUM_THREADS=20
export OMP_DYNAMIC=true

cd "$(dirname $0)"

### Set executable and data locations
execdir="../bin"

dataDir="../../mbir-demos"
dataName="shepp"
#dataName="xradia"

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

### Perform MBIR
#
### Form 1: Reconstruct with a single call
# $execdir/mbir_ct -i $parName -j $parName -k $parName -s $sinoName \
#    -w $wgtName -r $recName -m $matName -e $matName
# exit 0
#
### Form 2: First command call pre-computes system matrix and initial 
###   projection, the second call reconstructs. Advantage is the first call
###   only has to be done once for a given set of image/sinogram parameters.
# $execdir/mbir_ct -i $parName -j $parName -m $matName -f $matName
# $execdir/mbir_ct -i $parName -j $parName -k $parName -s $sinoName \
#    -w $wgtName -r $recName -m $matName -e $matName
# exit 0

### The following form is more interesting. It first checks to see if the
### matrix for the input problem dimensions has previously been computed and
### saved in the $matDir folder. If yes, skips ahead to the reconstruction.

### PRE-COMPUTE STAGE
# generate the hash value and check the exit status
HASH="$(./genMatrixHash.sh $parName)"
if [[ $? -eq 0 ]]; then
   matName="$matDir/$HASH"
else
   echo "Matrix hash generation failed. Can't read parameter files?"
   [[ -f "$matName.2Dsysmatrix" ]] && /bin/rm "$matName.2Dsysmatrix" 
fi

# check for matrix file, and compute if not present
if [[ -f "$matName.2Dsysmatrix" ]]; then
   echo "System matrix file found: $matName.2Dsysmatrix"
else
   echo "Generating system matrix file: $matName.2Dsysmatrix"
   echo "Generating projection file: $matName.2Dprojection"
   $execdir/mbir_ct -i $parName -j $parName -m $matName -f $matName
fi

### RECONSTRUCTION STAGE
# check for initial projection file existance before setting input -e option
if [[ -f "$matName.2Dprojection" ]]; then
   $execdir/mbir_ct -i $parName -j $parName -k $parName -s $sinoName \
      -w $wgtName -r $recName -m $matName -e $matName
else
   $execdir/mbir_ct -i $parName -j $parName -k $parName -s $sinoName \
      -w $wgtName -r $recName -m $matName 
fi

exit 0

### SEE THE FOLLOWING FOR MORE INFORMATION AND EXAMPLE FORMS

### the following prints a usage statement
# $execdir/mbir_ct -help

# Generate System Matrix and/or initial projection
#
# ./mbir_ct 
#	-i <basename>[.imgparams]
#	-j <basename>[.sinoparams]
#  (plus one or more of the following)
#	-m <basename>[.2Dsysmatrix]   (this is an output filename)
#	-f <basename>[.2Dprojection]  (this is an output filename)
#	OR
#	-f <basename>[_sliceNNN.2Dprojection] -t <basename>[_sliceNNN.2Dimgdata]

### Pre-processing Forms ###

### Pre-compute and write system matrix for given img/sino dimensions
# $execdir/mbir_ct -i $parName -j $parName -m $matName

### Pre-compute/write system matrix AND projection of default initial condition
# $execdir/mbir_ct -i $parName -j $parName -m $matName -f $matName

### Similar to above but initial projection is for supplied input image (-t)
# $execdir/mbir_ct -i $parName -j $parName -m $matName \
#    -f proj/$imgName -t init/$imgName

# In the latter two examples, can omit -m option to skip writing system matrix


# Compute MBIR Reconstruction
#
# ./mbir_ct 
#	-i <basename>[.imgparams]
#	-j <basename>[.sinoparams]
#	-k <basename>[.reconparams]
#	-s <basename>[_sliceNNN.2Dsinodata]
#	-w <basename>[_sliceNNN.2Dweightdata]
#	-r <basename>[_sliceNNN.2Dimgdata] 
#  (following are optional)
#	-m <basename>[.2Dsysmatrix]	: optional but encouraged
#	-t <basename>[_sliceNNN.2Dimgdata]
#	-e <basename>[.2Dprojection]	: optional but encouraged
#	OR
#	-e <basename>[_sliceNNN.2Dprojection] -t <basename>[_sliceNNN.2Dimgdata]
#	-f <basename>[_sliceNNN.2Dprojection]

### Reconstruction Examples ###

### Minimum arguments required..slowest by not using pre-computed matrix
# $execdir/mbir_ct -i $parName -j $parName -k $parName -s $sinoName \
#    -w $wgtNname -r $recName

### Use pre-computed matrix and initial projection (both aren't mandatory)
# $execdir/mbir_ct -i $parName -j $parName -k $parName -s $sinoName \
#    -w $wgtNname -r $recName -m $matName -e $projName

### Use provided initial condition for image
# $execdir/mbir_ct -i $parName -j $parName -k $parName -s $sinoName \
#    -w $wgtNname -r $recName -m $matName -t $initName


### Useful forms for Plug & Play mode ###

### Recon step of Plug & Play mode.. -t specifies initial image state
### This example writes out to the same file names as the initial condition
### but it doesn't have to.
#
# $execdir/mbir_ct -i $parName -j $parName -k $parName -s $sinoName \
#    -w $wgtName -m $matName -p $proxmapName -t $imgName -r $imgName

### This form also reads in the projection (-e) of the input image state,
### and writes out the projection of the final image state (-f)
#
# $execdir/mbir_ct -i $parName -j $parName -k $parName -s $sinoName \
#    -w $wgtName -m $matName -p $proxmapName -t $imgName -r $imgName \
#    -e $projName -f $projName


