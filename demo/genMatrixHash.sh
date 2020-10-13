#!/bin/bash
#
# This generates a truncated hash label for creating a unique tag for 
# a Sytem Matrix file. It pulls out the relevant fields from the image
# and sinogram sino parameter files including the list of view angles.
#
# usage 1:  ./genMatrixHash.sh <basename>
#    The argument <basename> includes the relative or full path to the files
#       <basename>.imgparams  and  <basename>.sinoparams
#
# usage 2:  ./genMatrixHash.sh <imgparams_filename> <sinoparams_filename>
#    The filenames include the path (relative or full) and extension.
#    This allows different basenames between the two files.
#
# examples: 
#    If there are parameter files "../data/shepp.{imgparams,sinoparams}"
#    Ex.1: ./genMatrixHash.sh ../data/shepp
#    Ex.2: ./genMatrixHash.sh ../data/shepp.imgparams ../data/shepp.sinoparams

# parse arguments
if [[ "$#" -lt 1 ]]; then
  echo "Not enough arguments"
  exit 1
elif [[ "$#" -eq 1 ]]; then
  imgfile="$1.imgparams"
  sinofile="$1.sinoparams"
else
  imgfile="$1"
  sinofile="$2"
fi
#echo $imgfile
#echo $sinofile

# check img/sino params file existance
if [[ ! -f "$imgfile" ]] || [[ ! -f "$sinofile" ]]; then
  echo "Can't read params files"
  exit 1
fi

# It's forgiving if sino/img arguments are swapped, HOWEVER the views 
# file path is pulled specifically from the last argument
viewsfile="$(cat $imgfile $sinofile |grep ViewAngle |cut -d : -f 2 |tr -d ' ')"
viewsfile="$(dirname $sinofile)/$viewsfile"
#echo "$viewsfile"

# check views list file existance 
if [[ ! -f "$viewsfile" ]]; then
  echo "Can't read view list file $viewsfile"
  exit 1
fi

# Pull the relevant parameter fields; remove spaces; sort to remove field order 
# depencence; add views file; generate and truncate sha1 hash
cat "$imgfile" "$sinofile" \
  | grep -e Nx -e Ny -e Deltaxy -e ROIRadius -e NChannels -e NViews -e DeltaChannel -e CenterOffset \
  | tr -d '[:blank:]' \
  | sort \
  | cat "$viewsfile" - \
  | shasum -a 1 \
  | cut -c 1-20

exit 0

