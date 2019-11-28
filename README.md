# sv-mbirct

### HIGH PERFORMANCE MODEL BASED IMAGE RECONSTRUCTION (MBIR) </br> FOR PARALLEL-BEAM COMPUTED TOMOGRAPHY
*Optimized for Intel multi-core processors*

Source code available at:  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
https://github.com/sjkisner/sv-mbirct

Demo scripts and data files for running this program are available
for download under a separate repository:  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
https://github.com/sjkisner/mbir-demos

Further references on MBIR and the technology in the high-performance implementation of this
code can be found at the bottom of this page, in the documentation accompanying the OpenMBIR
project and on Charles Bouman's website:  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
https://github.com/cabouman/OpenMBIR  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
http://engineering.purdue.edu/~bouman/publications/pub_tomography.html  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
http://engineering.purdue.edu/~bouman/publications/pub_security.html  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
http://engineering.purdue.edu/~bouman/publications/pdf/MBIP-book.pdf

## SYSTEM REQUIREMENTS

1. Intel-based CPU(s)
2. Intel "icc" compiler (included in "Parallel Studio XE", available from Intel for Linux, macOS)

## COMPILING

1. From a terminal prompt, enter the *src* folder and type *make*. If
compiling is successful the binary *mbir_ct* will be created and moved into
the *bin* folder. 
```
cd src  
make
```

Note: Initially after installing Parallel Studio XE, there may be complaints
of missing libraries when linking and running the code.
Most issues can be resolved by executing the following line, which should be
included in your .profile (or .bashrc, or whatever relevant resource file
is executed each time a shell is launched).
```
source /opt/intel/bin/compilervars.sh intel64
```
The path above assumes the install point for Parallel Studio XE was "/opt",
which is the default for Linux and macOS. If not, you may need to poke around
to find the location of "compilervars.sh".

## RUNNING

To print a usage statement:

     ./mbir_ct -help

The program is able to run completely with a single command call, but it's 
usually preferrable to run the reconstruction in two stages. In the 
first stage, the sytem matrix is precomputed and stored, and the second
stage is the reconstruction itself. Both stages use the executable *mbir_ct*.
The system matrix can take a significant time to compute,
however the matrix is fixed for a given geometry and data/image 
dimensions, so the matrix file can be reused for any scan that uses the 
same sinogram and image parameters.

(Note the accompanying demo scripts include a utility that detects whether
the necessary sytem matrix file has already been computed and is available, 
given the input image/sino parameters, and the script automatically reads
the file if available, or computes/stores it if not.)

### Stage 1: Compute and store the System Matrix (and initial projection)

    ./mbir_ct
       -i <basename>[.imgparams]     : Input image parameters
       -j <basename>[.sinoparams]    : Input sinogram parameters
    (plus one or more of the following)
       -m <basename>[.2Dsysmatrix]   : Output matrix file
       -f <basename>[.2Dprojection]  : Output projection of default or input IC
       -f <basename>[_sliceNNN.2Dprojection] -t <basename>[_sliceNNN.2Dimgdata]

In the above arguments, the exensions given in the '[]' symbols must be part
of the file names but should be omitted from the command line.
In the last line that includes both -f and -t arguments, the initial 
projection of the initial condition provided by -t is computed and 
saved to a file(s). Further description of data/image filenames is provided
further down.

Examples: (written as if file names have been assigned 
           to variables in a shell script)

To compute/write the system matrix and the projection of the default initial condition:  

    ./mbir_ct -i $parName -j $parName -m $matName -f $matName

To compute/write only the system matrix:  
 
    ./mbir_ct -i $parName -j $parName -m $matName


Similar to the above but the initial projection is computed for the supplied input image (-t):  

    ./mbir_ct -i $parName -j $parName -m $matName -f proj/$imgName -t init/$imgName

The -m option can be omitted if you only want to compute/store the
projection, however the system matrix will need to be computed in any case.


### Stage 2: Compute MBIR Reconstruction

    ./mbir_ct
       -i <basename>[.imgparams]           : Input image parameters
       -j <basename>[.sinoparams]          : Input sinogram parameters
       -k <basename>[.reconparams]         : Input reconstruction parameters
       -s <basename>[_sliceNNN.2Dsinodata] : Input sinogram projection file(s)
       -w <basename>[_sliceNNN.2Dweightdata] : Input sinogram weight file(s)
       -r <basename>[_sliceNNN.2Dimgdata]  : Output reconstructed image file(s)
    (following are optional)
       -m <basename>[.2Dsysmatrix]         : INPUT matrix (params must match!)
       -t <basename>[_sliceNNN.2Dimgdata]  : Input initial condition image(s)
       -e <basename>[_sliceNNN.2Dprojection] : Input projection of init. cond.
                                           : ** default IC if -t not specified
       -f <basename>[_sliceNNN.2Dprojection] : Output projection of final image
       -p <basename>[_sliceNNN.2Dimgdata]  : Proximal map for Plug & Play
                                           : * -p will apply proximal prior
                                           : * generally use with -t -e -f

Examples:

    ./mbir_ct -i $parName -j $parName -k $parName -s $sinoName \
       -w $wgtNname -r $recName -m $matName -e $projName

If either -m or -e are omitted, the corresponding entity (matrix or
projection) will be computed prior to starting the reconstruction.
The default prior model is a q-QGGMRF with a 10-pt 3D neighborhood, unless
the -p argument is included (Plug & Play).

### Useful forms for Plug & Play mode

The program can perform the inversion step for Plug & Play MBIR. 
The argument -t specifies the initial image state,
and the -r argument specifies the output.

The following example uses the same file names for the input 
(initial state) and the output. These file names can be different
if you want to save the intermediate image state in each ADMM iteration.

    ./mbir_ct -i $parName -j $parName -k $parName -s $sinoName \
       -w $wgtName -m $matName -p $proxmapName -t $imgName -r $imgName

The following form also reads in the projection (-e) of the input image state,
and writes out the projection of the final image state (-f) so that the
projection doesn't have to be re-computed in each *mbir_ct* call.

    ./mbir_ct -i $parName -j $parName -k $parName -s $sinoName \
       -w $wgtName -m $matName -p $proxmapName -t $imgName -r $imgName \
       -e $projName -f $projName

## DESCRIPTION OF DATA FILES

For a description of the file contents for all the data and parameter
files used in this program, see the documentation in the OpenMBIR project
referenced at the top of this readme, directly linked to 
[here](https://github.com/cabouman/OpenMBIR/raw/master/Documentation/MBIR-Modular-specification.docx).

In all the *mbir_ct* command arguments specifying a basename for sinogram or image data,
the associated 3D data is split across files, one file per slice.
The file naming convention is as follows, depending on the file type:

     <basename>_sliceNNN.2Dimgdata
     <basename>_sliceNNN.2Dsinodata
     <basename>_sliceNNN.2Dweightdata
     <basename>_sliceNNN.2Dprojection

where "NNN" is a slice index. The slice indices must be non-negative integers
that include leading zeros and no spaces (e.g. 0000 to 1023). The number of
digits used for the slice index is flexible (up to 4) but must be consistent
across all the files used in a given reconstruction call.

For example, one of the 
[demos](https://github.com/sjkisner/mbir-demos)
contains a subvolume of a microCT scan in following data files:

     xradia_slice0700.2Dsinodata
     xradia_slice0701.2Dsinodata
     xradia_slice0702.2Dsinodata
     xradia_slice0703.2Dsinodata

The reconstruction produces the following output files:

     xradia_slice0700.2Dimgdata
     xradia_slice0701.2Dimgdata
     xradia_slice0702.2Dimgdata
     xradia_slice0703.2Dimgdata

The names *xradia_slice700.2Dsinodata*, *xradia_slice701.2Dsinodata*, etc., would 
also have worked if that were used consistently with all the associated data and image
files.  


## References

##### Xiao Wang, Amit Sabne, Putt Sakdhnagool, Sherman J. Kisner, Charles A. Bouman, and Samuel P. Midkiff, "Massively Parallel 3D Image Reconstruction," *Proceedings of the International Conference for High Performance Computing, Networking, Storage and Analysis (SC'17)*, November 13-16, 2017. (One of three finalists for 2017 ACM Gordon Bell Prize.)

##### Amit Sabne, Xiao Wang, Sherman J. Kisner, Charles A. Bouman, Anand Raghunathan, and Samuel P. Midkiff, "Model-based Iterative CT Imaging Reconstruction on GPUs," *22nd ACM SIGPLAN Symposium on Principles and Practice of Parallel Programming (PPoPP '17)*, February 4-8, 2017.

##### Xiao Wang, K. Aditya Mohan, Sherman J. Kisner, Charles Bouman, and Samuel Midkiff, "Fast voxel line update for time-space image reconstruction," *Proceedings of the IEEE International Conference on Acoustics Speech and Signal Processing (ICASSP)*, pp. 1209-1213, March 20-25, 2016.

##### Xiao Wang, Amit Sabne, Sherman Kisner, Anand Raghunathan, Charles Bouman, and Samuel Midkiff, "High Performance Model Based Image Reconstruction," *21st ACM SIGPLAN Symposium on Principles and Practice of Parallel Programming (PPoPP '16)*, March 12-16, 2016. 

##### Suhas Sreehari, S. Venkat Venkatakrishnan, Brendt Wohlberg, Gregery T. Buzzard, Lawrence F. Drummy, Jeffrey P. Simmons, and Charles A. Bouman, "Plug-and-Play Priors for Bright Field Electron Tomography and Sparse Interpolation," *IEEE Transactions on Computational Imaging*, vol. 2, no. 4, Dec. 2016. 

##### Jean-Baptiste Thibault, Ken Sauer, Charles Bouman, and Jiang Hsieh, "A Three-Dimensional Statistical Approach to Improved Image Quality for Multi-Slice Helical CT," *Medical Physics*, pp. 4526-4544, vol. 34, no. 11, November 2007.
