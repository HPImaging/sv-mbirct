#ifndef MBIR_MODULAR_UTILS_2D_H
#define MBIR_MODULAR_UTILS_2D_H


/* Define constants that will be used in modular MBIR framework */
#define MBIR_MODULAR_UTIL_VERSION "0.0";

#define MBIR_MODULAR_SINOTYPE_2DPARALLEL 0;
#define MBIR_MODULAR_SINOTYPE_2DFAN 1; /* for future implementation */
#define MBIR_MODULAR_SINOTYPE_3DPARALLEL 2;

#define MBIR_MODULAR_IMAGETYPE_2D 0;
#define MBIR_MODULAR_IMAGETYPE_3D 1;
#define MBIR_MODULAR_IMAGETYPE_4D 2; /* for future implementation */

#define MBIR_MODULAR_RECONTYPE_QGGMRF_2D 0;
#define MBIR_MODULAR_RECONTYPE_QGGMRF_3D 1;
#define MBIR_MODULAR_RECONTYPE_PandP 2 /* for future implementation */

#define MBIR_MODULAR_YES 1;
#define MBIR_MODULAR_NO 0;
#define MBIR_MODULAR_MAX_NUMBER_OF_SLICE_DIGITS 4 /* allows up to 10,000 slices */

#define PI 3.1415926535897932384
#define MUWATER 0.0202527   /* mm-1 */
#define mu2hu(Mu, MuAir, MuWater) (1000.0*(Mu-MuAir)/(MuWater-MuAir)) /* (mm^-1) to HU units conversion */
#define hu2mu(HU, MuAir, MuWater) (HU*(MuWater-MuAir)/1000.0)+MuAir   /* (mm^-1) to HU units conversion */


/* The following utilities are used for managing data structures and files associated */
/* with the Modular MBIR Framework */

struct SinoParams2DParallel
{
   int NChannels;	/* Number of channels in detector */
   float DeltaChannel;	/* Detector spacing (mm) */
   float CenterOffset;	/* Offset of center-of-rotation ... */
			/* Computed from center of detector in increasing direction (no. of channels) */
			/* This can be fractional */
   int NViews;		/* Number of view angles */
   float *ViewAngles;	/* Array of NTheta view angle entries in degrees */
};

/* 2D Sinogram Data Structure */
struct Sino2DParallel
{
   struct SinoParams2DParallel sinoparams; /* Sinogram Parameters */
   float *sino;		/* Array of sinogram entries indexed by sino[NChannels*view + channel] */
   float *weight;	/* Weights for each measurement */
			/* If data arrays empty, then set the pointer = NULL */
};

/* 2D Image parameters*/
struct ImageParams2D
{
   int Nx;		/* Number of columns in image */
   int Ny;		/* Number of rows in image */
   float Deltaxy;	/* Spacing between pixels in x and y direction (mm) */
   float ROIRadius;	/* Radius of the reconstruction (mm) */
};

/* 2D Image Data Structure */
struct Image2D
{
   struct ImageParams2D imgparams;	/* Image parameters */
   float *image;	/* Array of image entries indexed by image[Nx*row + column] */
			/* If data array is empty, then set image = NULL */
};

/* Sparse Column Vector - Data Structure */
struct SparseColumn
{
   int Nnonzero;	/* Nnonzero is the number of nonzero entries in the column */
   int *RowIndex;	/* RowIndex[j] is the row index of the jth nonzero entry in the column */
   float *Value;	/* Value[j] is the value of the jth nonzero entry in the column of the matrix */
};

/* Sparse System Matrix Data Structure */
struct SysMatrix2D
{
   int Ncolumns;		/* Number of columns in sparse matrix */
   struct SparseColumn *column;	/* column[i] is the i-th column of the matrix in sparse format */
};


/**********************************************/
/*  Utilities for reading/writing 2D sinogram */
/**********************************************/

/* Utility for writing out 2D parallel beam sinogram data */
/* Returns 0 if no error occurs */
int WriteSino2DParallel(
	char *fname,				/* Input: Writes sinogram parameters to <fname>.2dsinodata */
	struct Sino2DParallel *sinogram);	/* Input: Writes out sinogram parameters and data */

int WriteWeights2D(
	char *fname,				/* Input: Writes sinogram measurement weights to <fname>.wght */
	struct Sino2DParallel *sinogram);	/* Input: Sinogram data structure */

/* Utility for reading 2D parallel beam sinogram data */
/* Warning: Memory must be allocated before use */
/* Returns 0 if no error occurs */
int ReadSinoData2DParallel(
	char *fname,				/* Input: Reads sinogram data from <fname>.2dsinodata */
	struct Sino2DParallel *sinogram);	/* Input/Output: Uses sinogram parameters and reads into sinogram data structure */

int ReadWeights2D(
	char *fname,				/* Input: Read sinogram weights from <fname>.wght */
	struct Sino2DParallel *sinogram);	/* Input: Stores weights into Sinogram Data Structure  */

/* Utility for allocating memory for Sino */
/* Returns 0 if no error occurs */
int AllocateSinoData2DParallel(
	struct Sino2DParallel *sinogram);	/* Input: Sinogram parameters data structure */

/* Utility for freeing memory allocated for ViewAngles and Sino */
/* Returns 0 if no error occurs */
int FreeSinoData2DParallel(
	struct Sino2DParallel *sinogram);	/* Input: Sinogram parameters data structure */

/*******************************************/
/* Utilities for reading/writing 2D images */
/*******************************************/

/* Utility for reading 2D image parameters and data */
/* Warning: Memory must be allocated before use */
/* Returns 0 if no error occurs */
int ReadImage2D(
	char *fname,		/* Input: Reads 2D image data from <fname>.2dimgdata */
	struct Image2D *Image);	/* Output: Reads into data structure */

/* Utility for writing 2D image parameters and data */
/* Returns 0 if no error occurs */
int WriteImage2D(
	char *fname,		/* Input: Writes to image data to <fname>.2dimgdata */
	struct Image2D *Image);	/* Input: Image data structure (both data and params) */

/* Utility for allocating memory for Image */
/* Returns 0 if no error occurs */
int AllocateImageData2D(
	struct Image2D *Image);	/* Input: Image data structure */

/* Utility for freeing memory for Image */
/* Returns 0 if no error occurs */
int FreeImageData2D(
	struct Image2D *Image);	/* Input: Image data structure */

/******************************************************/
/* Utilities for reading/writing sparse System matrix */
/******************************************************/

/* Utility for reading/allocating the Sparse System Matrix */
/* Returns 0 if no error occurs */
/* Warning: Memory is allocated for the data structure inside subroutine */
int ReadSysMatrix2D(
	char *fname,		/* Input: Basename of Sparse System Matrix file <fname>.2dsysmatrix */
	struct SysMatrix2D *A);	/* Ouput: Sparse system matrix structure */

/* Utility for writing the Sparse System Matrix */
/* Returns 0 if no error occurs */
int WriteSysMatrix2D(
	char *fname,		/* Input: Basename of output file <fname>.2dsysmatrix */
	struct SysMatrix2D *A);	/* Input: Sparse system matrix structure */

/* Utility for freeing memory from Sparse System Matrix */
/* Returns 0 if no error occurs */
int FreeSysMatrix2D(
	struct SysMatrix2D *A);	/* Input: Free all memory from data structure */


#endif /* MBIR_MODULAR_UTILS_2D_H*/

