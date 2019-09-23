
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
    int NChannels;         /* Number of channels in detector */
    float DeltaChannel;    /* Detector spacing (mm) */
    float CenterOffset;    /* Offset of center-of-rotation ... */
                           /* Computed from center of detector in increasing direction (no. of channels) */
                           /* This can be fractional though */
    int NViews;            /* Number of view angles */
    float *ViewAngles;     /* Array of NTheta view angle entries in degrees */
};


/* 2D Sinogram Data Structure */
struct Sino2DParallel
{
  struct SinoParams2DParallel sinoparams; /* Sinogram Parameters */
  float *sino;           /* Array of sinogram entries */
                         /* The array is indexed by 2DSino[ NChannels*view + channel ] */
                         /* If data array is empty, then set Sino = NULL */
  float *weight;         /* Weights for each measurement */
};

/* 2D Image parameters*/
struct ImageParams2D
{
    int Nx;                 /* Number of columns in image */
    int Ny;                 /* Number of rows in image */
    float Deltaxy;          /* Spacing between pixels in x and y direction (mm) */
    float ROIRadius;        /* Radius of the reconstruction (mm) */
};

/* 2D Image Data Structure */
struct Image2D
{
  struct ImageParams2D imgparams; /* Image parameters */
  float *image;                   /* Output: Array of image entries */
                                  /* The array is indexed by image[ Nx*row + column ] */
                                  /* If data array is empty, then set Image = NULL */
};


/* Reconstruction Parameters Data Structure */
struct ReconParamsQGGMRF2D
{
  float p;               /* q-GGMRF p parameter */
  float q;               /* q-GGMRF q parameter (q=2 is typical choice) */
  float T;               /* q-GGMRF T parameter */
  float SigmaX;          /* q-GGMRF sigma_x parameter (mm-1) */
  float SigmaY;          /* Scaling constant for weight matrix (W<-W/SigmaY^2); */
                          /* If SigmaY=0, then it is estimated */
  float b_nearest;       /* Relative nearest neighbor weight [default = 1] */
  float b_diag;          /* Relative diagonal neighbor weight in (x,y) plane [default = 1/sqrt(2)] */
  int Positivity;         /* Options: MBIR_MODULAR_YES or MBIR_MODULAR_NO */
  float StopThreshold;   /* Stopping threshold in percent */
  int MaxIterations;      /* Maximum number of iterations */
    
  float MuWater;         /* Attenuation coefficient of water (mm^-1) [default = 0.0202527 mm-1] */
  float MuAir ;          /* Attenuation coefficient of air [default = 0.0 mm-1] */
};


/* VS- Introduced this to make coding easier and modular (easier extension to 3D tomography applications) */
/* Sparse Column Vector - Data Structure */
/*
struct SparseColumn
{
    int Nnonzero;  
    int *RowIndex; 
    float *Value; 
};
*/

/* Sparse System Matrix Data Structure */
/*
struct SysMatrix2D
{
  int Ncolumns;                
  struct SparseColumn *column; 
};
*/


/**********************************************/
/*  Utilities for reading/writing 2D sinogram */
/**********************************************/


/* Utility for reading 2D parallel beam sinogram parameters */
/* Returns 0 if no error occurs */
int ReadSinoParams2DParallel(
  char *fname,                               /* Input: Reads sinogram parameters from <fname>.sinoparams */
  struct SinoParams2DParallel *sinoparams);  /* Output: Reads sinogram parameters into data structure */

/* Utility for writing out 2D parallel beam sinogram parameters and data */
/* Returns 0 if no error occurs */
int WriteSino2DParallel(
                        char *fname,             /* Input: Writes sinogram parameters to <fname>.sinoparams and data (if available) to <fname>.2dsinodata */
                        struct Sino2DParallel *sinogram);/* Input: Writes out sinogram parameters and data */

int WriteWeights2D(
                char *fname,             /* Input: Writes sinogram measurement weights <fname>.wght */
                struct Sino2DParallel *sinogram); /* Input: Sinogram data structure */

/* Utility for reading 2D parallel beam sinogram data */
/* Warning: Memory must be allocated before use */
/* Returns 0 if no error occurs */
int ReadSinoData2DParallel(
                           char *fname,               /* Input: Reads sinogram data from <fname>.2dsinodata */
                           struct Sino2DParallel *sinogram);  /* Input/Output: Uses sinogram parameters and reads sinogram data into data structure */

int ReadWeights2D(
                 char *fname,             /* Input: Read sinogram measurement weights from <fname>.wght */
                 struct Sino2DParallel *sinogram); /* Input: Stores weights into Sinogram Data Structure  */


/* Utility for allocating memory for Sino */
/* Returns 0 if no error occurs */
int AllocateSinoData2DParallel(
                               struct Sino2DParallel *sinogram);  /* Input: Sinogram parameters data structure */

/* Utility for freeing memory allocated for ViewAngles and Sino */
/* Returns 0 if no error occurs */
int FreeSinoData2DParallel(
                           struct Sino2DParallel *sinogram); /* Input: Sinogram parameters data structure */

/*******************************************/
/* Utilities for reading/writing 2D images */
/*******************************************/

/* VS : Utility for reading 2D Image parameters */
/* Returns 0 if no error occurs */
int ReadImageParams2D(
  char *fname,                         /* Input: Reads image type parameter from <fname>.imgparams */
  struct ImageParams2D *imgparams);    /* Output: Reads image parameters into data structure */

/* Utility for reading 2D image parameters and data */
/* Warning: Memory must be allocated before use */
/* Returns 0 if no error occurs */
int ReadImage2D(
                char *fname,              /* Input: Reads 2D image data from <fname>.2dimgdata */
                struct Image2D *Image);   /* Output: Reads Image parameters and data (if available) into data structure */


/* Utility for writing 2D image parameters and data */
/* Returns 0 if no error occurs */
int WriteImage2D(
                 char *fname,              /* Input: Writes to image parameters to <fname>.imgparams and data (if available) to <fname>.2dimgdata */
                 struct Image2D *Image);   /* Input: Image data structure (both data and params) */


/* Utility for allocating memory for Image */
/* Returns 0 if no error occurs */
int AllocateImageData2D(
                        struct Image2D *Image);  /* Input: Image data structure */


/* Utility for freeing memory for Image */
/* Returns 0 if no error occurs */
int FreeImageData2D(
                    struct Image2D *Image);    /* Input: Image data structure */

/******************************************************/
/* Utilities for reading/writing sparse System matrix */
/******************************************************/

/* Utility for reading/allocating the Sparse System Matrix */
/* Returns 0 if no error occurs */
/*
int ReadSysMatrix2D(
  char *fname,                
  struct SysMatrix2D *A);     
*/                              

/* Utility for writing the Sparse System Matrix */
/* Returns 0 if no error occurs */
/*
int WriteSysMatrix2D(
  char *fname,                
  struct SysMatrix2D *A);     
*/

/* Utility for freeing memory from Sparse System Matrix */
/* Returns 0 if no error occurs */
/*
int FreeSysMatrix2D(
  struct SysMatrix2D *A);      
*/
/**************************************************/
/* Utilities for reading in reconstruction params */
/**************************************************/

int ReadReconParamsQGGMRF2D(
                             char *fname,
                             struct ReconParamsQGGMRF2D *reconparams);

/***********************************/
/* Miscellanous Functions         */
/* Remove or shift them out later */
/***********************************/

void printReconParamsQGGMRF2D(struct ReconParamsQGGMRF2D *reconparams);
void printImageParams2D(struct ImageParams2D *imgparams);
void printSinoParams2DParallel(struct SinoParams2DParallel *sinoparams);


#endif /* MBIR_MODULAR_UTILS_2D_H*/


