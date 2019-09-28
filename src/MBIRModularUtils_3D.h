#ifndef MBIR_MODULAR_UTILS_3D_H
#define MBIR_MODULAR_UTILS_3D_H

/* The following utilities are used for managing data structures and files associated */
/* with the Modular MBIR Framework */

/* 3D Sinogram Parameters */
struct SinoParams3DParallel
{
    int NChannels;         /* Number of channels in detector */
    float DeltaChannel;    /* Detector spacing (mm) */
    float CenterOffset;    /* Offset of center-of-rotation ... */
                           /* Computed from center of detector in increasing direction (no. of channels) */
                           /* This can be fractional though */
    int NViews;            /* Number of view angles */
    float *ViewAngles;     /* Array of NTheta view angle entries in degrees */
    
    int NSlices;            /* Number of rows (slices) stored in Sino array */
    float DeltaSlice;      /* Spacing along row (slice) direction (mm) */
    int FirstSliceNumber;   /* Row (slice) index coresponding to first row (slice) stored in Sino array */
                            /* This is in absolute coordinates and is used if a partial set of slices is needed */
                            /* Otherwise, it is set to 0. */
};


/* 3D Sinogram Data Structure */
struct Sino3DParallel
{
  struct SinoParams3DParallel sinoparams; /* Sinogram Parameters */
  float **sino;           /* The array is indexed by sino[Slice][ View * NChannels + Channel ] */
                          /* If data array is empty, then set Sino = NULL */
  float **weight;         /* Weights for each measurement */
};

/* 3D Image parameters*/
struct ImageParams3D
{
    int Nx;                 /* Number of columns in image */
    int Ny;                 /* Number of rows in image */
    float Deltaxy;          /* Spacing between pixels in x and y direction (mm) */
    float ROIRadius;        /* Radius of the reconstruction (mm) */
    
    float DeltaZ;           /* Spacing between pixels in z direction (mm) [This should be equal to DeltaSlice */
    int Nz;                 /* Number of rows (slices) in image */
    int FirstSliceNumber;   /* Detector row (slice) index cooresponding to first row (slice) stored in Image array */
                            /* This is in absolute coordinates and is used if a partial set of slices is needed */
                            /* Otherwise, it is set to 0. */
};

/* 3D Image Data Structure */
struct Image3D
{
  struct ImageParams3D imgparams; /* Image parameters */
  float **image;                  /* The array is indexed by image[SliceIndex][ Row * Nx + Column ], Nx=NColumns */
                                  /* If data array is empty, then set Image = NULL */
};


/* Reconstruction Parameters Data Structure */
struct ReconParamsQGGMRF3D
{
  double p;               /* q-GGMRF p parameter */
  double q;               /* q-GGMRF q parameter (q=2 is typical choice) */
  double T;               /* q-GGMRF T parameter */
  double SigmaX;          /* q-GGMRF sigma_x parameter (mm-1) */
  double SigmaY;          /* Scaling constant for weight matrix (W<-W/SigmaY^2); */
                          /* If SigmaY=0, then it is estimated */
  double b_nearest;       /* Relative nearest neighbor weight [default = 1] */
  double b_diag;          /* Relative diagonal neighbor weight in (x,y) plane [default = 1/sqrt(2)] */
  double b_interslice;    /* Relative neighbor weight along z direction [default = 1] */
    
  int Positivity;         /* Positivity constraint: 1=yes, 0=no */
  double StopThreshold;   /* Stopping threshold in percent */
  int MaxIterations;      /* Maximum number of iterations */
    
  double InitImageValue;  /* Initial Condition pixel value. In our examples usually chosen as ... */
};


/**********************************************/
/*  Utilities for reading/writing 3D sinogram */
/**********************************************/

/* Utility for reading 3D parallel beam sinogram parameters */
/* Returns 0 if no error occurs */
int ReadSinoParams3DParallel(
  char *fname,                               /* Input: Reads sinogram parameters from <fname>.sinoparams */
  struct SinoParams3DParallel *sinoparams);  /* Output: Reads sinogram parameters into data structure */

/* Utility for writing out 3D parallel beam sinogram parameters and data */
/* Returns 0 if no error occurs */
int WriteSino3DParallel(
                        char *fname,             /* Input: Writes sinogram parameters to <fname>.sinoparams and data (if available) to ... */
                                                 /* <fname>_slice<InitialIndex>.2Dsinodata to <fname>_slice<FinalIndex>.2Dsinodata  */
                        struct Sino3DParallel *sinogram); /* Input: Writes out sinogram parameters and data */

/* Utility for writing out weights for 3D parallel beam sinogram data */
/* Returns 0 if no error occurs */
int WriteWeights3D(
                char *fname,             /* Input: Writes sinogram measurement weights <fname>_slice<InitialIndex>.2Dweightdata to <fname>_slice<FinalIndex>.2Dweightdata */
                struct Sino3DParallel *sinogram); /* Input: Sinogram data structure */

/* Utility for reading 3D parallel beam sinogram data */
/* Warning: Memory must be allocated before use */
/* Returns 0 if no error occurs */
int ReadSinoData3DParallel(
                           char *fname,               /* Input: Reads sinogram data from <fname>_slice<InitialIndex>.2Dsinodata to <fname>_slice<FinalIndex>.2Dsinodata */
                           struct Sino3DParallel *sinogram);  /* Input/Output: Uses sinogram parameters and reads sinogram data into data structure */

int ReadWeights3D(
                 char *fname,             /* Input: Read sinogram measurement weights from <fname>_slice<InitialIndex>.2Dweightdata to <fname>_slice<FinalIndex>.2Dweightdata */
                 struct Sino3DParallel *sinogram); /* Input: Stores weights into Sinogram Data Structure  */


/* Utility for allocating memory for Sino */
/* Returns 0 if no error occurs */
int AllocateSinoData3DParallel(
                               struct Sino3DParallel *sinogram);  /* Input: Sinogram parameters data structure */

/* Utility for freeing memory allocated for ViewAngles and Sino */
/* Returns 0 if no error occurs */
int FreeSinoData3DParallel(
                           struct Sino3DParallel *sinogram); /* Input: Sinogram parameters data structure */

/*******************************************/
/* Utilities for reading/writing 3D images */
/*******************************************/

/* VS : Utility for reading 3D Image parameters */
/* Returns 0 if no error occurs */
int ReadImageParams3D(
  char *fname,                         /* Input: Reads image type parameter from <fname>.imgparams */
  struct ImageParams3D *imgparams);    /* Output: Reads image parameters into data structure */

/* Utility for reading 3D image data */
/* Warning: Memory must be allocated before use */
/* Returns 0 if no error occurs */
int ReadImage3D(
                char *fname,              /* Input: Reads 2D image data from <fname>_slice<InitialIndex>.2Dimgdata to <fname>_slice<FinalIndex>.2Dimgdata */
                struct Image3D *Image);   /* Output: Reads Image parameters and data (if available) into data structure */


/* Utility for writing 3D image parameters and data */
/* Returns 0 if no error occurs */
int WriteImage3D(
                 char *fname,              /* Input: Writes to image parameters to <fname>.imgparams and data (if available) to .. */
                                           /* <fname>_slice<InitialIndex>.2Dimgdata to <fname>_slice<FinalIndex>.2Dimgdata */
                 struct Image3D *Image);   /* Input: Image data structure (both data and params) */


/* Utility for allocating memory for Image */
/* Returns 0 if no error occurs */
int AllocateImageData3D(
                        struct Image3D *Image);  /* Input: Image data structure */


/* Utility for freeing memory for Image */
/* Returns 0 if no error occurs */
int FreeImageData3D(
                    struct Image3D *Image);    /* Input: Image data structure */


/**************************************************/
/* Utilities for reading in reconstruction params */
/**************************************************/

int ReadReconParamsQGGMRF3D(
                             char *fname,
                             struct ReconParamsQGGMRF3D *reconparams);

/***********************************/
/* Miscellanous Functions         */
/* Remove or shift them out later */
/***********************************/

void printReconParamsQGGMRF3D(struct ReconParamsQGGMRF3D *reconparams);
void printImageParams3D(struct ImageParams3D *imgparams);
void printSinoParams3DParallel(struct SinoParams3DParallel *sinoparams);


#endif /* MBIR_MODULAR_UTILS_3D_H */

