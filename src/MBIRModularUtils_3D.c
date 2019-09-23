
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>	/* strcmp */
#include <getopt.h>	/* getopt */

#include "MBIRModularUtils_2D.h"
#include "MBIRModularUtils_3D.h"
#include "allocate.h"


/**********************************************/
/*  Utilities for reading/writing 3D sinogram */
/**********************************************/


void printSinoParams3DParallel(struct SinoParams3DParallel *sinoparams)
{
    fprintf(stdout, "\nSINOGRAM PARAMETERS:\n");
    fprintf(stdout, " - Number of sinogram views per slice    = %d\n", sinoparams->NViews);
    fprintf(stdout, " - Number of detector channels per slice = %d\n", sinoparams->NChannels);
    fprintf(stdout, " - Number of slices                      = %d\n", sinoparams->NSlices);
    fprintf(stdout, " - Spacing between Detector Channels     = %.7f (mm) \n", sinoparams->DeltaChannel);
    fprintf(stdout, " - center of rotation offset             = %.7f (channels)\n", sinoparams->CenterOffset);
    fprintf(stdout, " - Spacing between slices                = %.7f (mm)\n", sinoparams->DeltaSlice);
    fprintf(stdout, " - First Slice Index                     = %d \n", sinoparams->FirstSliceNumber);
}

/* Utility for reading 3D parallel beam sinogram parameters */
/* Returns 0 if no error occurs */
int ReadSinoParams3DParallel(
                             char *fname,                               /* Input: Reads sinogram parameters from <fname>.sinoparams */
                             struct SinoParams3DParallel *sinoparams)  /* Output: Reads sinogram parameters into data structure */
{
    FILE *fp;
    char tag[200];
    char AngleListFileName[200];
    int i;
    
    strcat(fname,".sinoparams"); /* append file extension */
    
    if ((fp = fopen(fname, "r")) == NULL)
    {
        fprintf(stderr, "ERROR in ReadSinoParams3DParallel: can't open file %s.\n", fname);
        exit(-1);
    }
    
    fgets(tag, 200, fp);
    fscanf(fp, "%d\n", &(sinoparams->NChannels));
    if(sinoparams->NChannels <= 0){
        fprintf(stderr,"ERROR in ReadSinoParams3DParallel: Number of channels must be a positive integer. And it must be specified.\n");
        exit(-1);
    }
    
    fgets(tag, 200, fp);
    fscanf(fp, "%d\n", &(sinoparams->NViews));
    if(sinoparams->NViews <= 0){
        fprintf(stderr,"ERROR in ReadSinoParams3DParallel: Number of views must be a positive integer. And it must be specified.\n");
        exit(-1);
    }
    
    fgets(tag, 200, fp);
    fscanf(fp, "%f\n", &(sinoparams->DeltaChannel));
    if(sinoparams->DeltaChannel <= 0){
        fprintf(stderr,"ERROR in ReadSinoParams3DParallel: Detector-channel spacing must be a positive floating point. And it must be specified.\n");
        exit(-1);
    }
    
    fgets(tag, 200, fp);
    fscanf(fp, "%f\n", &(sinoparams->CenterOffset)); /* NOTE : THIS IS IN UNITS OF NUMBER OF CHANNELS RATHER THAN ACTUAL DISPLACEMENT in mm */
    if(fabs(sinoparams->CenterOffset) >= sinoparams->NChannels){
        fprintf(stderr,"ERROR in ReadSinoParams3DParallel: Detector-center offset cannot be greater than number of channels \n");
        exit(-1);
    }

    fgets(tag, 200, fp);
    fscanf(fp, "%d\n", &(sinoparams->NSlices));
    if(sinoparams->NSlices <= 0){
        fprintf(stderr,"ERROR in ReadSinoParams3DParallel: Number of slices must be a positive integer. And it must be specified.\n");
        exit(-1);
    }
    
    fgets(tag, 200, fp);
    fscanf(fp, "%f\n", &(sinoparams->DeltaSlice));
    if(sinoparams->DeltaSlice <= 0){
        fprintf(stderr,"ERROR in ReadSinoParams3DParallel: Spacing between slices must be a positive floating point. And it must be specified.\n");
        exit(-1);
    }
    
    fgets(tag, 200, fp);
    fscanf(fp, "%d\n", &(sinoparams->FirstSliceNumber));
    if(sinoparams->FirstSliceNumber <= 0){
        fprintf(stderr,"ERROR in ReadSinoParams3DParallel: First slice must be a non-negative integer. And it must be specified.\n");
        exit(-1);
    }
    
    fgets(tag, 200, fp);
    fscanf(fp, "%s\n", AngleListFileName); /* List of View angles */
    
    fclose(fp);
    
    /* Read in list of View Angles */
    sinoparams->ViewAngles = (float *)get_spc(sinoparams->NViews, sizeof(float));
    
    if ((fp = fopen(AngleListFileName, "r")) == NULL)
    {
        fprintf(stderr, "ERROR in ReadSinoParams3DParallel: can't open file containing list of view angles %s.\n", AngleListFileName);
        exit(-1);
    }
    
    for(i=0;i<sinoparams->NViews;i++)
    {
        if(fscanf(fp,"%f\n",&(sinoparams->ViewAngles[i])) == 0)
       {
         fprintf(stderr, "ERROR in ReadSinoParams3DParallel: List of view angles in file %s terminated early.\n", AngleListFileName);
         exit(-1);
       }
    }
    
    fclose(fp);
    
    return 0;
}


/* Utility for reading 3D parallel beam sinogram data */
/* Warning: Memory must be allocated before use */
/* Returns 0 if no error occurs */
int ReadSinoData3DParallel(
                           char *fname,   /* Input: Reads sinogram data from <fname>_slice<InitialIndex>.2Dsinodata to <fname>_slice<FinalIndex>.2Dsinodata */
                           struct Sino3DParallel *sinogram)  /* Input/Output: Uses sinogram parameters and reads sinogram data into data structure */
{
    char slicefname[200];
    char *sliceindex;
    int i,NSlices,NChannels,NViews,FirstSliceNumber;
    struct Sino2DParallel SingleSliceSinogram;
    
    strcat(fname,"_slice"); /* <fname>_slice */
    sliceindex= (char *)malloc(MBIR_MODULAR_MAX_NUMBER_OF_SLICE_DIGITS);
    
    NSlices = sinogram->sinoparams.NSlices ;
    NChannels = sinogram->sinoparams.NChannels;
    NViews = sinogram->sinoparams.NViews;
    FirstSliceNumber=sinogram->sinoparams.FirstSliceNumber;
    
    /* Copy necessary slice information */
    SingleSliceSinogram.sinoparams.NChannels = NChannels;
    SingleSliceSinogram.sinoparams.NViews = NViews;
 
    // printf("\nReading 3-D Projection Data ... \n");
    
    for(i=0;i<NSlices;i++)
    {
        SingleSliceSinogram.sino = sinogram->sino[i];  /* pointer to beginning of data for i-th slice */
        
        /* slice index : integer to string conversion with fixed precision */
        sprintf(sliceindex,"%.*d",MBIR_MODULAR_MAX_NUMBER_OF_SLICE_DIGITS,i+FirstSliceNumber);
        
        /* Obtain file name for the given slice */
        strcpy(slicefname,fname);
        strcat(slicefname,sliceindex); /* append slice index */
        
        if(ReadSinoData2DParallel(slicefname, &SingleSliceSinogram))
        {   fprintf(stderr, "Error in ReadSinoData3DParallel : Unable to read sinogram data for slice %d from file %s \n",i,slicefname);
            exit(-1);
        }
    }
    
    free((void *)sliceindex);    
    
    return 0;
}


/* Utility for reading weights for 3D sinogram projections data */
/* Warning: Memory must be allocated before use */
/* Returns 0 if no error occurs */
int ReadWeights3D(
                    char *fname,       /* Input: Reads sinogram data from <fname>_slice<InitialIndex>.2Dweightdata to <fname>_slice<FinalIndex>.2Dweightdata */
                    struct Sino3DParallel *sinogram) /* Input/Output: Uses sinogram parameters and reads sinogram data into data structure */
{
    char slicefname[200];
    char *sliceindex;
    int i,NSlices,NChannels,NViews,FirstSliceNumber;
    struct Sino2DParallel SingleSliceSinogram;
    
    strcat(fname,"_slice"); /* <fname>_slice */
    sliceindex= (char *)malloc(MBIR_MODULAR_MAX_NUMBER_OF_SLICE_DIGITS);
    
    NSlices = sinogram->sinoparams.NSlices ;
    NChannels = sinogram->sinoparams.NChannels;
    NViews = sinogram->sinoparams.NViews;
    FirstSliceNumber = sinogram->sinoparams.FirstSliceNumber;
    
    /* Copy necessary slice information */
    SingleSliceSinogram.sinoparams.NChannels = NChannels;
    SingleSliceSinogram.sinoparams.NViews = NViews;
    
    // printf("\nReading 3-D Sinogram Weights Data ... \n");
    
    for(i=0;i<NSlices;i++)
    {
        SingleSliceSinogram.weight = sinogram->weight[i]; /* pointer to beginning of data for i-th slice */
        
        /* slice index : integer to string conversion with fixed precision */
        sprintf(sliceindex,"%.*d",MBIR_MODULAR_MAX_NUMBER_OF_SLICE_DIGITS,i+FirstSliceNumber);
        
        /* Obtain file name for the given slice */
        strcpy(slicefname,fname);
        strcat(slicefname,sliceindex); /* append slice index */
        
        if(ReadWeights2D(slicefname, &SingleSliceSinogram))
        {   fprintf(stderr, "Error in ReadWeights3D : Unable to read sinogram weight data for slice %d from file %s \n",i,slicefname);
            exit(-1);
        }
    }
    
    free((void *)sliceindex);    
    
    return 0;
}

/* Utility for writing out 3D parallel beam sinogram parameters and data */
/* Returns 0 if no error occurs */
int WriteSino3DParallel(
                        char *fname,             /* Input: Writes sinogram parameters to <fname>.sinoparams and data (if available) to ... */
                                                 /* <fname>_slice<InitialIndex>.2Dsinodata to <fname>_slice<FinalIndex>.2Dsinodata  */
                        struct Sino3DParallel *sinogram)  /* Input: Writes out sinogram parameters and data */
{
    char slicefname[200];
    char *sliceindex;
    int i,NSlices,NChannels,NViews,FirstSliceNumber;
    struct Sino2DParallel SingleSliceSinogram;
    
    strcat(fname,"_slice"); /* <fname>_slice */
    sliceindex= (char *)malloc(MBIR_MODULAR_MAX_NUMBER_OF_SLICE_DIGITS);
    
    NSlices = sinogram->sinoparams.NSlices ;
    NChannels = sinogram->sinoparams.NChannels;
    NViews = sinogram->sinoparams.NViews;
    FirstSliceNumber = sinogram->sinoparams.FirstSliceNumber;
    
    /* Copy necessary slice information */
    SingleSliceSinogram.sinoparams.NChannels = NChannels;
    SingleSliceSinogram.sinoparams.NViews = NViews;
    
    printf("\nWriting 3-D Projection Data ... \n");
    
    for(i=0;i<NSlices;i++)
    {
        SingleSliceSinogram.sino = sinogram->sino[i];  /* pointer to beginning of data for i-th slice */
        
        /* slice index : integer to string conversion with fixed precision */
        sprintf(sliceindex,"%.*d",MBIR_MODULAR_MAX_NUMBER_OF_SLICE_DIGITS,i+FirstSliceNumber);
        
        /* Obtain file name for the given slice */
        strcpy(slicefname,fname);
        strcat(slicefname,sliceindex); /* append slice index */
        
        if(WriteSino2DParallel(slicefname, &SingleSliceSinogram))
        {   fprintf(stderr, "Error in WriteSinoData3DParallel : Unable to write sinogram data for slice %d from file %s \n",i,slicefname);
            exit(-1);
        }
    }

    free((void *)sliceindex);    
    
    return 0;
}


/* Utility for writing out weights for 3D parallel beam sinogram data */
/* Returns 0 if no error occurs */
int WriteWeights3D(
                   char *fname,        /* Input: Writes sinogram measurement weights <fname>_slice<InitialIndex>.2Dweightdata to <fname>_slice<FinalIndex>.2Dweightdata */
                   struct Sino3DParallel *sinogram) /* Input: Sinogram data structure */
{
    char slicefname[200];
    char *sliceindex;
    int i,NSlices,NChannels,NViews,FirstSliceNumber;
    struct Sino2DParallel SingleSliceSinogram;
    
    strcat(fname,"_slice"); /* <fname>_slice */
    sliceindex= (char *)malloc(MBIR_MODULAR_MAX_NUMBER_OF_SLICE_DIGITS);
    
    NSlices = sinogram->sinoparams.NSlices ;
    NChannels = sinogram->sinoparams.NChannels;
    NViews = sinogram->sinoparams.NViews;
    FirstSliceNumber = sinogram->sinoparams.FirstSliceNumber;
    
    /* Copy necessary slice information */
    SingleSliceSinogram.sinoparams.NChannels = NChannels;
    SingleSliceSinogram.sinoparams.NViews = NViews;
    
    printf("\nWriting 3-D Sinogram Weights Data ... \n");
    
    for(i=0;i<NSlices;i++)
    {
        SingleSliceSinogram.weight = sinogram->weight[i];  /* pointer to beginning of data for i-th slice */
        
        /* slice index : integer to string conversion with fixed precision */
        sprintf(sliceindex,"%.*d",MBIR_MODULAR_MAX_NUMBER_OF_SLICE_DIGITS,i+FirstSliceNumber);
        
        /* Obtain file name for the given slice */
        strcpy(slicefname,fname);
        strcat(slicefname,sliceindex); /* append slice index */
        
        if(WriteWeights2D(slicefname, &SingleSliceSinogram))
        {   fprintf(stderr, "Error in WriteWeights3D: Unable to write sinogram weight data for slice %d from file %s \n",i,slicefname);
            exit(-1);
        }
    }
    
    free((void *)sliceindex);    
    
    return 0;
}

/* Utility for allocating memory for Sino */
/* Returns 0 if no error occurs */
int AllocateSinoData3DParallel(
                               struct Sino3DParallel *sinogram)  /* Input: Sinogram parameters data structure */
{
    // printf("\nAllocating Sinogram Memory ... \n");

	sinogram->weight = (float **)multialloc(sizeof(float), 2, sinogram->sinoparams.NSlices,sinogram->sinoparams.NViews * sinogram->sinoparams.NChannels);
	sinogram->sino = (float **)multialloc(sizeof(float), 2, sinogram->sinoparams.NSlices,sinogram->sinoparams.NViews * sinogram->sinoparams.NChannels);   
 
    return 0;
}


/* Utility for freeing memory allocated for ViewAngles and Sino */
/* Returns 0 if no error occurs */
int FreeSinoData3DParallel(
                           struct Sino3DParallel *sinogram)  /* Input: Sinogram parameters data structure */
{
    multifree(sinogram->sino,2);
    multifree(sinogram->weight,2);
    return 0;
}

/*******************************************/
/* Utilities for reading/writing 3D images */
/*******************************************/

/* VS : Utility for reading 2D Image parameters */
/* Returns 0 if no error occurs */
int ReadImageParams3D(
                      char *fname,                         /* Input: Reads image type parameter from <fname>.imgparams */
                      struct ImageParams3D *imgparams)     /* Output: Reads image parameters into data structure */
{
    FILE *fp;
    char tag[200];
    
    strcat(fname,".imgparams"); /* append file extension */
    
    if ((fp = fopen(fname, "r")) == NULL)
    {
        fprintf(stderr, "ERROR in ReadImageParams3D: can't open file %s.\n", fname);
        exit(-1);
    }
    
    fgets(tag, 200, fp);
    fscanf(fp, "%d\n", &(imgparams->Nx));
    if(imgparams->Nx <= 0){
        fprintf(stderr,"ERROR in ReadImageParams3D: No. of pixels along horizontal direction, Nx, must be a positive integer. And it must be specified.\n");
        exit(-1);
    }
    
    fgets(tag, 200, fp);
    fscanf(fp, "%d\n", &(imgparams->Ny));
    if(imgparams->Ny <= 0){
        fprintf(stderr,"ERROR in ReadImageParams3D: No. of pixels along vertical direction, Ny, must be a positive integer. And it must be specified. \n");
        exit(-1);
    }
    
    fgets(tag, 200, fp);
    fscanf(fp, "%f\n", &(imgparams->Deltaxy));
    if(imgparams->Deltaxy <= 0){
        fprintf(stderr,"ERROR in ReadImageParams3D: Pixel Dimension (mm) must be a positive floating point. And it must be specified.\n");
        exit(-1);
    }
    
    fgets(tag, 200, fp);
    fscanf(fp, "%f\n", &(imgparams->ROIRadius));
    if(imgparams->ROIRadius <= 0){
        fprintf(stderr,"ERROR in ReadImageParams3D: Region-of-Interest Radius (mm) must be a positive floating point. And it must be specified.\n");
        exit(-1);
    }
    
    fgets(tag, 200, fp);
    fscanf(fp, "%d\n", &(imgparams->Nz));
    if(imgparams->Nz <= 0){
        fprintf(stderr,"ERROR in ReadImageParams3D: Number of slices must be a positive integer. And it must be specified.\n");
        exit(-1);
    }
    
    fgets(tag, 200, fp);
    fscanf(fp, "%f\n", &(imgparams->DeltaZ));
    if(imgparams->DeltaZ <= 0){
        fprintf(stderr,"ERROR in ReadImageParams3D: spacing between slices (mm) must be a positive floating point. And it must be specified.\n");
        exit(-1);
    }
    
    fgets(tag, 200, fp);
    fscanf(fp, "%d\n", &(imgparams->FirstSliceNumber));
    if(imgparams->FirstSliceNumber <= 0){
        fprintf(stderr,"ERROR in ReadImageParams3D: First Slice Index must be a positive integer. And it must be specified.\n");
        exit(-1);
    }
    
    fclose(fp);
    return 0 ;
}


void printImageParams3D(struct ImageParams3D *imgparams)
{
    fprintf(stdout, "\nIMAGE PARAMETERS:\n");
    fprintf(stdout, " - Number of Pixels within a single slice in X direction = %d\n", imgparams->Nx);
    fprintf(stdout, " - Number of Pixels within a single slice in Y direction = %d\n", imgparams->Ny);
    fprintf(stdout, " - Number of Slices to reconstruct                       = %d \n", imgparams->Nz);
    fprintf(stdout, " - Pixel width  in XY plane                              = %.7f (mm)\n", imgparams->Deltaxy);
    fprintf(stdout, " - Spacing between slices                                = %.7f (mm)\n", imgparams->DeltaZ);
    fprintf(stdout, " - First Slice Index                                     = %d \n", imgparams->FirstSliceNumber);
    fprintf(stdout, " - ROIRadius                                             = %.7f (mm)\n", imgparams->ROIRadius);
}
           
/* Utility for reading 3D image data */
/* Warning: Memory must be allocated before use */
/* Returns 0 if no error occurs */
int ReadImage3D(
                char *fname,              /* Input: Reads 2D image data from <fname>_slice<InitialIndex>.2Dimgdata to <fname>_slice<FinalIndex>.2Dimgdata */
                struct Image3D *Image)    /* Output: Reads Image parameters and data (if available) into data structure */
{
    char slicefname[200];
    char *sliceindex;
    int i,Nx,Ny,Nz,FirstSliceNumber;
    struct Image2D SingleSliceImage;
    
    strcat(fname,"_slice"); /* <fname>_slice */
    sliceindex= (char *)malloc(MBIR_MODULAR_MAX_NUMBER_OF_SLICE_DIGITS);
    
    Nx = Image->imgparams.Nx ;
    Ny = Image->imgparams.Ny;
    Nz = Image->imgparams.Nz;
    FirstSliceNumber = Image->imgparams.FirstSliceNumber;
    
    /* Copy necessary slice information */
    SingleSliceImage.imgparams.Nx = Nx;
    SingleSliceImage.imgparams.Ny = Ny;
    
    for(i=0;i<Nz;i++)
    {
        SingleSliceImage.image = Image->image[i];  /* pointer to beginning of data for i-th slice */
        
        /* slice index : integer to string conversion with fixed precision */
        sprintf(sliceindex,"%.*d",MBIR_MODULAR_MAX_NUMBER_OF_SLICE_DIGITS,i+FirstSliceNumber);
        
        /* Obtain file name for the given slice */
        strcpy(slicefname,fname);
        strcat(slicefname,sliceindex); /* append slice index */
        
        if(ReadImage2D(slicefname, &SingleSliceImage))
        {   fprintf(stderr, "Error in ReadImage3D : Unable to read image data for slice %d from file %s \n",i,slicefname);
            exit(-1);
        }
    }

    free((void *)sliceindex);    
    
    return 0;
}
           

/* Utility for writing 3D image parameters and data */
/* Returns 0 if no error occurs */
int WriteImage3D(
                char *fname,            /* Input: Writes to image parameters to <fname>.imgparams and data (if available) to .. */
                                        /* <fname>_slice<InitialIndex>.2Dimgdata to <fname>_slice<FinalIndex>.2Dimgdata */
                struct Image3D *Image)  /* Input: Image data structure (both data and params) */
{
    char slicefname[200];
    char *sliceindex;
    int i,Nx,Ny,Nz,FirstSliceNumber;
    struct Image2D SingleSliceImage;
    
    strcat(fname,"_slice"); /* <fname>_slice */
    sliceindex= (char *)malloc(MBIR_MODULAR_MAX_NUMBER_OF_SLICE_DIGITS);
    
    Nx = Image->imgparams.Nx ;
    Ny = Image->imgparams.Ny;
    Nz = Image->imgparams.Nz;
    FirstSliceNumber = Image->imgparams.FirstSliceNumber;
    
    /* Copy necessary slice information */
    SingleSliceImage.imgparams.Nx = Nx;
    SingleSliceImage.imgparams.Ny = Ny;
    
    // printf("\nWriting 3-D Image ... \n");
    for(i=0;i<Nz;i++)
    {
        SingleSliceImage.image = Image->image[i];  /* pointer to beginning of data for i-th slice */
        
        /* slice index : integer to string conversion with fixed precision */
        sprintf(sliceindex,"%.*d",MBIR_MODULAR_MAX_NUMBER_OF_SLICE_DIGITS,i+FirstSliceNumber);
        
        /* Obtain file name for the given slice */
        strcpy(slicefname,fname);
        strcat(slicefname,sliceindex); /* append slice index */
        
        if(WriteImage2D(slicefname, &SingleSliceImage))
        {   fprintf(stderr, "Error in WriteImage3D : Unable to write image data for slice %d from file %s \n",i,slicefname);
            exit(-1);
        }
    }
    
    free((void *)sliceindex);
    
    return 0;
}
           

/* Utility for allocating memory for Image */
/* Returns 0 if no error occurs */
int AllocateImageData3D(
                        struct Image3D *Image)    /* Input: Image data structure */
{
    Image->image = (float **)multialloc(sizeof(float), 2, Image->imgparams.Nz,Image->imgparams.Ny*Image->imgparams.Nx); 
    return 0;
}

/* Utility for freeing memory for Image */
/* Returns 0 if no error occurs */
int FreeImageData3D(
                    struct Image3D *Image)    /* Input: Image data structure */
{
    multifree(Image->image,2);
    return 0;
}


/**************************************************/
/* Utilities for reading in reconstruction params */
/**************************************************/

/* Read prior model information */

int ReadReconParamsQGGMRF3D(char *fname, struct ReconParamsQGGMRF3D *reconparams)
{
    FILE *fp;
    char tag[200];
    
    strcat(fname,".reconparams"); /* append file extension */
    
    if ((fp = fopen(fname, "r")) == NULL)
    {
        fprintf(stderr, "ERROR in ReadReconParamsQGGMRF3D: can't open file %s.\n", fname);
        exit(-1);
    }
    
    fgets(tag, 200, fp);
    fscanf(fp, "%f\n", &(reconparams->MuWater));
    if(reconparams->MuWater<= 0){
        fprintf(stderr,"ERROR in ReadReconParamsQGGMRF3D: Attenuation coefficient of water (mm^-1) cannot be a negative number \n");
        exit(-1);
    }
    
    fgets(tag, 200, fp);
    fscanf(fp, "%f\n", &(reconparams->q));
    if(reconparams->q<= 1) {
        fprintf(stderr,"ERROR in ReadReconParamsQGGMRF3D: Parameter q must be greater than 1, typical choice is q=2 \n");
        exit(-1);
    }
    
    fgets(tag, 200, fp);
    fscanf(fp, "%f\n", &(reconparams->p));
    if(reconparams->p< 1 || reconparams->p >= reconparams->q){
        fprintf(stderr,"ERROR in ReadReconParamsQGGMRF3D: Parameter p must be in the range [1,q) \n");
        exit(-1);
    }
    
    fgets(tag, 200, fp);
    fscanf(fp, "%f\n", &(reconparams->T));
    if(reconparams->T <= 0){
        fprintf(stderr,"ERROR in ReadReconParamsQGGMRF3D: Parameter T is must be greater than zero \n");
        exit(-1);
    }
    
    fgets(tag, 200, fp);
    fscanf(fp, "%f\n", &(reconparams->SigmaX));
    if(reconparams->SigmaX <= 0){
        fprintf(stderr,"ERROR in ReadReconParamsQGGMRF3D: Parameter SigmaX must be greater than zero \n");
        exit(-1);
    }
    
    /* b_nearest and b_diag are neighborhood weights (2D) */
    /* default values : b_nearest=1 and b_diag=1 */
    
    fgets(tag, 200, fp);
    fscanf(fp, "%f\n", &(reconparams->b_nearest));
    if(reconparams->b_nearest <= 0){
        fprintf(stderr,"ERROR in ReadReconParamsQGGMRF3D: Parameter b_nearest must be greater than zero \n");
        exit(-1);
    }
    
    fgets(tag, 200, fp);
    fscanf(fp, "%f\n", &(reconparams->b_diag));
    if(reconparams->b_diag <= 0){
        fprintf(stderr,"ERROR in ReadReconParamsQGGMRF3D: Parameter b_diag must be greater than zero \n");
        exit(-1);
    }
    
    fgets(tag, 200, fp);
    fscanf(fp, "%f\n", &(reconparams->b_interslice));
    if(reconparams->b_interslice < 0){
        fprintf(stderr,"ERROR in ReadReconParamsQGGMRF3D: Parameter b_diag must be greater than zero \n");
        exit(-1);
    }
    
    fgets(tag, 200, fp);
    fscanf(fp, "%f\n", &(reconparams->StopThreshold));
    //SJK: allow 0 (or negative) StopThreshold to disable and "run maximum iterations"
    //if(reconparams->StopThreshold <= 0){
    //    fprintf(stderr,"ERROR in ReadReconParamsQGGMRF3D: Stop Threshold in %% must be greater than zero \n");
    //    exit(-1);
    //}
    
    fgets(tag, 200, fp);
    fscanf(fp, "%d\n",  &(reconparams->MaxIterations));
    if(reconparams->MaxIterations <= 0){
        fprintf(stderr,"ERROR in ReadReconParamsQGGMRF3D: Maximum no. of iterations must be a positive integer \n");
        exit(-1);
    }
    
    fgets(tag, 200, fp);
    fscanf(fp, "%d\n", &(reconparams->NSlices));
    if(reconparams->NSlices <= 0){
        fprintf(stderr,"ERROR in ReadReconParamsQGGMRF3D: Number of slices must be a positive integer. And it must be specified.\n");
        exit(-1);
    }
    
    fgets(tag, 200, fp);
    fscanf(fp, "%d\n", &(reconparams->FirstSliceNumber));
    if(reconparams->FirstSliceNumber <= 0){
        fprintf(stderr,"ERROR in ReadReconParamsQGGMRF3D: First Slice Index must be a positive integer. And it must be specified.\n");
        exit(-1);
    }

    
    reconparams->MuAir = 0.0 ; /* default value */
    
    fclose(fp);
    
    return 0;
}

/* Print prior model information */

void printReconParamsQGGMRF3D(struct ReconParamsQGGMRF3D *reconparams)
{
    fprintf(stdout, "\nPRIOR PARAMETERS:\n");
    fprintf(stdout, " - Q-GGMRF Prior Parameter, q                                    = %f\n", reconparams->p);
    fprintf(stdout, " - Q-GGMRF Prior Parameter, p                                    = %f\n", reconparams->q);
    fprintf(stdout, " - Q-GGMRF Prior Parameter, T                                    = %f \n", reconparams->T);
    fprintf(stdout, " - Prior Regularization parameter, sigmaX                        = %.7f (mm^-1)\n", reconparams->SigmaX);
    fprintf(stdout, " - Linear attenuation coefficient of water                       = %.7f (mm^-1)\n", reconparams->MuWater);
    fprintf(stdout, " - Linear attenuation coefficient of air                         = %.7f (mm^-1)\n", reconparams->MuAir);
    fprintf(stdout, " - Stop threshold for convergence (avg update)                   = %.7f \n", reconparams->StopThreshold);
    fprintf(stdout, " - Maximum number of ICD iterations                              = %d\n", reconparams->MaxIterations);
    fprintf(stdout, " - Prior weight for nearest voxel neighbors within same slice    = %.7f\n", reconparams->b_nearest);
    fprintf(stdout, " - Prior weight for diagonal voxel neighbors within same slice   = %.7f\n", reconparams->b_diag);
    fprintf(stdout, " - Prior weight for nearest voxel neighbors from adjacent slices = %.7f\n", reconparams->b_interslice);
    fprintf(stdout, " - Number of slices to reconstruct                               = %d\n", reconparams->NSlices);
    fprintf(stdout, " - First slice Index                                             = %d\n", reconparams->FirstSliceNumber);
}





