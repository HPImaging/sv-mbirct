
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>	/* strcmp */
#include <getopt.h>	/* getopt */

#include "MBIRModularUtils_2D.h"
#include "allocate.h"

/***************************************************/
/*  Utilities for reading/writing 2D System matrix */
/***************************************************/


/* write the System matrix to hard drive */
/* Utility for writing the Sparse System Matrix */
/* Returns 0 if no error occurs */
/*
int WriteSysMatrix2D(
                     char *fname,               
                     struct SysMatrix2D *A)      
{
    FILE *fp;
    int i, Nnonzero, Ncolumns;
    
    strcat(fname,".2Dsysmatrix"); 
    
    if ((fp = fopen(fname, "w")) == NULL)
    {
        fprintf(stderr, "ERROR in WriteSysMatrix2D: can't open file %s.\n", fname);
        exit(-1);
    }
    
    Ncolumns = A->Ncolumns;
    
    for (i = 0; i < Ncolumns; i++)
    {
        Nnonzero = A->column[i].Nnonzero;
        fwrite(&Nnonzero, sizeof(int), 1, fp);
        
        if(Nnonzero>0)
        {
            fwrite(A->column[i].RowIndex, sizeof(int), Nnonzero, fp);
            fwrite(A->column[i].Value, sizeof(float), Nnonzero, fp);
        }
    }
    
    fclose(fp);
    
    return 0;
}
*/

/* read the A matrix from hard drive */
/* Utility for reading/allocating the Sparse System Matrix */
/* Returns 0 if no error occurs */
/*
int ReadSysMatrix2D(
                    char *fname,              
                    struct SysMatrix2D *A)     
{                                              
    FILE *fp;
    int i, Ncolumns, Nnonzero;
    
    strcat(fname,".2Dsysmatrix"); 
    
    Ncolumns=A->Ncolumns;
    A->column = (struct SparseColumn *)get_spc(Ncolumns, sizeof(struct SparseColumn));
    
    printf("\nReading System-matrix ... \n");
    
    if ((fp = fopen(fname, "r")) == NULL)
    {
        fprintf(stderr, "ERROR in ReadSysMatrix2D: can't open file %s.\n", fname);
        exit(-1);
    }
    
    for (i = 0; i < Ncolumns; i++)
    {
        fread(&Nnonzero, sizeof(int), 1, fp);
        A->column[i].Nnonzero = Nnonzero;
        
        if(Nnonzero > 0)
        {
            A->column[i].RowIndex = (int *)get_spc(Nnonzero, sizeof(int));
            A->column[i].Value    = (float *)get_spc(Nnonzero, sizeof(float));
            
            if(fread(A->column[i].RowIndex, sizeof(int), Nnonzero, fp)!= Nnonzero)
            {
                fprintf(stderr, "ERROR in ReadSysMatrix2D: file terminated early %s.\n", fname);
                exit(-1);
            }
            
            if(fread(A->column[i].Value, sizeof(float), Nnonzero, fp) != Nnonzero)
            {
                fprintf(stderr, "ERROR in ReadSysMatrix2D: file terminated early %s.\n", fname);
                exit(-1);
            }
        }
    }
    
    fclose(fp);
    
    return 0;
    
}
*/
/* Utility for freeing memory from Sparse System Matrix */
/* Returns 0 if no error occurs */
/*
int FreeSysMatrix2D(
                    struct SysMatrix2D *A)       
{
    int i, Ncolumns;
    
    Ncolumns=A->Ncolumns;
    
    for (i = 0; i < Ncolumns; i++)
    {
        free((void *)A->column[i].RowIndex);
        free((void *)A->column[i].Value);
    }
    
    return 0;
}
*/

/**********************************************/
/*  Utilities for reading/writing 2D sinogram */
/**********************************************/


void printSinoParams2DParallel(struct SinoParams2DParallel *sinoparams)
{
    float DeltaViewAngle ;
    DeltaViewAngle = (sinoparams->ViewAngles[sinoparams->NViews-1]-sinoparams->ViewAngles[0])/sinoparams->NViews ;
    
    fprintf(stdout, "\nSINOGRAM PARAMETERS:\n");
    fprintf(stdout, " - NViews         = %-10d                (Number of view angles)\n", sinoparams->NViews);
    fprintf(stdout, " - NChannels      = %-10d                (Number of channels in detector)\n", sinoparams->NChannels);
    fprintf(stdout, " - DeltaChannel   = %-10f mm             (Detector spacing)\n", sinoparams->DeltaChannel);
    fprintf(stdout, " - CenterOffset   = %-10f channels       (Offset of center-of-rotation)\n", sinoparams->CenterOffset);
}


/* Utility for reading 2D parallel beam sinogram parameters */
/* Returns 0 if no error occurs */
int ReadSinoParams2DParallel(
                             char *fname,                               /* Input: Reads sinogram parameters from <fname>.sinoparams */
                             struct SinoParams2DParallel *sinoparams)  /* Output: Reads sinogram parameters into data structure */
{
    FILE *fp;
    char tag[200];
    char AngleListFileName[200];
    int i;
    
    strcat(fname,".sinoparams"); /* append file extension */
    
    if ((fp = fopen(fname, "r")) == NULL)
    {
        fprintf(stderr, "ERROR in ReadSinoParams2DParallel: can't open file %s.\n", fname);
        exit(-1);
    }
    
    fgets(tag, 200, fp);
    fscanf(fp, "%d\n", &(sinoparams->NChannels));
    if(sinoparams->NChannels <= 0){
        fprintf(stderr,"ERROR in ReadSinoParams2DParallel: Number of channels must be a positive integer. And it must be specified.\n");
        exit(-1);
    }
    
    fgets(tag, 200, fp);
    fscanf(fp, "%d\n", &(sinoparams->NViews));
    if(sinoparams->NViews <= 0){
        fprintf(stderr,"ERROR in ReadSinoParams2DParallel: Number of views must be a positive integer. And it must be specified.\n");
        exit(-1);
    }
    
    fgets(tag, 200, fp);
    fscanf(fp, "%f\n", &(sinoparams->DeltaChannel));
    if(sinoparams->DeltaChannel <= 0){
        fprintf(stderr,"ERROR in ReadSinoParams2DParallel: Detector-channel spacing must be a positive floating point. And it must be specified.\n");
        exit(-1);
    }
    
    fgets(tag, 200, fp);
    fscanf(fp, "%f\n", &(sinoparams->CenterOffset)); /* NOTE : THIS IS IN UNITS OF NUMBER OF CHANNELS RATHER THAN ACTUAL DISPLACEMENT in mm */
    if(fabs(sinoparams->CenterOffset) >= sinoparams->NChannels){
        fprintf(stderr,"ERROR in ReadSinoParams2DParallel: Detector-center offset cannot be greater than number of channels \n");
        exit(-1);
    }

    fgets(tag, 200, fp);
    fscanf(fp, "%s\n", AngleListFileName); /* List of View angles */
    
    fclose(fp);
    
    /* Read in list of View Angles */
    sinoparams->ViewAngles = (float *)get_spc(sinoparams->NViews, sizeof(float));
    
    if ((fp = fopen(AngleListFileName, "r")) == NULL)
    {
        fprintf(stderr, "ERROR in ReadSinoParams2DParallel: can't open file containing list of view angles %s.\n", AngleListFileName);
        exit(-1);
    }
    
    for(i=0;i<sinoparams->NViews;i++)
    {
        if(fscanf(fp,"%f\n",&(sinoparams->ViewAngles[i])) == 0)
       {
         fprintf(stderr, "ERROR in ReadSinoParams2DParallel: List of view angles in file %s terminated early.\n", AngleListFileName);
         exit(-1);
       }
    }
    
    fclose(fp);
    
    return 0;
}


/* Utility for reading 2D parallel beam sinogram data */
/* Warning: Memory must be allocated before use */
/* Returns 0 if no error occurs */
int ReadSinoData2DParallel(
                           char *fname,               /* Input: Reads sinogram data from <fname>.2dsinodata */
                           struct Sino2DParallel *sinogram)  /* Input/Output: Uses sinogram parameters and reads sinogram data into data structure */
{
    FILE *fp;
    int M;
    
    strcat(fname,".2Dsinodata"); /* append file extension */
    
    /* NOTE: As of now only write out sinogram data, no parameters */
    /* Fix that in next version of code */
    
    if ((fp = fopen(fname, "r")) == NULL)
    {
        fprintf(stderr, "ERROR in ReadSino2DParallel: can't open file %s.\n", fname);
        exit(-1);
    }
    
    M = sinogram->sinoparams.NViews * sinogram->sinoparams.NChannels;
    
    if(fread(sinogram->sino,sizeof(float),M,fp)!=M)
    {
        fprintf(stderr, "ERROR in ReadSino2DParallel: file terminated early \n");
        fclose(fp);
        exit(1);
    }
    
    fclose(fp);
    
    return 0;
    
}


/* Utility for reading weights for 2D sinogram projections data */
/* Warning: Memory must be allocated before use */
/* Returns 0 if no error occurs */
int ReadWeights2D(
                  char *fname,             /* Input: Read sinogram measurement weights from <fname>.2Dweightdata */
                  struct Sino2DParallel *sinogram) /* Input: Stores weights into Sinogram Data Structure  */
{
    FILE *fp;
    int M;
    
    strcat(fname,".2Dweightdata"); /* append file extension */
    
    /* NOTE: As of now only write out sinogram data, no parameters */
    /* Fix that in next version of code */
    
    if ((fp = fopen(fname, "r")) == NULL)
    {
        fprintf(stderr, "ERROR in ReadWeights2D: can't open file %s.\n", fname);
        exit(-1);
    }
    
    M = sinogram->sinoparams.NViews * sinogram->sinoparams.NChannels;
    
    if(fread(sinogram->weight,sizeof(float),M,fp)!=M)
    {
        fprintf(stderr, "ERROR in ReadWeights2D: file terminated early \n");
        fclose(fp);
        exit(1);
    }
    
    fclose(fp);
    
    return 0;
}


/* Utility for writing out 2D parallel beam sinogram parameters and data */
/* Returns 0 if no error occurs */
int WriteSino2DParallel(
                        char *fname,            /* Input: Writes sinogram parameters to <fname>.sinoparams and data (if available) to <fname>.2Dsinodata */
                        struct Sino2DParallel *sinogram) /* Input: Writes out sinogram parameters and data */
{
    FILE *fp;
    int M;
    
    strcat(fname,".2Dsinodata"); /* append file extension */
    
    /* NOTE: As of now only write out sinogram data, no parameters */
    /* Fix that in next version of code */
    
    if ((fp = fopen(fname, "w")) == NULL)
    {
        fprintf(stderr, "ERROR in WriteSino2DParallel: can't open file %s.\n", fname);
        exit(-1);
    }
    
    M = sinogram->sinoparams.NViews * sinogram->sinoparams.NChannels;
    
    if(fwrite(sinogram->sino,sizeof(float),M,fp)!=M)
    {
        fprintf(stderr, "ERROR in WriteSino2DParallel: file terminated early \n");
        fclose(fp);
        exit(1);
    }
    
    fclose(fp);
    
    return 0;
    
}

/* Utility for writing weights for 2D sinogram projections data */
/* Returns 0 if no error occurs */
int WriteWeights2D(
                   char *fname,              /* Input: Writes sinogram measurement weights to <fname>.2Dweightdata */
                   struct Sino2DParallel *sinogram)   /* Input: Sinogram data structure */
{
    FILE *fp;
    int M;
    
    strcat(fname,".2Dweightdata"); /* append file extension */
    
    /* NOTE: As of now only write out sinogram data, no parameters */
    /* Fix that in next version of code */
    
    if ((fp = fopen(fname, "w")) == NULL)
    {
        fprintf(stderr, "ERROR in WriteWeights: can't open file %s.\n", fname);
        exit(-1);
    }
    
    M = sinogram->sinoparams.NViews * sinogram->sinoparams.NChannels;
    
    if(fwrite(sinogram->weight,sizeof(float),M,fp)!=M)
    {
        fprintf(stderr, "ERROR in WriteWeights: file terminated early \n");
        fclose(fp);
        exit(1);
    }
    
    fclose(fp);
    
    return 0;
}


/* Utility for allocating memory for Sino */
/* Returns 0 if no error occurs */
int AllocateSinoData2DParallel(
                               struct Sino2DParallel *sinogram)  /* Input: Sinogram parameters data structure */
{
    sinogram->sino   = (float *)get_spc(sinogram->sinoparams.NViews * sinogram->sinoparams.NChannels, sizeof(float));
    sinogram->weight = (float *)get_spc(sinogram->sinoparams.NViews * sinogram->sinoparams.NChannels, sizeof(float));
 
    return 0;
}


/* Utility for freeing memory allocated for ViewAngles and Sino */
/* Returns 0 if no error occurs */
int FreeSinoData2DParallel(
                           struct Sino2DParallel *sinogram)  /* Input: Sinogram parameters data structure */
{
    free((void *)sinogram->sino);
    free((void *)sinogram->weight);
    return 0;
}

/*******************************************/
/* Utilities for reading/writing 2D images */
/*******************************************/

/* VS : Utility for reading 2D Image parameters */
/* Returns 0 if no error occurs */
int ReadImageParams2D(
                      char *fname,                         /* Input: Reads image type parameter from <fname>.imgparams */
                      struct ImageParams2D *imgparams)     /* Output: Reads image parameters into data structure */
{
    FILE *fp;
    char tag[200];
    
    strcat(fname,".imgparams"); /* append file extension */
    
    if ((fp = fopen(fname, "r")) == NULL)
    {
        fprintf(stderr, "ERROR in ReadImageParams2D: can't open file %s.\n", fname);
        exit(-1);
    }
    
    fgets(tag, 200, fp);
    fscanf(fp, "%d\n", &(imgparams->Nx));
    if(imgparams->Nx <= 0){
        fprintf(stderr,"ERROR in ReadImageParams2D: No. of pixels along horizontal direction, Nx, must be a positive integer. And it must be specified.\n");
        exit(-1);
    }
    
    fgets(tag, 200, fp);
    fscanf(fp, "%d\n", &(imgparams->Ny));
    if(imgparams->Ny <= 0){
        fprintf(stderr,"ERROR in ReadImageParams2D: No. of pixels along vertical direction, Ny, must be a positive integer. And it must be specified. \n");
        exit(-1);
    }
    
    fgets(tag, 200, fp);
    fscanf(fp, "%f\n", &(imgparams->Deltaxy));
    if(imgparams->Deltaxy <= 0){
        fprintf(stderr,"ERROR in ReadImageParams2D: Pixel Dimension (mm) must be a positive floating point. And it must be specified.\n");
        exit(-1);
    }
    
    fgets(tag, 200, fp);
    fscanf(fp, "%f\n", &(imgparams->ROIRadius));
    if(imgparams->ROIRadius <= 0){
        fprintf(stderr,"ERROR in ReadImageParams2D: Region-of-Interest Radius (mm) must be a positive floating point. And it must be specified.\n");
        exit(-1);
    }
    
    /* Note : imgparams->Deltaxy needs to be read in. For some reason it is now calculated bade on field of view */
    /* Remove this dependency later once Parameter file is changed */
    
    fclose(fp);
    return 0 ;
}


void printImageParams2D(struct ImageParams2D *imgparams)
{
    fprintf(stdout, "\nIMAGE PARAMETERS:\n");
    fprintf(stdout, " - Nx        = %-10d                     (Number of pixels along x axis)\n", imgparams->Nx);
    fprintf(stdout, " - Ny        = %-10d                     (Number of pixels along y axis)\n", imgparams->Ny);
    fprintf(stdout, " - Deltaxy   = %-10f mm                  (spacing between pixels in x and y direction)\n", imgparams->Deltaxy);
    fprintf(stdout, " - ROIRadius = %-10f mm                  (radius of the reconstruction)\n", imgparams->ROIRadius);
}


/* Here image data read in units of (mm^-1) */
int ReadImage2D(char *fname, struct Image2D *Image)
{
    FILE *fp;
    int N;
    
    strcat(fname,".2Dimgdata"); /* append file extension */
    
    if ((fp = fopen(fname, "r")) == NULL)
    {
        fprintf(stderr, "ERROR in ReadImage2D: can't open file %s.\n", fname);
        exit(-1);
    }
    
    N=Image->imgparams.Nx * Image->imgparams.Ny ;
    
    if(fread(Image->image,sizeof(float),N,fp)!=N)
    {
        fprintf(stderr, "ERROR in ReadImage2D: file terminated early \n");
        fclose(fp);
        exit(1);
    }
    
    
    fclose(fp);
    return 0;
}

/* Utility for writing 2D image parameters and data */
/* Returns 0 if no error occurs */
int WriteImage2D(
                 char *fname,              /* Input: Writes to image parameters to <fname>.imgparams and data (if available) to <fname>.2dimgdata */
                 struct Image2D *Image)    /* Input: Image data structure (both data and params) */
{
    FILE *fp;
    int N;
    
    strcat(fname,".2Dimgdata"); /* append file extension */
    
    if ((fp = fopen(fname, "w")) == NULL)
    {
        fprintf(stderr, "ERROR in WriteImage2D: can't open file %s.\n", fname);
        exit(-1);
    }
    
    N=Image->imgparams.Nx * Image->imgparams.Ny ;
    
    if(fwrite(Image->image,sizeof(float),N,fp)!=N)
    {
        fprintf(stderr, "ERROR in ReadImage2D: file terminated early \n");
        fclose(fp);
        exit(1);
    }
    
    fclose(fp);
    return 0;
}


/* Utility for allocating memory for Image */
/* Returns 0 if no error occurs */
int AllocateImageData2D(
                        struct Image2D *Image)    /* Input: Image data structure */
{
    Image->image = (float *)get_spc(Image->imgparams.Nx * Image->imgparams.Ny, sizeof(float));
    
    return 0;
}

/* Utility for freeing memory for Image */
/* Returns 0 if no error occurs */
int FreeImageData2D(
                    struct Image2D *Image)    /* Input: Image data structure */
{
    free((void *)Image->image);
    return 0;
}


/**************************************************/
/* Utilities for reading in reconstruction params */
/**************************************************/

/* Read prior model information */

int ReadReconParamsQGGMRF2D(char *fname, struct ReconParamsQGGMRF2D *reconparams)
{
    FILE *fp;
    char tag[200];
    
    strcat(fname,".reconparams"); /* append file extension */
    
    if ((fp = fopen(fname, "r")) == NULL)
    {
        fprintf(stderr, "ERROR in ReadReconParamsQGGMRF2D: can't open file %s.\n", fname);
        exit(-1);
    }
    
    fgets(tag, 200, fp);
    fscanf(fp, "%f\n", &(reconparams->MuWater));
    if(reconparams->MuWater<= 0){
        fprintf(stderr,"ERROR in ReadReconParamsQGGMRF2D: Attenuation coefficient of water (mm^-1) cannot be a negative number \n");
        exit(-1);
    }
    
    fgets(tag, 200, fp);
    fscanf(fp, "%f\n", &(reconparams->q));
    if(reconparams->q<= 1) {
        fprintf(stderr,"ERROR in ReadReconParamsQGGMRF2D: Parameter q must be greater than 1, typical choice is q=2 \n");
        exit(-1);
    }
    
    fgets(tag, 200, fp);
    fscanf(fp, "%f\n", &(reconparams->p));
    if(reconparams->p< 1 || reconparams->p >= reconparams->q){
        fprintf(stderr,"ERROR in ReadReconParamsQGGMRF2D: Parameter p must be in the range [1,q) \n");
        exit(-1);
    }
    
    fgets(tag, 200, fp);
    fscanf(fp, "%f\n", &(reconparams->T));
    if(reconparams->T <= 0){
        fprintf(stderr,"ERROR in ReadReconParamsQGGMRF2D: Parameter T is must be greater than zero \n");
        exit(-1);
    }
    
    fgets(tag, 200, fp);
    fscanf(fp, "%f\n", &(reconparams->SigmaX));
    if(reconparams->SigmaX <= 0){
        fprintf(stderr,"ERROR in ReadReconParamsQGGMRF2D: Parameter SigmaX must be greater than zero \n");
        exit(-1);
    }
    
    /* b_nearest and b_diag are neighborhood weights (2D) */
    /* default values : b_nearest=1 and b_diag=1 */
    
    fgets(tag, 200, fp);
    fscanf(fp, "%f\n", &(reconparams->b_nearest));
    if(reconparams->b_nearest <= 0){
        fprintf(stderr,"ERROR in ReadReconParamsQGGMRF2D: Parameter b_nearest must be greater than zero \n");
        exit(-1);
    }
    
    fgets(tag, 200, fp);
    fscanf(fp, "%f\n", &(reconparams->b_diag));
    if(reconparams->b_nearest <= 0){
        fprintf(stderr,"ERROR in ReadReconParamsQGGMRF2D: Parameter b_diag must be greater than zero \n");
        exit(-1);
    }
    
    fgets(tag, 200, fp);
    fscanf(fp, "%f\n", &(reconparams->StopThreshold));
    //SJK: allow 0 (or negative) StopThreshold to disable and "run maximum iterations"
    //if(reconparams->StopThreshold <= 0){
    //    fprintf(stderr,"ERROR in ReadReconParamsQGGMRF2D: Stop Threshold in %% must be greater than zero \n");
    //    exit(-1);
    //}
    
    fgets(tag, 200, fp);
    fscanf(fp, "%d\n",  &(reconparams->MaxIterations));
    if(reconparams->MaxIterations <= 0){
        fprintf(stderr,"ERROR in ReadReconParamsQGGMRF2D: Maximum no. of iterations must be a positive integer \n");
        exit(-1);
    }
    
    reconparams->MuAir = 0.0 ; /* default value */
    
    fclose(fp);
    
    return 0;
}

/* Print prior model information */

void printReconParamsQGGMRF2D(struct ReconParamsQGGMRF2D *reconparams)
{
    fprintf(stdout, "\nPRIOR PARAMETERS:\n");
    fprintf(stdout, " - p              = %-10f                (q-GGMRF p parameter)\n", reconparams->p);
    fprintf(stdout, " - q              = %-10f                (q-GGMRF q parameter)\n", reconparams->q);
    fprintf(stdout, " - T              = %-10f                (q-GGMRF T parameter)\n", reconparams->T);
    fprintf(stdout, " - SigmaX         = %-10f (mm-1)         (q-GGMRF sigma_x parameter)\n", reconparams->SigmaX);
    fprintf(stdout, " - MuWater        = %-10f (mm-1)         (the water attenuation coefficient)\n", reconparams->MuWater);
    fprintf(stdout, " - MuAir          = %-10f (mm-1)         (the air attenuation coefficient)\n", reconparams->MuAir);
    fprintf(stdout, " - StopThreshold  = %-10f %%              (Stopping threshold in percent)\n", reconparams->StopThreshold);
    fprintf(stdout, " - MaxIterations  = %-10d                (Maximum number of iterations)\n", reconparams->MaxIterations);
    fprintf(stdout, " - b_nearest      = %-10f                (Relative nearest neighbor weight)\n", reconparams->b_nearest);
    fprintf(stdout, " - b_diag         = %-10f                (Relative diagonal neighbor weight in (x,y) plane)\n", reconparams->b_diag);
    
}




