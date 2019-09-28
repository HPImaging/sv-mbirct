
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "allocate.h"
#include "MBIRModularUtils_2D.h"


/***************************************************/
/*  Utilities for reading/writing 2D System matrix */
/***************************************************/

/* write the System matrix to hard drive */
/* Utility for writing the Sparse System Matrix */
/* Returns 0 if no error occurs */
int WriteSysMatrix2D(
                     char *fname,                /* Input: Basename of output file <fname>.2dsysmatrix */
                     struct SysMatrix2D *A)      /* Input: Sparse system matrix structure */
{
    FILE *fp;
    int i, Nnonzero, Ncolumns;
    
    strcat(fname,".2Dsysmatrix"); /* append file extension */
    
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

/* read the A matrix from hard drive */
/* Utility for reading/allocating the Sparse System Matrix */
/* Returns 0 if no error occurs */
int ReadSysMatrix2D(
                    char *fname,               /* Input: Basename of Sparse System Matrix file <fname>.2dsysmatrix */
                    struct SysMatrix2D *A)     /* Ouput: Sparse system matrix structure */
{                                              /* Warning: Memory is allocated for the data structure inside subroutine */
    FILE *fp;
    int i, Ncolumns, Nnonzero;
    
    strcat(fname,".2Dsysmatrix"); /* append file extension */
    
    /* Allocate memory */
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

/* Utility for freeing memory from Sparse System Matrix */
/* Returns 0 if no error occurs */
int FreeSysMatrix2D(
                    struct SysMatrix2D *A)       /* Input: Free all memory from data structure */
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


/**********************************************/
/*  Utilities for reading/writing 2D sinogram */
/**********************************************/

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

