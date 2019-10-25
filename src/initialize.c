
#include <stdio.h>
#include "mbir_ct.h"
#include "MBIRModularDefs.h"
#include "MBIRModularUtils.h"
#include "allocate.h"
#include "initialize.h"


void initSVParams(struct SVParams *svpar,struct ImageParams3D imgparams,struct SinoParams3DParallel sinoparams)
{
	int i,j;
	svpar->SVLength=SVLENGTH;
	svpar->overlap=OVERLAPPINGDISTANCE;
	svpar->SVDepth=SVDEPTH;
	svpar->Nsv=0;
	svpar->pieceLength=computePieceLength(sinoparams.NViews);

        for(i=0;i<imgparams.Ny;i+=(svpar->SVLength*2-svpar->overlap))
        for(j=0;j<imgparams.Nx;j+=(svpar->SVLength*2-svpar->overlap))
                svpar->Nsv++;

	#if 0
	fprintf(stdout,"SVlength %d\n",svpar->SVLength);
	fprintf(stdout,"overlap %d\n",svpar->overlap);
	fprintf(stdout,"SVDepth %d\n",svpar->SVDepth);
	fprintf(stdout,"Nsv %d\n",svpar->Nsv);
	fprintf(stdout,"pieceLength %d\n",svpar->pieceLength);
	#endif
}


/* The pieceLength is the block size in the super-voxel buffer. From past 
 * experiments a good block size is about 1/16 of the views but it has to divide
 * evenly into Nviews. For example, if we have 900 views, and 900/16 = 56.25, 
 * check integers from 56 to 1 in a descending order. The first integer divisible 
 * is a viable block size.
 */
int computePieceLength(int NViews)
{
        int pieceLength=NViews/16;
        if(pieceLength<1)
                pieceLength=1;

        while( pieceLength>1 && (NViews%pieceLength!=0) )
        {
                pieceLength--;
        }
        //fprintf(stderr, "Nviews %d, pieceLength %d\n",NViews,pieceLength);

        return(pieceLength);
}


/* Initialize image state */
void Initialize_Image(
	struct Image3D *Image,
	struct CmdLine *cmdline,
	char *ImageReconMask,
	float InitValue,
	float OutsideROIValue)
{
    int j,jz;
    int Nxy = Image->imgparams.Nx * Image->imgparams.Ny;
    int Nz = Image->imgparams.Nz;

    //fprintf(stdout, "\nInitializing Image ... \n");
    
    if(cmdline->InitImageDataFileFlag == 0) /* Image file not available */
    {
        //GenConstImage(Image, InitValue);               /* generate image with uniform pixel value */
    }
    else {
        //ReadImage3D(cmdline->InitImageDataFile, Image); /* read image which has values in HU units */
        fprintf(stdout, "Note initial image feature not implemented--using constant initial condition\n");
    }

    /* Generate constant image */
    // ***move this up when projector is fixed
    for(jz=0; jz<Nz; jz++)
    for(j=0; j<Nxy; j++)
    if(ImageReconMask[j]==0)
        Image->image[jz][j] = OutsideROIValue;
    else
        Image->image[jz][j] = InitValue;

}


/* Allocate and generate Image Reconstruction mask */
char *GenImageReconMask(struct ImageParams3D *imgparams)
{
    int jx, jy, jz, Nx, Ny, Nz;
    float x_0, y_0, Deltaxy, x, y, yy, ROIRadius, R_sq, R_sq_max;
    char *ImageReconMask;
    
    Nx = imgparams->Nx;
    Ny = imgparams->Ny;
    Nz = imgparams->Nz;
    Deltaxy = imgparams->Deltaxy;
    ROIRadius = imgparams->ROIRadius;
    x_0 = -(Nx-1)*Deltaxy/2;
    y_0 = -(Ny-1)*Deltaxy/2;
    
    ImageReconMask = (char *)get_spc(Ny*Nx,sizeof(char));
    
    if (ROIRadius < 0.0)
    {
        printf("Error in GenImageReconMask : Invalid Value for Radius of Reconstruction \n");
        exit(-1);
    }
    else
    {
        R_sq_max = ROIRadius*ROIRadius;
        for (jy = 0; jy < Ny; jy++)
        {
            y = y_0 + jy*Deltaxy;
            yy = y*y;
            for (jx = 0; jx < Nx; jx++)
            {
                x = x_0 + jx*Deltaxy;
                R_sq = x*x + yy;
                if (R_sq > R_sq_max)
                    ImageReconMask[jy*Nx+jx] = 0;
                else
                    ImageReconMask[jy*Nx+jx] = 1;
            }
        }
    }
    return(ImageReconMask);
}


/* Normalize weights to sum to 1 */
/* Only neighborhood specific */
void NormalizePriorWeights3D(struct ReconParamsQGGMRF3D *reconparams)
{
    double sum;
    
    /* All neighbor weights must sum to one, assuming 8 pt neighborhood */
    sum = 4.0*reconparams->b_nearest + 4.0*reconparams->b_diag + 2.0*reconparams->b_interslice;
    
    reconparams->b_nearest /= sum;
    reconparams->b_diag /= sum;
    reconparams->b_interslice /= sum;
}



