
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "MBIRModularDefs.h"
#include "MBIRModularUtils.h"
#include "allocate.h"
#include "A_comp.h"
#include "initialize.h"
#include "recon3d.h"

/* #ifdef STORE_A_MATRIX - This option set by default in A_comp.h */
/* This option is to precompute and store the forward matrix rather than compute it on the fly */

int main(int argc, char *argv[])
{

	struct Image3D Image;
    	struct Sino3DParallel sinogram;
	struct ReconParamsQGGMRF3D reconparams;
	/*struct SysMatrix2D A;*/
	struct minStruct *bandMinMap;
	struct maxStruct *bandMaxMap;
	struct AValues_char ** A_Padded_Map; 
	float* max_num_pointer;	
	unsigned int sum=0;
	int i,j;
	int myChunk=0;	
        struct timeval tm1,tm2;    	   	      
	int numprocs=1;	 		    
    
    	struct CmdLineMBIR cmdline;
    
    	char *ImageReconMask; /* Image reconstruction mask (determined by ROI) */
    	float InitValue;     /* Image data initial condition is read in from a file if available ... */
                          /* else intialize it to a uniform image with value InitValue */
    	float OutsideROIValue;/* Image pixel value outside ROI Radius */

    
	/* read command line */
	readCmdLineMBIR(argc, argv, &cmdline);
    
	/* read parameters */
    	readSystemParams_MBIR (&cmdline, &Image.imgparams, &sinogram.sinoparams, &reconparams);
    	
       	if((sinogram.sinoparams.NSlices%numprocs)!=0)
       	{	
       		fprintf(stderr, "the number of slices must be a multiple of the number of nodes.\n");
                exit(-1);
       	} 
	myChunk=sinogram.sinoparams.NSlices/numprocs;
	int firstSliceOfVolume=sinogram.sinoparams.FirstSliceNumber;	
	Image.imgparams.Nz=myChunk;
	Image.imgparams.FirstSliceNumber=firstSliceOfVolume;  
	sinogram.sinoparams.NSlices=myChunk;
	sinogram.sinoparams.FirstSliceNumber=firstSliceOfVolume;

	/* Read System Matrix */
	int pieceLength=computePieceLength(sinogram.sinoparams.NViews);
        if(sinogram.sinoparams.NViews%pieceLength!=0){
                fprintf(stderr, "Error: NViews mod pieceLength must be 0\n");
                fprintf(stderr, "Exiting %s\n",argv[0]);
                exit(-1);
        }
	for(i=0;i<Image.imgparams.Ny;i+=(SVLength*2-overlappingDistance1)){
	  	for(j=0;j<Image.imgparams.Nx;j+=(SVLength*2-overlappingDistance2)){
	    		sum++;
	  	}
	}     
    	bandMinMap = (struct minStruct *)get_spc(sum,sizeof(struct minStruct));
    	bandMaxMap = (struct maxStruct *)get_spc(sum,sizeof(struct maxStruct));
    	A_Padded_Map = (struct AValues_char **)multialloc(sizeof(struct AValues_char), 2, sum, (SVLength*2+1)*(SVLength*2+1));	

	max_num_pointer = (float *)malloc(Image.imgparams.Ny*Image.imgparams.Nx*sizeof(float));	
	char fname[200];
    	sprintf(fname,"%s.2Dsysmatrix",cmdline.SysMatrixFile);		
	readAmatrix(fname,A_Padded_Map, max_num_pointer,&Image.imgparams,&sinogram.sinoparams, sum,bandMinMap,bandMaxMap,pieceLength);


    /* Read Sinogram and Weights */
    	if(AllocateSinoData3DParallel(&sinogram))
    	{   
    		fprintf(stderr, "Error in allocating sinogram data (and weights) memory through function AllocateSinoData3DParallel \n");
        	exit(-1);
    	}
    	if(ReadSinoData3DParallel(cmdline.SinoDataFile, &sinogram))
    	{   fprintf(stderr, "Error in reading sinogram data from file %s through function ReadSinoData3DParallel \n",cmdline.SinoDataFile);
        	exit(-1);
    	}
    	if(ReadWeights3D(cmdline.SinoWeightsFile, &sinogram))
    	{   fprintf(stderr, "Error in reading sinogram weights from file %s through function ReadWeights3D \n", cmdline.SinoWeightsFile);
        	exit(-1);
    	}
    	
    	/* Allocate memory for image */
    	if(AllocateImageData3D(&Image))
    	{   fprintf(stderr, "Error in allocating memory for image through function AllocateImageData3D \n");
       		exit(-1);
    	}

    	/* Initialize image and reconstruction mask */
	//InitValue = reconparams.MuWater;
	//OutsideROIValue = reconparams.MuAir;
	InitValue = MUWATER;  /* careful..the initial forward projection is written out in GenSysMatrix, so the inital image has to match */
	OutsideROIValue = 0;
	Initialize_Image(&Image, &cmdline, InitValue);
	ImageReconMask = GenImageReconMask(&Image,OutsideROIValue);

	//gettimeofday(&tm1,NULL);

    	/* MBIR - Reconstruction */
    	MBIRReconstruct3D(&Image,&sinogram,reconparams,ImageReconMask,bandMinMap,bandMaxMap,A_Padded_Map,max_num_pointer,&cmdline,sum,pieceLength);
    		
        //gettimeofday(&tm2,NULL);
        //unsigned long long tt = 1000 * (tm2.tv_sec - tm1.tv_sec) + (tm2.tv_usec - tm1.tv_usec) / 1000;
        //printf("run time %llu ms \n", tt);
	
    
    	/* Write out reconstructed image */
    	if(WriteImage3D(cmdline.ReconImageDataFile, &Image))
    	{
        	fprintf(stderr, "Error in writing out reconstructed image file through function WriteImage3D \n");
        	exit(-1);
    	}
    	
	/* free image, sinogram and system matrix memory allocation */
    	if(FreeImageData3D(&Image))
    	{  fprintf(stderr, "Error image memory could not be freed through function FreeImageData3D \n");
        	exit(-1);
    	}
    		
    	if(FreeSinoData3DParallel(&sinogram))
    	{  fprintf(stderr, "Error sinogram memory could not be freed through function FreeSinoData3DParallel \n");
        	exit(-1);
    	}

	free((char *)ImageReconMask);	    
    
    return 0;
}


