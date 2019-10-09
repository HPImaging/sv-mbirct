
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>

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
	struct CmdLineMBIR cmdline;
	/*struct SysMatrix2D A;*/
	struct minStruct *bandMinMap;
	struct maxStruct *bandMaxMap;
	struct AValues_char ** A_Padded_Map; 
	float* max_num_pointer;	
	unsigned int sum=0;
	int i,j;
	struct timeval tm1,tm2;

	char *ImageReconMask;	/* Image reconstruction mask (determined by ROI) */
	float InitValue;	/* Image data initial condition is read in from a file if available ... */
				/* else intialize it to a uniform image with value InitValue */
	float OutsideROIValue;	/* Image pixel value outside ROI Radius */

	fprintf(stdout, "Starting Reconstruction...\n\n");

	/* read command line */
	readCmdLineMBIR(argc, argv, &cmdline);

	/* read parameters */
	readSystemParams_MBIR(&cmdline, &Image.imgparams, &sinogram.sinoparams, &reconparams);
	printSinoParams3DParallel(&sinogram.sinoparams);
	printImageParams3D(&Image.imgparams);
	printReconParamsQGGMRF3D(&reconparams);
	fprintf(stdout,"\n");

	/* Looks like this is a partial step in splitting volume across nodes */
	/* Note it reassigns image indices, which should be used to specificy the recon volume */
	#if 1
	int numprocs=1;
	if((sinogram.sinoparams.NSlices%numprocs)!=0) {
		fprintf(stderr, "the number of slices must be a multiple of the number of nodes.\n");
		exit(-1);
	} 
	int myChunk=sinogram.sinoparams.NSlices/numprocs;
	int firstSliceOfVolume=sinogram.sinoparams.FirstSliceNumber;	
	/********* THIS!!! *******/
	Image.imgparams.Nz=myChunk;
	Image.imgparams.FirstSliceNumber=firstSliceOfVolume;  
	/*************************/
	sinogram.sinoparams.NSlices=myChunk;
	sinogram.sinoparams.FirstSliceNumber=firstSliceOfVolume;
	#endif

	/* Read System Matrix */
	int pieceLength=computePieceLength(sinogram.sinoparams.NViews);
	if(sinogram.sinoparams.NViews%pieceLength!=0){
		fprintf(stderr, "Error: NViews mod pieceLength must be 0\n");
		fprintf(stderr, "Exiting %s\n",argv[0]);
		exit(-1);

	}
	for(i=0;i<Image.imgparams.Ny;i+=(SVLength*2-overlappingDistance1))
	for(j=0;j<Image.imgparams.Nx;j+=(SVLength*2-overlappingDistance2))
		sum++;

	bandMinMap = (struct minStruct *)get_spc(sum,sizeof(struct minStruct));
	bandMaxMap = (struct maxStruct *)get_spc(sum,sizeof(struct maxStruct));
	A_Padded_Map = (struct AValues_char **)multialloc(sizeof(struct AValues_char), 2, sum, (SVLength*2+1)*(SVLength*2+1));	

	max_num_pointer = (float *)malloc(Image.imgparams.Ny*Image.imgparams.Nx*sizeof(float));	
	char fname[200];
	sprintf(fname,"%s.2Dsysmatrix",cmdline.SysMatrixFile);
	readAmatrix(fname,A_Padded_Map,max_num_pointer,&Image.imgparams,&sinogram.sinoparams,sum,bandMinMap,bandMaxMap,pieceLength);

	/* Read Sinogram and Weights */
	AllocateSinoData3DParallel(&sinogram);
	ReadSinoData3DParallel(cmdline.SinoDataFile, &sinogram);
	ReadWeights3D(cmdline.SinoWeightsFile, &sinogram);

	/* Allocate memory for image */
	AllocateImageData3D(&Image);

	/* allocate and generate recon mask based on ROIRadius--do this before image initialization */
	ImageReconMask = GenImageReconMask(&(Image.imgparams));

	/* Initialize image and reconstruction mask */
	//InitValue = reconparams.MuWater;
	//OutsideROIValue = reconparams.MuAir;
	InitValue = MUWATER;  /* careful..the initial forward projection is written out in GenSysMatrix, so the inital image has to match */
	OutsideROIValue = 0;
	Initialize_Image(&Image, &cmdline, ImageReconMask, InitValue, OutsideROIValue);

	gettimeofday(&tm1,NULL);

	MBIRReconstruct3D(&Image,&sinogram,reconparams,ImageReconMask,bandMinMap,bandMaxMap,A_Padded_Map,max_num_pointer,&cmdline,sum,pieceLength);

	gettimeofday(&tm2,NULL);
	unsigned long long tt = 1000 * (tm2.tv_sec - tm1.tv_sec) + (tm2.tv_usec - tm1.tv_usec) / 1000;
	printf("run time %llu ms \n", tt);

	fprintf(stdout, "Done with reconstruction. Writing output files...\n");

	/* Write out reconstructed image */
	WriteImage3D(cmdline.ReconImageDataFile, &Image);

	/* free image, sinogram and system matrix memory allocation */
	FreeImageData3D(&Image);
	FreeSinoData3DParallel(&sinogram);
	free((char *)ImageReconMask);	    
    
	return(0);
}


