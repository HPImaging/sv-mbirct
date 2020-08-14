
#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <string.h>
//#include <time.h>
#include <sys/time.h>

#include "mbir_ct.h"
#include "MBIRModularDefs.h"
#include "MBIRModularUtils.h"
#include "allocate.h"
#include "A_comp.h"
#include "initialize.h"
#include "recon3d.h"

/* Internal Functions */
void readCmdLine(int argc, char *argv[], struct CmdLine *cmdline);
void printCmdLineUsage(char *ExecFileName);
int CmdLineHelpOption(char *string);
void setNumSliceDigits(char *basename, char *ext, int slice, struct SinoParams3DParallel *sinoparams, struct ImageParams3D *imgparams);


int main(int argc, char *argv[])
{
	struct CmdLine cmdline;
	struct Image3D Image;
	struct Image3D ProxMap;
	struct Sino3DParallel sinogram;
	struct ReconParams reconparams;
	struct SVParams svpar;
	struct AValues_char **A_Padded_Map; 
	float *max_num_pointer;	
	char *ImageReconMask;	/* Image reconstruction mask (determined by ROI) */
	char fname[1024];
	struct timeval tm1,tm2;
	unsigned long long tdiff;
	int i,j,jz;
	float **e;

	fprintf(stdout,"SUPER-VOXEL MBIR RECONSTRUCTION FOR 3D PARALLEL-BEAM CT\n");
	fprintf(stdout,"build time: %s, %s\n\n", __DATE__,  __TIME__);

	readCmdLine(argc, argv, &cmdline);

	/* Read image/sino parameter files */
	ReadSinoParams3DParallel(cmdline.SinoParamsFile,&sinogram.sinoparams);
	ReadImageParams3D(cmdline.ImageParamsFile,&Image.imgparams);
	printSinoParams3DParallel(&sinogram.sinoparams);
	printImageParams3D(&Image.imgparams);
	if(cmdline.reconFlag)
	{
		ReadReconParams(cmdline.ReconParamsFile,&reconparams);
		NormalizePriorWeights3D(&reconparams);

		if(cmdline.reconFlag == MBIR_MODULAR_RECONTYPE_QGGMRF_3D)
			printReconParamsQGGMRF3D(&reconparams);
		if(cmdline.reconFlag == MBIR_MODULAR_RECONTYPE_PandP)
			printReconParamsPandP(&reconparams);

		if(reconparams.ReconType != cmdline.reconFlag) {
			fprintf(stdout,"**\nWarning: \"PriorModel\" field in reconparams file doesn't agree with\n");
			fprintf(stdout,"Warning: what the command line is doing. Proceeding anyway.\n**\n");
			reconparams.ReconType = cmdline.reconFlag;
		}
	}

	initSVParams(&svpar, Image.imgparams, sinogram.sinoparams);  /* Initialize/allocate SV parameters */
	fprintf(stdout,"\n");

	/* The image parameters specify the relevant slice range to reconstruct, so re-set the  */
	/* relevant sinogram parameters so it pulls the correct slices and indexes them consistently */
	sinogram.sinoparams.NSlices = Image.imgparams.Nz;
	sinogram.sinoparams.FirstSliceNumber = Image.imgparams.FirstSliceNumber;

	int NvNc = sinogram.sinoparams.NViews * sinogram.sinoparams.NChannels;
	int Nxy = Image.imgparams.Nx * Image.imgparams.Ny;
	int Nz = Image.imgparams.Nz;
	int FirstSliceNumber = Image.imgparams.FirstSliceNumber;

	/* Detect the number of slice number digits in input file names */
	if(cmdline.reconFlag)
		setNumSliceDigits(cmdline.SinoDataFile,"2Dsinodata",FirstSliceNumber,&sinogram.sinoparams,&Image.imgparams);
	else if(cmdline.readInitImageFlag)
		setNumSliceDigits(cmdline.InitImageFile,"2Dimgdata",FirstSliceNumber,&sinogram.sinoparams,&Image.imgparams);

	int NumSliceDigits = Image.imgparams.NumSliceDigits;

	/* Allocate and generate recon mask based on ROIRadius */
	ImageReconMask = GenImageReconMask(&Image.imgparams);

	/* Read/compute/write System Matrix */
	A_Padded_Map = (struct AValues_char **)multialloc(sizeof(struct AValues_char),2,svpar.Nsv,(2*svpar.SVLength+1)*(2*svpar.SVLength+1));
	max_num_pointer = (float *) get_spc(Nxy,sizeof(float));
	if(cmdline.readAmatrixFlag)
	{
		fprintf(stdout,"Reading system matrix...\n");
		sprintf(fname,"%s.2Dsvmatrix",cmdline.SysMatrixFile);
		readAmatrix(fname, A_Padded_Map, max_num_pointer, &Image.imgparams, &sinogram.sinoparams, svpar);
	}
	else
	{
		fprintf(stdout,"Computing system matrix...\n");
	        A_comp(A_Padded_Map,max_num_pointer,svpar,&sinogram.sinoparams,ImageReconMask,&Image.imgparams);
        }
	if(cmdline.writeAmatrixFlag)
	{
		fprintf(stdout,"Writing system matrix...\n");
		sprintf(fname,"%s.2Dsvmatrix",cmdline.SysMatrixFile);
		writeAmatrix(fname,A_Padded_Map,max_num_pointer,&Image.imgparams,&sinogram.sinoparams,svpar);
	}

	/* Initialize image and forward project, if necessary */
	if(cmdline.reconFlag || cmdline.writeProjectionFlag)
	{
		/* Initialize image */
		AllocateImageData3D(&Image);
		if(cmdline.readInitImageFlag)
		{
			fprintf(stdout,"Reading initial image...\n");
			ReadImage3D(cmdline.InitImageFile,&Image);
		}
		else
			initConstImage(&Image, ImageReconMask, MUWATER, 0);

		/* Initialize Forward Projection of initial image */
		e = (float **)multialloc(sizeof(float),2,sinogram.sinoparams.NSlices,NvNc);

		if(cmdline.readInitProjectionFlag)
		{
			if(cmdline.readInitImageFlag)  /* If input image specified, initial projection is multi-slice */
			{
				for(jz=0;jz<Nz;jz++)
				{
					sprintf(fname,"%s_slice%.*d.2Dprojection",cmdline.inputProjectionFile,NumSliceDigits,jz+FirstSliceNumber);
					if(ReadFloatArray(fname,e[jz],NvNc)) {
						fprintf(stderr,"Error: can't read %s\n",fname);
						exit(-1);
					}
				}
			}
			else    /* default initial condition, so all slices have same projection */
			{
				sprintf(fname,"%s.2Dprojection",cmdline.inputProjectionFile);
				if(ReadFloatArray(fname,e[0],NvNc)) {
					fprintf(stdout,"Warning: can't read %s\n",fname);
					fprintf(stdout,"Projecting initial image...\n");
					forwardProject2D(e[0],Image.image[0],A_Padded_Map,max_num_pointer,&sinogram.sinoparams,&Image.imgparams,svpar);
				}
				for(jz=1;jz<Nz;jz++)
					memcpy(e[jz],e[0],NvNc*sizeof(float));
			}
		}
		else	/* Compute initial projection */
		{
			fprintf(stdout,"Projecting initial image...\n");
			gettimeofday(&tm1,NULL);

			//compProjectionError(e,&Image,&sinogram,A_Padded_Map,max_num_pointer,svpar);
			if(cmdline.readInitImageFlag)
			{
				/* here we need to project each slice */
				#pragma omp parallel for schedule(dynamic)
				for(jz=0;jz<Nz;jz++)
					forwardProject2D(e[jz],Image.image[jz],A_Padded_Map,max_num_pointer,&sinogram.sinoparams,&Image.imgparams,svpar);
			}
			else	/* if using default initial image, only need to project 1st slice and copy */
			{
				forwardProject2D(e[0],Image.image[0],A_Padded_Map,max_num_pointer,&sinogram.sinoparams,&Image.imgparams,svpar);
				for(jz=1;jz<Nz;jz++)
					memcpy(e[jz],e[0],NvNc*sizeof(float));
			}

			gettimeofday(&tm2,NULL);
			tdiff = 1000 * (tm2.tv_sec - tm1.tv_sec) + (tm2.tv_usec - tm1.tv_usec) / 1000;
			fprintf(stdout,"\tProjection time = %llu ms\n",tdiff);
		}
	}

	/***** Reconstruction mode *****/
	if(cmdline.reconFlag)
	{
		/* Allocate and Read sinogram data */
		AllocateSinoData3DParallel(&sinogram);
		ReadSinoData3DParallel(cmdline.SinoDataFile, &sinogram);
		if(cmdline.SinoWeightsFileFlag)
			ReadWeights3D(cmdline.SinoWeightsFile, &sinogram);
		else
			ComputeSinoWeights(sinogram,reconparams);

		/* Read Proximal map if necessary */
		if(reconparams.ReconType == MBIR_MODULAR_RECONTYPE_PandP)
		{
			ProxMap.imgparams.Nx = Image.imgparams.Nx;
			ProxMap.imgparams.Ny = Image.imgparams.Ny;
			ProxMap.imgparams.Nz = Image.imgparams.Nz;
			ProxMap.imgparams.FirstSliceNumber = Image.imgparams.FirstSliceNumber;
			ProxMap.imgparams.NumSliceDigits = Image.imgparams.NumSliceDigits;
			AllocateImageData3D(&ProxMap);
			ReadImage3D(cmdline.ProxMapImageFile,&ProxMap);
			reconparams.proximalmap = ProxMap.image;  // **ptr to proximal map image
		}

		/* Start Reconstruction */
		fprintf(stdout,"Reconstructing...\n");
		gettimeofday(&tm1,NULL);

		/* "e" will hold the sinogram error (y-Ax) during reconstruction */
		for(jz=0; jz<Nz; jz++)
		for(i=0; i<NvNc; i++)
			e[jz][i] = sinogram.sino[jz][i]-e[jz][i];

		MBIRReconstruct3D(&Image,&sinogram,e,reconparams,svpar,A_Padded_Map,max_num_pointer,ImageReconMask,&cmdline);

		gettimeofday(&tm2,NULL);
		tdiff = 1000 * (tm2.tv_sec - tm1.tv_sec) + (tm2.tv_usec - tm1.tv_usec) / 1000;
		fprintf(stdout,"\tReconstruction time = %llu ms\n",tdiff);

		/* Write out reconstructed image(s) */
		fprintf(stdout,"Writing image files...\n");
		WriteImage3D(cmdline.ReconImageFile, &Image);

		if(cmdline.writeProjectionFlag)  /* flip it back to get projection Ax */
		{
			for(jz=0; jz<Nz; jz++)
			for(i=0; i<NvNc; i++)
				e[jz][i] = sinogram.sino[jz][i]-e[jz][i];
		}

		FreeSinoData3DParallel(&sinogram);
		if(reconparams.ReconType == MBIR_MODULAR_RECONTYPE_PandP)
			FreeImageData3D(&ProxMap);

	}

	/* Write Projection of image state if requested */
	if(cmdline.writeProjectionFlag)
	{
		fprintf(stdout,"Writing projection of image state...\n");

		if(cmdline.writeProjectionFlag == 1)  // single slice case
		{
			sprintf(fname,"%s.2Dprojection",cmdline.outputProjectionFile);
			if( WriteFloatArray(fname,e[0],NvNc) ) {
				fprintf(stderr,"Error: can't open file %s for writing\n",fname);
				exit(-1);
			}
		}
		if(cmdline.writeProjectionFlag == 2)  // multi-slice case
		{
			for(jz=0; jz<Nz; jz++)
			{
				sprintf(fname,"%s_slice%.*d.2Dprojection",cmdline.outputProjectionFile,NumSliceDigits,jz+FirstSliceNumber);
				if( WriteFloatArray(fname,e[jz],NvNc) ) {
					fprintf(stderr,"Error: can't open file %s for writing\n",fname);
					exit(-1);
				}
			}
		}
	}

	if(cmdline.reconFlag || cmdline.writeProjectionFlag)
	{
		multifree(e,2);
		FreeImageData3D(&Image);
	}

	/* Free SV memory */
	for(j=0;j<svpar.Nsv;j++)  free((void *)svpar.bandMinMap[j].bandMin);
	for(j=0;j<svpar.Nsv;j++)  free((void *)svpar.bandMaxMap[j].bandMax);
	free((void *)svpar.bandMinMap);
	free((void *)svpar.bandMaxMap);

	/* Free system matrix */
	for(i=0;i<svpar.Nsv;i++)
	for(j=0;j<(2*svpar.SVLength+1)*(2*svpar.SVLength+1);j++)
	if(A_Padded_Map[i][j].length>0)
	{
		free((void *)A_Padded_Map[i][j].val);
		free((void *)A_Padded_Map[i][j].pieceWiseMin);
		free((void *)A_Padded_Map[i][j].pieceWiseWidth);
	}
	multifree(A_Padded_Map,2);
	free((void *)max_num_pointer);
	free((void *)ImageReconMask);

	fprintf(stdout,"Done.\n\n");
	return(0);
}


/* Read Command-line */
void readCmdLine(int argc, char *argv[], struct CmdLine *cmdline)
{
    char ch;
    
    /* set defaults */
    cmdline->SinoParamsFileFlag=0;
    cmdline->ImageParamsFileFlag=0;
    cmdline->ReconParamsFileFlag=0;
    cmdline->SinoDataFileFlag=0;
    cmdline->SinoWeightsFileFlag=0;
    cmdline->ReconImageFileFlag=0;
    cmdline->SysMatrixFileFlag=0;

    cmdline->reconFlag = MBIR_MODULAR_RECONTYPE_QGGMRF_3D;
    cmdline->readInitImageFlag=0;
    cmdline->readInitProjectionFlag=0;
    cmdline->writeProjectionFlag=0;

    /* Print usage statement if no arguments, or help argument given */
    if(argc==1 || CmdLineHelpOption(argv[1]))
    {
        //fprintf(stdout,"Printing usage statement for %s\n",argv[0]);
        printCmdLineUsage(argv[0]);
        exit(0);
    }
    
    /* get options */
    while ((ch = getopt(argc, argv, "i:j:k:s:w:r:m:t:e:f:p:v")) != EOF)
    {
        switch (ch)
        {
            case 'i':
            {
                cmdline->ImageParamsFileFlag=1;
                sprintf(cmdline->ImageParamsFile, "%s", optarg);
                break;
            }
            case 'j':
            {
                cmdline->SinoParamsFileFlag=1;
                sprintf(cmdline->SinoParamsFile, "%s", optarg);
                break;
            }
            case 'k':
            {
                cmdline->ReconParamsFileFlag=1;
                sprintf(cmdline->ReconParamsFile, "%s", optarg);
                break;
            }
            case 's':
            {
                cmdline->SinoDataFileFlag=1;
                sprintf(cmdline->SinoDataFile, "%s", optarg);
                break;
            }
            case 'w':
            {
                cmdline->SinoWeightsFileFlag=1;
                sprintf(cmdline->SinoWeightsFile, "%s", optarg);
                break;
            }
            case 'r':
            {
                cmdline->ReconImageFileFlag=1;
                sprintf(cmdline->ReconImageFile, "%s", optarg);
                break;
            }
            case 'm':
            {
                cmdline->SysMatrixFileFlag=1;
                sprintf(cmdline->SysMatrixFile, "%s", optarg);
                break;
            }
            case 't':
            {
                cmdline->readInitImageFlag=1;
                sprintf(cmdline->InitImageFile, "%s", optarg);
                break;
            }
            case 'e':
            {
                cmdline->readInitProjectionFlag=1;
                sprintf(cmdline->inputProjectionFile, "%s", optarg);
                break;
            }
            case 'f':
            {
                cmdline->writeProjectionFlag=1;
                sprintf(cmdline->outputProjectionFile, "%s", optarg);
                break;
            }
            case 'p':
            {
                cmdline->reconFlag = MBIR_MODULAR_RECONTYPE_PandP;
                sprintf(cmdline->ProxMapImageFile, "%s", optarg);
                break;
            }
            // Reserve this for verbose-mode flag
            case 'v':
            {
                break;
            }
            default:
            {
                //fprintf(stderr,"%s: invalid option '%c'\n",argv[0],ch);  //getopt does this already
                fprintf(stderr,"Try '%s -help' for more information.\n",argv[0]);
                exit(-1);
                break;
            }
        }
    }

    fprintf(stdout,"Parsing command line...\n");

    /* Check for mandatory arguments */
    if(!cmdline->SinoParamsFileFlag || !cmdline->ImageParamsFileFlag){
        fprintf(stderr,"Error: Either sinoparams or imgparams (-i,-j) file wasn't specified\n");
        fprintf(stderr,"Try '%s -help' for more information.\n",argv[0]);
        exit(-1);
    }

    /* Determine what to do based on supplied options */
    cmdline->readAmatrixFlag=0;
    cmdline->writeAmatrixFlag=0;

    if(cmdline->ReconImageFileFlag)  /* reconstruction mode */
    {
        //cmdline->reconFlag already set above
        if(cmdline->SysMatrixFileFlag)
            cmdline->readAmatrixFlag=1;

        if(cmdline->writeProjectionFlag)
            cmdline->writeProjectionFlag=2;  /* this means write full multi-slice projection */

    }
    else  /* precompute mode */
    {
        cmdline->reconFlag=0;
        cmdline->readInitProjectionFlag=0;

        if(cmdline->SysMatrixFileFlag)
            cmdline->writeAmatrixFlag=1;

        /* Case where we just want to apply projector to an input image */
        /* Here, if matrix file is given take it as an input (so it must be pre-computed in this case) */
        if(cmdline->writeProjectionFlag && cmdline->readInitImageFlag) {
            cmdline->writeProjectionFlag=2;  /* this means write full multi-slice projection */
            if(cmdline->SysMatrixFileFlag) {
                cmdline->writeAmatrixFlag=0;
                cmdline->readAmatrixFlag=1;
            }
        }
    }

    /* Print output and check errors of above parsing sequence  */
    if(cmdline->reconFlag)
    {
        fprintf(stdout,"-> will perform reconstruction ");
        if(cmdline->reconFlag == MBIR_MODULAR_RECONTYPE_QGGMRF_3D)
            fprintf(stdout,"(QGGMRF)\n");
        if(cmdline->reconFlag == MBIR_MODULAR_RECONTYPE_PandP)
            fprintf(stdout,"(Plug & Play)\n");

        if(cmdline->readAmatrixFlag)
            fprintf(stdout,"-> will read system matrix from file\n");
        else
        {
            fprintf(stdout,"-> will compute system matrix\n");
            fprintf(stdout," *** NOTE if you precompute/store the system matrix, any further reconstruction\n");
            fprintf(stdout," *** with the same image/sinogram in-slice dimensions will execute MUCH faster.\n");
            fprintf(stdout," *** See help (-m option)\n");
        //  fprintf(stdout,"***80 columns*******************************************************************\n\n");
        }
        if(!cmdline->SinoWeightsFileFlag)
            fprintf(stdout,"-> will compute sinogram weights internally (no file provided)\n");
        if(cmdline->readInitImageFlag)
            fprintf(stdout,"-> will read initial condition from file(s)\n");
        if(cmdline->readInitProjectionFlag)
            fprintf(stdout,"-> will read projection of initial condition\n");
        else
        {
            fprintf(stdout,"-> will compute forward projection of initial condition\n");
            fprintf(stdout," *** NOTE you may save run time by pre-computing the initial projection\n");
            fprintf(stdout," *** See help (-e option)\n");
        }

        if(cmdline->writeProjectionFlag)
            fprintf(stdout,"-> will save projection of output image state to file(s)\n");

        if(!cmdline->ReconParamsFileFlag || !cmdline->SinoDataFileFlag)
        {
            fprintf(stderr,"Error: Either input data or reconstruction parameters weren't specified\n");
            fprintf(stderr,"Try '%s -help' for more information.\n",argv[0]);
            exit(-1);
        }

    }
    else if(cmdline->writeAmatrixFlag || cmdline->writeProjectionFlag)
    {
        fprintf(stdout,"-> pre-compute mode..no reconstruction\n");
        if(cmdline->writeAmatrixFlag)
            fprintf(stdout,"-> will compute system matrix and write to file\n");
        if(cmdline->writeProjectionFlag)
            fprintf(stdout,"-> will compute projection of initial condition and write to file(s)\n");
        if(cmdline->ReconParamsFileFlag || cmdline->SinoDataFileFlag || cmdline->SinoWeightsFileFlag || cmdline->readInitProjectionFlag)
            fprintf(stdout,"Note some command line options are being ignored.\n");
    }
    else
    {
        fprintf(stderr,"Error: From the given command options, not sure what you want to do.\n");
        fprintf(stderr,"Try '%s -help' for more information.\n",argv[0]);
        exit(-1);
    }

    fprintf(stdout,"\n");

    fprintf(stdout,"Filenames provided:\n");

    if(cmdline->SinoParamsFileFlag)
        fprintf(stdout,"   Sino params   = %s.sinoparams\n",cmdline->SinoParamsFile);
    if(cmdline->ImageParamsFileFlag)
        fprintf(stdout,"   Image params  = %s.imgparams\n",cmdline->ImageParamsFile);
    if(cmdline->ReconParamsFileFlag)
        fprintf(stdout,"   Recon params  = %s.reconparams\n",cmdline->ReconParamsFile);
    if(cmdline->SinoDataFileFlag)
        fprintf(stdout,"   Sinogram data = %s_sliceNNN.2Dsinodata\n",cmdline->SinoDataFile);
    if(cmdline->SinoWeightsFileFlag)
        fprintf(stdout,"   Weight data   = %s_sliceNNN.2Dweightdata\n",cmdline->SinoWeightsFile);
    if(cmdline->ReconImageFileFlag)
        fprintf(stdout,"   Output images = %s_sliceNNN.2Dimgdata\n",cmdline->ReconImageFile);
    if(cmdline->readInitImageFlag)
        fprintf(stdout,"   Initial image = %s_sliceNNN.2Dimgdata\n",cmdline->InitImageFile);
    if(cmdline->SysMatrixFileFlag)
        fprintf(stdout,"   System matrix = %s.2Dsvmatrix\n",cmdline->SysMatrixFile);
    if(cmdline->readInitProjectionFlag)
        fprintf(stdout,"   Initial projection = %s.2Dprojection\n",cmdline->inputProjectionFile);
    if(cmdline->writeProjectionFlag)
    {
        if(cmdline->readInitImageFlag)
            fprintf(stdout,"   Output projection = %s_sliceNNN.2Dprojection\n",cmdline->outputProjectionFile);
        else
            fprintf(stdout,"   Output projection = %s.2Dprojection\n",cmdline->outputProjectionFile);
    }

}



int NumSliceDigits(char *basename, char *ext, int slice)
{
    FILE *fp;
    char fname[1024];
    int Ndigits = MBIR_MODULAR_MAX_NUMBER_OF_SLICE_DIGITS;

    while(Ndigits > 0)
    {
        sprintf(fname,"%s_slice%.*d.%s",basename, Ndigits, slice, ext);
        //printf("%s\n",fname);
        if( (fp=fopen(fname,"r")) ) {
            fclose(fp);
            break;
        }
        else
            Ndigits--;
    }
    return(Ndigits);
}


void setNumSliceDigits(
	char *basename,
	char *ext,
	int slice,
	struct SinoParams3DParallel *sinoparams,
	struct ImageParams3D *imgparams)
{
	int Ndigits;
	if( (Ndigits = NumSliceDigits(basename,ext,slice)) > 0 )
	{
		sinoparams->NumSliceDigits = Ndigits;
		imgparams->NumSliceDigits = Ndigits;
	}
	else
	{
		fprintf(stderr,"Error: Can't determine number of slice digits from given input file.\n");
		fprintf(stderr,"* Looking for file with this format: %s_slice%d.%s\n",basename,slice,ext);
		fprintf(stderr,"* where the slice number can contain leading zeros but no spaces.\n");
		exit(-1);
	}
}


void printBanner(void)
{
    fprintf(stdout,"MBIR RECONSTRUCTION FOR 3D PARALLEL-BEAM CT\n");
    fprintf(stdout,"build time: %s, %s\n\n", __DATE__,  __TIME__);
}


void printCmdLineUsage(char *ExecFileName)
{
//  fprintf(stdout,"***80 columns*******************************************************************\n\n");
    fprintf(stdout,"Command Line Help\n\n");
    fprintf(stdout,"There are two forms for the command line. The first computes the system matrix,\n");
    fprintf(stdout,"writes it to a file and exits. This matrix is required for the reconstruction\n");
    fprintf(stdout,"but is fixed for a given set of sinogram and image parameters so it saves a lot\n");
    fprintf(stdout,"of time in the reconstruction phase to precompute this. You can also choose to\n");
    fprintf(stdout,"pre-compute the initial projection of the default or input initial condition\n");
    fprintf(stdout,"which saves having to run the projection during the reconstruction phase.\n");
    fprintf(stdout,"Command Line Syntax: (printed on multiple lines for clarity)\n");
    fprintf(stdout,"\n");
    fprintf(stdout,"  %s\n",ExecFileName);
    fprintf(stdout,"\t-i <filename>[.imgparams]    : Input image parameters\n");
    fprintf(stdout,"\t-j <filename>[.sinoparams]   : Input sinogram parameters\n");
    fprintf(stdout,"\t-m <filename>[.2Dsvmatrix]  : Output matrix file\n");
    fprintf(stdout,"    (following are optional)\n");
//  fprintf(stdout,"***80 columns*******************************************************************\n\n");
    fprintf(stdout,"\t-f <baseFilename>            : Output projection of default or input IC\n");
    fprintf(stdout,"\t-t <baseFilename> (+ \"-f\")   : Input initial image(s)\n");
    fprintf(stdout,"\n");
    fprintf(stdout,"In the above arguments, the exensions given in the '[]' symbols must be part of\n");
    fprintf(stdout,"the file names but should be omitted from the command line.\n");
    fprintf(stdout,"\n");
    fprintf(stdout,"The second form is for the actual reconstruction:\n");
    fprintf(stdout,"\n");
    fprintf(stdout,"  %s\n",ExecFileName);
    fprintf(stdout,"\t-i <filename>[.imgparams]    : Input image parameters\n");
    fprintf(stdout,"\t-j <filename>[.sinoparams]   : Input sinogram parameters\n");
    fprintf(stdout,"\t-k <filename>[.reconparams]  : Reconstruction parameters\n");
    fprintf(stdout,"\t-s <baseFilename>            : Input sinogram projection file(s)\n");
    fprintf(stdout,"\t-w <baseFilename>            : Input sinogram weight file(s)\n");
    fprintf(stdout,"\t-r <baseFilename>            : Output reconstruced image file(s)\n");
    fprintf(stdout,"    (following are optional)\n");
    fprintf(stdout,"\t-p <baseFilename>            : Proximal map image(s) for Plug & Play\n");
    fprintf(stdout,"\t                             : ** -p specifies to use proximal prior\n");
    fprintf(stdout,"\t                             : ** generally use with -t -e -f\n");
    fprintf(stdout,"\t-m <filename>[.2Dsvmatrix]  : INPUT matrix file (params must match!)\n");
    fprintf(stdout,"\t-t <baseFilename>            : Input initial condition image(s)\n");
    fprintf(stdout,"\t-e <baseFilename>            : Input projection of initial condition\n");
    fprintf(stdout,"\t                             : ** default IC if -t not specified\n");
    fprintf(stdout,"\t-f <baseFilename>            : Output projection of final image state\n");
    fprintf(stdout,"\t-v                           : verbose mode (TBD)\n");
    fprintf(stdout,"\n");
    fprintf(stdout,"For all the arguments specifying <baseFilename>, the relevant 3D data is split\n");
    fprintf(stdout,"across files, one file per slice. The file naming convention is as follows,\n");
    fprintf(stdout,"depending on the data contents:\n");
    fprintf(stdout,"\n");
    fprintf(stdout,"\t<baseFilename>_slice<sliceIndex>.2Dimgdata\n");
    fprintf(stdout,"\t<baseFilename>_slice<sliceIndex>.2Dsinodata\n");
    fprintf(stdout,"\t<baseFilename>_slice<sliceIndex>.2Dweightdata\n");
    fprintf(stdout,"\t<baseFilename>_slice<sliceIndex>.2Dprojection\n");
    fprintf(stdout,"\n");
    fprintf(stdout,"where <sliceIndex> (skip '<>' symbols) is a non-negative integer including\n");
    fprintf(stdout,"leading zeros and no spaces (e.g. 0000 to 1023). The number of digits\n");
    fprintf(stdout,"is flexible (up to %d) but must be consistent.\n",MBIR_MODULAR_MAX_NUMBER_OF_SLICE_DIGITS);
    fprintf(stdout,"\n");
}


int CmdLineHelpOption(char *string)
{
    if( (strcmp(string,"-h")==0) || (strcmp(string,"-help")==0) || (strcmp(string,"--help")==0) || (strcmp(string,"help")==0) )
        return 1;
    else
        return 0;
}



