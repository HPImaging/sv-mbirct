
#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <string.h>
#include <time.h>
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
void setNumSliceDigits(struct SinoParams3DParallel *sinoparams, struct ImageParams3D *imgparams, struct CmdLine *cmdline, char *cmd);


int main(int argc, char *argv[])
{
	struct CmdLine cmdline;
	struct Image3D Image;
	struct Sino3DParallel sinogram;
	struct ReconParamsQGGMRF3D reconparams;
	struct SVParams svpar;

	struct AValues_char **A_Padded_Map; 
	float *max_num_pointer;	

	char *ImageReconMask;	/* Image reconstruction mask (determined by ROI) */
	float InitValue=MUWATER;/* If initial image not specified, set to uniform image w/ */
	float OutsideROIValue=0;/* value InitValue inside ROI radius and OutsideROIValue outside */
	char fname[200];
	struct timeval tm1,tm2;

	fprintf(stdout,"MBIR RECONSTRUCTION FOR 3D PARALLEL-BEAM CT\n");
	fprintf(stdout,"build time: %s, %s\n\n", __DATE__,  __TIME__);

	readCmdLine(argc, argv, &cmdline);

	/* Read image/sino parameter files */
	ReadImageParams3D(cmdline.ImageParamsFile,&Image.imgparams);
	ReadSinoParams3DParallel(cmdline.SinoParamsFile,&sinogram.sinoparams);
	fprintf(stdout,"INPUT ");
	printSinoParams3DParallel(&sinogram.sinoparams);
	fprintf(stdout,"OUTPUT ");
	printImageParams3D(&Image.imgparams);
	if(cmdline.reconFlag)
	{
		ReadReconParamsQGGMRF3D(cmdline.ReconParamsFile,&reconparams);
		printReconParamsQGGMRF3D(&reconparams);
		NormalizePriorWeights3D(&reconparams);
	}
	fprintf(stdout,"\n");

	int NvNc = sinogram.sinoparams.NViews * sinogram.sinoparams.NChannels;
	int NxNy = Image.imgparams.Nx * Image.imgparams.Ny;
	int Nz = Image.imgparams.Nz;

	/* Initialize SV parameters */
	initSVParams(&svpar, Image.imgparams, sinogram.sinoparams);
	svpar.bandMinMap = (struct minStruct *)get_spc(svpar.Nsv,sizeof(struct minStruct));
	svpar.bandMaxMap = (struct maxStruct *)get_spc(svpar.Nsv,sizeof(struct maxStruct));

	/* Allocate and generate recon mask based on ROIRadius */
	ImageReconMask = GenImageReconMask(&Image.imgparams);

	/* Set up System Matrix */
	A_Padded_Map = (struct AValues_char **)multialloc(sizeof(struct AValues_char), 2, svpar.Nsv, (svpar.SVLength*2+1)*(svpar.SVLength*2+1));
	max_num_pointer = (float *) get_spc(NxNy,sizeof(float));

	if(cmdline.readAmatrixFlag)
	{
		fprintf(stdout,"Reading system matrix...\n");
		sprintf(fname,"%s.2Dsysmatrix",cmdline.SysMatrixFile);
		readAmatrix(fname, A_Padded_Map, max_num_pointer, &Image.imgparams, &sinogram.sinoparams, svpar);
	}
	else
	{
		fprintf(stdout,"Computing system matrix...\n");
	        A_comp(A_Padded_Map,max_num_pointer,svpar,&sinogram.sinoparams,ImageReconMask,&Image.imgparams);
        }

	/***** Precompute mode *****/
	if(cmdline.precompAmatrixFlag)
	{
		fprintf(stdout,"Writing system matrix...\n");
		sprintf(fname,"%s.2Dsysmatrix",cmdline.SysMatrixFile);
		remove(fname);
		writeAmatrix(fname,A_Padded_Map,max_num_pointer,&Image.imgparams,&sinogram.sinoparams,svpar);

		if(cmdline.precompInitProjFlag)
		{
			/* compute initial projection, assuming constant initial image */
			float *initProjection = (float *) get_spc(NvNc,sizeof(float));
			float *x = (float *) get_spc(NxNy,sizeof(float));
			int i;
			for(i=0; i<NxNy; i++) 
				x[i]=InitValue;
			forwardProject2D(initProjection, x, A_Padded_Map, max_num_pointer, &sinogram.sinoparams, &Image.imgparams, svpar);

			fprintf(stdout,"Writing initial projection...\n");
			sprintf(fname,"%s.initProj",cmdline.InitProjFile);
			if( WriteFloatArray(fname,initProjection,NvNc) ) {
				fprintf(stderr,"ERROR: can't open file %s for writing\n",fname);
				exit(-1);
			}
			free((void *)initProjection);
			free((void *)x);
		}
	}

	/***** Reconstruction mode *****/
	if(cmdline.reconFlag)
	{
		/* The image parameters specify the relevant slice range to reconstruct, so re-set the  */
		/* relevant sinogram parameters so it pulls the correct slices and indexes consistently */
		sinogram.sinoparams.NSlices = Image.imgparams.Nz;
		sinogram.sinoparams.FirstSliceNumber = Image.imgparams.FirstSliceNumber;

		/* Determine and SET number of slice index digits for input data and output images */
		setNumSliceDigits(&sinogram.sinoparams, &Image.imgparams, &cmdline, argv[0]);

		/* Allocate memory for data and image volume */
		AllocateSinoData3DParallel(&sinogram);
		AllocateImageData3D(&Image);
		float **e = (float **)multialloc(sizeof(float), 2, sinogram.sinoparams.NSlices,NvNc);     /* projection error */

		/* Read sinogram data */
		ReadSinoData3DParallel(cmdline.SinoDataFile, &sinogram);
		ReadWeights3D(cmdline.SinoWeightsFile, &sinogram);

		/* Initialize image */
		/** Note the pixels outside the ROI radius never get touched by the projector or updates */
		/** so these values will pass through to the output */
		initImage(&Image, &cmdline, ImageReconMask, InitValue, OutsideROIValue);
		initProjectionError(e,&Image,&sinogram,svpar,A_Padded_Map,max_num_pointer,&cmdline);

		fprintf(stdout,"Reconstructing...\n");
		gettimeofday(&tm1,NULL);
		MBIRReconstruct3D(&Image,&sinogram,e,reconparams,svpar,A_Padded_Map,max_num_pointer,&cmdline);
		gettimeofday(&tm2,NULL);

		unsigned long long tt = 1000 * (tm2.tv_sec - tm1.tv_sec) + (tm2.tv_usec - tm1.tv_usec) / 1000;
		printf("\trun time %llu ms \n", tt);

		/* Write out reconstructed image */
		fprintf(stdout,"Writing image files...\n");
		WriteImage3D(cmdline.ReconImageDataFile, &Image);
		multifree(e,2);
		FreeImageData3D(&Image);
		FreeSinoData3DParallel(&sinogram);
	}

	free((void *)ImageReconMask);
	free((void *)svpar.bandMinMap);
	free((void *)svpar.bandMaxMap);
    
	fprintf(stdout,"Done.\n");
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
    cmdline->ReconImageDataFileFlag=0;
    cmdline->SysMatrixFileFlag=0;
    cmdline->InitImageDataFileFlag=0;
    cmdline->InitProjFileFlag=0;

    /* Print usage statement if no arguments, or help argument given */
    if(argc==1 || CmdLineHelpOption(argv[1]))
    {
        //fprintf(stdout,"Printing usage statement for %s\n",argv[0]);
        printCmdLineUsage(argv[0]);
        exit(0);
    }
    
    /* get options */
    while ((ch = getopt(argc, argv, "i:j:k:m:e:s:w:r:t:v")) != EOF)
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
            case 'm':
            {
                cmdline->SysMatrixFileFlag=1;
                sprintf(cmdline->SysMatrixFile, "%s", optarg);
                break;
            }
            case 'e':
            {
                cmdline->InitProjFileFlag=1;
                sprintf(cmdline->InitProjFile, "%s", optarg);
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
                cmdline->ReconImageDataFileFlag=1;
                sprintf(cmdline->ReconImageDataFile, "%s", optarg);
                break;
            }
            case 't':
            {
                cmdline->InitImageDataFileFlag=1;
                sprintf(cmdline->InitImageDataFile, "%s", optarg);
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
    if(!cmdline->SinoParamsFileFlag){
        fprintf(stderr,"OOPS: No [].sinoparams file specified\n");
        fprintf(stderr,"Try '%s -help' for more information.\n",argv[0]);
        exit(-1);
    }
    if(!cmdline->ImageParamsFileFlag){
        fprintf(stderr,"OOPS: No [].imgparams file specified\n");
        fprintf(stderr,"Try '%s -help' for more information.\n",argv[0]);
        exit(-1);
    }

    /* Determine what to do based on supplied options */
    cmdline->precompAmatrixFlag=0;
    cmdline->precompInitProjFlag=0;
    cmdline->reconFlag=0;
    cmdline->readAmatrixFlag=0;
    cmdline->readInitProjFlag=0;

    //If output images specified assume we want to reconstruct
    if(cmdline->ReconImageDataFileFlag)
    {
        cmdline->reconFlag=1;
        if(cmdline->SysMatrixFileFlag)
            cmdline->readAmatrixFlag=1;

        if(cmdline->InitProjFileFlag) 
            cmdline->readInitProjFlag=1;

        if(cmdline->InitImageDataFileFlag)  /* if initial condition supplied, ignore initial projection file */
            cmdline->readInitProjFlag=0;
    }
    //If no output images specified assume we just want to compute/store forward matrix
    else if(cmdline->SysMatrixFileFlag)
    {
        cmdline->precompAmatrixFlag=1;
        if(cmdline->InitProjFileFlag)
            cmdline->precompInitProjFlag=1;
    }
    else
    {
        fprintf(stderr,"%s: From the given command options, not sure what you want to do.\n",argv[0]);
        fprintf(stderr,"Try '%s -help' for more information.\n",argv[0]);
        exit(-1);
    }

    //Print output and check errors of above parsing sequence 
    if(cmdline->reconFlag)
    {
        fprintf(stdout,"-> will perform reconstruction\n");
        if(!cmdline->ReconParamsFileFlag || !cmdline->SinoDataFileFlag || !cmdline->SinoWeightsFileFlag)
        {
            fprintf(stderr,"OOPS: Output images specified, but input data or parameters are missing.\n");
            fprintf(stderr,"Try '%s -help' for more information.\n",argv[0]);
            exit(-1);
        }
        if(cmdline->readAmatrixFlag)
            fprintf(stdout,"-> will read precomputed system matrix\n");
        else
        {
            fprintf(stdout,"-> will compute system matrix\n");
            fprintf(stdout," * NOTE if you precompute/store the forward matrix, any further reconstruction\n");
            fprintf(stdout," * with the same image/sinogram in-slice dimensions will execute MUCH faster.\n");
            fprintf(stdout," * See help (-m option)\n");
        }
        if(cmdline->readInitProjFlag)
            fprintf(stdout,"-> will read projection of default initial condition\n");
        else
        {
            fprintf(stdout,"-> will forward project the default initial condition\n");
            fprintf(stdout," * NOTE if using the default initial condition (no -t option), run time can be\n");
            fprintf(stdout," * reduced by precomputing the initial projection. See help (-e option)\n");
        }
    }
    else if(cmdline->precompAmatrixFlag)
    {
        fprintf(stdout,"-> will compute/write system matrix only\n");
        if(cmdline->precompInitProjFlag)
            fprintf(stdout,"-> will compute/write projection of default initial condition\n");
        if(cmdline->ReconParamsFileFlag || cmdline->SinoDataFileFlag || cmdline->SinoWeightsFileFlag)
            fprintf(stdout,"Note some command line options are being ignored.\n");
    }

    fprintf(stdout,"\n");

}


void setNumSliceDigits(
	struct SinoParams3DParallel *sinoparams,
	struct ImageParams3D *imgparams,
	struct CmdLine *cmdline,
	char *cmd)
{
	int i,Ndigits;
	if( (Ndigits = NumSinoSliceDigits(cmdline->SinoDataFile, sinoparams->FirstSliceNumber))==0 )
	{
		fprintf(stderr,"Error: Can't read the first data file. Looking for any one of the following:\n");
		for(i=MBIR_MODULAR_MAX_NUMBER_OF_SLICE_DIGITS; i>0; i--)
			fprintf(stderr,"\t%s_slice%.*d.2Dsinodata\n",cmdline->SinoDataFile, i, sinoparams->FirstSliceNumber);
		fprintf(stderr,"Try '%s -help' for more information.\n",cmd);
		exit(-1);
	}
	sinoparams->NumSliceDigits = Ndigits;
	imgparams->NumSliceDigits = Ndigits;
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
    fprintf(stdout,"pre-compute the initial projection of the default initial condition, which saves\n");
    fprintf(stdout,"having to run the projection during the reconstruction phase.\n");
    fprintf(stdout,"Usage (printed on multiple lines for clarity):\n");
    fprintf(stdout,"\n");
    fprintf(stdout,"  %s\n",ExecFileName);
    fprintf(stdout,"\t-i <filename>[.imgparams]    \\\\Input image parameters\n");
    fprintf(stdout,"\t-j <filename>[.sinoparams]   \\\\Input sinogram parameters\n");
    fprintf(stdout,"\t-m <filename>[.2Dsysmatrix]  \\\\Output matrix file\n");
    fprintf(stdout,"    (following is optional)\n");
    fprintf(stdout,"\t-e <filename>[.initProj]     \\\\Projection of constant initial condition\n");
    fprintf(stdout,"\n");
    fprintf(stdout,"In the above arguments, the exensions given in the '[]' symbols must be part of\n");
    fprintf(stdout,"the file names but should be omitted from the command line.\n");
    fprintf(stdout,"\n");
    fprintf(stdout,"The second form is for the actual reconstruction:\n");
    fprintf(stdout,"\n");
    fprintf(stdout,"  %s\n",ExecFileName);
    fprintf(stdout,"\t-i <filename>[.imgparams]    \\\\Input image parameters\n");
    fprintf(stdout,"\t-j <filename>[.sinoparams]   \\\\Input sinogram parameters\n");
    fprintf(stdout,"\t-k <filename>[.reconparams]  \\\\Reconstruction parameters\n");
    fprintf(stdout,"\t-s <baseFilename>            \\\\Input sinogram projection file(s)\n");
    fprintf(stdout,"\t-w <baseFilename>            \\\\Input sinogram weight file(s)\n");
    fprintf(stdout,"\t-r <baseFilename>            \\\\Output reconstruced image file(s)\n");
    fprintf(stdout,"    (following are optional)\n");
    fprintf(stdout,"\t-m <filename>[.2Dsysmatrix]  \\\\INPUT matrix file (params must match!)\n");
    fprintf(stdout,"\t-e <filename>[.initProj]     \\\\Projection of constant initial condition\n");
    fprintf(stdout,"\t-t <baseFilename>            \\\\Input initial image(s)\n");
    fprintf(stdout,"\t-v                           \\\\verbose mode (TBD)\n");
    fprintf(stdout,"\n");
    fprintf(stdout,"For all the arguments specifying <baseFilename>, the relevant 3D data is split\n");
    fprintf(stdout,"across files, one file per slice. The file naming convention is as follows,\n");
    fprintf(stdout,"depending on the data contents:\n");
    fprintf(stdout,"\n");
    fprintf(stdout,"\t<baseFilename>_slice<sliceIndex>.2Dimgdata\n");
    fprintf(stdout,"\t<baseFilename>_slice<sliceIndex>.2Dsinodata\n");
    fprintf(stdout,"\t<baseFilename>_slice<sliceIndex>.2Dweightdata\n");
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



