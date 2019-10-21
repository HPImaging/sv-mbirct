
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
void PrintCmdLineUsage(char *ExecFileName);
int CmdLineHelp(char *string);


int main(int argc, char *argv[])
{
	struct CmdLine cmdline;

	struct Image3D Image;
	struct Sino3DParallel sinogram;
	struct ReconParamsQGGMRF3D reconparams;
	struct SVParams svpar;
	struct AValues_char **A_Padded_Map; 
	float *max_num_pointer;	
	struct timeval tm1,tm2;

	char *ImageReconMask;	/* Image reconstruction mask (determined by ROI) */
	float InitValue=MUWATER;/* If initial image not specified, set to uniform image w/ */
	float OutsideROIValue=0;/* value InitValue inside ROI radius and OutsideROIValue outside */
	char fname[200];

	fprintf(stdout, "Starting Reconstruction...\n\n");

	/* Parse command line and parameter files */
	readCmdLine(argc, argv, &cmdline);
	readSystemParams(&cmdline, &Image.imgparams, &sinogram.sinoparams, &reconparams);
	fprintf(stdout,"INPUT ");
	printSinoParams3DParallel(&sinogram.sinoparams);
	fprintf(stdout,"OUTPUT ");
	printImageParams3D(&Image.imgparams);
	printReconParamsQGGMRF3D(&reconparams);
	fprintf(stdout,"\n");

	/* The image parameters specify the relevant slice range to reconstruct, so re-set the  */
	/* relevant sinogram parameters so it pulls the correct slices and indexes consistently */
	sinogram.sinoparams.NSlices = Image.imgparams.Nz;
	sinogram.sinoparams.FirstSliceNumber = Image.imgparams.FirstSliceNumber;

	/* Read Sinogram and Weights */
	AllocateSinoData3DParallel(&sinogram);
	ReadSinoData3DParallel(cmdline.SinoDataFile, &sinogram);
	ReadWeights3D(cmdline.SinoWeightsFile, &sinogram);

	/* Allocate memory for image */
	AllocateImageData3D(&Image);

	/* Allocate and generate recon mask based on ROIRadius */
	/* Initialize image and reconstruction mask */
	/* --Careful..the initial forward projection is written out in GenSysMatrix, so the inital image has to match */
	/* --Pixels outside ROI radius never get touched so the these values will pass */
	/*     through to the output and will NOT be included in the projection */
	ImageReconMask = GenImageReconMask(&(Image.imgparams));
	Initialize_Image(&Image, &cmdline, ImageReconMask, InitValue, OutsideROIValue);

	/* Initialize SV parameters */
	initSVParams(&svpar, Image.imgparams, sinogram.sinoparams);
	svpar.bandMinMap = (struct minStruct *)get_spc(svpar.Nsv,sizeof(struct minStruct));
	svpar.bandMaxMap = (struct maxStruct *)get_spc(svpar.Nsv,sizeof(struct maxStruct));

	/* Read System Matrix */
	A_Padded_Map = (struct AValues_char **)multialloc(sizeof(struct AValues_char), 2, svpar.Nsv, (svpar.SVLength*2+1)*(svpar.SVLength*2+1));
	max_num_pointer = (float *)malloc(Image.imgparams.Ny*Image.imgparams.Nx*sizeof(float));	
	sprintf(fname,"%s.2Dsysmatrix",cmdline.SysMatrixFile);
	//readAmatrix(fname,A_Padded_Map,max_num_pointer,&Image.imgparams,&sinogram.sinoparams,svpar.Nsv,svpar.bandMinMap,svpar.bandMaxMap,svpar.pieceLength);
	readAmatrix(fname, A_Padded_Map, max_num_pointer, &Image.imgparams, &sinogram.sinoparams, svpar);

	gettimeofday(&tm1,NULL);

	MBIRReconstruct3D(&Image,&sinogram,reconparams,svpar,A_Padded_Map,max_num_pointer,&cmdline);

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
	free((void *)svpar.bandMinMap);
	free((void *)svpar.bandMaxMap);
    
	return(0);
}





/* Read Command-line */
void readCmdLine(int argc, char *argv[], struct CmdLine *cmdline)
{
    char ch;
    
    strcpy(cmdline->InitImageDataFile, "NA"); /* default */
    
    if(argc<15)
    {
        if(argc==2 && CmdLineHelp(argv[1]))
        {
            fprintf(stdout,"\n=========HELP==========\n");
            PrintCmdLineUsage(argv[0]);
            exit(-1);
        }
        else
        {
         fprintf(stderr, "\nError : Improper Command line for exec-program %s, Number of arguments lower than needed \n",argv[0]);
         PrintCmdLineUsage(argv[0]);
         exit(-1);
        }
    }
    
    /* get options */
    while ((ch = getopt(argc, argv, "i:j:k:m:s:w:r:t:v")) != EOF)
    {
        switch (ch)
        {
            case 'i':
            {
                sprintf(cmdline->ImageParamsFile, "%s", optarg);
                break;
            }
            case 'j':
            {
                sprintf(cmdline->SinoParamsFile, "%s", optarg);
                break;
            }
            case 'k':
            {
                sprintf(cmdline->ReconParamsFile, "%s", optarg);
                break;
            }
            case 'm':
            {
                sprintf(cmdline->SysMatrixFile, "%s", optarg);
                break;
            }
            case 's':
            {
                sprintf(cmdline->SinoDataFile, "%s", optarg);
                break;
            }
            case 'w':
            {
                sprintf(cmdline->SinoWeightsFile, "%s", optarg);
                break;
            }
            case 'r':
            {
                sprintf(cmdline->ReconImageDataFile, "%s", optarg);
                break;
            }
            case 't':
            {
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
                printf("\nError : Command line Symbol not recongized\n");
                PrintCmdLineUsage(argv[0]);
                exit(-1);
                break;
            }
        }
    }

}

void PrintCmdLineUsage(char *ExecFileName)
{
    fprintf(stdout, "\nBASELINE MBIR RECONSTRUCTION SOFTWARE FOR 3D PARALLEL-BEAM  CT \n");
    fprintf(stdout, "build time: %s, %s\n", __DATE__,  __TIME__);
    fprintf(stdout, "\nCommand line Format for Executable File %s : \n", ExecFileName);
    fprintf(stdout, "%s -i <InputFileName>[.imgparams] -j <InputFileName>[.sinoparams]\n",ExecFileName);
    fprintf(stdout, "   -k <InputFileName>[.reconparams] -m <InputFileName>[.2Dsysmatrix]\n");
    fprintf(stdout, "   -s <InputProjectionsBaseFileName> -w <InputWeightsBaseFileName>\n");
    fprintf(stdout, "   -r <OutputImageBaseFileName>\n\n");
    fprintf(stdout, "Additional option to read in initial image: -t <InitialImageBaseFileName> \n\n");
    fprintf(stdout, "Note : The necessary extensions for certain input files are mentioned above within\n");
    fprintf(stdout, "a \"[]\" symbol above, however the extensions should be OMITTED in the command line\n\n");
    fprintf(stdout, "The following instructions pertain to the -s, -w and -r options: \n");
    fprintf(stdout, "A) The Sinogram Projection data files should be stored slice by slice in a single\n");
    fprintf(stdout, "   directory. All files within this directory must share a common BaseFileName and\n");
    fprintf(stdout, "   adhere to the following format :\n");
    fprintf(stdout, "      <ProjectionsBaseFileName>_slice<SliceIndex>.2Dsinodata \n");
    fprintf(stdout, "   where \"SliceIndex\" is a non-negative integer indexing each slice and is printed\n");
    fprintf(stdout, "   with 4 digits. Eg : 0000 to 9999 is a valid descriptor for \"SliceIndex\" \n");
    fprintf(stdout, "B) Similarly, the format for the Sinogram Weights files is :\n");
    fprintf(stdout, "      <WeightsBaseFileName>_slice<SliceIndex>.2Dweightdata \n");
    fprintf(stdout, "C) Similarly, the Reconstructed (Output) Image is organized slice by slice :\n");
    fprintf(stdout, "      <ImageBaseFileName>_slice<SliceIndex>.2Dimgdata\n\n");
}


int CmdLineHelp(char *string)
{
    if( (strcmp(string,"-h")==0) || (strcmp(string,"-help")==0) || (strcmp(string,"--help")==0) || (strcmp(string,"help")==0) )
        return 1;
    else
        return 0;
}







