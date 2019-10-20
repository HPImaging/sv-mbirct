
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <string.h>
#include <getopt.h>

#include "mbir_ct.h"
#include "MBIRModularDefs.h"
#include "MBIRModularUtils.h"
#include "allocate.h"
#include "recon3d.h"
#include "A_comp.h"


/* Command Line structure for Generating System matrix */
struct CmdLineSysGen
{
    char imgparamsFileName[256];  /* input file */
    char sinoparamsFileName[256]; /* input file */
    char SysMatrixFileName[256]; /* output file */
};

/* Internal fns */
void readCmdLineSysGen(int argc, char *argv[], struct CmdLineSysGen *cmdline);
void PrintCmdLineUsage(char *ExecFileName);
int CmdLineHelpCheck(char *string);


int main(int argc, char *argv[])
{
	struct CmdLineSysGen cmdline;
	struct ImageParams3D imgparams;
	struct SinoParams3DParallel sinoparams;
	//struct ReconParamsQGGMRF3D reconparams;
	float **PixelDetector_profile;
	struct minStruct *bandMinMap;
	struct maxStruct *bandMaxMap;
	struct AValues_char ** A_Padded_Map;
	int *order;
	float x_0, y_0, Deltaxy, x, y, yy, ROIRadius, R_sq, R_sq_max;
	int jx, jy,Nxy;
	/*struct SysMatrix2D *A ;*/

	/* read command line and parameter files */
	readCmdLineSysGen(argc, argv, &cmdline);
	ReadSinoParams3DParallel(cmdline.sinoparamsFileName, &sinoparams);
	ReadImageParams3D(cmdline.imgparamsFileName, &imgparams);
	printSinoParams3DParallel(&sinoparams);
	printImageParams3D(&imgparams);
	fprintf(stdout, "\n");

	unsigned int sum=0;
	int i,j,p,t;
	int Ny=imgparams.Ny;
	int Nx=imgparams.Nx;
	int NViews=sinoparams.NViews;
	int NChannels=sinoparams.NChannels;
	int NvNc = NViews*NChannels;
	int SVLength=SVLENGTH;
	int overlappingDistance=OVERLAPPINGDISTANCE;

	fprintf(stdout, "Generating System Matrix...\n\n");

	int pieceLength=computePieceLength(NViews);

	if(NViews%pieceLength!=0){
		fprintf(stderr, "Error: NViews mod pieceLength must be 0\n");
		fprintf(stderr, "Exiting %s\n",argv[0]);
		exit(-1);
        }        

	for(i=0;i<Ny;i+=(SVLength*2-overlappingDistance))
	for(j=0;j<Nx;j+=(SVLength*2-overlappingDistance))
		sum++;

	//fprintf(stdout, "Ny is %d Nx %d sum %d channels %d views %d\n",Ny,Nx,sum,sinoparams.NChannels,sinoparams.NViews);

	order = (int *)_mm_malloc(sum*sizeof(int),64);
	
	t=0;	
	
	for(i=0;i<Ny;i+=(SVLength*2-overlappingDistance))
	for(j=0;j<Nx;j+=(SVLength*2-overlappingDistance)){
		order[t]=i*Nx+j;  /* order is the first voxel coordinate, not the center */
		t++;
	}	

    	bandMinMap = (struct minStruct *)get_spc(sum,sizeof(struct minStruct));
    	bandMaxMap = (struct maxStruct *)get_spc(sum,sizeof(struct maxStruct));
    	A_Padded_Map = (struct AValues_char **)multialloc(sizeof(struct AValues_char), 2, sum, (SVLength*2+1)*(SVLength*2+1));	

	float* max_num_pointer = (float *)malloc(Ny*Nx*sizeof(float));	
	char **ImageReconMask = (char **)multialloc(sizeof(char), 2, Ny, Nx); 
    	PixelDetector_profile = ComputePixelProfile3DParallel(&sinoparams, &imgparams);  /* pixel-detector profile function */

	//fprintf(stdout, "after allocation \n");

    	Deltaxy = imgparams.Deltaxy;
    	ROIRadius = imgparams.ROIRadius;
    
    	x_0 = -(Nx-1)*Deltaxy/2;
    	y_0 = -(Ny-1)*Deltaxy/2;
    	Nxy = Nx*Ny;
    
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
                		{
                    			ImageReconMask[jy][jx] = 0;
                    		}	
                		else
                		{
                    			ImageReconMask[jy][jx] = 1;
                		}
            		}
        	}
    	}
	
	A_comp(bandMinMap,bandMaxMap,A_Padded_Map,max_num_pointer,&sinoparams,sum,ImageReconMask,order,&imgparams,PixelDetector_profile,cmdline.SysMatrixFileName,pieceLength);
	    
	fprintf(stdout, "Done generating system matrix\n");


	/***************************************************************/
	/* COMPUTE INTIAL ERROR IMAGE, ASSUMING CONSTANT INITIAL IMAGE */
	/***************************************************************/
	fprintf(stdout, "Computing projection of initial image\n");

	//float InitValue = reconparams.InitImageValue;
	float InitValue = MUWATER;   // careful here..this has to match initial value in reconstruction
	float *initialError = (float *)malloc(sizeof(float)*NvNc);
	forwardProject2D(initialError, InitValue, max_num_pointer,A_Padded_Map,bandMinMap, &sinoparams, &imgparams, pieceLength);
    
	char fname[200];
	sprintf(fname,"%s.initialError",cmdline.SysMatrixFileName);
	int exitcode;
	if( (exitcode=WriteFloatArray(fname,initialError,NvNc)) ) {
		if(exitcode==1) fprintf(stderr, "ERROR in Gen_SysMatrix: can't open file %s for writing\n",fname);
		if(exitcode==2) fprintf(stderr, "ERROR in Gen_SysMatrix: write to file %s terminated early\n",fname);
		exit(-1);
	}
	free((void *)initialError);

	fprintf(stdout, "Done computing initial projection\n");


	_mm_free(order);
	for(j=0;j<sum;j++){
		free((void *)bandMinMap[j].bandMin);
		free((void *)bandMaxMap[j].bandMax);
	}
	//fprintf(stdout, "after bandMin and bandMax \n");

	free((void *)bandMinMap);
	free((void *)bandMaxMap);
	for(i=0;i<sum;i++){
		for(j=0;j<((2*SVLength+1)*(2*SVLength+1));j++){
			if(A_Padded_Map[i][j].length>0){
				free((void *)A_Padded_Map[i][j].val);
				free((void *)A_Padded_Map[i][j].pieceWiseMin);
				free((void *)A_Padded_Map[i][j].pieceWiseWidth);
			}
		}
	}
	//fprintf(stdout, "A_Padded_Map\n");
	multifree(A_Padded_Map,2);
	multifree(ImageReconMask,2);
	free((void *)max_num_pointer);

	return 0;
}



void readCmdLineSysGen(
    int argc,
    char *argv[],
    struct CmdLineSysGen *cmdline)
{
    char ch;
    
    if(argc<7)
    {
        if(argc==2 && CmdLineHelpCheck(argv[1]))
        {
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

    //fprintf(stdout,"print argc is %d\n",argc);
 
    /* get options */
    while ((ch = getopt(argc, argv, "i:j:m:")) != EOF)
    {
        switch (ch)
        {
            case 'i':
            {
                sprintf(cmdline->imgparamsFileName, "%s", optarg);
		//fprintf(stdout,"image param file %s optarg %s\n",cmdline->imgparamsFileName,optarg);
                break;
            }
            case 'j':
            {
                sprintf(cmdline->sinoparamsFileName, "%s", optarg);
		//fprintf(stdout,"sino param file %s optarg %s\n",cmdline->sinoparamsFileName,optarg);
                break;
            }
            case 'm':
            {
                sprintf(cmdline->SysMatrixFileName, "%s", optarg);
		//fprintf(stdout,"Sys param file %s optarg %s\n",cmdline->SysMatrixFileName,optarg);
                break;
            }
            default:
            {
                fprintf(stderr, "\nError : Unrecognized Command-line Symbol for exec-program %s \n",argv[0]);
                PrintCmdLineUsage(argv[0]);
                exit(-1);
                break;
            }
        }
    }
}



void PrintCmdLineUsage(char *ExecFileName)
{
    fprintf(stdout, "\nBASELINE MBIR RECONSTRUCTION SOFTWARE FOR 3D PARALLEL-BEAM CT \n");
    fprintf(stdout, "build time: %s, %s\n", __DATE__,  __TIME__);
    fprintf(stdout, "\nCommand line Format for Executable File %s : \n", ExecFileName);
    fprintf(stdout, "%s ./<Executable File Name>  -i <InputFileName>[.imgparams] -j <InputFileName>[.sinoparams] -m <OutputFileName>[.2Dsysmatrix] \n\n",ExecFileName);
    fprintf(stdout, "Note : The file extensions above enclosed in \"[ ]\" symbols are necessary \n");
    fprintf(stdout, "but should be omitted from the command line arguments\n");
}


int CmdLineHelpCheck(char *string)
{
    if( (strcmp(string,"-h")==0) || (strcmp(string,"-help")==0) || (strcmp(string,"--help")==0) || (strcmp(string,"help")==0) )
        return 1;
    else
        return 0;
}


