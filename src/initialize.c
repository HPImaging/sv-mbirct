
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>	/* strcmp */
#include <getopt.h>	/* getopt */

#include "MBIRModularUtils_2D.h"
#include "MBIRModularUtils_3D.h"
#include "allocate.h"
#include "initialize.h"

/* Initialize constant image */
void Initialize_Image(struct Image3D *Image, struct CmdLineMBIR *cmdline, float InitValue)
{
    //fprintf(stdout, "\nInitializing Image ... \n");
    
    if(strcmp(cmdline->InitImageDataFile,"NA") == 0) /* Image file not available */
        GenConstImage(Image, InitValue);               /* generate image with uniform pixel value */
    else
        ReadImage3D(cmdline->InitImageDataFile, Image); /* read image which has values in HU units */
}

/* create constant image. Each pixel value is the intial condition. */
void GenConstImage(struct Image3D *Image, float value)
{
    int i,j, N;
    
    N = Image->imgparams.Nx * Image->imgparams.Ny;
    
    for (i = 0; i < Image->imgparams.Nz;i++){
    for (j = 0; j < N; j++)
    {
        Image->image[i][j] = value;
    }
    }
}

/* Generate Image Reconstruction mask */
char *GenImageReconMask (struct Image3D *Image, float OutsideROIValue)
{
    int jx, jy, jz, Nx, Ny, Nz, Nxy;
    float x_0, y_0, Deltaxy, x, y, yy, ROIRadius, R_sq, R_sq_max;
    char *ImageReconMask;
    
    Nx = Image->imgparams.Nx;
    Ny = Image->imgparams.Ny;
    Nz = Image->imgparams.Nz;
    Deltaxy = Image->imgparams.Deltaxy;
    ROIRadius = Image->imgparams.ROIRadius;
    
    x_0 = -(Nx-1)*Deltaxy/2;
    y_0 = -(Ny-1)*Deltaxy/2;
    Nxy = Nx*Ny;
    
    /* Reconstruction Mask same for each slice, hence 2-D */
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
                {
                    ImageReconMask[jy*Nx+jx] = 0;
                    for(jz=0;jz<Nz;jz++)
                        Image->image[jz][jy*Nx+jx] = OutsideROIValue;
                }
                else
                {
                    ImageReconMask[jy*Nx+jx] = 1;
                }
            }
        }
    }
    
    return ImageReconMask;
}


/* Normalize weights to sum to 1 */
/* Only neighborhood specific */
void NormalizePriorWeights3D(
                         struct ReconParamsQGGMRF3D *reconparams)
{
    double sum;
    
    /* All neighbor weights must sum to one, assuming 8 pt neighborhood */
    sum = 4.0*reconparams->b_nearest + 4.0*reconparams->b_diag + 2.0*reconparams->b_interslice;
    
    reconparams->b_nearest /= sum;
    reconparams->b_diag /= sum;
    reconparams->b_interslice /= sum;
}

/* Wrapper to read in Image, sinogram and reconstruction parameters */
void readSystemParams_MBIR  (
                         struct CmdLineMBIR *cmdline,
                         struct ImageParams3D *imgparams,
                         struct SinoParams3DParallel *sinoparams,
                         struct ReconParamsQGGMRF3D *reconparams)
{
    //printf("\nReading Image, Sinogram and Reconstruction Parameters ... \n");
    
    if(ReadImageParams3D(cmdline->ImageParamsFile, imgparams))
      {
        fprintf(stdout,"Error in reading image parameters \n");
        exit(-1);
      }
    if(ReadSinoParams3DParallel(cmdline->SinoParamsFile, sinoparams))
      {
        fprintf(stdout,"Error in reading sinogram parameters \n");
        exit(-1);
      }
    if(ReadReconParamsQGGMRF3D(cmdline->ReconParamsFile ,reconparams))
      {
        fprintf(stdout,"Error in reading reconstruction parameters \n");
        exit(-1);
      }
          
    /* Tentatively initialize weights. Remove once this is read in directly from params file */
    NormalizePriorWeights3D(reconparams);
    
    /* Print paramters */
    //printImageParams3D(imgparams);
    //printSinoParams3DParallel(sinoparams);
    //printReconParamsQGGMRF3D(reconparams);
}

/* Read Command-line */
void readCmdLineMBIR(int argc, char *argv[], struct CmdLineMBIR *cmdline)
{
    char ch;
    
    strcpy(cmdline->InitImageDataFile, "NA"); /* default */
    
    if(argc<15)
    {
        if(argc==2 && CmdLineHelp_MBIR(argv[1]))
        {
            fprintf(stdout,"\n=========HELP==========\n");
            PrintCmdLineUsage_MBIR(argv[0]);
            exit(-1);
        }
        else
        {
         fprintf(stderr, "\nError : Improper Command line for exec-program %s, Number of arguments lower than needed \n",argv[0]);
         PrintCmdLineUsage_MBIR(argv[0]);
         exit(-1);
        }
    }
    
    /* get options */
    while ((ch = getopt(argc, argv, "i:j:k:m:s:w:r:v:")) != EOF)
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
            case 'v':
            {
                sprintf(cmdline->InitImageDataFile, "%s", optarg);
                break;
            }
            default:
            {
                printf("\nError : Command line Symbol not recongized by readParamsSysMatrix \n");
                PrintCmdLineUsage_MBIR(argv[0]);
                exit(-1);
                break;
            }
        }
    }

}

void PrintCmdLineUsage_MBIR(char *ExecFileName)
{
    fprintf(stdout, "\nBASELINE MBIR RECONSTRUCTION SOFTWARE FOR 3D PARALLEL-BEAM  CT \n");
    fprintf(stdout, "build time: %s, %s\n", __DATE__,  __TIME__);
    fprintf(stdout, "\nCommand line Format for Executable File %s : \n", ExecFileName);
    fprintf(stdout, "%s -i <InputFileName>[.imgparams] -j <InputFileName>[.sinoparams]  -k <InputFileName>[.reconparams] \
-m <InputFileName>[.2Dsysmatrix] -s <InputProjectionsBaseFileName> -w <InputWeightsBaseFileName> -r <OutputImageBaseFileName> \n",ExecFileName);
    fprintf(stdout, "\nAdditional option (to read in initial image from hard disk): -v <InitialImageBaseFileName> \n\n");
    fprintf(stdout, "Note : The necessary extensions for certain input files are mentioned above within a \"[]\" symbol above \n");
    fprintf(stdout, "However, they are NOT to be included as part of the file name in the command line arguments \n\n");
    fprintf(stdout, "The below instructions pertain to the -s, -w and -r options in the above Command line format: \n");
    fprintf(stdout, "A) \n");
    fprintf(stdout, "The Sinogram Projections measured for the 3-D volume are organized slice by slice in a single directory \n");
    fprintf(stdout, "All files within this directory must share a common BaseFileName and follow the below format :\n");
    fprintf(stdout, "<ProjectionsBaseFileName>_slice<SliceIndex>.2Dsinodata \n");
    fprintf(stdout, ", where \"SliceIndex\" is a non-negative integer indexing each slice and is printed with 4 digits. Eg : 0000 to 9999 is a valid descriptor for \"SliceIndex\" \n");
    fprintf(stdout, "B) \n");
    fprintf(stdout, "Similarly, in the case of the Sinogram Weights, the format for the file-names are :\n");
    fprintf(stdout, "<WeightsBaseFileName>_slice<SliceIndex>.2Dweightdata \n");
    fprintf(stdout, "C) \n");
    fprintf(stdout, "Similarly, the Reconstructed (Output) Image is organized slice by slice. The format for the file-names are :\n");
    fprintf(stdout, "<ImageBaseFileName>_slice<SliceIndex>.2Dimgdata \n");
    fprintf(stdout, "\nAlso see sample files run.sh under Data/Demo/ for the correct format \n\n");
}

int CmdLineHelp_MBIR(char *string)
{
    if( (strcmp(string,"-h")==0) || (strcmp(string,"-help")==0) || (strcmp(string,"--help")==0) || (strcmp(string,"help")==0) )
        return 1;
    else
        return 0;
}


