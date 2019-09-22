
#ifndef _INITIALIZE_H_
#define _INITIALIZE_H_

struct CmdLineMBIR{
    
    char ImageParamsFile[200];
    char InitImageDataFile[200]; /* optional input */
    char ReconImageDataFile[200]; /* output */
    char SinoParamsFile[200];
    char SinoDataFile[200];
    char SinoWeightsFile[200];
    char ReconParamsFile[200];
    char SysMatrixFile[200];
};


void Initialize_Image(struct Image3D *Image, struct CmdLineMBIR *cmdline, float InitValue);
void GenConstImage(struct Image3D *Image, float value);
char *GenImageReconMask (struct Image3D *Image, float OutsideROIValue);
void readSystemParams_MBIR ( struct CmdLineMBIR *cmdline, struct ImageParams3D *imgparams, struct SinoParams3DParallel *sinoparams, struct ReconParamsQGGMRF3D *reconparams);
void NormalizePriorWeights3D(struct ReconParamsQGGMRF3D *reconparams);
void readCmdLineMBIR(int argc, char *argv[], struct CmdLineMBIR *cmdline);
void PrintCmdLineUsage_MBIR(char *ExecFileName);
int CmdLineHelp_MBIR(char *string);

#endif
