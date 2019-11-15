#ifndef _MBIRCT_H_
#define _MBIRCT_H_


//#define USE_INTEL_MEMCPY
//#define find_RMSE 

#define SVLENGTH 9
#define OVERLAPPINGDISTANCE 2
#define SVDEPTH 4

struct SVParams
{
	struct minStruct *bandMinMap;
	struct maxStruct *bandMaxMap;
	int SVLength;
	int overlap;
	int SVDepth;
	int SV_per_Z;
	int SVsPerLine;
	int Nsv;
	int pieceLength;
};

#if 0
// Experimental
struct SysMatrixSV
{
	struct AValues_char **A_Padded_Map;
	float *max_num_pointer;
};
#endif


struct CmdLine
{
    char ReconType;            /* 1:QGGMRF, 2:PandP */
    char SinoParamsFile[256], SinoParamsFileFlag;
    char ImageParamsFile[256], ImageParamsFileFlag;
    char ReconParamsFile[256], ReconParamsFileFlag;
    char SinoDataFile[256], SinoDataFileFlag;
    char SinoWeightsFile[256], SinoWeightsFileFlag;
    char ReconImageDataFile[256], ReconImageDataFileFlag;
    char SysMatrixFile[256], SysMatrixFileFlag;
    char InitImageDataFile[256], InitImageDataFileFlag;
    char ProxMapImageDataFile[256], ProxMapImageDataFileFlag;
    char InitProjFile[256], InitProjFileFlag;
    char precompAmatrixFlag;
    char precompInitProjFlag;
    char reconFlag;
    char readAmatrixFlag;
    char readInitProjFlag;
};



#endif
