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
    char SinoParamsFile[256];
    char ImageParamsFile[256];
    char ReconParamsFile[256];
    char SinoDataFile[256];
    char SinoWeightsFile[256];
    char ReconImageDataFile[256]; /* output */
    char SysMatrixFile[256];
    char InitImageDataFile[256]; /* optional input */
};



#endif
