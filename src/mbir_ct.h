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
    char SinoParamsFile[256], SinoParamsFileFlag;
    char ImageParamsFile[256], ImageParamsFileFlag;
    char ReconParamsFile[256], ReconParamsFileFlag;
    char SinoDataFile[256], SinoDataFileFlag;
    char SinoWeightsFile[256], SinoWeightsFileFlag;
    char ReconImageFile[256], ReconImageFileFlag;
    char SysMatrixFile[256], SysMatrixFileFlag;
    char InitImageFile[256];
    char inputProjectionFile[256];
    char outputProjectionFile[256];
    char ProxMapImageFile[256];
    /* operation flags */
    char reconFlag;              /* 0=pre-compute mode; 1=reconstruct (QGGMRF), 2=reconstruct (PandP) */
    char readInitImageFlag;
    char readInitProjectionFlag;
    char writeProjectionFlag;    /* 0=don't; 1=write Projection default IC (single slice); 2=write multi-slice */
    char readAmatrixFlag;        /* 0=compute A; 1=read A */
    char writeAmatrixFlag;
};



#endif
