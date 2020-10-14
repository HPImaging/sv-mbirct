#ifndef _MBIRCT_H_
#define _MBIRCT_H_


//This directive is set with the compiler option -DICC
//#define ICC

#ifdef ICC
    /* Can't find an Intel header that prototypes this, so adding it here to
     * get rid of the compiler "implicit declaration" warnings */
    void *_intel_fast_memcpy(void *dest, const void * src, size_t n);
#endif

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
	int SVsPerRow;
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
    char SinoParamsFile[1024], SinoParamsFileFlag;
    char ImageParamsFile[1024], ImageParamsFileFlag;
    char ReconParamsFile[1024], ReconParamsFileFlag;
    char SinoDataFile[1024], SinoDataFileFlag;
    char SinoWeightsFile[1024], SinoWeightsFileFlag;
    char ReconImageFile[1024], ReconImageFileFlag;
    char SysMatrixFile[1024], SysMatrixFileFlag;
    char InitImageFile[1024];
    char inputProjectionFile[1024];
    char outputProjectionFile[1024];
    char ProxMapImageFile[1024];
    /* operation flags */
    char reconFlag;              /* 0=pre-compute mode; 1=reconstruct (QGGMRF), 2=reconstruct (PandP) */
    char readInitImageFlag;
    char readInitProjectionFlag;
    char writeProjectionFlag;    /* 0=don't; 1=write Projection default IC (single slice); 2=write multi-slice */
    char readAmatrixFlag;        /* 0=compute A; 1=read A */
    char writeAmatrixFlag;
};



#endif
