#ifndef _ACOMP_H_
#define _ACOMP_H_

#include "MBIRModularDefs.h"


struct ACol
{
	int n_index;
	unsigned char *countTheta;
	int *minIndex;
};

struct AValues_char{
	unsigned char *val;
	int *pieceWiseMin;
	int *pieceWiseWidth;
	int length;
};

struct pointerAddress{
	struct ACol** addressA;
	struct AValues_char** addressB;
};


struct minStruct
{
	int *bandMin;
};
struct maxStruct
{
	int *bandMax;
};


float **ComputePixelProfile3DParallel(struct SinoParams3DParallel *sinoparams, struct ImageParams3D *imgparams);
struct pointerAddress A_comp(struct minStruct *bandMinMap,struct maxStruct *bandMaxMap,struct AValues_char ** A_Padded_Map,float * max_num_pointer,struct SinoParams3DParallel *sinoparams,int sum,char **recon_mask,int *order,struct ImageParams3D *imgparams, float** pix_prof,char* sysMatrixPath,int pieceLength);
void readAmatrix(char *fname,struct AValues_char ** A_Padded_Map, float * max_num_pointer,struct ImageParams3D *imgparams, struct SinoParams3DParallel *sinoparams, int sum,struct minStruct *bandMinMap,struct maxStruct *bandMaxMap, int pieceLength);
void writeAmatrix(char *fname,struct AValues_char ** A_Padded_Map, float * max_num_pointer,struct ImageParams3D *imgparams, struct SinoParams3DParallel *sinoparams, int sum,struct minStruct *bandMinMap,struct maxStruct *bandMaxMap, int pieceLength);
void A_comp_ij(int im_row,int im_col,struct SinoParams3DParallel *sinoparams,struct ImageParams3D *imgparams,float **pix_prof,struct ACol *A_col,float *A_Values);
int computePieceLength(int NViews);

#endif
