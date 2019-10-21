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


/* Functions */
float **ComputePixelProfile3DParallel(struct SinoParams3DParallel *sinoparams, struct ImageParams3D *imgparams);
struct pointerAddress A_comp(
	struct AValues_char ** A_Padded_Map,
	float * max_num_pointer,
	struct SVParams svpar,
	struct SinoParams3DParallel *sinoparams,
	char **recon_mask,
	int *order,
	struct ImageParams3D *imgparams,
	float** pix_prof,
	char* sysMatrixPath);
void A_comp_ij(int im_row,int im_col,struct SinoParams3DParallel *sinoparams,struct ImageParams3D *imgparams,float **pix_prof,struct ACol *A_col,float *A_Values);

void readAmatrix(
	char *fname,
	struct AValues_char ** A_Padded_Map,
	float * max_num_pointer,
	struct ImageParams3D *imgparams,
	struct SinoParams3DParallel *sinoparams,
	struct SVParams svpar);
void writeAmatrix(
	char *fname,
	struct AValues_char ** A_Padded_Map,
	float * max_num_pointer,
	struct ImageParams3D *imgparams,
	struct SinoParams3DParallel *sinoparams,
	struct SVParams svpar);


#endif
