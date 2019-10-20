#ifndef _ACOMP_H_
#define _ACOMP_H_

#include "MBIRModularDefs.h"

/* Computation options */
#define WIDE_BEAM   /* Finite element analysis of detector channel, accounts for sensitivity variation across its aperture */
#define LEN_PIX 511 /* determines the spatial resolution for Detector-Pixel computation. Higher LEN_PIX, higher resolution */
                    /* In this implementation, spatial resolution is : [2*PixelDimension/LEN_PIX]^(-1) */
#define LEN_DET 101 /* No. of Detector Elements */
                    /* Each detector channel is "split" into LEN_DET smaller elements ... */
                    /* to account for detector sensitivity variation across its aperture */

int compA_ij_first;
int Ntheta, NChannels, N_x, N_y;
float DeltaChannel, DeltaPix, t_0,x_0,y_0;
float dprof[101];


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
}address_arr;


struct minStruct
{
	int *bandMin;
};
struct maxStruct
{
	int *bandMax;
};



/* The System matrix does not vary with slice for 3-D Parallel Geometry */
/* So, the method of compuatation is same as that of 2-D Parallel Geometry */

/* Compute Pixel-Detector Profile for 3D Parallel Beam Geometry */
float **ComputePixelProfile3DParallel(struct SinoParams3DParallel *sinoparams, struct ImageParams3D *imgparams);
struct pointerAddress A_comp(struct minStruct *bandMinMap,struct maxStruct *bandMaxMap,struct AValues_char ** A_Padded_Map,float * max_num_pointer,struct SinoParams3DParallel *sinoparams,int sum,char **recon_mask,int *order,struct ImageParams3D *imgparams, float** pix_prof,char* sysMatrixPath,int pieceLength);
void readAmatrix(char *fname,struct AValues_char ** A_Padded_Map, float * max_num_pointer,struct ImageParams3D *imgparams, struct SinoParams3DParallel *sinoparams, int sum,struct minStruct *bandMinMap,struct maxStruct *bandMaxMap, int pieceLength);
void writeAmatrix(char *fname,struct AValues_char ** A_Padded_Map, float * max_num_pointer,struct ImageParams3D *imgparams, struct SinoParams3DParallel *sinoparams, int sum,struct minStruct *bandMinMap,struct maxStruct *bandMaxMap, int pieceLength);
void A_comp_ij(int im_row,int im_col,struct SinoParams3DParallel *sinoparams,struct ImageParams3D *imgparams,float **pix_prof,struct ACol *A_col,float *A_Values);
int computePieceLength(int NViews);

#endif
