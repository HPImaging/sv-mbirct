
#ifndef _RECON3D_H_
#define _RECON3D_H_

#include "MBIRModularDefs.h"
#include "A_comp.h"
#include "initialize.h"

#define  c_ratio 0.07
	
#define convergence_rho 0.7

// #define find_RMSE 

void MBIRReconstruct3D(struct Image3D *Image,struct Sino3DParallel *sinogram,struct ReconParamsQGGMRF3D reconparams,char *ImageReconMask,struct minStruct *bandMinMap,struct maxStruct *bandMaxMap,struct AValues_char ** A_Padded_Map,float *max_num_pointer,struct CmdLineMBIR * cmdLine,int sum, int pieceLength);

float MAPCostFunction3D(float **e, struct Image3D *Image, struct Sino3DParallel *sinogram, struct ReconParamsQGGMRF3D *reconparams);

/*void forwardProject3D(float *AX, struct Image3D *X, struct SysMatrix2D *A, struct SinoParams3DParallel sinoparams); */

void shuffle(int *order, int len);

#endif
