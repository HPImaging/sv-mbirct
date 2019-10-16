#ifndef _RECON3D_H_
#define _RECON3D_H_

#include "MBIRModularDefs.h"
#include "A_comp.h"
#include "initialize.h"

// #define find_RMSE 

void MBIRReconstruct3D(
	struct Image3D *Image,
	struct Sino3DParallel *sinogram,
	struct ReconParamsQGGMRF3D reconparams,
	struct minStruct *bandMinMap,
	struct maxStruct *bandMaxMap,
	struct AValues_char ** A_Padded_Map,
	float *max_num_pointer,
	struct CmdLineMBIR * cmdLine,
	int sum,
	int pieceLength);

void forwardProject2D(
	float *e,
	float InitValue,
	float *max_num_pointer,
	struct AValues_char ** A_Padded_Map,
	struct minStruct *bandMinMap,
	struct SinoParams3DParallel *sinoparams,
	struct ImageParams3D *imgparams,
	int pieceLength);


#endif
