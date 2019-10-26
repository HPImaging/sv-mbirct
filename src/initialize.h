#ifndef _INITIALIZE_H_
#define _INITIALIZE_H_

#include "mbir_ct.h"
#include "MBIRModularDefs.h"

void initSVParams(struct SVParams *svpar,struct ImageParams3D imgparams,struct SinoParams3DParallel sinoparams);
int computePieceLength(int NViews);
char *GenImageReconMask(struct ImageParams3D *imgparams);
void NormalizePriorWeights3D(struct ReconParamsQGGMRF3D *reconparams);

void initImage(
	struct Image3D *Image,
	struct CmdLine *cmdline,
	char *ImageReconMask,
	float InitValue,
	float OutsideROIValue);

void initProjectionError(
	float **e,
	struct Image3D *Image,
	struct Sino3DParallel *sinogram,
	struct SVParams svpar,
	struct AValues_char **A_Padded_Map,
	float *max_num_pointer,
	struct CmdLine *cmdline);


#endif
