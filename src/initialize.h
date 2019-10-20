#ifndef _INITIALIZE_H_
#define _INITIALIZE_H_

#include "mbir_ct.h"
#include "MBIRModularDefs.h"

void Initialize_Image(struct Image3D *Image, struct CmdLine *cmdline, char *ImageReconMask, float InitValue, float OutsideROIValue);
char *GenImageReconMask(struct ImageParams3D *imgparams);
void readSystemParams(struct CmdLine *cmdline, struct ImageParams3D *imgparams, struct SinoParams3DParallel *sinoparams, struct ReconParamsQGGMRF3D *reconparams);
void NormalizePriorWeights3D(struct ReconParamsQGGMRF3D *reconparams);

#endif
