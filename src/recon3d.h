#ifndef _RECON3D_H_
#define _RECON3D_H_

#include "mbir_ct.h"
#include "MBIRModularDefs.h"
#include "A_comp.h"


void MBIRReconstruct3D(
	struct Image3D *Image,
	struct Sino3DParallel *sinogram,
	float **e,
	struct ReconParams reconparams,
	struct SVParams svpar,
	struct AValues_char **A_Padded_Map,
	float *Aval_max_ptr,
	char *ImageReconMask,
	char verboseLevel);

void MBIRReconstruct(
    float *image,
    float *sino,
    float *weight,
    float *proj_init,
    struct ImageParams3D imgparams,
    struct SinoParams3DParallel sinoparams,
    struct ReconParams reconparams,
    char *Amatrix_fname,
    char verboseLevel);

void forwardProject2D(
	float *e,
	float *x,
	struct AValues_char **A_Padded_Map,
	float *Aval_max_ptr,
	struct SinoParams3DParallel *sinoparams,
	struct ImageParams3D *imgparams,
	struct SVParams svpar);

void forwardProject(
    float *image,
    float *proj,
    struct ImageParams3D imgparams,
    struct SinoParams3DParallel sinoparams,
    char *Amatrix_fname,
    char verboseLevel);

#endif
