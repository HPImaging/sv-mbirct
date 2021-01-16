#ifndef _RECON3D_H_
#define _RECON3D_H_

#include "MBIRModularDefs.h"
#include "A_comp.h"

#define SVLENGTH 9
#define OVERLAPPINGDISTANCE 2
#define SVDEPTH 4

/* Can't find an Intel header that prototypes this, so adding it here to
 * get rid of the compiler "implicit declaration" warnings */
#ifdef ICC
    void *_intel_fast_memcpy(void *dest, const void * src, size_t n);
#endif

void MBIRReconstruct(
    float *image,
    float *sino,
    float *weight,
    float *proj_init,
    float *proximalmap,
    struct ImageParams3D imgparams,
    struct SinoParams3DParallel sinoparams,
    struct ReconParams reconparams,
    char *Amatrix_fname,
    char verboseLevel);

void forwardProject(
    float *proj,
    float *image,
    struct ImageParams3D imgparams,
    struct SinoParams3DParallel sinoparams,
    char *Amatrix_fname,
    char verboseLevel);

#endif
