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
struct pointerAddress A_comp(
    struct AValues_char ** A_Padded_Map,
    float * max_num_pointer,
    struct SVParams svpar,
    struct SinoParams3DParallel *sinoparams,
    char *recon_mask,
    struct ImageParams3D *imgparams);

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
