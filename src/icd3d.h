
#ifndef _ICD3D_H_
#define _ICD3D_H_

float ICDStep3D(
	struct ReconParams reconparams,
	float THETA1,
	float THETA2,
	float tempV,
	float *neighbors);

float QGGMRF3D_Update(
	struct ReconParams reconparams,
	float tempV,
	float *neighbors,
	float THETA1,
	float THETA2);

float QGGMRF_SurrogateCoeff(
	float delta,
	struct ReconParams reconparams);

float QGGMRF_Potential(float delta, struct ReconParams *Rparams);
void ExtractNeighbors3D(float *neighbors,int jx,int jy,float *image,struct ImageParams3D imgparams);

#endif
