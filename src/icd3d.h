
#ifndef _ICD3D_H_
#define _ICD3D_H_

float ICDStep3D(
	struct ReconParamsQGGMRF3D reconparams,
	float THETA1,
	float THETA2,
	float tempV,
	float *neighbors,
	float pow_sigmaX_p,
	float pow_sigmaX_q,
	float pow_T_qmp);

float QGGMRF3D_Update(
	struct ReconParamsQGGMRF3D reconparams,
	float tempV,
	float *neighbors,
	float THETA1,
	float THETA2,
	float pow_sigmaX_p,
	float pow_sigmaX_q,
	float pow_T_qmp);

float QGGMRF_SurrogateCoeff(
	float delta,
	struct ReconParamsQGGMRF3D reconparams,
	float pow_sigmaX_p,
	float pow_sigmaX_q,
	float pow_T_qmp);

float QGGMRF_Potential(float delta, struct ReconParamsQGGMRF3D *Rparams);
void ExtractNeighbors3D(float *neighbors,int jx,int jy,float *image,struct ImageParams3D imgparams);

#endif
