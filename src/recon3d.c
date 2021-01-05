
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <string.h>
#include <omp.h>

#include "mbir_ct.h"
#include "MBIRModularDefs.h"
#include "MBIRModularUtils.h"
#include "allocate.h"
#include "icd3d.h"
#include "heap.h"
#include "A_comp.h"
#include "initialize.h"
#include "recon3d.h"

//#define COMP_RMSE

/* Internal functions */
void super_voxel_recon(int jj,struct SVParams svpar,unsigned long *NumUpdates,float *totalValue,float *totalChange,int iter,
	char *phaseMap,long *order,int *indexList,float **w,float **e,
	struct AValues_char **A_Padded_Map,float *Aval_max_ptr,struct heap_node *headNodeArray,
	struct SinoParams3DParallel sinoparams,struct ReconParams reconparams,struct Image3D *Image,
	float *voxelsBuffer1,float *voxelsBuffer2,char *group_array,int group_id);
void coordinateShuffle(int *order1, int *order2,int len);
void three_way_shuffle(long *order1, char *order2, struct heap_node *headNodeArray,int len);
float MAPCostFunction3D(float **e,struct Image3D *Image,struct Sino3DParallel *sinogram,struct ReconParams *reconparams);


void MBIRReconstruct3D(
	struct Image3D *Image,
	struct Sino3DParallel *sinogram,
	float **e,  /* e=y-Ax, error */
	struct ReconParams reconparams,
	struct SVParams svpar,
	struct AValues_char ** A_Padded_Map,
	float *Aval_max_ptr,
	char *ImageReconMask,
	char verboseLevel)
{
	int i,j,jj,p,t,iter,it_print=1;
	float **w;  /* projections weights data */
	float *voxelsBuffer1;  /* the first N entries are the voxel values.  */
	float *voxelsBuffer2;
	unsigned long NumUpdates=0;
	float totalValue=0,totalChange=0,equits=0;
	float avg_update,avg_update_rel;

	struct heap priorityheap;
	initialize_heap(&priorityheap);
	long *order;
	char *phaseMap;
	struct timeval tm1,tm2;

	w = sinogram->weight; /* vector of weights for each sinogram measurement */
	int Nx = Image->imgparams.Nx;
	int Ny = Image->imgparams.Ny;
	int Nxy = Nx*Ny;
	int Nz = Image->imgparams.Nz;
	int MaxIterations = reconparams.MaxIterations;
	float StopThreshold = reconparams.StopThreshold;
	int SVLength = svpar.SVLength;
	int overlappingDistance = svpar.overlap;
	int SV_depth = svpar.SVDepth;
	int SV_per_Z = svpar.SV_per_Z;
	int SVsPerRow = svpar.SVsPerRow;
	int sum = svpar.Nsv;

	float c_ratio=0.07;
	float convergence_rho=0.7;
	int rep_num=(int)ceil(1/(4*c_ratio*convergence_rho));

	int NumMaskVoxels=0;
	for(j=0;j<Nxy;j++)
	if(ImageReconMask[j])
		NumMaskVoxels++;

	#ifdef COMP_RMSE
		struct Image3D Image_ref;
		Image_ref.imgparams.Nx = Image->imgparams.Nx;
		Image_ref.imgparams.Ny = Image->imgparams.Ny;
		Image_ref.imgparams.Nz = Image->imgparams.Nz;
		Image_ref.imgparams.FirstSliceNumber = Image->imgparams.FirstSliceNumber;
		Image_ref.imgparams.NumSliceDigits = Image->imgparams.NumSliceDigits;
		AllocateImageData3D(&Image_ref);
		ReadImage3D("ref/ref",&Image_ref);

		float ** image_ref = Image_ref.image;
		double rms_err=0,rms_val=0;
		int Nz0=0, Nz1=Nz;
		for(i=Nz0; i<Nz1; i++)
		for(j=0; j<Nxy; j++)
		if(ImageReconMask[j]) {
			rms_val += image_ref[i][j]*image_ref[i][j];
			rms_err += (Image->image[i][j]-image_ref[i][j])*(Image->image[i][j]-image_ref[i][j]);
		}
		rms_val = sqrt(rms_val/((float)NumMaskVoxels*(Nz1-Nz0)));
		rms_err = sqrt(rms_err/((float)NumMaskVoxels*(Nz1-Nz0)));
		FILE *fp_mse=fopen("rmse.txt","w");
		fprintf(fp_mse,"equits|rms_err|rms_val|rms_err/rms_val\n");
		fprintf(fp_mse,"%.2f %g %g %g\n",equits,rms_err,rms_val,rms_err/rms_val);
	#endif

	order = (long *) mget_spc(sum*SV_per_Z,sizeof(long));

	/* Order of pixel updates need NOT be raster order, just initialize */
	t=0;
	for(p=0;p<Nz;p+=SV_depth)
	for(i=0;i<Ny;i+=(SVLength*2-overlappingDistance))
	for(j=0;j<Nx;j+=(SVLength*2-overlappingDistance))
	{
		order[t]=(long)p*Nxy+i*Nx+j;  /* order is the first voxel coordinate, not the center */
		t++;
	}

	phaseMap = (char *) mget_spc(sum*SV_per_Z,sizeof(char));

	#pragma omp parallel for private(jj) schedule(dynamic)
	for(i=0;i<SV_per_Z;i++)
	for(jj=0;jj<sum;jj++)
	{
		if((jj/SVsPerRow)%2==0)
		{
			if((jj%SVsPerRow)%2==0)
				phaseMap[i*sum+jj]=0;
			else
				phaseMap[i*sum+jj]=1;			
		}
		else
		{
			if((jj%SVsPerRow)%2==0)
				phaseMap[i*sum+jj]=2;
			else
				phaseMap[i*sum+jj]=3;			
		}
	}

	char group_id_list[SV_per_Z][4];

	for(i=0;i<SV_per_Z;i++){
		if(i%4==0){
			group_id_list[i][0]=0;
			group_id_list[i][1]=3;
			group_id_list[i][2]=1;
			group_id_list[i][3]=2;										
		}
		else if(i%4==1){	
			group_id_list[i][0]=3;
			group_id_list[i][1]=0;
			group_id_list[i][2]=2;
			group_id_list[i][3]=1;										
		}
		else if(i%4==2){			
			group_id_list[i][0]=1;
			group_id_list[i][1]=2;
			group_id_list[i][2]=3;
			group_id_list[i][3]=0;										
		}
		else{
			group_id_list[i][0]=2;
			group_id_list[i][1]=1;
			group_id_list[i][2]=0;
			group_id_list[i][3]=3;										
		}				
	}
	srand(time(NULL));

	struct heap_node *headNodeArray;
	headNodeArray = (struct heap_node *) mget_spc(sum*SV_per_Z,sizeof(struct heap_node));

	for(i=0;i<SV_per_Z;i++)
	for(jj=0;jj<sum;jj++)
	{
		headNodeArray[i*sum+jj].pt=i*sum+jj;
		headNodeArray[i*sum+jj].x=0.0;
	}
	int indexList_size=(int) sum*SV_per_Z*4*c_ratio*(1-convergence_rho);	
	int indexList[indexList_size];   	             	    
    
	#ifdef ICC
		voxelsBuffer1 = (float *)_mm_malloc(Nxy*sizeof(float),64);
		voxelsBuffer2 = (float *)_mm_malloc(Nxy*sizeof(float),64);
	#else
		voxelsBuffer1 = (float *) malloc(Nxy*sizeof(float));
		voxelsBuffer2 = (float *) malloc(Nxy*sizeof(float));
		//voxelsBuffer1 = (float *) aligned_alloc(64,Nxy*sizeof(float));
		//voxelsBuffer2 = (float *) aligned_alloc(64,Nxy*sizeof(float));
		//voxelsBuffer1 = (float *) _aligned_malloc(Nxy*sizeof(float),64);
		//voxelsBuffer2 = (float *) _aligned_malloc(Nxy*sizeof(float),64);
	#endif

	for(i=0;i<Nxy;i++) {
		voxelsBuffer1[i]=0;
		voxelsBuffer2[i]=0;
	}

	//coordinateShuffle(&order[0],&phaseMap[0],sum*SV_per_Z);
	long tmp_long;
	char tmp_char;
	for(i=0; i<sum*SV_per_Z-1; i++)
	{
		j = i + (rand() % (sum*SV_per_Z-i));
		tmp_long = order[j];
		order[j] = order[i];
		order[i] = tmp_long;
		tmp_char = phaseMap[j];
		phaseMap[j] = phaseMap[i];
		phaseMap[i] = tmp_char;
	}

	gettimeofday(&tm1,NULL);
         
	iter=0;
	char stop_FLAG=0;
	int startIndex=0;
	int endIndex=0;

	#pragma omp parallel
	{
		while(stop_FLAG==0 && equits<MaxIterations && iter<100*MaxIterations)
		{
			#pragma omp single
			{		
				if(iter==0)
				{
					startIndex=0;
					endIndex=sum*SV_per_Z;
				}	
				else
				{
					if((iter-1)%(2*rep_num)==0 && iter!=1)
						three_way_shuffle(&order[0],&phaseMap[0],&headNodeArray[0],sum*SV_per_Z);
			
					if(iter%2==1)
					{
						priorityheap.size=0;
						for(jj=0;jj<sum*SV_per_Z;jj++){
							heap_insert(&priorityheap, &(headNodeArray[jj]));
						}						
						startIndex=0;					
						endIndex=indexList_size;					

						for(i=0;i<endIndex;i++){
							struct heap_node tempNode;
							get_heap_max(&priorityheap, &tempNode);
							indexList[i]=tempNode.pt;
						}	
					}				
					else{					
						startIndex=((iter-2)/2)%rep_num*sum*SV_per_Z/rep_num;
						endIndex=(((iter-2)/2)%rep_num+1)*sum*SV_per_Z/rep_num;
					}
				}
			}

			int group=0;

			for (group = 0; group < 4; group++)
			{
				#pragma omp for schedule(dynamic)  reduction(+:NumUpdates) reduction(+:totalValue) reduction(+:totalChange)
				for (jj = startIndex; jj < endIndex; jj+=1)
					super_voxel_recon(jj,svpar,&NumUpdates,&totalValue,&totalChange,iter,&phaseMap[0],order,&indexList[0],w,e,A_Padded_Map,&Aval_max_ptr[0],&headNodeArray[0],sinogram->sinoparams,reconparams,Image,voxelsBuffer1,voxelsBuffer2,&group_id_list[0][0],group);

			}

			#pragma omp single
			{
 
			if(NumUpdates>0)
			{
				avg_update = totalChange/NumUpdates;
				float avg_value = totalValue/NumUpdates;
				avg_update_rel = avg_update/avg_value * 100;
				//printf("avg_update %f, avg_value %f, avg_update_rel %f\n",avg_update,avg_value,avg_update_rel);
			}
			
			/*
			float cost = MAPCostFunction3D(e, Image, sinogram, &reconparams);
			fprintf(stdout, "it %d cost = %-15f, avg_update %f \n", iter, cost, avg_update);
			*/

			if (avg_update_rel < StopThreshold && (endIndex!=0))
				stop_FLAG = 1;

			iter++;
			equits += (float)NumUpdates/((float)NumMaskVoxels*Nz);

			if(verboseLevel)
			if(equits > it_print)
			{
				fprintf(stdout,"\titeration %d, average change %.4f %%\n",it_print,avg_update_rel);
				it_print++;
			}

			#ifdef COMP_RMSE
				rms_err=0;
				for(i=Nz0; i<Nz1; i++)
				for(j=0; j<Nxy; j++)
				if(ImageReconMask[j])
					rms_err += (Image->image[i][j]-image_ref[i][j])*(Image->image[i][j]-image_ref[i][j]);
				rms_err = sqrt(rms_err/((float)NumMaskVoxels*(Nz1-Nz0)));
				fprintf(fp_mse,"%.2f %g %g %g\n",equits,rms_err,rms_val,rms_err/rms_val);
			#endif

			NumUpdates=0;
			totalValue=0;
			totalChange=0;

			}

		}
	}
	
	if(verboseLevel)
	{
		if(StopThreshold <= 0)
			fprintf(stdout,"\tNo stopping condition--running fixed iterations\n");
		else if(stop_FLAG == 1)
			fprintf(stdout,"\tReached stopping condition\n");
		else
			fprintf(stdout,"\tWarning: Didn't reach stopping condition\n");
	}

	if(verboseLevel>1)
	{
		fprintf(stdout,"\tEquivalent iterations = %.1f, (non-homogeneous iterations = %d)\n",equits,iter);
		fprintf(stdout,"\tAverage update in last iteration (relative) = %f %%\n",avg_update_rel);
		fprintf(stdout,"\tAverage update in last iteration (magnitude) = %.4g\n",avg_update);
		gettimeofday(&tm2,NULL);
		//unsigned long long tt = 1000 * (tm2.tv_sec - tm1.tv_sec) + (tm2.tv_usec - tm1.tv_usec) / 1000;
		//printf("\tReconstruction time = %llu ms (iterations only)\n", tt);
	}

	#ifdef ICC
		_mm_free((void *)voxelsBuffer1);
		_mm_free((void *)voxelsBuffer2);
	#else
		free((void *)voxelsBuffer1);
		free((void *)voxelsBuffer2);
		//_aligned_free((void *)voxelsBuffer1);
		//_aligned_free((void *)voxelsBuffer2);
	#endif

	free((void *)headNodeArray);
	free_heap((void *)&priorityheap);
	free((void *)phaseMap);
	free((void *)order);

	#ifdef COMP_RMSE
		FreeImageData3D(&Image_ref);
		fclose(fp_mse);
	#endif

}   /*  END MBIRReconstruct3D()  */


				

void forwardProject2D(
	float *e,
	float *x,
	struct AValues_char ** A_Padded_Map,
	float *Aval_max_ptr,
	struct SinoParams3DParallel *sinoparams,
	struct ImageParams3D *imgparams,
	struct SVParams svpar)
{
	int jx,jy,i,r,j,p;
	int SVLength = svpar.SVLength;
	int overlappingDistance = svpar.overlap;
	struct minStruct * bandMinMap = svpar.bandMinMap;
	int pieceLength = svpar.pieceLength;
	int SVsPerRow = svpar.SVsPerRow;
	int NViewSets = sinoparams->NViews/pieceLength;
	int Nx = imgparams->Nx;
	int Ny = imgparams->Ny;
	int NChannels = sinoparams->NChannels;
	int M = sinoparams->NViews*sinoparams->NChannels;

	for (i = 0; i < M; i++)
		e[i] = 0.0;

	for (jy = 0; jy < Ny; jy++)
	for (jx = 0; jx < Nx; jx++)
	{
		int SV_ind_y = jy/(2*SVLength-overlappingDistance);
		int SV_ind_x = jx/(2*SVLength-overlappingDistance);
		int SVPosition = SV_ind_y*SVsPerRow + SV_ind_x;

		int SV_jy = SV_ind_y*(2*SVLength-overlappingDistance);
		int SV_jx = SV_ind_x*(2*SVLength-overlappingDistance);
		int VoxelPosition = (jy-SV_jy)*(2*SVLength+1)+(jx-SV_jx);

		// fprintf(stdout,"jy %d jx %d SVPosition %d SV_jy %d SV_jx %d VoxelPosition %d \n",jy,jx,SVPosition,SV_jy,SV_jx,VoxelPosition);
		// I think the second condition will always be true
		if (A_Padded_Map[SVPosition][VoxelPosition].length > 0 && VoxelPosition < ((2*SVLength+1)*(2*SVLength+1)))
		{
			unsigned char* A_padd_Tr_ptr = &A_Padded_Map[SVPosition][VoxelPosition].val[0];
			float rescale = Aval_max_ptr[jy*Nx+jx]*(1.0/255);
			float xval = x[jy*Nx+jx];

			for(p=0;p<NViewSets;p++)
			{
				int myCount = A_Padded_Map[SVPosition][VoxelPosition].pieceWiseWidth[p];
				int pieceWiseMin = A_Padded_Map[SVPosition][VoxelPosition].pieceWiseMin[p];
				int position = p*pieceLength*NChannels + pieceWiseMin;

				for(r=0;r<myCount;r++)
				for(j=0;j<pieceLength;j++)
				{
					int bandMin = bandMinMap[SVPosition].bandMin[p*pieceLength+j];
					int proj_idx = position + j*NChannels + bandMin + r;

					if((pieceWiseMin + bandMin + r) >= NChannels || proj_idx >= M ) {
						fprintf(stderr,"forwardProject() out of bounds: p %d r %d j %d\n",p,r,j);
						fprintf(stderr,"forwardProject() out of bounds: total_1 %d total_2 %d\n",pieceWiseMin+bandMin+r,proj_idx);
						exit(-1);
					}
					else
						e[proj_idx] += A_padd_Tr_ptr[r*pieceLength+j]*rescale * xval;
				}
				A_padd_Tr_ptr += myCount*pieceLength;
			}
		}
	}

}   /* END forwardProject2D() */


void super_voxel_recon(
	int jj,
	struct SVParams svpar,
	unsigned long *NumUpdates,
	float *totalValue,
	float *totalChange,
	int iter,
	char *phaseMap,
	long *order,
	int *indexList,
	float **w,
	float **e,
	struct AValues_char ** A_Padded_Map,
	float *Aval_max_ptr,
	struct heap_node *headNodeArray,
	struct SinoParams3DParallel sinoparams,
	struct ReconParams reconparams,
	struct Image3D *Image,
	float *voxelsBuffer1,
	float *voxelsBuffer2,
	char *group_array,
	int group_id)
{
	int jy,jx,p,i,q,t,j,currentSlice,startSlice;
	int SV_depth_modified;
	int NumUpdates_loc=0;
	float totalValue_loc=0,totalChange_loc=0;

	float ** image = Image->image;
	float ** proximalmap = reconparams.proximalmap;
	struct ImageParams3D imgparams = Image->imgparams;
	int Nx = imgparams.Nx;
	int Ny = imgparams.Ny;
	int Nz = imgparams.Nz;
	char PositivityFlag = reconparams.Positivity;

	int SVLength = svpar.SVLength;
	int overlappingDistance = svpar.overlap;
	int SV_depth = svpar.SVDepth;
	int SVsPerRow = svpar.SVsPerRow;
	struct minStruct * bandMinMap = svpar.bandMinMap;
	struct maxStruct * bandMaxMap = svpar.bandMaxMap;
	int pieceLength = svpar.pieceLength;
	int NViewSets = sinoparams.NViews/pieceLength;

	int jj_new;
	if(iter%2==0)
		jj_new=jj;
	else
		jj_new=indexList[jj];

	startSlice = order[jj_new] / Nx / Ny;
	jy = (order[jj_new] - startSlice* Nx * Ny) / Nx;
	jx = (order[jj_new] - startSlice* Nx * Ny) % Nx;

	if(phaseMap[jj_new]!=group_array[startSlice/SV_depth*4+group_id])
		return;

	if((startSlice+SV_depth)>Nz)
		SV_depth_modified=Nz-startSlice;
	else
		SV_depth_modified=SV_depth;

	int SVPosition=jy/(2*SVLength-overlappingDistance)*SVsPerRow+jx/(2*SVLength-overlappingDistance);

	int countNumber=0;	/*XW: the number of voxels chosen for a certain radius of circle*/
	int radius =SVLength;	/*XW: choose the circle radius*/
	int coordinateSize=1;	/*XW: CoordinateSize is the size of a minimum square enclosing the circle. For the baseline code, coordinateSize=1 because only 1 voxel*/
	/*XW: imagine that this is a minimum square enclosing the circle. The square 
	 * "touches" the circle on 4 coordinate points. CoordianteSize is the possible 
	 * maximum number of voxels for a certain radius of circle */
	if(radius!=0)
		coordinateSize=(2*radius+1)*(2*radius+1);
	int k_newCoordinate[coordinateSize];
	int j_newCoordinate[coordinateSize];
	int j_newAA=0;
	int k_newAA=0;
	int voxelIncrement=0;

	/*XW: choosing the voxels locations in a circle*/
	for(j_newAA=jy;j_newAA<=(jy+2*radius);j_newAA++)
	for(k_newAA=jx;k_newAA<=(jx+2*radius);k_newAA++)
	{
		if(j_newAA>=0 && k_newAA >=0 && j_newAA <Ny && k_newAA < Nx)
		{
			if(A_Padded_Map[SVPosition][voxelIncrement].length >0) {
				j_newCoordinate[countNumber]=j_newAA;
				k_newCoordinate[countNumber]=k_newAA;
				countNumber++;
			} 
		}
		voxelIncrement++;
	}

	/*XW: if no voxel chosen, we skip this loop iteration*/
	if(countNumber==0)
		return;

	coordinateShuffle(&j_newCoordinate[0],&k_newCoordinate[0],countNumber);

	/*XW: for a supervoxel, bandMin records the starting position of the sinogram band at each view*/
	/*XW: for a supervoxel, bandMax records the end position of the sinogram band at each view */
	channel_t bandMin[sinoparams.NViews]__attribute__((aligned(32)));
	channel_t bandMax[sinoparams.NViews]__attribute__((aligned(32)));
	channel_t bandWidthTemp[sinoparams.NViews]__attribute__((aligned(32)));
	channel_t bandWidth[NViewSets]__attribute__((aligned(32)));

	#ifdef ICC
	_intel_fast_memcpy(&bandMin[0],&bandMinMap[SVPosition].bandMin[0],sizeof(channel_t)*(sinoparams.NViews));
	_intel_fast_memcpy(&bandMax[0],&bandMaxMap[SVPosition].bandMax[0],sizeof(channel_t)*(sinoparams.NViews));
	#else
	memcpy(&bandMin[0],&bandMinMap[SVPosition].bandMin[0],sizeof(channel_t)*(sinoparams.NViews));
	memcpy(&bandMax[0],&bandMaxMap[SVPosition].bandMax[0],sizeof(channel_t)*(sinoparams.NViews));
	#endif

	#pragma vector aligned 
	for(p=0;p< sinoparams.NViews;p++)
		bandWidthTemp[p]=bandMax[p]-bandMin[p];

	for (p = 0; p < NViewSets; p++)
	{
		int bandWidthMax=bandWidthTemp[p*pieceLength];
		for(t=0;t<pieceLength;t++){
			if(bandWidthTemp[p*pieceLength+t]>bandWidthMax)
				bandWidthMax=bandWidthTemp[p*pieceLength+t];
		}
		bandWidth[p]=bandWidthMax;
	}

	float ** newWArray = (float **)malloc(sizeof(float *) * NViewSets);
	float ** newEArray = (float **)malloc(sizeof(float *) * NViewSets);
	float ** CopyNewEArray = (float **)malloc(sizeof(float *) * NViewSets);

	for (p = 0; p < NViewSets; p++)
	{
		newWArray[p] = (float *)malloc(sizeof(float)*bandWidth[p]*pieceLength*SV_depth_modified);
		newEArray[p] = (float *)malloc(sizeof(float)*bandWidth[p]*pieceLength*SV_depth_modified);
		CopyNewEArray[p] = (float *)malloc(sizeof(float)*bandWidth[p]*pieceLength*SV_depth_modified);
	}

	float *newWArrayPointer=&newWArray[0][0];
	float *newEArrayPointer=&newEArray[0][0];

	/*XW: copy the interlaced we into the memory buffer*/ 
	for (p = 0; p < NViewSets; p++)
	{
		newWArrayPointer=&newWArray[p][0];
		newEArrayPointer=&newEArray[p][0];
		for(i=0;i<SV_depth_modified;i++)
		for(q=0;q<pieceLength;q++) 
		{
			#ifdef ICC
			_intel_fast_memcpy(newWArrayPointer,&w[startSlice+i][p*pieceLength*sinoparams.NChannels+q*sinoparams.NChannels+bandMin[p*pieceLength+q]],sizeof(float)*(bandWidth[p]));
			_intel_fast_memcpy(newEArrayPointer,&e[startSlice+i][p*pieceLength*sinoparams.NChannels+q*sinoparams.NChannels+bandMin[p*pieceLength+q]],sizeof(float)*(bandWidth[p]));
			#else
			memcpy(newWArrayPointer,&w[startSlice+i][p*pieceLength*sinoparams.NChannels+q*sinoparams.NChannels+bandMin[p*pieceLength+q]],sizeof(float)*(bandWidth[p]));
			memcpy(newEArrayPointer,&e[startSlice+i][p*pieceLength*sinoparams.NChannels+q*sinoparams.NChannels+bandMin[p*pieceLength+q]],sizeof(float)*(bandWidth[p]));
			#endif
			newWArrayPointer+=bandWidth[p];
			newEArrayPointer+=bandWidth[p];
		}
	}

	for (p = 0; p < NViewSets; p++)
	{
		#ifdef ICC
		_intel_fast_memcpy(&CopyNewEArray[p][0],&newEArray[p][0],sizeof(float)*bandWidth[p]*pieceLength*SV_depth_modified);
		#else
		memcpy(&CopyNewEArray[p][0],&newEArray[p][0],sizeof(float)*bandWidth[p]*pieceLength*SV_depth_modified);
		#endif
	}

	float ** newWArrayTransposed = (float **)malloc(sizeof(float *) * NViewSets);
	float ** newEArrayTransposed = (float **)malloc(sizeof(float *) * NViewSets);

	for (p = 0; p < NViewSets; p++)
	{
		newWArrayTransposed[p] = (float *)malloc(sizeof(float)*bandWidth[p]*pieceLength*SV_depth_modified);
		newEArrayTransposed[p] = (float *)malloc(sizeof(float)*bandWidth[p]*pieceLength*SV_depth_modified);
	}

	float *WTransposeArrayPointer=&newWArrayTransposed[0][0];
	float *ETransposeArrayPointer=&newEArrayTransposed[0][0];
	
	for (p = 0; p < NViewSets; p++)
	for(currentSlice=0;currentSlice<(SV_depth_modified);currentSlice++) 
	{
		WTransposeArrayPointer=&newWArrayTransposed[p][currentSlice*bandWidth[p]*pieceLength];
		ETransposeArrayPointer=&newEArrayTransposed[p][currentSlice*bandWidth[p]*pieceLength];
		newEArrayPointer=&newEArray[p][currentSlice*bandWidth[p]*pieceLength];
		newWArrayPointer=&newWArray[p][currentSlice*bandWidth[p]*pieceLength];
		for(q=0;q<bandWidth[p];q++)
		{
			#pragma vector aligned 
			for(t=0;t<pieceLength;t++)
			{
				ETransposeArrayPointer[q*pieceLength+t]=newEArrayPointer[bandWidth[p]*t+q];
				WTransposeArrayPointer[q*pieceLength+t]=newWArrayPointer[bandWidth[p]*t+q];
			}
		}
	}

	WTransposeArrayPointer=&newWArrayTransposed[0][0];
	ETransposeArrayPointer=&newEArrayTransposed[0][0];
	newEArrayPointer=&newEArray[0][0];

	for (p = 0; p < NViewSets; p++)
		free((void *)newWArray[p]);

	free((void **)newWArray);

	/* Turn off zero-skipping for 1st iteration */
	char zero_skip_enable=0;  // 1: enable, 0: disable
	if(iter>0 && PositivityFlag)
		zero_skip_enable=1;

	/*XW: the start of the loop to compute theta1, theta2*/
	for(i=0;i<countNumber;i++)
	{
		const short j_new = j_newCoordinate[i];   /*XW: get the voxel's x,y location*/
		const short k_new = k_newCoordinate[i];
		float Aval_max = Aval_max_ptr[j_new*Nx+k_new];
		float tempV[SV_depth_modified];
		float tempProxMap[SV_depth_modified];
		float neighbors[SV_depth_modified][10];
		char zero_skip_FLAG[SV_depth_modified];
		float diff[SV_depth_modified];

		float THETA1[SV_depth_modified];
		float THETA2[SV_depth_modified];
		memset(&THETA1[0],0.0, sizeof(THETA1));
		memset(&THETA2[0],0.0, sizeof(THETA2));	

		int theVoxelPosition=(j_new-jy)*(2*SVLength+1)+(k_new-jx); 
		unsigned char * A_padd_Tranpose_pointer = &A_Padded_Map[SVPosition][theVoxelPosition].val[0];

		for(currentSlice=0;currentSlice<SV_depth_modified;currentSlice++)
		{
			tempV[currentSlice] = (float)(image[startSlice+currentSlice][j_new*Nx+k_new]); /*XW: current voxel's value*/

			zero_skip_FLAG[currentSlice] = 0;

			if(reconparams.ReconType == MBIR_MODULAR_RECONTYPE_QGGMRF_3D)
			{
				ExtractNeighbors3D(&neighbors[currentSlice][0],k_new,j_new,&image[startSlice+currentSlice][0],imgparams);

				if((startSlice+currentSlice)==0)
					neighbors[currentSlice][8]=voxelsBuffer1[j_new*Nx+k_new];
				else
					neighbors[currentSlice][8]=image[startSlice+currentSlice-1][j_new*Nx+k_new];

				if((startSlice+currentSlice)<(Nz-1))
					neighbors[currentSlice][9]=image[startSlice+currentSlice+1][j_new*Nx+k_new];
				else
					neighbors[currentSlice][9]=voxelsBuffer2[j_new*Nx+k_new];

				if(zero_skip_enable)
				if(tempV[currentSlice] == 0.0)
				{
					zero_skip_FLAG[currentSlice] = 1;
					for (j = 0; j < 10; j++)
					{
						if (neighbors[currentSlice][j] != 0.0)
						{
							zero_skip_FLAG[currentSlice] = 0;
							break; 
						}
					}
				}
			}
			if(reconparams.ReconType == MBIR_MODULAR_RECONTYPE_PandP)
				tempProxMap[currentSlice] = proximalmap[startSlice+currentSlice][j_new*Nx+k_new];
		}

		A_padd_Tranpose_pointer = &A_Padded_Map[SVPosition][theVoxelPosition].val[0];
		for(p=0;p<NViewSets;p++)
		{
			const int myCount=A_Padded_Map[SVPosition][theVoxelPosition].pieceWiseWidth[p];
			const int pieceMin=A_Padded_Map[SVPosition][theVoxelPosition].pieceWiseMin[p];
			#pragma vector aligned
			for(currentSlice=0;currentSlice<SV_depth_modified;currentSlice++)
			if(zero_skip_FLAG[currentSlice] == 0)
			{
				WTransposeArrayPointer=&newWArrayTransposed[p][currentSlice*bandWidth[p]*pieceLength];
				ETransposeArrayPointer=&newEArrayTransposed[p][currentSlice*bandWidth[p]*pieceLength];
				WTransposeArrayPointer+=pieceMin*pieceLength;
				ETransposeArrayPointer+=pieceMin*pieceLength;
				float tempTHETA1=0.0;
				float tempTHETA2=0.0;
				//Not finding evidence this makes a difference --SJK
				//Deprecated by Intel anyway
				//#pragma vector aligned
				//#pragma simd reduction(+:tempTHETA2,tempTHETA1)
				for(t=0;t<myCount*pieceLength;t++)
				{	/* summing over voxels which are not skipped or masked*/     
					tempTHETA1 += A_padd_Tranpose_pointer[t]*WTransposeArrayPointer[t]*ETransposeArrayPointer[t];
					tempTHETA2 += A_padd_Tranpose_pointer[t]*WTransposeArrayPointer[t]*A_padd_Tranpose_pointer[t];
				}
				THETA1[currentSlice]+=tempTHETA1;
				THETA2[currentSlice]+=tempTHETA2;
			}
			A_padd_Tranpose_pointer += myCount*pieceLength;
		}
		for(currentSlice=0;currentSlice<SV_depth_modified;currentSlice++)
		{
			THETA1[currentSlice]=-THETA1[currentSlice]*Aval_max*(1.0/255);
			THETA2[currentSlice]=THETA2[currentSlice]*Aval_max*(1.0/255)*Aval_max*(1.0/255);
		}

		ETransposeArrayPointer=&newEArrayTransposed[0][0];

		A_padd_Tranpose_pointer = &A_Padded_Map[SVPosition][theVoxelPosition].val[0];
	
		for(currentSlice=0;currentSlice<SV_depth_modified;currentSlice++)
		if(zero_skip_FLAG[currentSlice] == 0)
		{
			float pixel,step;
			if(reconparams.ReconType == MBIR_MODULAR_RECONTYPE_QGGMRF_3D)
			{
				step = QGGMRF3D_Update(reconparams,tempV[currentSlice],&neighbors[currentSlice][0],THETA1[currentSlice],THETA2[currentSlice]);
			}
			else if(reconparams.ReconType == MBIR_MODULAR_RECONTYPE_PandP)
			{
				step = PandP_Update(reconparams,tempV[currentSlice],tempProxMap[currentSlice],THETA1[currentSlice],THETA2[currentSlice]);
			}
			else
			{
				fprintf(stderr,"Error** Unrecognized ReconType in ICD update\n");
				exit(-1);
			}

			pixel = tempV[currentSlice] + step;  /* can apply over-relaxation to the step size here */

			if(PositivityFlag)
				image[startSlice+currentSlice][j_new*Nx+k_new] = ((pixel < 0.0) ? 0.0 : pixel);
			else
				image[startSlice+currentSlice][j_new*Nx+k_new] = pixel;

			diff[currentSlice] = image[startSlice+currentSlice][j_new*Nx+k_new] - tempV[currentSlice];

			totalChange_loc += fabs(diff[currentSlice]);
			totalValue_loc += fabs(tempV[currentSlice]);
			NumUpdates_loc++;

			diff[currentSlice]=diff[currentSlice]*Aval_max*(1.0/255);
		}

		for(p=0;p<NViewSets;p++)
		{
			const int myCount=A_Padded_Map[SVPosition][theVoxelPosition].pieceWiseWidth[p];
			const int pieceMin=A_Padded_Map[SVPosition][theVoxelPosition].pieceWiseMin[p]; 
			#pragma vector aligned
			for(currentSlice=0;currentSlice<SV_depth_modified;currentSlice++)
			if(diff[currentSlice]!=0 && zero_skip_FLAG[currentSlice] == 0)
			{
				ETransposeArrayPointer=&newEArrayTransposed[p][currentSlice*bandWidth[p]*pieceLength];
				ETransposeArrayPointer+=pieceMin*pieceLength;

				#pragma vector aligned
				for(t=0;t<(myCount*pieceLength);t++)
					ETransposeArrayPointer[t]= ETransposeArrayPointer[t]-A_padd_Tranpose_pointer[t]*diff[currentSlice];
			}
			A_padd_Tranpose_pointer+=myCount*pieceLength;
		}
	}

	for (p = 0; p < NViewSets; p++)
		free((void *)newWArrayTransposed[p]);

	free((void **)newWArrayTransposed);

	headNodeArray[jj_new].x=totalChange_loc;

	for (p = 0; p < NViewSets; p++)
	for(currentSlice=0;currentSlice<SV_depth_modified;currentSlice++)
	{
		ETransposeArrayPointer=&newEArrayTransposed[p][currentSlice*bandWidth[p]*pieceLength];
		newEArrayPointer=&newEArray[p][currentSlice*bandWidth[p]*pieceLength]; 
		for(q=0;q<bandWidth[p];q++)
		{
			#pragma vector aligned
			for(t=0;t<pieceLength;t++)
				newEArrayPointer[bandWidth[p]*t+q]=ETransposeArrayPointer[q*pieceLength+t];
		}
	}

	for (p = 0; p < NViewSets; p++)
		free((void *)newEArrayTransposed[p]);

	free((void **)newEArrayTransposed);

	newEArrayPointer=&newEArray[0][0];
	float* CopyNewEArrayPointer=&CopyNewEArray[0][0];
	float* eArrayPointer=&e[0][0];

	for (p = 0; p < NViewSets; p++)      /*XW: update the error term in the memory buffer*/
	{
		newEArrayPointer=&newEArray[p][0];
		CopyNewEArrayPointer=&CopyNewEArray[p][0];
		for (currentSlice=0; currentSlice< SV_depth_modified;currentSlice++)
		{
			#pragma vector aligned
			for(q=0;q<pieceLength;q++)
			{
				eArrayPointer=&e[startSlice+currentSlice][p*pieceLength*sinoparams.NChannels+q*sinoparams.NChannels+bandMin[p*pieceLength+q]];
				for(t=0;t<bandWidth[p];t++)
				{
					#pragma omp atomic
					*eArrayPointer += (*newEArrayPointer)-(*CopyNewEArrayPointer); 
					newEArrayPointer++;
					CopyNewEArrayPointer++;
					eArrayPointer++;
				}
			}
		}
	}

	for (p = 0; p < NViewSets; p++)
	{
		free((void *)newEArray[p]);
		free((void *)CopyNewEArray[p]);
	}
	free((void **)newEArray);
	free((void **)CopyNewEArray);

	*NumUpdates += NumUpdates_loc;
	*totalValue += totalValue_loc;
	*totalChange += totalChange_loc;

}   /* END super_voxel_recon() */




void coordinateShuffle(int *order1, int *order2,int len)
{
	int i, j, tmp1,tmp2;

	for (i = 0; i < len-1; i++)
	{
		j = i + (rand() % (len-i));
		tmp1 = order1[j];
		tmp2 = order2[j];
		order1[j] = order1[i];
		order2[j] = order2[i];
		order1[i] = tmp1;
		order2[i] = tmp2;
	}
}

void three_way_shuffle(long *order1, char *order2, struct heap_node *headNodeArray, int len)
{
	int i,j;
	long tmp_long;
	char tmp_char;
	float temp_x;

	for (i = 0; i < len-1; i++)
	{
		j = i + (rand() % (len-i));
		tmp_long = order1[j];
		order1[j] = order1[i];
		order1[i] = tmp_long;
		tmp_char = order2[j];
		order2[j] = order2[i];
		order2[i] = tmp_char;
		temp_x=headNodeArray[j].x;
		headNodeArray[j].x=headNodeArray[i].x;
		headNodeArray[i].x=temp_x;
	}
}



float MAPCostFunction3D(float **e,struct Image3D *Image,struct Sino3DParallel *sinogram,struct ReconParams *reconparams)
{
	int i, M, j, jx, jy, jz, Nx, Ny, Nz, plusx, minusx, plusy, plusz ;
    	float **x ;
    	float **w ;
	float nloglike, nlogprior_nearest, nlogprior_diag, nlogprior_interslice ;
    
    	x = Image->image;
    	w = sinogram->weight;
    
	M = sinogram->sinoparams.NViews * sinogram->sinoparams.NChannels ;
	Nx = Image->imgparams.Nx;
	Ny = Image->imgparams.Ny;
    	Nz = Image->imgparams.Nz;
    
	nloglike = 0.0;
	for (i = 0; i <sinogram->sinoparams.NSlices; i++){
		for (j = 0; j < M; j++)
    			nloglike += e[i][j]*w[i][j]*e[i][j];
    	}

	nloglike /= 2.0;
	nlogprior_nearest = 0.0;
	nlogprior_diag = 0.0;
    	nlogprior_interslice = 0.0;
    
	for (jz = 0; jz < Nz; jz++)
	for (jy = 0; jy < Ny; jy++)
	for (jx = 0; jx < Nx; jx++)
	{
		plusx = jx + 1;
		plusx = ((plusx < Nx) ? plusx : 0);
		minusx = jx - 1;
		minusx = ((minusx < 0) ? Nx-1 : minusx);
		plusy = jy + 1;
		plusy = ((plusy < Ny) ? plusy : 0);
		plusz = jz + 1;
		plusz = ((plusz < Nz) ? plusz : 0);

		j = jy*Nx + jx; 

		nlogprior_nearest += QGGMRF_Potential((x[jz][j] - x[jz][jy*Nx+plusx]), reconparams);
		nlogprior_nearest += QGGMRF_Potential((x[jz][j] - x[jz][plusy*Nx+jx]),reconparams);

		nlogprior_diag += QGGMRF_Potential((x[jz][j] - x[jz][plusy*Nx+minusx]),reconparams);
		nlogprior_diag += QGGMRF_Potential((x[jz][j] - x[jz][plusy*Nx+plusx]),reconparams);

		nlogprior_interslice += QGGMRF_Potential((x[jz][j] - x[plusz][jy*Nx+jx]),reconparams);
	}

	return (nloglike + reconparams->b_nearest * nlogprior_nearest + reconparams->b_diag * nlogprior_diag + reconparams->b_interslice * nlogprior_interslice) ;
}




void forwardProject(
    float *image,
    float *proj,
    struct ImageParams3D imgparams,
    struct SinoParams3DParallel sinoparams,
    char *Amatrix_fname,
    char verboseLevel)
{
    struct AValues_char **A_Padded_Map;
    float *Aval_max_ptr;
    struct SVParams svpar;

    size_t i;
    int j,jz;
    int Nx = imgparams.Nx;
    int Ny = imgparams.Ny;
    int Nz = imgparams.Nz;
    int NChannels = sinoparams.NChannels;
    int M = sinoparams.NViews * sinoparams.NChannels;

    /* Initialize/allocate SV parameters */
    initSVParams(&svpar, imgparams, sinoparams);
    int Nsv = svpar.Nsv;
    int SVLength = svpar.SVLength;
    int pieceLength = svpar.pieceLength;
    int SVsPerRow = svpar.SVsPerRow;
    int NViewSets = sinoparams.NViews/pieceLength;
    struct minStruct * bandMinMap = svpar.bandMinMap;

    /* Allocate and generate recon mask based on ROIRadius */
    char * ImageReconMask = GenImageReconMask(&imgparams);

    /* Read/compute/write System Matrix */
    A_Padded_Map = (struct AValues_char **)multialloc(sizeof(struct AValues_char),2,Nsv,(2*SVLength+1)*(2*SVLength+1));
    Aval_max_ptr = (float *) mget_spc(Nx*Ny,sizeof(float));
    if(Amatrix_fname != NULL)
    {
        if(verboseLevel)
            fprintf(stdout,"Reading system matrix...\n");
        readAmatrix(Amatrix_fname, A_Padded_Map, Aval_max_ptr, &imgparams, &sinoparams, svpar);
    }
    else
    {
        if(verboseLevel)
            fprintf(stdout,"Computing system matrix...\n");
        A_comp(A_Padded_Map,Aval_max_ptr,svpar,&sinoparams,ImageReconMask,&imgparams);
    }

    /* initialize projection */
    for (i = 0; i < M*Nz; i++)
        proj[i] = 0.0;

    #pragma omp parallel for schedule(dynamic)
    for(jz=0;jz<Nz;jz++)
    {
	    int jx,jy,k,r,p;

        for (jy = 0; jy < Ny; jy++)
        for (jx = 0; jx < Nx; jx++)
        {
            int SV_ind_y = jy/(2*SVLength-svpar.overlap);
            int SV_ind_x = jx/(2*SVLength-svpar.overlap);
            int SVPosition = SV_ind_y*SVsPerRow + SV_ind_x;

            int SV_jy = SV_ind_y*(2*SVLength-svpar.overlap);
            int SV_jx = SV_ind_x*(2*SVLength-svpar.overlap);
            int VoxelPosition = (jy-SV_jy)*(2*SVLength+1)+(jx-SV_jx);

            // The second condition should always be true
            if (A_Padded_Map[SVPosition][VoxelPosition].length > 0 && VoxelPosition < ((2*SVLength+1)*(2*SVLength+1)))
            {
                unsigned char* A_padd_Tr_ptr = &A_Padded_Map[SVPosition][VoxelPosition].val[0];
                float rescale = Aval_max_ptr[jy*Nx+jx]*(1.0/255);
                float xval = image[jz*Nx*Ny + jy*Nx + jx];

                for(p=0;p<NViewSets;p++)
                {
                    int myCount = A_Padded_Map[SVPosition][VoxelPosition].pieceWiseWidth[p];
                    int pieceWiseMin = A_Padded_Map[SVPosition][VoxelPosition].pieceWiseMin[p];
                    int position = p*pieceLength*NChannels + pieceWiseMin;

                    for(r=0;r<myCount;r++)
                    for(k=0;k<pieceLength;k++)
                    {
                        int bandMin = bandMinMap[SVPosition].bandMin[p*pieceLength+k];
                        int proj_idx = jz*M + position + k*NChannels + bandMin + r;

                        if((pieceWiseMin + bandMin + r) >= NChannels || (position + k*NChannels + bandMin + r) >= M ) {
                            fprintf(stderr,"forwardProject() out of bounds: p %d r %d k %d\n",p,r,k);
                            fprintf(stderr,"forwardProject() out of bounds: total_1 %d total_2 %d\n",pieceWiseMin+bandMin+r,position+k*NChannels+bandMin+r);
                            exit(-1);
                        }
                        else
                            proj[proj_idx] += A_padd_Tr_ptr[r*pieceLength+k]*rescale * xval;
                    }
                    A_padd_Tr_ptr += myCount*pieceLength;
                }
            }
        }
    }

    /* Free SV memory */
    for(i=0;i<Nsv;i++) {
        free((void *)svpar.bandMinMap[i].bandMin);
        free((void *)svpar.bandMaxMap[i].bandMax);
    }
    free((void *)svpar.bandMinMap);
    free((void *)svpar.bandMaxMap);

    /* Free system matrix */
    for(i=0;i<Nsv;i++)
    for(j=0;j<(2*SVLength+1)*(2*SVLength+1);j++)
    if(A_Padded_Map[i][j].length>0)
    {
        free((void *)A_Padded_Map[i][j].val);
        free((void *)A_Padded_Map[i][j].pieceWiseMin);
        free((void *)A_Padded_Map[i][j].pieceWiseWidth);
    }
    multifree(A_Padded_Map,2);
    free((void *)Aval_max_ptr);
    free((void *)ImageReconMask);

}   /* END forwardProject() */


