
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

#define  c_ratio 0.07
#define convergence_rho 0.7

/* Internal functions */
void super_voxel_recon(int jj,struct SVParams svpar,unsigned long *NumUpdates,float *totalValue,float *totalChange,int iter,
	char *phaseMap,long *order,int *indexList,float **w,float **e,
	struct AValues_char ** A_Padded_Map,float *max_num_pointer,struct heap_node *headNodeArray,
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
	float *max_num_pointer,
	char *ImageReconMask,
	struct CmdLine *cmdline)
{
	int i,j,jj,p,t,iter,it_print=1;
	int NumMaskVoxels=0;
	//float **x;  /* image data */
	//float **y;  /* sinogram projections data */
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
	//struct tm1,tm2;

	//x = Image->image;   /* x is the image vector */
	//y = sinogram->sino;  /* y is the sinogram vector */
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
	//int pieceLength = svpar.pieceLength;
	//struct minStruct * bandMinMap = svpar.bandMinMap;
	//struct maxStruct * bandMaxMap = svpar.bandMaxMap;

	int rep_num=(int)ceil(1/(4*c_ratio*convergence_rho));

        for(j=0;j<Nxy;j++)
        if(ImageReconMask[j])
                NumMaskVoxels++;

	#if 0
	//#ifdef find_RMSE
	float **golden;
	golden = (float **)multialloc(sizeof(float), 2, Nz,Nxy);
	// note this needs to be updated
	read_golden(cmdline->ReconImageDataFile,golden,Nz,Nxy,Image);
	float updatedVoxelsList[300];

	float sumOfSE=0;
	for(i=0;i<Nz;i++)
	for(j=0;j<Nxy;j++)
		sumOfSE+=(Image->image[i][j]-golden[i][j])*(Image->image[i][j]-golden[i][j]);

	float MSE=sumOfSE/Nxy/Nz;
	float RMSE=sqrt(MSE);

	fprintf(stdout,"Rho: %f initial_RMSE: %f \n",convergence_rho,RMSE);
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

	//struct heap_node headNodeArray[sum*SV_per_Z];  //This allocation crashed silently for a large problem
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

	for(i=0;i<Nxy;i++) voxelsBuffer1[i]=0;
	for(i=0;i<Nxy;i++) voxelsBuffer2[i]=0;

	iter=0;

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

	int startIndex=0;
	int endIndex=0;        		

	//gettimeofday(&tm1,NULL);
         
	char stop_FLAG=0;

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
						initialize_heap(&priorityheap);						
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
					super_voxel_recon(jj,svpar,&NumUpdates,&totalValue,&totalChange,iter,&phaseMap[0],order,&indexList[0],w,e,A_Padded_Map,&max_num_pointer[0],&headNodeArray[0],sinogram->sinoparams,reconparams,Image,voxelsBuffer1,voxelsBuffer2,&group_id_list[0][0],group);

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

			#if 0
			//#ifdef find_RMSE
			float sumOfSE=0;
			for(i=0;i<Nz;i++)
			for(j=0;j<Nxy;j++)
				sumOfSE+=(Image->image[i][j]-golden[i][j])*(Image->image[i][j]-golden[i][j]);
			float MSE=sumOfSE/Nxy/Nz;
			float RMSE=sqrt(MSE);
			if(iter<300) {
				updatedVoxelsList[iter]=NumUpdates*1.0/Nxy/Nz;
				fprintf(stdout,"Rho: %f Equits: %f RMSE: %f \n",convergence_rho,updatedVoxelsList[iter],RMSE);
			}
			#endif		

			if (avg_update_rel < StopThreshold && (endIndex!=0))
				stop_FLAG = 1;

			iter++;
			equits += (float)NumUpdates/((float)NumMaskVoxels*Nz);

			if(cmdline->verboseLevel)
			if(equits > it_print)
			{
				fprintf(stdout,"\titeration %d, average change %.4f %%\n",it_print,avg_update_rel);
				it_print++;
			}

			NumUpdates=0;
			totalValue=0;
			totalChange=0;

			}

		}
	}
	
        //gettimeofday(&tm2,NULL);
        //unsigned long long tt = 1000 * (tm2.tv_sec - tm1.tv_sec) + (tm2.tv_usec - tm1.tv_usec) / 1000;
        //printf("\trun time %llu ms (iterations only)\n", tt);

	if(cmdline->verboseLevel)
	{
		if(StopThreshold <= 0)
			fprintf(stdout,"\tNo stopping condition--running fixed iterations\n");
		else if(stop_FLAG == 1)
			fprintf(stdout,"\tReached stopping condition\n");
		else
			fprintf(stdout,"\tWARNING: Didn't reach stopping condition\n");

		fprintf(stdout,"\tEquivalent iterations = %.1f, (non-homogeneous iterations = %d)\n",equits,iter);
		fprintf(stdout,"\tAverage update in last iteration (relative) = %f %%\n",avg_update_rel);
		fprintf(stdout,"\tAverage update in last iteration (magnitude) = %.4g\n",avg_update);
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
	if(priorityheap.size>0)
		free_heap((void *)&priorityheap); 
	free((void *)phaseMap);
	free((void *)order);

	#ifdef find_RMSE
	multifree(golden,2);
	#endif

}   /*  END MBIRReconstruct3D()  */


				

void forwardProject2D(
	float *e,
	float *x,
	struct AValues_char ** A_Padded_Map,
	float *max_num_pointer,
	struct SinoParams3DParallel *sinoparams,
	struct ImageParams3D *imgparams,
	struct SVParams svpar)
{
	int jx,jy,Nx,Ny,i,M,r,j,p;
	float inverseNumber=1.0/255;
	int SVLength = svpar.SVLength;
	int overlappingDistance = svpar.overlap;
	struct minStruct * bandMinMap = svpar.bandMinMap;
	int pieceLength = svpar.pieceLength;
	int SVsPerRow = svpar.SVsPerRow;

	const int NViewsdivided=(sinoparams->NViews)/pieceLength;

	Nx = imgparams->Nx;
	Ny = imgparams->Ny;
	M = sinoparams->NViews*sinoparams->NChannels;

	for (i = 0; i < M; i++)
		e[i] = 0.0;

	for (jy = 0; jy < Ny; jy++)
	for (jx = 0; jx < Nx; jx++)
	{
		int SV_ind_x = jx/(2*SVLength-overlappingDistance);
		int SV_ind_y = jy/(2*SVLength-overlappingDistance);
		int SVPosition = SV_ind_y*SVsPerRow + SV_ind_x;

		int SV_jy = SV_ind_y*(2*SVLength-overlappingDistance);
		int SV_jx = SV_ind_x*(2*SVLength-overlappingDistance);
		int VoxelPosition = (jy-SV_jy)*(2*SVLength+1)+(jx-SV_jx);
		/*
		fprintf(stdout,"jy %d jx %d SVPosition %d SV_jy %d SV_jx %d VoxelPosition %d \n",jy,jx,SVPosition,SV_jy,SV_jx,VoxelPosition);
		*/
		// I think the second condition will always be true
		if (A_Padded_Map[SVPosition][VoxelPosition].length > 0 && VoxelPosition < ((2*SVLength+1)*(2*SVLength+1)))
		{
			/*XW: remove the index field in struct ACol and exploit the spatial locality */
			unsigned char* A_padd_Tranpose_pointer = &A_Padded_Map[SVPosition][VoxelPosition].val[0];
			for(p=0;p<NViewsdivided;p++) 
			{
				const int myCount=A_Padded_Map[SVPosition][VoxelPosition].pieceWiseWidth[p];
				int position=p*pieceLength*sinoparams->NChannels+A_Padded_Map[SVPosition][VoxelPosition].pieceWiseMin[p];

				for(r=0;r<myCount;r++)
				for(j=0;j< pieceLength;j++)
				{
					if((A_Padded_Map[SVPosition][VoxelPosition].pieceWiseMin[p]+bandMinMap[SVPosition].bandMin[p*pieceLength+j]+r)>=sinoparams->NChannels)
						fprintf(stdout, "p %d r %d j %d total_1 %d \n",p,r,j,A_Padded_Map[SVPosition][VoxelPosition].pieceWiseMin[p]+bandMinMap[SVPosition].bandMin[p*pieceLength+j]+r);

					if((position+j*sinoparams->NChannels+bandMinMap[SVPosition].bandMin[p*pieceLength+j]+r)>= M)
						fprintf(stdout, "p %d r %d j %d total_2 %d \n",p,r,j,position+j*sinoparams->NChannels+bandMinMap[SVPosition].bandMin[p*pieceLength+j]+r);

					if((position+j*sinoparams->NChannels+bandMinMap[SVPosition].bandMin[p*pieceLength+j]+r)< M)
						e[position+j*sinoparams->NChannels+bandMinMap[SVPosition].bandMin[p*pieceLength+j]+r] += A_padd_Tranpose_pointer[r*pieceLength+j]*max_num_pointer[jy*Nx+jx]*inverseNumber*x[jy*Nx+jx];

				}
				A_padd_Tranpose_pointer+=myCount*pieceLength;
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
	float *max_num_pointer,
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
	int NViewsdivided = sinoparams.NViews/pieceLength;

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
	int bandMin[sinoparams.NViews]__attribute__((aligned(32)));
	int bandMax[sinoparams.NViews]__attribute__((aligned(32)));
	int bandWidthTemp[sinoparams.NViews]__attribute__((aligned(32)));
	int bandWidth[NViewsdivided]__attribute__((aligned(32)));

	#ifdef ICC
	_intel_fast_memcpy(&bandMin[0],&bandMinMap[SVPosition].bandMin[0],sizeof(int)*(sinoparams.NViews));
	_intel_fast_memcpy(&bandMax[0],&bandMaxMap[SVPosition].bandMax[0],sizeof(int)*(sinoparams.NViews)); 
	#else
	memcpy(&bandMin[0],&bandMinMap[SVPosition].bandMin[0],sizeof(int)*(sinoparams.NViews));
	memcpy(&bandMax[0],&bandMaxMap[SVPosition].bandMax[0],sizeof(int)*(sinoparams.NViews));
	#endif

	#pragma vector aligned 
	for(p=0;p< sinoparams.NViews;p++)
		bandWidthTemp[p]=bandMax[p]-bandMin[p];

	for (p = 0; p < NViewsdivided; p++)
	{
		int bandWidthMax=bandWidthTemp[p*pieceLength];
		for(t=0;t<pieceLength;t++){
			if(bandWidthTemp[p*pieceLength+t]>bandWidthMax)
				bandWidthMax=bandWidthTemp[p*pieceLength+t];
		}
		bandWidth[p]=bandWidthMax;
	}

	float ** newWArray = (float **)malloc(sizeof(float *) * NViewsdivided);
	float ** newEArray = (float **)malloc(sizeof(float *) * NViewsdivided);
	float ** CopyNewEArray = (float **)malloc(sizeof(float *) * NViewsdivided); 

	for (p = 0; p < NViewsdivided; p++)
	{
		newWArray[p] = (float *)malloc(sizeof(float)*bandWidth[p]*pieceLength*SV_depth_modified);
		newEArray[p] = (float *)malloc(sizeof(float)*bandWidth[p]*pieceLength*SV_depth_modified);
		CopyNewEArray[p] = (float *)malloc(sizeof(float)*bandWidth[p]*pieceLength*SV_depth_modified);
	}

	float *newWArrayPointer=&newWArray[0][0];
	float *newEArrayPointer=&newEArray[0][0];

	/*XW: copy the interlaced we into the memory buffer*/ 
	for (p = 0; p < NViewsdivided; p++)
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

	for (p = 0; p < NViewsdivided; p++)
	{
		#ifdef ICC
		_intel_fast_memcpy(&CopyNewEArray[p][0],&newEArray[p][0],sizeof(float)*bandWidth[p]*pieceLength*SV_depth_modified);
		#else
		memcpy(&CopyNewEArray[p][0],&newEArray[p][0],sizeof(float)*bandWidth[p]*pieceLength*SV_depth_modified);
		#endif
	}

	float ** newWArrayTransposed = (float **)malloc(sizeof(float *) * NViewsdivided);
	float ** newEArrayTransposed = (float **)malloc(sizeof(float *) * NViewsdivided);

	for (p = 0; p < NViewsdivided; p++)
	{
		newWArrayTransposed[p] = (float *)malloc(sizeof(float)*bandWidth[p]*pieceLength*SV_depth_modified);
		newEArrayTransposed[p] = (float *)malloc(sizeof(float)*bandWidth[p]*pieceLength*SV_depth_modified);
	}

	float *WTransposeArrayPointer=&newWArrayTransposed[0][0];
	float *ETransposeArrayPointer=&newEArrayTransposed[0][0];
	
	for (p = 0; p < NViewsdivided; p++)
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
	float inverseNumber=1.0/255;

	for (p = 0; p < NViewsdivided; p++)
		free((void *)newWArray[p]);

	free((void **)newWArray);

	/* Turn off zero-skipping for 1st iteration */
	char zero_skip_enable=0;  // 1: enable, 0: disable
	if(iter>0 && PositivityFlag)
		zero_skip_enable=1;

	/*XW: the start of the loop to compute theta1, theta2*/
	for(i=0;i<countNumber;i++)
	{
		const short j_new=j_newCoordinate[i];   /*XW: get the voxel's x,y location*/
		const short k_new=k_newCoordinate[i];
		float tempV[SV_depth_modified];
		float tempProxMap[SV_depth_modified];
		float neighbors[SV_depth_modified][10];
		char zero_skip_FLAG[SV_depth_modified];
		float max=max_num_pointer[j_new*Nx+k_new];
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
		for(p=0;p<NViewsdivided;p++)
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
			A_padd_Tranpose_pointer+=myCount*pieceLength;
		}
		for(currentSlice=0;currentSlice<SV_depth_modified;currentSlice++)
		{
			THETA1[currentSlice]=-THETA1[currentSlice]*max*inverseNumber;
			THETA2[currentSlice]=THETA2[currentSlice]*max*inverseNumber*max*inverseNumber;
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

			diff[currentSlice]=diff[currentSlice]*max*inverseNumber;
		}

		for(p=0;p<NViewsdivided;p++)
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

	for (p = 0; p < NViewsdivided; p++)
		free((void *)newWArrayTransposed[p]);

	free((void **)newWArrayTransposed);

	headNodeArray[jj_new].x=totalChange_loc;

	for (p = 0; p < NViewsdivided; p++)
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

	for (p = 0; p < NViewsdivided; p++)
		free((void *)newEArrayTransposed[p]);

	free((void **)newEArrayTransposed);

	newEArrayPointer=&newEArray[0][0];
	float* CopyNewEArrayPointer=&CopyNewEArray[0][0];
	float* eArrayPointer=&e[0][0];

	for (p = 0; p < NViewsdivided; p++)      /*XW: update the error term in the memory buffer*/
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

	for (p = 0; p < NViewsdivided; p++)
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




#if 0
// NOTE this needs to be updated
void read_golden(char *fname,float **golden,int Nz,int N, struct Image3D *Image)
{
	FILE *fp;
	int i;
        char slicefname[1024];
        char *sliceindex;
	sliceindex= (char *)malloc(MBIR_MODULAR_MAX_NUMBER_OF_SLICE_DIGITS);
	for(i=0;i<Nz;i++){
		sprintf(sliceindex,"%.*d",MBIR_MODULAR_MAX_NUMBER_OF_SLICE_DIGITS,Image->imgparams.FirstSliceNumber+i);

		/* Obtain file name for the given slice */
		strcpy(slicefname,fname);
		strcat(slicefname,"_slice"); 
		strcat(slicefname,sliceindex); /* append slice index */
		strcat(slicefname,".2Dimgdata");
		if ((fp = fopen(slicefname, "r")) == NULL)
		{
			fprintf(stderr, "ERROR in read golden: can't open file %s.\n", slicefname);
			exit(-1);
		}

		fread(&golden[i][0],sizeof(float),N,fp);
		fclose(fp);
	}

}
#endif

#if 0
static __inline__ unsigned long long rdtsc()
{
   unsigned hi,lo;
   __asm __volatile__ ("rdtsc" : "=a"(lo),"=d"(hi));
   return ( (unsigned long long)lo)|( ((unsigned long long)hi)<<32 );
}
#endif




