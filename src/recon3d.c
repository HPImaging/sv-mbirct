
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
void super_voxel_recon(int jj,struct SVParams svpar,float *total_updates,int it,int *phaseMap,int *order,int *indexList,int Nx,int Ny,
	float **w,float **e, struct AValues_char ** A_Padded_Map,const int Nslices,struct heap_node *headNodeArray,const int NViewsdivided,
	struct SinoParams3DParallel sinoparams,struct ReconParamsQGGMRF3D reconparams,struct ImageParams3D imgparams,
	float *max_num_pointer,float **image,float *voxelsBuffer1,float *voxelsBuffer2,int* group_array,int group_id,
	int SV_per_Z,int SVsPerLine,long *updatedVoxels,float pow_sigmaX_p,float pow_sigmaX_q,float pow_T_qmp);
void coordinateShuffle(int *order1, int *order2,int len);
void three_way_shuffle(int *order1, int *order2,struct heap_node *headNodeArray,int len);
float MAPCostFunction3D(float **e,struct Image3D *Image,struct Sino3DParallel *sinogram,struct ReconParamsQGGMRF3D *reconparams);


void MBIRReconstruct3D(
	struct Image3D *Image,
	struct Sino3DParallel *sinogram,
	float **e,  /* e=y-Ax, error */
	struct ReconParamsQGGMRF3D reconparams,
	struct SVParams svpar,
	struct AValues_char ** A_Padded_Map,
	float *max_num_pointer,
	struct CmdLine *cmdline)
{
	int it,i,j,jj,p,t;
	float **x;  /* image data */
	float **y;  /* sinogram projections data */
	float **w;  /* projections weights data */
	float *voxelsBuffer1;  /* the first N entries are the voxel values.  */
	float *voxelsBuffer2;
	float cost, avg_update, total_updates=0, equits=0;

	struct heap priorityheap;
	initialize_heap(&priorityheap);
	int *order;
	struct timeval tm1,tm2;

	x = Image->image;   /* x is the image vector */
	y = sinogram->sino;  /* y is the sinogram vector */
	w = sinogram->weight; /* vector of weights for each sinogram measurement */
	int Nx = Image->imgparams.Nx;
	int Ny = Image->imgparams.Ny;
	int Nxy = Nx*Ny;
	int Nz = Image->imgparams.Nz;
	int NvNc = sinogram->sinoparams.NViews * sinogram->sinoparams.NChannels;
	int NViews = sinogram->sinoparams.NViews;
	int MaxIterations = reconparams.MaxIterations;
	float StopThreshold = reconparams.StopThreshold;
	int SVLength = svpar.SVLength;
	int overlappingDistance = svpar.overlap;
	int SV_depth = svpar.SVDepth;
	int sum = svpar.Nsv;
	int pieceLength = svpar.pieceLength;
	struct minStruct * bandMinMap = svpar.bandMinMap;
	struct maxStruct * bandMaxMap = svpar.bandMaxMap;

	int SV_per_Z=0;
	int rep_num=(int)ceil(1/(4*c_ratio*convergence_rho));

	float pow_sigmaX_p = pow(reconparams.SigmaX,reconparams.p);
	float pow_sigmaX_q = pow(reconparams.SigmaX,reconparams.q);
	float pow_T_qmp = pow(reconparams.T,reconparams.q - reconparams.p);

	if((Nz%SV_depth)==0)
		SV_per_Z=Nz/SV_depth;
	else
		SV_per_Z=Nz/SV_depth+1;

	if(NViews%pieceLength !=0)
	{
		fprintf(stderr, "ERROR in MBIRReconstruct3D: pieceLength should divide evenly into NViews\n");
		exit(-1);
	}

	const int NViewsdivided=NViews/pieceLength;

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

	/* Order of pixel updates need NOT be raster order, just initialize */
	order = (int *)_mm_malloc(sum*SV_per_Z*sizeof(int),64);

	t=0;

	for(p=0;p<Nz;p+=SV_depth)
	for(i=0;i<Ny;i+=(SVLength*2-overlappingDistance))
	for(j=0;j<Nx;j+=(SVLength*2-overlappingDistance))
	{
		order[t]=p*Nxy+i*Nx+j;  /* order is the first voxel coordinate, not the center */
		t++;
	}

	int phaseMap[sum*SV_per_Z];
	int SVsPerLine=0;
	if((Nx%(2*SVLength-overlappingDistance))==0)
		SVsPerLine=Nx/(2*SVLength-overlappingDistance);
	else	
		SVsPerLine=Nx/(2*SVLength-overlappingDistance)+1; 

	#pragma omp parallel for private(jj) schedule(dynamic)
	for(i=0;i<SV_per_Z;i++)
	for(jj=0;jj<sum;jj++)
	{
		if((jj/SVsPerLine)%2==0)
		{
			if((jj%SVsPerLine)%2==0)
				phaseMap[i*sum+jj]=0;
			else
				phaseMap[i*sum+jj]=1;			
		}
		else
		{
			if((jj%SVsPerLine)%2==0)
				phaseMap[i*sum+jj]=2;
			else
				phaseMap[i*sum+jj]=3;			
		}
	}

	int group_id_list[SV_per_Z][4];

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
	struct heap_node headNodeArray[sum*SV_per_Z];

	for(i=0;i<SV_per_Z;i++)
	for(jj=0;jj<sum;jj++)
	{
		headNodeArray[i*sum+jj].pt=i*sum+jj;
		headNodeArray[i*sum+jj].x=0.0;
	}
	int indexList_size=(int) sum*SV_per_Z*4*c_ratio*(1-convergence_rho);	
	int indexList[indexList_size];   	             	    
    

	voxelsBuffer1 = (float *)_mm_malloc(Nxy*sizeof(float),64);
	voxelsBuffer2 = (float *)_mm_malloc(Nxy*sizeof(float),64);

	for(i=0;i<Nxy;i++) voxelsBuffer1[i]=0;
	for(i=0;i<Nxy;i++) voxelsBuffer2[i]=0;

	it=0;

	long updatedVoxels=0;
	
	coordinateShuffle(&order[0],&phaseMap[0],sum*SV_per_Z);
	
	int startIndex=0;
	int endIndex=0;        		

	//gettimeofday(&tm1,NULL);
         
	char stop_FLAG=0;

	#pragma omp parallel
	{
		while(stop_FLAG==0 && equits<MaxIterations && it<20*MaxIterations)
		{
			#pragma omp single
			{		
				if(it==0)
				{
					startIndex=0;
					endIndex=sum*SV_per_Z;
				}	
				else
				{
					if((it-1)%(2*rep_num)==0 && it!=1)
						three_way_shuffle(&order[0],&phaseMap[0],&headNodeArray[0],sum*SV_per_Z);
			
					if(it%2==1)
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
						startIndex=((it-2)/2)%rep_num*sum*SV_per_Z/rep_num;
						endIndex=(((it-2)/2)%rep_num+1)*sum*SV_per_Z/rep_num;
					}
				}
			}

			int group=0;

			for (group = 0; group < 4; group++)
			{
			
				#pragma omp for schedule(dynamic)  reduction(+:total_updates)
				for (jj = startIndex; jj < endIndex; jj+=1)
					super_voxel_recon(jj,svpar,&total_updates,it, &phaseMap[0],order,&indexList[0],Nx,Ny, w,e,A_Padded_Map,Nz,&headNodeArray[0],NViewsdivided,sinogram->sinoparams,reconparams,Image->imgparams,&max_num_pointer[0],Image->image,voxelsBuffer1,voxelsBuffer2,&group_id_list[0][0],group,SV_per_Z,SVsPerLine,&updatedVoxels,pow_sigmaX_p,pow_sigmaX_q,pow_T_qmp);
			}

			#pragma omp single
			{
 
			if(it==0)
				avg_update = total_updates/(float)Nxy/Nz;
			else
			{
				if(it%2==1)
					avg_update = (total_updates/(float)Nxy)*(4*convergence_rho/(1-convergence_rho))/Nz;
				else
					avg_update = (total_updates/(float)Nxy)*4/Nz;
			}           		
			
			/*
			cost = MAPCostFunction3D(e, Image, sinogram, &reconparams);
			fprintf(stdout, "it %d cost = %-15f, avg_update %f \n", it, cost, avg_update);           	
			*/

			#if 0
			//#ifdef find_RMSE
			if(it<300)
				updatedVoxelsList[it]=updatedVoxels*1.0/Nxy/Nz;
	
			float sumOfSE=0;
			for(i=0;i<Nz;i++)
			for(j=0;j<Nxy;j++)
				sumOfSE+=(Image->image[i][j]-golden[i][j])*(Image->image[i][j]-golden[i][j]);

			float MSE=sumOfSE/Nxy/Nz;
			float RMSE=sqrt(MSE);

			fprintf(stdout,"Rho: %f Equits: %f RMSE: %f \n",convergence_rho,updatedVoxelsList[it],RMSE);

			#endif		

			//if (avg_update < StopThreshold && (endIndex!=0))
			//	stop_FLAG = 1;

			it++;
			equits += total_updates/(Nxy*Nz);
			//fprintf(stdout,"iteration %d, equits %f\n",it,equits);

			total_updates = 0.0;

			}

		}
	}
	
        //gettimeofday(&tm2,NULL);
        //unsigned long long tt = 1000 * (tm2.tv_sec - tm1.tv_sec) + (tm2.tv_usec - tm1.tv_usec) / 1000;
        //printf("\trun time %llu ms (iterations only)\n", tt);

	if (stop_FLAG == 0 && StopThreshold > 0)
	{
		fprintf(stdout, "\tWARNING: Didn't reach stopping condition.\n");
		fprintf(stdout, "\tAverage update magnitude = %f\n", avg_update);
	}
	else if (stop_FLAG == 0 && StopThreshold <= 0)
	{
		fprintf(stdout,"\tNo stopping condition.\n");
		fprintf(stdout,"\titerations %d, equits %.1f\n",it,equits);
		//fprintf(stdout,"\tAverage update magnitude = %f\n", avg_update);
	}
	//if(AvgVoxelValue>0)
	//fprintf(stdout, "Average Update to Average Voxel-Value Ratio = %f %% \n", ratio);
	
	_mm_free((void *)order);
	_mm_free((void *)voxelsBuffer1);
	_mm_free((void *)voxelsBuffer2);
	if(priorityheap.size>0)
		free_heap((void *)&priorityheap); 

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
	int jx,jy,Nx,Ny,i,M,r,j,p,SVNumPerRow;
	float inverseNumber=1.0/255;
	int SVLength = svpar.SVLength;
	int overlappingDistance = svpar.overlap;
	struct minStruct * bandMinMap = svpar.bandMinMap;
	int pieceLength = svpar.pieceLength;

	const int NViewsdivided=(sinoparams->NViews)/pieceLength;

	Nx = imgparams->Nx;
	Ny = imgparams->Ny;
	M = sinoparams->NViews*sinoparams->NChannels;

	for (i = 0; i < M; i++)
		e[i] = 0.0;

	if((Nx%(2*SVLength-overlappingDistance))==0)
		SVNumPerRow=Nx/(2*SVLength-overlappingDistance);
	else
		SVNumPerRow=Nx/(2*SVLength-overlappingDistance)+1;

	for (jy = 0; jy < Ny; jy++)
	for (jx = 0; jx < Nx; jx++)
	{
		int temp1=jy/(2*SVLength-overlappingDistance);
		if(temp1==SVNumPerRow)  // I don't think this will happen
			temp1=SVNumPerRow-1;

		int temp2=jx/(2*SVLength-overlappingDistance);
		if(temp2==SVNumPerRow)  // I don't think this will happen
			temp2=SVNumPerRow-1;

		int SVPosition=temp1*SVNumPerRow+temp2;
 
		int SV_jy=temp1*(2*SVLength-overlappingDistance);
		int SV_jx=temp2*(2*SVLength-overlappingDistance);
		int VoxelPosition=(jy-SV_jy)*(2*SVLength+1)+(jx-SV_jx);
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
	float *total_updates,
	int it,
	int *phaseMap,
	int *order,
	int *indexList,
	int Nx,
	int Ny,
	float **w,
	float **e,
	struct AValues_char ** A_Padded_Map,
	const int Nslices,
	struct heap_node *headNodeArray,
	const int NViewsdivided,
	struct SinoParams3DParallel sinoparams,
	struct ReconParamsQGGMRF3D reconparams,
	struct ImageParams3D imgparams,
	float *max_num_pointer,
	float **image,
	float *voxelsBuffer1,
	float *voxelsBuffer2,
	int* group_array,
	int group_id,
	int SV_per_Z,
	int SVsPerLine,
	long *updatedVoxels,
	float pow_sigmaX_p,
	float pow_sigmaX_q,
	float pow_T_qmp)
{

	int jy,jx,p,i,q,t,j,currentSlice;
	int startSlice;	
	int SV_depth_modified;	
	int SVLength = svpar.SVLength;
	int overlappingDistance = svpar.overlap;
	int SV_depth = svpar.SVDepth;
	struct minStruct * bandMinMap = svpar.bandMinMap;
	struct maxStruct * bandMaxMap = svpar.bandMaxMap;
	int pieceLength = svpar.pieceLength;

	if(it%2==0)
	{
		startSlice = order[jj] / Nx / Ny;
		jy = (order[jj] - startSlice* Nx * Ny) / Nx;  
		jx = (order[jj] - startSlice* Nx * Ny) % Nx;
	}
	else
	{
		startSlice = order[indexList[jj]] / Nx / Ny;
		jy=(order[indexList[jj]] - startSlice* Nx * Ny) /Nx;
		jx=(order[indexList[jj]] - startSlice* Nx * Ny) %Nx;	
	}

	if((startSlice+SV_depth)>Nslices)
		SV_depth_modified=Nslices-startSlice;
	else
		SV_depth_modified=SV_depth;

	int theSVPosition=jy/(2*SVLength-overlappingDistance)*SVsPerLine+jx/(2*SVLength-overlappingDistance);
	if(it%2==0)
	{
		if(phaseMap[jj]!=group_array[startSlice/SV_depth*4+group_id])
			return;
	}
	else
	{
		if(phaseMap[indexList[jj]]!=group_array[startSlice/SV_depth*4+group_id])
			return;
	}

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
			if(A_Padded_Map[theSVPosition][voxelIncrement].length >0) {
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
	int bandWidth[(sinoparams.NViews)/pieceLength]__attribute__((aligned(32)));	

	#ifdef USE_INTEL_MEMCPY
	_intel_fast_memcpy(&bandMin[0],&bandMinMap[theSVPosition].bandMin[0],sizeof(int)*(sinoparams.NViews));
	_intel_fast_memcpy(&bandMax[0],&bandMaxMap[theSVPosition].bandMax[0],sizeof(int)*(sinoparams.NViews)); 
	#else
	memcpy(&bandMin[0],&bandMinMap[theSVPosition].bandMin[0],sizeof(int)*(sinoparams.NViews));
	memcpy(&bandMax[0],&bandMaxMap[theSVPosition].bandMax[0],sizeof(int)*(sinoparams.NViews));
	#endif

	#pragma vector aligned 
	for(p=0;p< sinoparams.NViews;p++)
		bandWidthTemp[p]=bandMax[p]-bandMin[p];

	for (p = 0; p < (sinoparams.NViews)/pieceLength; p++)
	{
		int bandWidthMax=bandWidthTemp[p*pieceLength];
		for(t=0;t<pieceLength;t++){
			if(bandWidthTemp[p*pieceLength+t]>bandWidthMax)
				bandWidthMax=bandWidthTemp[p*pieceLength+t];
		}
		bandWidth[p]=bandWidthMax;
	}

	int tempCount=0;

	#pragma vector aligned
	#pragma simd reduction(+:tempCount) 
	for (p = 0; p < (sinoparams.NViews)/pieceLength; p++)
		tempCount+=bandWidth[p]*pieceLength;

	float ** newWArray = (float **)malloc(sizeof(float *)*(sinoparams.NViews)/pieceLength);
	float ** newEArray = (float **)malloc(sizeof(float *)*(sinoparams.NViews)/pieceLength);
	float ** CopyNewEArray = (float **)malloc(sizeof(float *)*(sinoparams.NViews)/pieceLength); 

	for (p = 0; p < (sinoparams.NViews)/pieceLength; p++)
	{
		newWArray[p] = (float *)malloc(sizeof(float)*bandWidth[p]*pieceLength*SV_depth_modified);
		newEArray[p] = (float *)malloc(sizeof(float)*bandWidth[p]*pieceLength*SV_depth_modified);
		CopyNewEArray[p] = (float *)malloc(sizeof(float)*bandWidth[p]*pieceLength*SV_depth_modified);
	}

	float *newWArrayPointer=&newWArray[0][0];
	float *newEArrayPointer=&newEArray[0][0];

	const int n_theta=sinoparams.NViews;

	/*XW: copy the interlaced we into the memory buffer*/ 
	for (p = 0; p < (sinoparams.NViews)/pieceLength; p++)
	{
		newWArrayPointer=&newWArray[p][0];
		newEArrayPointer=&newEArray[p][0];
		for(i=0;i<SV_depth_modified;i++)
		for(q=0;q<pieceLength;q++) 
		{
			#ifdef USE_INTEL_MEMCPY
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

	for (p = 0; p < (sinoparams.NViews)/pieceLength; p++)
	{
		#ifdef USE_INTEL_MEMCPY
		_intel_fast_memcpy(&CopyNewEArray[p][0],&newEArray[p][0],sizeof(float)*bandWidth[p]*pieceLength*SV_depth_modified);
		#else
		memcpy(&CopyNewEArray[p][0],&newEArray[p][0],sizeof(float)*bandWidth[p]*pieceLength*SV_depth_modified);
		#endif
	}

	float ** newWArrayTransposed = (float **)malloc(sizeof(float *)*(sinoparams.NViews)/pieceLength);
	float ** newEArrayTransposed = (float **)malloc(sizeof(float *)*(sinoparams.NViews)/pieceLength);

	for (p = 0; p < (sinoparams.NViews)/pieceLength; p++)
	{
		newWArrayTransposed[p] = (float *)malloc(sizeof(float)*bandWidth[p]*pieceLength*SV_depth_modified);
		newEArrayTransposed[p] = (float *)malloc(sizeof(float)*bandWidth[p]*pieceLength*SV_depth_modified);
	}

	float *WTransposeArrayPointer=&newWArrayTransposed[0][0];
	float *ETransposeArrayPointer=&newEArrayTransposed[0][0];
	
	for (p = 0; p < (sinoparams.NViews)/pieceLength; p++)
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
	float updateChange=0;
	float inverseNumber=1.0/255;

	for (p = 0; p < (sinoparams.NViews)/pieceLength; p++)
		free((void *)newWArray[p]);

	free((void **)newWArray);

	/*XW: the start of the loop to compute theta1, theta2*/
	for(i=0;i<countNumber;i++)
	{
		const short j_new=j_newCoordinate[i];   /*XW: get the voxel's x,y location*/
		const short k_new=k_newCoordinate[i];
		float tempV[SV_depth_modified];
		float neighbors[SV_depth_modified][10];
		char zero_skip_FLAG[SV_depth_modified];
		float max=max_num_pointer[j_new*Nx+k_new];
		float THETA1[SV_depth_modified];
		float THETA2[SV_depth_modified];
		memset(&THETA1[0],0.0, sizeof(THETA1));
		memset(&THETA2[0],0.0, sizeof(THETA2));	
		float diff[SV_depth_modified];

		int theVoxelPosition=(j_new-jy)*(2*SVLength+1)+(k_new-jx); 
		unsigned char * A_padd_Tranpose_pointer = &A_Padded_Map[theSVPosition][theVoxelPosition].val[0];

		for(currentSlice=0;currentSlice<SV_depth_modified;currentSlice++)
		{
			tempV[currentSlice] = (float)(image[startSlice+currentSlice][j_new*Nx+k_new]); /*XW: current voxel's value*/

			ExtractNeighbors3D(&neighbors[currentSlice][0],k_new,j_new,&image[startSlice+currentSlice][0],imgparams);

			if((startSlice+currentSlice)==0)
				neighbors[currentSlice][8]=voxelsBuffer1[j_new*Nx+k_new];
			else
				neighbors[currentSlice][8]=image[startSlice+currentSlice-1][j_new*Nx+k_new];

			if((startSlice+currentSlice)<(Nslices-1))
				neighbors[currentSlice][9]=image[startSlice+currentSlice+1][j_new*Nx+k_new];
			else
				neighbors[currentSlice][9]=voxelsBuffer2[j_new*Nx+k_new];

			zero_skip_FLAG[currentSlice] = 0;

			if (tempV[currentSlice] == 0.0)
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

		#ifdef find_RMSE
		for(currentSlice=0;currentSlice<SV_depth_modified;currentSlice++)
		{
			if(zero_skip_FLAG[currentSlice] == 0 ){
				#pragma omp atomic
				(*updatedVoxels)=(*updatedVoxels)+1;
			}
		}
		#endif

		A_padd_Tranpose_pointer = &A_Padded_Map[theSVPosition][theVoxelPosition].val[0];
		for(p=0;p<NViewsdivided;p++)
		{
			const int myCount=A_Padded_Map[theSVPosition][theVoxelPosition].pieceWiseWidth[p];
			const int pieceMin=A_Padded_Map[theSVPosition][theVoxelPosition].pieceWiseMin[p];
			#pragma vector aligned
			for(currentSlice=0;currentSlice<SV_depth_modified;currentSlice++)
			if(zero_skip_FLAG[currentSlice] == 0 )
			{
				WTransposeArrayPointer=&newWArrayTransposed[p][currentSlice*bandWidth[p]*pieceLength];
				ETransposeArrayPointer=&newEArrayTransposed[p][currentSlice*bandWidth[p]*pieceLength];
				WTransposeArrayPointer+=pieceMin*pieceLength;
				ETransposeArrayPointer+=pieceMin*pieceLength;
				float tempTHETA1=0.0;
				float tempTHETA2=0.0;
				#pragma vector aligned
				#pragma simd reduction(+:tempTHETA2,tempTHETA1)
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

		A_padd_Tranpose_pointer = &A_Padded_Map[theSVPosition][theVoxelPosition].val[0];
	
		/* float previousVoxel=0.0; */
		for(currentSlice=0;currentSlice<SV_depth_modified;currentSlice++)
		if(zero_skip_FLAG[currentSlice] == 0)
		{
			//if(currentSlice!=0)
			//	neighbors[currentSlice][8]=previousVoxel;
			float pixel = ICDStep3D(reconparams,THETA1[currentSlice],THETA2[currentSlice],tempV[currentSlice],&neighbors[currentSlice][0],pow_sigmaX_p,pow_sigmaX_q,pow_T_qmp); 
			/* For Synchrotron dataset.  Don't clip it */
			/* SJK: put positivity condition here */
			image[startSlice+currentSlice][j_new*Nx+k_new]= ((pixel < 0.0) ? 0.0 : pixel);  

			/* previousVoxel=image->img[currentSlice][j_new*Nx+k_new]; */
			diff[currentSlice] = image[startSlice+currentSlice][j_new*Nx+k_new] - tempV[currentSlice];
			updateChange += fabs(diff[currentSlice]);
			diff[currentSlice]=diff[currentSlice]*max*inverseNumber;
		}

		for(p=0;p<NViewsdivided;p++)
		{
			const int myCount=A_Padded_Map[theSVPosition][theVoxelPosition].pieceWiseWidth[p];
			const int pieceMin=A_Padded_Map[theSVPosition][theVoxelPosition].pieceWiseMin[p]; 
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

	for (p = 0; p < (sinoparams.NViews)/pieceLength; p++)
		free((void *)newWArrayTransposed[p]);

	free((void **)newWArrayTransposed);

	if((it%2)==0)
		headNodeArray[jj].x=updateChange;
	else
		headNodeArray[indexList[jj]].x=updateChange;

	for (p = 0; p < (sinoparams.NViews)/pieceLength; p++)
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

	for (p = 0; p < (sinoparams.NViews)/pieceLength; p++)
		free((void *)newEArrayTransposed[p]);

	free((void **)newEArrayTransposed);

	newEArrayPointer=&newEArray[0][0];
	float* CopyNewEArrayPointer=&CopyNewEArray[0][0];
	float* eArrayPointer=&e[0][0];

	for (p = 0; p < (sinoparams.NViews)/pieceLength; p++)      /*XW: update the error term in the memory buffer*/
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

	for (p = 0; p < (sinoparams.NViews)/pieceLength; p++)
	{
		free((void *)newEArray[p]);
		free((void *)CopyNewEArray[p]);
	}
	free((void **)newEArray);
	free((void **)CopyNewEArray);

	/* SJK: The stopping condition needs to be changed..for now I'm using total_updates to track number of voxel updates */
	//(*total_updates)+=updateChange;	
	(*total_updates)+= (float)(countNumber*SV_depth_modified);

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

void three_way_shuffle(int *order1, int *order2,struct heap_node *headNodeArray,int len)
{
	int i, j, tmp1,tmp2;

	float temp_x;

	for (i = 0; i < len-1; i++)
	{
		j = i + (rand() % (len-i));
		tmp1 = order1[j];
		tmp2 = order2[j];
		temp_x=headNodeArray[j].x;
		order1[j] = order1[i];
		order2[j] = order2[i];
		headNodeArray[j].x=headNodeArray[i].x;
		order1[i] = tmp1;
		order2[i] = tmp2;
		headNodeArray[i].x=temp_x;
	}
}



float MAPCostFunction3D(float **e,struct Image3D *Image,struct Sino3DParallel *sinogram,struct ReconParamsQGGMRF3D *reconparams)
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
void read_golden(char *fname,float **golden,int Nslices,int N, struct Image3D *Image)
{
	FILE *fp;
	int i;
        char slicefname[200];
        char *sliceindex;
	sliceindex= (char *)malloc(MBIR_MODULAR_MAX_NUMBER_OF_SLICE_DIGITS);
	for(i=0;i<Nslices;i++){
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




