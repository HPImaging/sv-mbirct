
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <string.h>
#include <omp.h>

#include "MBIRModularDefs.h"
#include "MBIRModularUtils.h"
#include "allocate.h"
#include "icd3d.h"
#include "heap.h"
#include "A_comp.h"
#include "initialize.h"
#include "recon3d.h"

/* The MBIR algorithm  */
/* Note : */
/* 1) Image must be intialized before this function is called */
/* 2) Image reconstruction Mask must be generated before this call */
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


void readErrorSinogram(char *fname,float *e,int M)
{
	FILE *fp;
	int i;
	if ((fp = fopen(fname, "r")) == NULL)
	{
		fprintf(stderr, "ERROR in writeErrorSinogram: can't open file %s.\n", fname);
		exit(-1);
	}
	fread(&e[0],sizeof(float),M,fp);
	fclose(fp);		
}

void read_golden(char *fname,float **golden,int mySize,int N, struct Image3D *Image)
{
	FILE *fp;
	int i;
        char slicefname[200];
        char *sliceindex;
    	sliceindex= (char *)malloc(MBIR_MODULAR_MAX_NUMBER_OF_SLICE_DIGITS);
	for(i=0;i<mySize;i++){
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

static __inline__ unsigned long long rdtsc()
{
   unsigned hi,lo;
   __asm __volatile__ ("rdtsc" : "=a"(lo),"=d"(hi));
   return ( (unsigned long long)lo)|( ((unsigned long long)hi)<<32 );
}




void super_voxel_recon(int jj,float *total_updates,int it, int *phaseMap,int *order,int *indexList,int Nx,int Ny,struct minStruct *bandMinMap, struct maxStruct *bandMaxMap,float **w,float **e, struct AValues_char ** A_Padded_Map,const int mySize,struct heap_node *headNodeArray,const int N_thetadivided,struct SinoParams3DParallel sinoparams,struct ReconParamsQGGMRF3D reconparams,struct ImageParams3D imgparams,float *max_num_pointer,float **image,float *voxelsBuffer1,float *voxelsBuffer2,int* group_array,int group_id,int SV_per_Z,int SVsPerLine,long *updatedVoxels,float pow_sigmaX_p,float pow_sigmaX_q,float pow_T_qmp,int pieceLength){

			int jy,jx,p,i,q,t,j,currentSlice;
			int startSlice;	
			int SV_depth_modified;	

			if(it%2==0){
				startSlice = order[jj] / Nx / Ny;
				jy = (order[jj] - startSlice* Nx * Ny) / Nx;  
				jx = (order[jj] - startSlice* Nx * Ny) % Nx;
			}
			else{
				startSlice = order[indexList[jj]] / Nx / Ny;
				jy=(order[indexList[jj]] - startSlice* Nx * Ny) /Nx;
				jx=(order[indexList[jj]] - startSlice* Nx * Ny) %Nx;	
							
			}

			
			if((startSlice+SV_depth)>mySize){
				SV_depth_modified=mySize-startSlice;
			}
			else{
				SV_depth_modified=SV_depth;
			}
			
			
			
            		int theSVPosition=jy/(2*SVLength-overlappingDistance1)*SVsPerLine+jx/(2*SVLength-overlappingDistance2);				
		
			if(it%2==0){
				if(phaseMap[jj]!=group_array[startSlice/SV_depth*4+group_id]){
					return;
				}
			}

			else{
				if(phaseMap[indexList[jj]]!=group_array[startSlice/SV_depth*4+group_id]){
					return;
				}
			}                                                 			

			int countNumber=0;  /*XW: the number of voxels chosen for a certain radius of circle*/
			int radius =SVLength;  /*XW: choose the circle radius*/
			int coordinateSize=1;  /*XW: CoordinateSize is the size of a minimum square enclosing the circle. For the baseline code, coordinateSize=1 because only 1 voxel*/
			if(radius!=0)
			  coordinateSize=(2*radius+1)*(2*radius+1);/*XW: imagine that this is a minimum square enclosing the circle. The square "touches" the circle on 4 coordinate points. CoordianteSize is the possible maximum number of voxels for a certain radius of circle */
			int k_newCoordinate[coordinateSize];
			int j_newCoordinate[coordinateSize];
			int j_newAA=0;
			int k_newAA=0;			
			
			int voxelIncrement=0;
			
			/*XW: choosing the voxels locations in a circle*/
			for(j_newAA=jy;j_newAA<=(jy+2*radius);j_newAA++){
			  for(k_newAA=jx;k_newAA<=(jx+2*radius);k_newAA++){
                              if(j_newAA>=0 && k_newAA >=0 && j_newAA <Ny && k_newAA < Nx){
				if(A_Padded_Map[theSVPosition][voxelIncrement].length >0){
  				  	j_newCoordinate[countNumber]=j_newAA;
  				  	k_newCoordinate[countNumber]=k_newAA;
  	   			  	countNumber++;
				} 
  	   		      }
  	   		      voxelIncrement++;	
			  }
			}
					
			
			/*XW: if no voxel chosen, we skip this loop iteration*/
				
			if(countNumber==0)
			   return;



			coordinateShuffle(&j_newCoordinate[0],&k_newCoordinate[0],countNumber);
			

                        int bandMin[sinoparams.NViews]__attribute__((aligned(32)));  /*XW: for a supervoxel, bandMin records the starting position of the sinogram band at each view*/
			int bandMax[sinoparams.NViews]__attribute__((aligned(32)));/*XW: for a supervoxel, bandMax records the end position of the sinogram band at each view */
            		int bandWidthTemp[sinoparams.NViews]__attribute__((aligned(32)));
            		int bandWidth[(sinoparams.NViews)/pieceLength]__attribute__((aligned(32)));	

            		_intel_fast_memcpy(&bandMin[0],&bandMinMap[theSVPosition].bandMin[0],sizeof(int)*(sinoparams.NViews));
            		_intel_fast_memcpy(&bandMax[0],&bandMaxMap[theSVPosition].bandMax[0],sizeof(int)*(sinoparams.NViews)); 
	    		
	    			
            		#pragma vector aligned 
            		for(p=0;p< sinoparams.NViews;p++){		
                		bandWidthTemp[p]=bandMax[p]-bandMin[p];
            		}


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
            		for (p = 0; p < (sinoparams.NViews)/pieceLength; p++){
                		tempCount+=bandWidth[p]*pieceLength;
            		}
         		           		           		                     
    			
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
     			
	
      			for (p = 0; p < (sinoparams.NViews)/pieceLength; p++){   /*XW: copy the interlaced we into the memory buffer*/ 
				newWArrayPointer=&newWArray[p][0];
				newEArrayPointer=&newEArray[p][0];
				for(i=0;i<SV_depth_modified;i++){  				 
					for(q=0;q<pieceLength;q++){      			      						
          		   			_intel_fast_memcpy(newWArrayPointer,&w[startSlice+i][p*pieceLength*sinoparams.NChannels+q*sinoparams.NChannels+bandMin[p*pieceLength+q]],sizeof(float)*(bandWidth[p]));
          		   			_intel_fast_memcpy(newEArrayPointer,&e[startSlice+i][p*pieceLength*sinoparams.NChannels+q*sinoparams.NChannels+bandMin[p*pieceLength+q]],sizeof(float)*(bandWidth[p]));
			   			newWArrayPointer+=bandWidth[p];
			   			newEArrayPointer+=bandWidth[p];              		   			 		
      					}  					
      				}
      			}
      			

			for (p = 0; p < (sinoparams.NViews)/pieceLength; p++){
				_intel_fast_memcpy(&CopyNewEArray[p][0],&newEArray[p][0],sizeof(float)*bandWidth[p]*pieceLength*SV_depth_modified);
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
            		{
				for(currentSlice=0;currentSlice<(SV_depth_modified);currentSlice++){
				WTransposeArrayPointer=&newWArrayTransposed[p][currentSlice*bandWidth[p]*pieceLength];
				ETransposeArrayPointer=&newEArrayTransposed[p][currentSlice*bandWidth[p]*pieceLength];
				newEArrayPointer=&newEArray[p][currentSlice*bandWidth[p]*pieceLength];
				newWArrayPointer=&newWArray[p][currentSlice*bandWidth[p]*pieceLength];					            	
                		for(q=0;q<bandWidth[p];q++){  			
						#pragma vector aligned 
                    				for(t=0;t<pieceLength;t++){
                        				ETransposeArrayPointer[q*pieceLength+t]=newEArrayPointer[bandWidth[p]*t+q];
                        				WTransposeArrayPointer[q*pieceLength+t]=newWArrayPointer[bandWidth[p]*t+q];
                    				}
                			}
            			}
            		}
   
            		WTransposeArrayPointer=&newWArrayTransposed[0][0];
            		ETransposeArrayPointer=&newEArrayTransposed[0][0];
            		newEArrayPointer=&newEArray[0][0];	
			float updateChange=0;
			float inverseNumber=1.0/255;
			 					
			
			
            		for (p = 0; p < (sinoparams.NViews)/pieceLength; p++)
            		{
            			free((void *)newWArray[p]);
        		}
    			free((void **)newWArray);
									
											
		        for(i=0;i<countNumber;i++){  /*XW: the start of the loop to compute theta1, theta2*/

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
							
				for(currentSlice=0;currentSlice<SV_depth_modified;currentSlice++){ 			 
			        	tempV[currentSlice] = (float)(image[startSlice+currentSlice][j_new*Nx+k_new]); /*XW: find the current voxel's value*/
			                
			      		ExtractNeighbors3D(&neighbors[currentSlice][0], k_new, j_new, &image[startSlice+currentSlice][0], imgparams);  /*XW: extract neighbors*/ 	      		
			      		
			      		if((startSlice+currentSlice)==0){
			      			neighbors[currentSlice][8]=voxelsBuffer1[j_new*Nx+k_new];
			      		}			      				
			      		else{
			      		        neighbors[currentSlice][8]=image[startSlice+currentSlice-1][j_new*Nx+k_new];
			      		}        
			      		        			      		        
			      		        
			      		if((startSlice+currentSlice)<(mySize-1)){		
			      			neighbors[currentSlice][9]=image[startSlice+currentSlice+1][j_new*Nx+k_new];
			      		}		
			      		else{
			      			neighbors[currentSlice][9]=voxelsBuffer2[j_new*Nx+k_new];
			      		}				
			      			
	                						
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
			      			
				for(currentSlice=0;currentSlice<SV_depth_modified;currentSlice++){
					if(zero_skip_FLAG[currentSlice] == 0 ){				
			      			#pragma omp atomic									
			        		(*updatedVoxels)=(*updatedVoxels)+1;
			        	}
			        }			     
			      	#endif	
			        			     

				A_padd_Tranpose_pointer = &A_Padded_Map[theSVPosition][theVoxelPosition].val[0];				
                    		for(p=0;p<N_thetadivided;p++){
                        		const int myCount=A_Padded_Map[theSVPosition][theVoxelPosition].pieceWiseWidth[p];
                        		const int pieceMin=A_Padded_Map[theSVPosition][theVoxelPosition].pieceWiseMin[p];
                        		#pragma vector aligned
					for(currentSlice=0;currentSlice<SV_depth_modified;currentSlice++){					
						if(zero_skip_FLAG[currentSlice] == 0 ){	
							WTransposeArrayPointer=&newWArrayTransposed[p][currentSlice*bandWidth[p]*pieceLength];
							ETransposeArrayPointer=&newEArrayTransposed[p][currentSlice*bandWidth[p]*pieceLength];
                        				WTransposeArrayPointer+=pieceMin*pieceLength;
                       					ETransposeArrayPointer+=pieceMin*pieceLength;
                       					float tempTHETA1=0.0;
                       					float tempTHETA2=0.0;
                        				#pragma vector aligned
                            				#pragma simd reduction(+:tempTHETA2,tempTHETA1)                       			
                        				for(t=0;t<myCount*pieceLength;t++){                    				      
  						/*XW: if not zero skipped, compute the theta1, theta2*/	
						/*XW: adding the number of updated voxels which are not skipped or masked*/     
                                    				tempTHETA1 += A_padd_Tranpose_pointer[t]*WTransposeArrayPointer[t]*ETransposeArrayPointer[t];
                                    				tempTHETA2 += A_padd_Tranpose_pointer[t]*WTransposeArrayPointer[t]*A_padd_Tranpose_pointer[t];
											
                        				}
                        				THETA1[currentSlice]+=tempTHETA1;
                        				THETA2[currentSlice]+=tempTHETA2;
                        			}
                    			}
                        		A_padd_Tranpose_pointer+=myCount*pieceLength;                    			
                    		}		
				for(currentSlice=0;currentSlice<SV_depth_modified;currentSlice++){				
					THETA1[currentSlice]=-THETA1[currentSlice]*max*inverseNumber;
                    			THETA2[currentSlice]=THETA2[currentSlice]*max*inverseNumber*max*inverseNumber;
                    		}		 				
	
						
				ETransposeArrayPointer=&newEArrayTransposed[0][0];

                    		A_padd_Tranpose_pointer = &A_Padded_Map[theSVPosition][theVoxelPosition].val[0];
		
	
				/* float previousVoxel=0.0; */
				for(currentSlice=0;currentSlice<SV_depth_modified;currentSlice++)		
				{
					if(zero_skip_FLAG[currentSlice] == 0){
						/*
						if(currentSlice!=0){
							neighbors[currentSlice][8]=previousVoxel;
						}
						*/					
						float pixel = ICDStep3D(reconparams,THETA1[currentSlice],THETA2[currentSlice],tempV[currentSlice],&neighbors[currentSlice][0],pow_sigmaX_p,pow_sigmaX_q,pow_T_qmp); /*XW: do the functional substitution to get updated voxel value*/
					/* For Synchrotron dataset.  Don't clip it */				
			      			image[startSlice+currentSlice][j_new*Nx+k_new]= ((pixel < 0.0) ? 0.0 : pixel);  

						/* previousVoxel=image->img[currentSlice][j_new*Nx+k_new]; */
						diff[currentSlice] = image[startSlice+currentSlice][j_new*Nx+k_new] - tempV[currentSlice];
						updateChange += fabs(diff[currentSlice]);
						diff[currentSlice]=diff[currentSlice]*max*inverseNumber;
	                    		}			
				}	
	
                        	for(p=0;p<N_thetadivided;p++){
                            		const int myCount=A_Padded_Map[theSVPosition][theVoxelPosition].pieceWiseWidth[p];
                                        const int pieceMin=A_Padded_Map[theSVPosition][theVoxelPosition].pieceWiseMin[p]; 
                      			#pragma vector aligned    		
                            		for(currentSlice=0;currentSlice<SV_depth_modified;currentSlice++){
                                         	if(diff[currentSlice]!=0 && zero_skip_FLAG[currentSlice] == 0){
							ETransposeArrayPointer=&newEArrayTransposed[p][currentSlice*bandWidth[p]*pieceLength];
                      					ETransposeArrayPointer+=pieceMin*pieceLength;
                      					
                                			#pragma vector aligned            					
                      					
                            				for(t=0;t<(myCount*pieceLength);t++){
                            				     ETransposeArrayPointer[t]= ETransposeArrayPointer[t]-A_padd_Tranpose_pointer[t]*diff[currentSlice];
                            				}
                        			}
                    			}
                                        A_padd_Tranpose_pointer+=myCount*pieceLength;			
				}
											
			}
						
			
            		for (p = 0; p < (sinoparams.NViews)/pieceLength; p++)
            		{
            			free((void *)newWArrayTransposed[p]);
        		}
    			free((void **)newWArrayTransposed);			
			
			
			

			if((it%2)==0){
				headNodeArray[jj].x=updateChange;
				
			}
			else{
				headNodeArray[indexList[jj]].x=updateChange;
			}
		
            		
            		
            		for (p = 0; p < (sinoparams.NViews)/pieceLength; p++)
            		{
				for(currentSlice=0;currentSlice<SV_depth_modified;currentSlice++){
				ETransposeArrayPointer=&newEArrayTransposed[p][currentSlice*bandWidth[p]*pieceLength];
				newEArrayPointer=&newEArray[p][currentSlice*bandWidth[p]*pieceLength]; 
                		for(q=0;q<bandWidth[p];q++){
					#pragma vector aligned 		                		  			
                    				for(t=0;t<pieceLength;t++){
                        				newEArrayPointer[bandWidth[p]*t+q]=ETransposeArrayPointer[q*pieceLength+t];
                    				}
                			}
            			}
            		}            		 		
            		
            		
            		for (p = 0; p < (sinoparams.NViews)/pieceLength; p++)
            		{
            			free((void *)newEArrayTransposed[p]);
        		}
    			free((void **)newEArrayTransposed);            		
            		 										

			newEArrayPointer=&newEArray[0][0];
			float* CopyNewEArrayPointer=&CopyNewEArray[0][0];
			float* eArrayPointer=&e[0][0];


			for (p = 0; p < (sinoparams.NViews)/pieceLength; p++)      /*XW: update the error term in the memory buffer*/
			{
				newEArrayPointer=&newEArray[p][0];
				CopyNewEArrayPointer=&CopyNewEArray[p][0];
				for (currentSlice=0; currentSlice< SV_depth_modified;currentSlice++){	
				#pragma vector aligned											
                                for(q=0;q<pieceLength;q++){
						eArrayPointer=&e[startSlice+currentSlice][p*pieceLength*sinoparams.NChannels+q*sinoparams.NChannels+bandMin[p*pieceLength+q]];
						for(t=0;t<bandWidth[p];t++){
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
			
			
			(*total_updates)+=updateChange;	

}



void MBIRReconstruct3D(
                       struct Image3D *Image,
                       struct Sino3DParallel *sinogram,
                       struct ReconParamsQGGMRF3D reconparams,
                       char *ImageReconMask,struct minStruct *bandMinMap,struct maxStruct *bandMaxMap,
	struct AValues_char ** A_Padded_Map,float *max_num_pointer,struct CmdLineMBIR * cmdline,int sum, int pieceLength)
{
    	int it, MaxIterations, j, k, l, Nx, Ny, Nz, Nxy, N, i, XYPixelIndex, M, currentSlice,jj,p,t,N_theta;
    	float **x;  /* image data */
    	float **y;  /* sinogram projections data */
    	float **e;  /* e=y-Ax, error */
    	float **w;  /* projections weights data */
    	float *voxelsBuffer1;  /* the first N entries are the voxel values.  */
    	float *voxelsBuffer2;   
    	struct heap priorityheap;
        struct timeval tm1,tm2;    	   	          	
        
    	initialize_heap(&priorityheap);
    	int *order;    	      
  
    	float voxel, diff;
    	float cost, avg_update, AvgVoxelValue, total_updates, StopThreshold, ratio;
    	char zero_skip_FLAG;
    	char stop_FLAG;
	char fname[200];
	float pow_sigmaX_p,pow_sigmaX_q,pow_T_qmp;
    
    	/*struct ICDInfo icd_info; */ /* Local Cost Function Information */
    
   	x = Image->image;   /* x is the image vector */
    	y = sinogram->sino;  /* y is the sinogram vector */
   	w = sinogram->weight; /* vector of weights for each sinogram measurement */
   	Nx = Image->imgparams.Nx;
    	Ny = Image->imgparams.Ny;
    	N= Nx*Ny;
    	M = sinogram->sinoparams.NViews * sinogram->sinoparams.NChannels;
    	N_theta=sinogram->sinoparams.NViews;
    	stop_FLAG=0;
    	MaxIterations = reconparams.MaxIterations;
    	StopThreshold = reconparams.StopThreshold;    	
    
    	int SV_per_Z=0;
    	int rep_num=(int)ceil(1/(4*c_ratio*convergence_rho));
	const int mySize=sinogram->sinoparams.NSlices;
	
	pow_sigmaX_p=pow(reconparams.SigmaX,reconparams.p);
	pow_sigmaX_q=pow(reconparams.SigmaX,reconparams.q);
	pow_T_qmp=pow(reconparams.T,reconparams.q - reconparams.p);	
	
    	if((mySize%SV_depth)==0)
		SV_per_Z=mySize/SV_depth;
    	else
		SV_per_Z=mySize/SV_depth+1;
	
    	if((N_theta)%pieceLength !=0)
    	{
		fprintf(stderr, "ERROR in pieceLength\n");
		exit(-1);
    	}
    
    	const int N_thetadivided=N_theta/pieceLength;


    	#ifdef find_RMSE

    	float **golden;
    	golden = (float **)multialloc(sizeof(float), 2, mySize,N);	

	read_golden(cmdline->ReconImageDataFile,golden,mySize,N,Image);
    	float updatedVoxelsList[300];


	float sumOfSE=0;
	for(i=0;i<mySize;i++){
		for(j=0;j<N;j++){
			sumOfSE+=(Image->image[i][j]-golden[i][j])*(Image->image[i][j]-golden[i][j]);
		}
	}
			
	float MSE=sumOfSE/N/(mySize);
	float RMSE=sqrt(MSE);


	fprintf(stdout,"Rho: %f initial_RMSE: %f \n",convergence_rho,RMSE);
    	#endif     

    	/* Order of pixel updates need NOT be raster order, just initialize */
    	order = (int *)_mm_malloc(sum*SV_per_Z*sizeof(int),64);
	
    	t=0;	
	
    	for(p=0;p<mySize;p+=SV_depth){
		for(i=0;i<Ny;i+=(SVLength*2-overlappingDistance1)){
			for(j=0;j<Nx;j+=(SVLength*2-overlappingDistance2)){
				order[t]=p*Nx*Ny+i*Nx+j;  /* order is the first voxel coordinate, not the center */
		    		t++;
			}
		}	
     	}	
     
     	int phaseMap[sum*SV_per_Z];
     	int SVsPerLine=0; 
	if((Nx%(2*SVLength-overlappingDistance2))==0)
		SVsPerLine=Nx/(2*SVLength-overlappingDistance2);
	else	
		SVsPerLine=Nx/(2*SVLength-overlappingDistance2)+1; 

	#pragma omp parallel for private(jj) schedule(dynamic)
	for(i=0;i<SV_per_Z;i++){
		for(jj=0;jj<sum;jj++)
		{
			if((jj/SVsPerLine)%2==0){
				if((jj%SVsPerLine)%2==0){
					phaseMap[i*sum+jj]=0;
				}
				else{
					phaseMap[i*sum+jj]=1;			
				}
			}	
			else{
				if((jj%SVsPerLine)%2==0){
					phaseMap[i*sum+jj]=2;
				}
				else{
					phaseMap[i*sum+jj]=3;			
				}		
			}
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

	for(i=0;i<SV_per_Z;i++){
		for(jj=0;jj<sum;jj++)
		{
			headNodeArray[i*sum+jj].pt=i*sum+jj;		
			headNodeArray[i*sum+jj].x=0.0;
		}
	}
	int indexList_size=(int) sum*SV_per_Z*4*c_ratio*(1-convergence_rho);	
	int indexList[indexList_size];   	             	    

    
    	/********************************************/
    	/* Forward Projection and Error Calculation */
    	/********************************************/
	e = (float **)multialloc(sizeof(float), 2, mySize,M);  	 /* error term memory allocation */
    	/* Initialize error to zero, since it is first computed as forward-projection Ax */
    	/* compute Ax (store it in e as of now) */
    	sprintf(fname,"%s.initialError",cmdline->SysMatrixFile);		    	    	   	
    	
	#pragma omp parallel for private(i) schedule(dynamic)
	for(currentSlice=0;currentSlice<mySize;currentSlice++){    	
    		readErrorSinogram(fname,e[currentSlice],M);    		
    		/* Compute the initial error e=y-Ax */    		
    		for (i = 0; i < M; i++){
			e[currentSlice][i]=y[currentSlice][i]-e[currentSlice][i];
		}	
    	}   
	voxelsBuffer1 = (float *)_mm_malloc(N*sizeof(float),64);
	voxelsBuffer2 = (float *)_mm_malloc(N*sizeof(float),64);

        for(i=0;i<N;i++){
        	voxelsBuffer1[i]=0;
        }
        
        for(i=0;i<N;i++){
        	voxelsBuffer2[i]=0;
        }
               
	it=0;

	long updatedVoxels=0;

	
	coordinateShuffle(&order[0],&phaseMap[0],sum*SV_per_Z);
	
	int startIndex=0;
	int endIndex=0;        		
    	 			
        gettimeofday(&tm1,NULL);
         
	#pragma omp parallel
	{
		while(stop_FLAG==0 && it <MaxIterations)
		{
		
			#pragma omp single
			{		
				if(it==0){
					startIndex=0;
					endIndex=sum*SV_per_Z;
				}	
				else{			
					if((it-1)%(2*rep_num)==0 && it!=1){
						three_way_shuffle(&order[0],&phaseMap[0],&headNodeArray[0],sum*SV_per_Z);
					}
			
					if(it%2==1){
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
			
			for (group = 0; group < 4; group++){
			
				#pragma omp for schedule(dynamic)  reduction(+:total_updates)
				for (jj = startIndex; jj < endIndex; jj+=1)
				{				
					super_voxel_recon(jj,&total_updates,it, &phaseMap[0],order,&indexList[0],Nx,Ny, bandMinMap, bandMaxMap,w,e,A_Padded_Map,mySize,&headNodeArray[0],N_thetadivided,sinogram->sinoparams,reconparams,Image->imgparams,&max_num_pointer[0],Image->image,voxelsBuffer1,voxelsBuffer2,&group_id_list[0][0],group,SV_per_Z,SVsPerLine,&updatedVoxels,pow_sigmaX_p,pow_sigmaX_q,pow_T_qmp,pieceLength);
				}
			}
			
			
			#pragma omp single
			{
           		
			if(it==0){
				avg_update = total_updates/(float)N/mySize;
			}		
			else{
				if(it%2==1)
					avg_update = (total_updates/(float)N)*(4*convergence_rho/(1-convergence_rho))/mySize;				
				else
					avg_update = (total_updates/(float)N)*4/mySize;
			}           		
			
			/*
			cost = MAPCostFunction3D(e, Image, sinogram, &reconparams);
        		fprintf(stdout, "it %d cost = %-15f, avg_update %f \n", it, cost, avg_update);           	
			*/
			
			#ifdef find_RMSE

			if(it<300){	
			
				updatedVoxelsList[it]=updatedVoxels*1.0/N/(mySize);
			}
			
	
			float sumOfSE=0;
			for(i=0;i<mySize;i++){
				for(j=0;j<N;j++){
					sumOfSE+=(Image->image[i][j]-golden[i][j])*(Image->image[i][j]-golden[i][j]);
				}
			}
			
			float MSE=sumOfSE/N/(mySize);
			float RMSE=sqrt(MSE);

			fprintf(stdout,"Rho: %f Equits: %f RMSE: %f \n",convergence_rho,updatedVoxelsList[it],RMSE);
			
			#endif		


        		
        		if ((avg_update) < StopThreshold && (endIndex!=0))
        		{
            			stop_FLAG = 1;
        		}
        		
		

			it++;
			total_updates = 0.0;
				
			}				
				
		}		
	}
	
        gettimeofday(&tm2,NULL);
        unsigned long long tt = 1000 * (tm2.tv_sec - tm1.tv_sec) + (tm2.tv_usec - tm1.tv_usec) / 1000;
        printf("run time %llu ms (iterations only)\n", tt);

/*	sprintf(fname, "%s.err", prior_info->oname);
	writeErrorSinogram(fname, e, sino_info);
	sprintf(fname, "%s.nmlerr", prior_info->oname);
	writeNormalizedErrorSinogram(fname, e, w, sino_info);*/


	if (stop_FLAG == 0 && StopThreshold > 0)
	{
		fprintf(stdout, "WARNING: Didn't reach stopping condition.\n");
		fprintf(stdout, "Average update magnitude = %f\n", avg_update);
	}
	else if (stop_FLAG == 0 && StopThreshold <= 0)
	{
		fprintf(stdout, "No stopping condition.\n");
		fprintf(stdout, "Average update magnitude = %f\n", avg_update);
	}
    	if(AvgVoxelValue>0)
    	fprintf(stdout, "Average Update to Average Voxel-Value Ratio = %f %% \n", ratio);	
	
	multifree(e,2);
	_mm_free((void *)order);
	_mm_free((void *)voxelsBuffer1);
	_mm_free((void *)voxelsBuffer2);
	if(priorityheap.size>0)
		free_heap((void *)&priorityheap); 

        for(jj=0;jj<sum;jj++){
            	free((void *)bandMinMap[jj].bandMin);
            	free((void *)bandMaxMap[jj].bandMax);
        }

    	free((void *)bandMinMap);
    	free((void *)bandMaxMap);
    	
	#ifdef find_RMSE
    	multifree(golden,2);        	
        #endif
            	
    	for(i=0;i<sum;i++){
        	for(jj=0;jj<((2*SVLength+1)*(2*SVLength+1));jj++){
            		if(A_Padded_Map[i][jj].length>0){
                		free((void *)A_Padded_Map[i][jj].val);
                		free((void *)A_Padded_Map[i][jj].pieceWiseMin);
                		free((void *)A_Padded_Map[i][jj].pieceWiseWidth);
            		}
        	}
    	}
    	multifree(A_Padded_Map,2);
    	free((void *)max_num_pointer);
}
				
						     



/* The function to compute cost function */

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
  {
	for (jy = 0; jy < Ny; jy++)
	{
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
	}
  }

	return (nloglike + reconparams->b_nearest * nlogprior_nearest + reconparams->b_diag * nlogprior_diag + reconparams->b_interslice * nlogprior_interslice) ;
}


/* compute A times X, A-matrix is pre-computed */
/*
void forwardProject3D(
                      float *AX,           
                      struct Image3D *X,
                      struct SysMatrix2D *A,
                      struct SinoParams3DParallel sinoparams)
{
    int j,k,n, jz, Nxy, NViewsTimesNChannels, NSlices ;
    struct SparseColumn A_column;
    
    printf("\nComputing Forward Projection ... \n");
    
    Nxy = X->imgparams.Nx * X->imgparams.Ny; 
    NViewsTimesNChannels = sinoparams.NViews * sinoparams.NChannels; 
    NSlices = sinoparams.NSlices;
    
    if(A->Ncolumns != Nxy)
    {
      fprintf(stderr,"Error in forwardProject3D : dimensions of System Matrix and Image are not compatible \n");
      exit(-1);
    }
    
    for (j = 0; j < A->Ncolumns; j++) 
    {
        A_column = A->column[j]; 
        if (A_column.Nnonzero > 0)
        {
            for (n = 0; n < A_column.Nnonzero; n++)
            {
                for(jz=0;jz<NSlices;jz++)                          
                {
                  k = A_column.RowIndex[n] + jz*NViewsTimesNChannels ; 
                  AX[k] += A_column.Value[n]*X->image[j+jz*Nxy] ;      
                }
            }
        }
    }
}
*/

/* shuffle the coordinate to enable random update */
void shuffle(int *order, int len)
{
    int i, j, tmp;
    
    srand(time(NULL));
    
    for (i = 0; i < len-1; i++)
    {
        j = i + (rand() % (len-i));
        tmp = order[j];
        order[j] = order[i];
        order[i] = tmp;
    }
}
