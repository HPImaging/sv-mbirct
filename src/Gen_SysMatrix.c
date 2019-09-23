
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "MBIRModularUtils_3D.h"
#include "MBIRModularUtils_2D.h"
#include "allocate.h"
#include "A_comp.h"

/* #ifdef STORE_A_MATRIX - This option set by default in A_comp.h */
/* This option is to precompute and store the forward matrix rather than compute it on the fly */

void writeErrorSinogram(
	char *fname,
	float *e,
	struct SinoParams3DParallel *sinoparams)
{
	FILE *fp;
	int i, M;
	M = sinoparams->NViews*sinoparams->NChannels;
	if ((fp = fopen(fname, "w")) == NULL)
	{
		fprintf(stderr, "ERROR in writeErrorSinogram: can't open file %s.\n", fname);
		exit(-1);
	}
	fwrite(&e[0],sizeof(float),M,fp);
	fclose(fp);	
}



void forwardProject2D(
	float *e,
	float InitValue,
	float *max_num_pointer,struct AValues_char ** A_Padded_Map,struct minStruct *bandMinMap,struct SinoParams3DParallel *sinoparams,struct ImageParams3D *imgparams)
{
	int jx=0; 
	int jy=0; 
	int Nx, Ny, i, M, r,j,p, SVNumPerRow;
	float inverseNumber=1.0/255;
    	const int NViewsdivided=(sinoparams->NViews)/pieceLength;		

	Nx = imgparams->Nx;
	Ny = imgparams->Ny;
	M = sinoparams->NViews*sinoparams->NChannels;

	for (i = 0; i < M; i++)
	{
		e[i] = 0.0;
	}
		

	if((Nx%(2*SVLength-overlappingDistance2))==0)
		SVNumPerRow=Nx/(2*SVLength-overlappingDistance2);
	else	
		SVNumPerRow=Nx/(2*SVLength-overlappingDistance2)+1; 


	for (jy = 0; jy < Ny; jy++)
	{
		for (jx = 0; jx < Nx; jx++)
		{
		
			int temp1=jy/(2*SVLength-overlappingDistance1);
			
			if(temp1==SVNumPerRow){
				temp1=SVNumPerRow-1;
			}
			
			int temp2=jx/(2*SVLength-overlappingDistance2);
			if(temp2==SVNumPerRow){
				temp2=SVNumPerRow-1;
			}
		
            		int theSVPosition=temp1*SVNumPerRow+temp2;	
            			
			int SV_jy=temp1*(2*SVLength-overlappingDistance1);
			int SV_jx=temp2*(2*SVLength-overlappingDistance2);
			int theVoxelPosition=(jy-SV_jy)*(2*SVLength+1)+(jx-SV_jx); 
			/*
			fprintf(stdout,"jy %d jx %d theSVPosition %d SV_jy %d SV_jx %d theVoxelPosition %d \n",jy,jx,theSVPosition,SV_jy,SV_jx,theVoxelPosition);		
			*/
			if (A_Padded_Map[theSVPosition][theVoxelPosition].length > 0 && theVoxelPosition < ((2*SVLength+1)*(2*SVLength+1)))
			{
				/*XW: remove the index field in struct ACol and exploit the spatial locality */
				unsigned char* A_padd_Tranpose_pointer = &A_Padded_Map[theSVPosition][theVoxelPosition].val[0];					
                    		for(p=0;p<NViewsdivided;p++){        		
                        		const int myCount=A_Padded_Map[theSVPosition][theVoxelPosition].pieceWiseWidth[p];
                   		
		   			int position=p*pieceLength*sinoparams->NChannels+A_Padded_Map[theSVPosition][theVoxelPosition].pieceWiseMin[p];
					                        		
                        		for(r=0;r<myCount;r++){                    		
						for(j=0;j< pieceLength;j++){
							if((A_Padded_Map[theSVPosition][theVoxelPosition].pieceWiseMin[p]+bandMinMap[theSVPosition].bandMin[p*pieceLength+j]+r)>=sinoparams->NChannels)
							fprintf(stdout, "p %d r %d j %d total_1 %d \n",p,r,j,A_Padded_Map[theSVPosition][theVoxelPosition].pieceWiseMin[p]+bandMinMap[theSVPosition].bandMin[p*pieceLength+j]+r);						
		
							if((position+j*sinoparams->NChannels+bandMinMap[theSVPosition].bandMin[p*pieceLength+j]+r)>= M)
							fprintf(stdout, "p %d r %d j %d total_2 %d \n",p,r,j,position+j*sinoparams->NChannels+bandMinMap[theSVPosition].bandMin[p*pieceLength+j]+r);													
						
							if((position+j*sinoparams->NChannels+bandMinMap[theSVPosition].bandMin[p*pieceLength+j]+r)< M)							{	
								e[position+j*sinoparams->NChannels+bandMinMap[theSVPosition].bandMin[p*pieceLength+j]+r] += A_padd_Tranpose_pointer[r*pieceLength+j]*max_num_pointer[jy*Nx+jx]*inverseNumber*InitValue;
							}
																
						}				
					}
                        		A_padd_Tranpose_pointer+=myCount*pieceLength;					
				}
			}
		}	
	}
}
	



int main(int argc, char *argv[])
{

    struct CmdLineSysGen cmdline;
    struct ImageParams3D imgparams;
    struct SinoParams3DParallel sinoparams;
    //struct ReconParamsQGGMRF3D reconparams;    
    float **PixelDetector_profile;
    struct minStruct *bandMinMap;
    struct maxStruct *bandMaxMap;
    struct AValues_char ** A_Padded_Map;    
    int *order;   
    float x_0, y_0, Deltaxy, x, y, yy, ROIRadius, R_sq, R_sq_max;
    int jx, jy,Nxy;         
    
    /*struct SysMatrix2D *A ;*/
    
    /* read Command Line */
    	fprintf(stdout,"argc %d argv[0] %s argv[1] %s argv[2] %s argv[3] %s argv[4] %s argv[5] %s argv[6] %s argv[7] %s argv[8] %s\n",argc,argv[0],argv[1],argv[2],argv[3],argv[4],argv[5],argv[6],argv[7],argv[8]);      
    
    
    	readCmdLineSysGen(argc, argv, &cmdline);
    	fprintf(stdout,"reach 1 system matrix directory %s\n",cmdline.SysMatrixFileName);  

	/* read input arguments and parameters */
	readParamsSysMatrix(&cmdline, &imgparams, &sinoparams);
        
    	if(sinoparams.NViews%pieceLength!=0){
		fprintf(stderr, "NViews mod pieceLength must be 0.\n");
                exit(-1);
        }        
    
    	fprintf(stdout, "\nGenerating System Matrix...\n");

	unsigned int sum=0;
	int i,j,p,t;
	int Ny=imgparams.Ny;
	int Nx=imgparams.Nx;
	
	for(i=0;i<Ny;i+=(SVLength*2-overlappingDistance1)){
	  	for(j=0;j<Nx;j+=(SVLength*2-overlappingDistance2)){
	    		sum++;
	  	}
	}
	
    	fprintf(stdout, "Ny is %d Nx %d sum %d channels %d views %d\n",Ny,Nx,sum,sinoparams.NChannels,sinoparams.NViews);	

	order = (int *)_mm_malloc(sum*sizeof(int),64);
	
	t=0;	
	
	for(i=0;i<Ny;i+=(SVLength*2-overlappingDistance1)){
		  for(j=0;j<Nx;j+=(SVLength*2-overlappingDistance2)){
		    	order[t]=i*Nx+j;  /* order is the first voxel coordinate, not the center */
		    	t++;
		  }
	}	


    	bandMinMap = (struct minStruct *)get_spc(sum,sizeof(struct minStruct));
    	bandMaxMap = (struct maxStruct *)get_spc(sum,sizeof(struct maxStruct));
    	A_Padded_Map = (struct AValues_char **)multialloc(sizeof(struct AValues_char), 2, sum, (SVLength*2+1)*(SVLength*2+1));	

	float* max_num_pointer = (float *)malloc(Ny*Nx*sizeof(float));	
	char **ImageReconMask = (char **)multialloc(sizeof(char), 2, Ny, Nx); 
	float *initialError = (float *)malloc(sizeof(float)*sinoparams.NViews*sinoparams.NChannels);    	   
    /* Compute Pixel-Detector Profile */
    	PixelDetector_profile = ComputePixelProfile3DParallel(&sinoparams, &imgparams);  /* pixel-detector profile function */
    /* Compute Forward Matrix */
    		

    	fprintf(stdout, "after allocation \n");	

    	Deltaxy = imgparams.Deltaxy;
    	ROIRadius = imgparams.ROIRadius;
    
    	x_0 = -(Nx-1)*Deltaxy/2;
    	y_0 = -(Ny-1)*Deltaxy/2;
    	Nxy = Nx*Ny;
    
    	if (ROIRadius < 0.0)
    	{
        	printf("Error in GenImageReconMask : Invalid Value for Radius of Reconstruction \n");
        	exit(-1);
    	}
    	else
   	{
        	R_sq_max = ROIRadius*ROIRadius;
        
        	for (jy = 0; jy < Ny; jy++)
        	{
            		y = y_0 + jy*Deltaxy;
            		yy = y*y;
            		for (jx = 0; jx < Nx; jx++)
            		{
                		x = x_0 + jx*Deltaxy;
                		R_sq = x*x + yy;
                		if (R_sq > R_sq_max)
                		{
                    			ImageReconMask[jy][jx] = 0;
                    		}	
                		else
                		{
                    			ImageReconMask[jy][jx] = 1;
                		}
            		}
        	}
    	}



	
	A_comp(bandMinMap,bandMaxMap,A_Padded_Map,max_num_pointer,&sinoparams,sum,ImageReconMask,order,&imgparams,PixelDetector_profile,cmdline.SysMatrixFileName);
	    
    
    	fprintf(stdout, "after A comp \n");	

	//float InitValue = reconparams.MuWater;    
	float InitValue = MUWATER;    
    
	forwardProject2D(initialError, InitValue, max_num_pointer,A_Padded_Map,bandMinMap, &sinoparams, &imgparams);	    
    
	char fname[200];
    	sprintf(fname,"%s.initialError",cmdline.SysMatrixFileName);
	writeErrorSinogram(fname,initialError,&sinoparams);   

    	fprintf(stdout, "after E writing \n");	 
    
    /*A = ComputeSysMatrix3DParallel(&sinoparams, &imgparams, PixelDetector_profile);*/
    
    
    /* Write out System Matrix */
    /*
    if(WriteSysMatrix2D(cmdline.SysMatrixFileName, A))
    {  fprintf(stderr, "Error in writing out System Matrix to file %s through function WriteSysMatrix2D \n", cmdline.SysMatrixFileName);
       exit(-1);
    }
    */
    /* Free System Matrix */
    /*
    if(FreeSysMatrix2D(A))
    {  fprintf(stderr, "Error System Matrix memory could not be freed through function FreeSysMatrix2D \n");
       exit(-1);
    }
    */
    	_mm_free(order);
        for(j=0;j<sum;j++){
            	free((void *)bandMinMap[j].bandMin);
            	free((void *)bandMaxMap[j].bandMax);
        }

    	fprintf(stdout, "after bandMin and bandMax \n");	 

    	free((void *)bandMinMap);
    	free((void *)bandMaxMap);    	
    	for(i=0;i<sum;i++){
        	for(j=0;j<((2*SVLength+1)*(2*SVLength+1));j++){
            		if(A_Padded_Map[i][j].length>0){
                		free((void *)A_Padded_Map[i][j].val);
                		free((void *)A_Padded_Map[i][j].pieceWiseMin);
                		free((void *)A_Padded_Map[i][j].pieceWiseWidth);
            		}
        	}
    	}
    	fprintf(stdout, "A_Padded_Map\n");	     	
    	multifree(A_Padded_Map,2);
    	multifree(ImageReconMask,2);
    	free((void *)initialError);    	
    	free((void *)max_num_pointer);    	
    
	return 0;
}
