
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "mbir_ct.h"
#include "MBIRModularDefs.h"
#include "MBIRModularUtils.h"
#include "allocate.h"
#include "A_comp.h"

/* Pixel profile params */
#define WIDE_BEAM   /* Finite element analysis of detector channel, accounts for sensitivity variation across its aperture */
#define LEN_PIX 511 /* determines the spatial resolution for Detector-Pixel computation. Higher LEN_PIX, higher resolution */
                    /* In this implementation, spatial resolution is : [2*PixelDimension/LEN_PIX]^(-1) */
#define LEN_DET 101 /* Each detector channel is "split" into LEN_DET sub-elements ... */
                    /* to account for detector sensitivity variation across its aperture */

/******************************************************************/
/* Compute line segment lengths through a pixel for the given set */
/* of sinogram angles and for multiple displacements (LEN_PIX)    */
/******************************************************************/
/* Fill profiles of pixels from all angles
   ONLY FOR SQUARE PIXELS NOW   
   Each profile assumed 2 pixels wide
   This is used to speed calculation (too many to store) of the
   entries of the projection matrix.
   Values are scaled by "scale"
 */


/* The System matrix does not vary with slice for 3-D Parallel Geometry */
/* So, the method of compuatation is same as that of 2-D Parallel Geometry */

float **ComputePixelProfile3DParallel(
        struct SinoParams3DParallel *sinoparams,
        struct ImageParams3D *imgparams)
{
	int i, j;
	float pi, ang, d1, d2, t, t_1, t_2, t_3, t_4, maxval, rc, DeltaPix;
	float **pix_prof ; /* Detector-pixel profile, indexed by view angle and detector-pixel displacement */

	DeltaPix = imgparams->Deltaxy;

	pix_prof = (float **)get_img(LEN_PIX, sinoparams->NViews, sizeof(float));

	pi = PI; /* defined in MBIRModularUtils_2D.h */
	rc = sin(pi/4.0); /* Constant sin(pi/4) */
    
	/* Compute 3 parameters of the profile function */
	/* Here the corresponding parameters are : maxval, d1 and d2 */
    
	for (i = 0; i < sinoparams->NViews; i++)
	{
		ang = sinoparams->ViewAngles[i];

		while(ang >= pi/2.0) ang -= pi/2.0;
		while(ang < 0.0) ang += pi/2.0;

		if (ang <= pi/4.0)
			maxval = DeltaPix/cos(ang);
		else
			maxval = DeltaPix/cos(pi/2.0-ang);

		d1 = rc*cos(pi/4.0-ang);
		d2 = rc*fabs(sin(pi/4.0-ang));

		t_1 = 1.0 - d1;
		t_2 = 1.0 - d2;
		t_3 = 1.0 + d2;
		t_4 = 1.0 + d1;
        
		/* Profile is a trapezoidal function of detector-pixel displacement*/
		for (j = 0; j < LEN_PIX; j++)
		{
			t = 2.0*j/(float)LEN_PIX;
			if(t <= t_1 || t > t_4) {
				pix_prof[i][j] = 0.0;
			}
			else if(t <= t_2) {
				pix_prof[i][j] = maxval*(t-t_1)/(t_2-t_1);
			}
			else if(t <= t_3) {
				pix_prof[i][j] = maxval;
			}
			else {
				pix_prof[i][j] = maxval*(t_4-t)/(t_4-t_3);
			}
		}
	}
    
	return pix_prof;
}

/* Compute the System Matrix column for a given pixel */

void A_comp_ij(
	int im_row,
	int im_col,
	struct SinoParams3DParallel *sinoparams,
	struct ImageParams3D *imgparams,
	float **pix_prof,
	struct ACol *A_col,float *A_Values)
{
	static int first_call=1;
	static int Ntheta, NChannels, N_x, N_y;
	static float DeltaChannel, DeltaPix, t_0, x_0, y_0, dprof[LEN_DET];
	int ind_min, ind_max, pr;
	int pix_prof_ind, i, proj_count;
	float Aval, t_min, t_max, ang, x, y;
	int k;
	float t, const1, const2, const3, const4;

	t = 0;
	pix_prof_ind = 0;

	if (first_call == 1)
	{
		first_call = 0;

		Ntheta = sinoparams->NViews;
		NChannels = sinoparams->NChannels;
		DeltaChannel = sinoparams->DeltaChannel;
		t_0 = -(sinoparams->NChannels-1)*sinoparams->DeltaChannel/2.0 - sinoparams->CenterOffset * sinoparams->DeltaChannel;

		N_x = imgparams->Nx;
		N_y = imgparams->Ny;
		DeltaPix = imgparams->Deltaxy;
		x_0 = -(N_x-1)*DeltaPix/2.0;
		y_0 = -(N_y-1)*DeltaPix/2.0;		

		/* delta profile */
		/*
		   for(k=0;k<LEN_DET;k++) 
		   dprof[k] =0;
		   dprof[(LEN_DET-1)/2]=1.0;
		 */

		/* square profile */
		for (k = 0; k < LEN_DET; k++)
		{
			dprof[k] = 1.0/(LEN_DET);
		}

		/* triangular profile */
		/*
		   float sum=0;
		   for(k=0;k<LEN_DET;k++) {
		   if(k<=(LEN_DET-1)/2)
		   dprof[k]=2.0*k/(LEN_DET-1);
		   else
		   dprof[k]=2.0*(LEN_DET-1-k)/(LEN_DET-1);
		   sum+=dprof[k];
		   }
		   for(k=0;k<LEN_DET;k++) 
		   dprof[k] /= sum;
		 */
	} 

	/* WATCH THIS; ONLY FOR SQUARE PIXELS NOW   */

	y = y_0 + im_row*DeltaPix;
	x = x_0 + im_col*DeltaPix;

	proj_count = 0;
	for (pr = 0; pr < Ntheta; pr++)
	{
		int countTemp=proj_count;
		int write=1;
	        int minCount=0;
		//pind = pr*NChannels;
		ang = sinoparams->ViewAngles[pr];

		//z = x*cos(ang) + y*sin(ang);
		//d1 = d_source + z;
		//d2 = d_detector - z;

		/* range for pixel profile.  Need profile to contain 2 pixel widths */
		t_min = y*cos(ang) - x*sin(ang) - DeltaPix;
		t_max = t_min + 2.0*DeltaPix;
		/* This also prevents over-reach (with rounding of negative numbers)  */
		if (t_max < t_0) 
		{
			A_col->countTheta[pr]=0;
			A_col->minIndex[pr]=0;
			continue;
		}

		/* Changed 10/2000 (TF): more general, old implementation assumed
		   that geom->dia=sino->NChannels*sino->delta_t; */ 

		/* Relevant detector indices */
		ind_min = ceil((t_min-t_0)/DeltaChannel - 0.5);
		ind_max = (t_max-t_0)/DeltaChannel + 0.5;

		/* Fix this 4/91 to prevent over-reach at ends  */
		ind_min = (ind_min<0) ? 0 : ind_min;
		ind_max = (ind_max>=NChannels) ? NChannels-1 : ind_max;

		const1 = t_0 - DeltaChannel/2.0;
		const2 = DeltaChannel/(float)(LEN_DET-1);
		const3 = DeltaPix - (y*cos(ang) - x*sin(ang));
		const4 = (float)(LEN_PIX-1)/(2.0*DeltaPix);

		for (i = ind_min; i <= ind_max; i++)
		{
			//ind = pind + i;

		#ifdef WIDE_BEAM
			/* step through values of detector profile, inner product with PIX prof */
			Aval = 0;
			for (k = 0; k < LEN_DET; k++)
			{
				t = const1 + (float)i*DeltaChannel + (float)k*const2;
				pix_prof_ind = (t+const3)*const4 +0.5;   /* +0.5 for rounding */
				if (pix_prof_ind >= 0 && pix_prof_ind < LEN_PIX)
				{
					Aval+= dprof[k]*pix_prof[pr][pix_prof_ind];
				}
			}
		#else
			/*** this block computes zero-beam-width projection model ****/
			int prof_ind = LEN_PIX*(t_0+i*DeltaChannel+const3)/(2.0*DeltaPix);

			if (prof_ind >= LEN_PIX || prof_ind < 0)
			{
				if (prof_ind == LEN_PIX)
				{
					prof_ind = LEN_PIX-1;
				}
				else if (prof_ind == -1)
				{
					prof_ind = 0;
				}
				else
				{
					fprintf(stderr,"\nExiting Program: input parameters inconsistant\n");
					exit(-1);
				}
			}
			Aval = pix_prof[pr][prof_ind];
		#endif

			if (Aval > 0.0)  
			{
				/*XW: record the starting position for each view. */

				if(write==1){
				  minCount=i;
				  write=0;
				}
				A_Values[proj_count] = Aval;
				proj_count++;
			}
		}
		A_col->countTheta[pr]=proj_count-countTemp;
		A_col->minIndex[pr]=minCount;
	} 

	A_col->n_index = proj_count;
}


void A_piecewise(
	struct pointerAddress twoAddresses,
	struct AValues_char ** A_Padded_Map,
	float *max_num_pointer,
	struct SVParams svpar,
	struct SinoParams3DParallel *sinoparams,
	char *recon_mask,
	struct ImageParams3D *imgparams)
{
	struct ACol ** A_Col_pointer=twoAddresses.addressA;
	struct AValues_char ** A_Values_pointer=twoAddresses.addressB;
	
	int j,jj,i,jy,jx,p,q,t;
	int Nx = imgparams->Nx;
	int Ny = imgparams->Ny;	
	int SVLength = svpar.SVLength; 
	struct minStruct * bandMinMap = svpar.bandMinMap;
	struct maxStruct * bandMaxMap = svpar.bandMaxMap;
	int sum = svpar.Nsv;
	int pieceLength = svpar.pieceLength;

	#ifdef USING_ICC
	int *order = (int *)_mm_malloc(svpar.Nsv*sizeof(int),64);
	#else
	int *order = (int *) aligned_alloc(64,svpar.Nsv*sizeof(int));
	#endif

	t=0;
	for(i=0;i<Ny;i+=(svpar.SVLength*2-svpar.overlap))
	for(j=0;j<Nx;j+=(svpar.SVLength*2-svpar.overlap)){
		order[t]=i*Nx+j;  /* order is the first voxel coordinate, not the center */
		t++;
	}

	for(jj=0;jj<sum;jj++)
	for(i=0;i<(SVLength*2+1)*(SVLength*2+1);i++) {
		A_Padded_Map[jj][i].val=NULL;
		A_Padded_Map[jj][i].length=0;
	}

        for (jj = 0; jj < sum; jj++)
        {
            jy = order[jj] / Nx;
            jx = order[jj] % Nx;
            int countNumber=0;
            int radius =SVLength;
            int coordinateSize=1;
            if(radius!=0)
            coordinateSize=(2*radius+1)*(2*radius+1);
            int k_newCoordinate[coordinateSize];
            int j_newCoordinate[coordinateSize];
            int j_newAA=0;
            int k_newAA=0;

            for(j_newAA=jy;j_newAA<=(jy+2*radius);j_newAA++)
            for(k_newAA=jx;k_newAA<=(jx+2*radius);k_newAA++)
            if(j_newAA>=0 && k_newAA >=0 && j_newAA <Ny && k_newAA < Nx){
                if(recon_mask[j_newAA*Nx + k_newAA]){
                    if(A_Col_pointer[j_newAA][k_newAA].n_index >0){
                        j_newCoordinate[countNumber]=j_newAA;
                        k_newCoordinate[countNumber]=k_newAA;
                        countNumber++;
                    }
                } 
            }

            int bandMin[sinoparams->NViews]__attribute__((aligned(64)));
            int bandMax[sinoparams->NViews]__attribute__((aligned(64)));

            for(p=0;p< sinoparams->NViews;p++){
                bandMin[p]=sinoparams->NChannels;
            }

            for(i=0;i<countNumber;i++)
            {
                int j_new= j_newCoordinate[i];
                int k_new= k_newCoordinate[i];
                for(p=0;p< sinoparams->NViews;p++)
                {
                    if(A_Col_pointer[j_new][k_new].minIndex[p]==0 && A_Col_pointer[j_new][k_new].countTheta[p]==0)
                    {
			if(p!=0){
                    		A_Col_pointer[j_new][k_new].minIndex[p]=A_Col_pointer[j_new][k_new].minIndex[p-1];
			}
			else{
				int k=0;
				while(A_Col_pointer[j_new][k_new].minIndex[k]==0){
					k++;
				}
				A_Col_pointer[j_new][k_new].minIndex[p]=A_Col_pointer[j_new][k_new].minIndex[k];
			}
                    }
                    else if(A_Col_pointer[j_new][k_new].minIndex[p]==(sinoparams->NChannels-1) && A_Col_pointer[j_new][k_new].countTheta[p]==0)
                    {
			if(p!=0){
                    		A_Col_pointer[j_new][k_new].minIndex[p]=A_Col_pointer[j_new][k_new].minIndex[p-1];
			}
			else{
				int k=0;
				while(A_Col_pointer[j_new][k_new].minIndex[k]==(sinoparams->NChannels-1)){
					k++;
				}
				A_Col_pointer[j_new][k_new].minIndex[p]=A_Col_pointer[j_new][k_new].minIndex[k];
			}
                    }	                    
                    	
                    if(A_Col_pointer[j_new][k_new].minIndex[p]<bandMin[p])
                    	bandMin[p]=A_Col_pointer[j_new][k_new].minIndex[p];
                }
            }

            for(p=0;p< sinoparams->NViews;p++){
                bandMax[p]=bandMin[p];
            }

            for(i=0;i<countNumber;i++){
                int j_new= j_newCoordinate[i];
                int k_new= k_newCoordinate[i];
                for(p=0;p< sinoparams->NViews;p++){
                    if((A_Col_pointer[j_new][k_new].minIndex[p]+A_Col_pointer[j_new][k_new].countTheta[p])>bandMax[p])
                    bandMax[p]=A_Col_pointer[j_new][k_new].minIndex[p]+A_Col_pointer[j_new][k_new].countTheta[p];
                }
            }
            
            int bandWidthTemp[sinoparams->NViews]__attribute__((aligned(64)));
            int bandWidth[(sinoparams->NViews)/pieceLength]__attribute__((aligned(64)));	            
            
            #pragma vector aligned
            for(p=0;p< sinoparams->NViews;p++){   		
                bandWidthTemp[p]=bandMax[p]-bandMin[p];
            }

            for (p = 0; p < (sinoparams->NViews)/pieceLength; p++)
            {
                int bandWidthMax=bandWidthTemp[p*pieceLength];
                for(t=0;t<pieceLength;t++){
                    	if(bandWidthTemp[p*pieceLength+t]>bandWidthMax){
                    		bandWidthMax=bandWidthTemp[p*pieceLength+t];
                	}
                	bandWidth[p]=bandWidthMax;
            	}
            }	
            		
            #pragma vector aligned            		
            for(p=0;p< sinoparams->NViews;p++){
		if((bandMin[p]+bandWidth[p/pieceLength])>=sinoparams->NChannels){
                	bandMin[p]=sinoparams->NChannels-bandWidth[p/pieceLength];
                }
            }            

            #ifdef USING_ICC
            _intel_fast_memcpy(&bandMinMap[jj].bandMin[0],&bandMin[0],sizeof(int)*(sinoparams->NViews));
            _intel_fast_memcpy(&bandMaxMap[jj].bandMax[0],&bandMax[0],sizeof(int)*(sinoparams->NViews));
            #else
            memcpy(&bandMinMap[jj].bandMin[0],&bandMin[0],sizeof(int)*(sinoparams->NViews));
            memcpy(&bandMaxMap[jj].bandMax[0],&bandMax[0],sizeof(int)*(sinoparams->NViews));
            #endif
            
            int piecewiseMinArray[countNumber][(sinoparams->NViews)/pieceLength]__attribute__((aligned(64)));
            int piecewiseMaxArray[countNumber][(sinoparams->NViews)/pieceLength]__attribute__((aligned(64)));
            int piecewiseWidth[countNumber][(sinoparams->NViews)/pieceLength]__attribute__((aligned(64)));
            int totalSumArray[countNumber];

	    for(i=0;i<countNumber;i++){
                for (p = 0; p < (sinoparams->NViews)/pieceLength; p++){
                	piecewiseMinArray[i][p]=0;
                	piecewiseMaxArray[i][p]=0;
                }
            }    

            for(i=0;i<countNumber;i++){
                int j_new= j_newCoordinate[i];
                int k_new= k_newCoordinate[i];
                totalSumArray[i]=0;
                for (p = 0; p < (sinoparams->NViews)/pieceLength; p++)
                {
                    piecewiseMinArray[i][p]=A_Col_pointer[j_new][k_new].minIndex[p*pieceLength]-bandMin[p*pieceLength];
                    piecewiseMaxArray[i][p]=A_Col_pointer[j_new][k_new].minIndex[p*pieceLength]-bandMin[p*pieceLength]+A_Col_pointer[j_new][k_new].countTheta[p*pieceLength];
                    for(t=0;t<pieceLength;t++)
                    {
                        if(piecewiseMinArray[i][p]>(A_Col_pointer[j_new][k_new].minIndex[p*pieceLength+t]-bandMin[p*pieceLength+t]))
                            piecewiseMinArray[i][p]=A_Col_pointer[j_new][k_new].minIndex[p*pieceLength+t]-bandMin[p*pieceLength+t];
                        if((A_Col_pointer[j_new][k_new].minIndex[p*pieceLength+t]-bandMin[p*pieceLength+t]+A_Col_pointer[j_new][k_new].countTheta[p*pieceLength+t])>piecewiseMaxArray[i][p])
                            piecewiseMaxArray[i][p]=A_Col_pointer[j_new][k_new].minIndex[p*pieceLength+t]-bandMin[p*pieceLength+t]+A_Col_pointer[j_new][k_new].countTheta[p*pieceLength+t];
                    }
                }
            }

            for(i=0;i<countNumber;i++)
            for(p = 0; p < (sinoparams->NViews)/pieceLength; p++)
                piecewiseWidth[i][p]=piecewiseMaxArray[i][p]-piecewiseMinArray[i][p];

            for(i=0;i<countNumber;i++)
            {
                #pragma vector aligned
                for (p = 0; p < (sinoparams->NViews)/pieceLength; p++)
                {
                    totalSumArray[i]+=(piecewiseWidth[i][p])*pieceLength;
                }
            }

            unsigned char ** AMatrixPadded= (unsigned char **)malloc(countNumber*sizeof(unsigned char *));
            unsigned char ** AMatrixPaddedTranspose=(unsigned char **)malloc(countNumber*sizeof(unsigned char *));

            for(i=0;i<countNumber;i++){
                AMatrixPadded[i]=(unsigned char *)malloc(totalSumArray[i]*sizeof(unsigned char));
                AMatrixPaddedTranspose[i]=(unsigned char *)malloc(totalSumArray[i]*sizeof(unsigned char));
            }

            unsigned char*    newProjectionValueArrayPointer=    &A_Values_pointer[0][0].val[0];

            for(i=0;i<countNumber;i++){
                int j_new= j_newCoordinate[i];
                int k_new= k_newCoordinate[i];
                unsigned char * A_padded_pointer=&AMatrixPadded[i][0];
                newProjectionValueArrayPointer=    &A_Values_pointer[j_new][k_new].val[0];
                for (p = 0; p < sinoparams->NViews; p++)
                {
                    #pragma vector aligned
                    for(t=0;t<(A_Col_pointer[j_new][k_new].minIndex[p]-piecewiseMinArray[i][p/pieceLength]-bandMin[p]);t++){
			*A_padded_pointer=0;
                        A_padded_pointer++;
                    }
                    #pragma vector aligned
                    for(t=0;t<A_Col_pointer[j_new][k_new].countTheta[p];t++){
                        *A_padded_pointer = *newProjectionValueArrayPointer;
                        A_padded_pointer++;
                        newProjectionValueArrayPointer++;
                    }
                    #pragma vector aligned
                    for(t=0;t<(piecewiseMaxArray[i][p/pieceLength]-A_Col_pointer[j_new][k_new].minIndex[p]-A_Col_pointer[j_new][k_new].countTheta[p]+bandMin[p]);t++){
			*A_padded_pointer=0;
                        A_padded_pointer++;
                    }
                }
            }

            for(i=0;i<countNumber;i++){
                unsigned char * A_padded_pointer=&AMatrixPadded[i][0];
                unsigned char * A_padd_Tranpose_pointer =&AMatrixPaddedTranspose[i][0];
                for (p = 0; p < (sinoparams->NViews)/pieceLength; p++)
                {
                    for(q=0;q<piecewiseWidth[i][p];q++){
                        for(t=0;t<pieceLength;t++){
                            A_padd_Tranpose_pointer[q*pieceLength+t]=A_padded_pointer[t*piecewiseWidth[i][p]+q];
                        }
                    }
                    A_padded_pointer+=piecewiseWidth[i][p]*pieceLength;
                    A_padd_Tranpose_pointer+=piecewiseWidth[i][p]*pieceLength;
                }
            }

            for(i=0;i<countNumber;i++){
                int j_new= j_newCoordinate[i];
                int k_new= k_newCoordinate[i];                
                int theVoxelPosition=(j_new-jy)*(2*SVLength+1)+(k_new-jx);
                	    
                A_Padded_Map[jj][theVoxelPosition].val= (unsigned char *)get_spc(totalSumArray[i], sizeof(unsigned char));
                A_Padded_Map[jj][theVoxelPosition].pieceWiseMin = (int *)get_spc((sinoparams->NViews)/pieceLength,sizeof(int));
                A_Padded_Map[jj][theVoxelPosition].pieceWiseWidth = (int *)get_spc((sinoparams->NViews)/pieceLength,sizeof(int));
                A_Padded_Map[jj][theVoxelPosition].length=totalSumArray[i];
                #ifdef USING_ICC
                _intel_fast_memcpy(&A_Padded_Map[jj][theVoxelPosition].val[0],&AMatrixPaddedTranspose[i][0],sizeof(unsigned char)*totalSumArray[i]);
                _intel_fast_memcpy(&A_Padded_Map[jj][theVoxelPosition].pieceWiseMin[0],&piecewiseMinArray[i][0],sizeof(int)*(sinoparams->NViews)/pieceLength);
                _intel_fast_memcpy(&A_Padded_Map[jj][theVoxelPosition].pieceWiseWidth[0],&piecewiseWidth[i][0],sizeof(int)*(sinoparams->NViews)/pieceLength);
                #else
                memcpy(&A_Padded_Map[jj][theVoxelPosition].val[0],&AMatrixPaddedTranspose[i][0],sizeof(unsigned char)*totalSumArray[i]);
                memcpy(&A_Padded_Map[jj][theVoxelPosition].pieceWiseMin[0],&piecewiseMinArray[i][0],sizeof(int)*(sinoparams->NViews)/pieceLength);
                memcpy(&A_Padded_Map[jj][theVoxelPosition].pieceWiseWidth[0],&piecewiseWidth[i][0],sizeof(int)*(sinoparams->NViews)/pieceLength);
                #endif
            }

            for(i=0;i<countNumber;i++){
                free((void *)AMatrixPadded[i]);
                free((void *)AMatrixPaddedTranspose[i]);
            }

            free((void *)AMatrixPadded);
            free((void *)AMatrixPaddedTranspose);
        }
    
	#ifdef USING_ICC
	_mm_free(order);
	#else
	free(order);
	#endif
}




/* Compute Entire System Matrix */
/* The System matrix does not vary with slice for 3-D Parallel Geometry */
/* So, the method of compuatation is same as that of 2-D Parallel Geometry */
       
struct pointerAddress A_comp(
	struct AValues_char ** A_Padded_Map,
	float * max_num_pointer,
	struct SVParams svpar,
	struct SinoParams3DParallel *sinoparams,
	char *recon_mask,
	struct ImageParams3D *imgparams)
{
	int i, j, r;
	int col_length, n_x, n_y;

	struct ACol A_col_sgl;
	float* A_Values_sgl;	
	struct pointerAddress address_arr;	

	col_length = sinoparams->NChannels*sinoparams->NViews;

	A_Values_sgl = (float *)get_spc(col_length, sizeof(float));
	A_col_sgl.countTheta=(unsigned char *)get_spc(sinoparams->NViews,sizeof(unsigned char));
	A_col_sgl.minIndex=(int *)get_spc(sinoparams->NViews,sizeof(int));

	n_x = imgparams->Nx;
	n_y = imgparams->Ny;
	address_arr.addressA = (struct ACol **)multialloc(sizeof(struct ACol), 2, n_y, n_x);
	address_arr.addressB = (struct AValues_char **)multialloc(sizeof(struct AValues_char), 2, n_y, n_x);

	float **pix_prof = ComputePixelProfile3DParallel(sinoparams,imgparams);

	for (i = 0; i < n_y; i++)
	for (j = 0; j < n_x; j++)
	{
		A_comp_ij(i, j, sinoparams, imgparams, pix_prof,  &A_col_sgl,A_Values_sgl);  
		address_arr.addressA[i][j].n_index = A_col_sgl.n_index;
		address_arr.addressA[i][j].countTheta=(unsigned char *)get_spc(sinoparams->NViews,sizeof(unsigned char));
		address_arr.addressA[i][j].minIndex=(int *)get_spc(sinoparams->NViews,sizeof(int));
		address_arr.addressB[i][j].val = (unsigned char *)get_spc(A_col_sgl.n_index, sizeof(unsigned char));

		float max=A_Values_sgl[0];
		for (r = 0; r < A_col_sgl.n_index; r++)
		{
			if(A_Values_sgl[r]>max)
				max = A_Values_sgl[r];
		}
		max_num_pointer[i*n_x+j]=max;

		for (r = 0; r < A_col_sgl.n_index; r++)
			address_arr.addressB[i][j].val[r] = (unsigned char)((A_Values_sgl[r])/max*255+0.5);

		for (r = 0; r < sinoparams->NViews; r++) {
			address_arr.addressA[i][j].countTheta[r]=A_col_sgl.countTheta[r];
			address_arr.addressA[i][j].minIndex[r]=A_col_sgl.minIndex[r];
		}

	}
	free((void *)A_Values_sgl);
	free((void *)A_col_sgl.countTheta);
	free((void *)A_col_sgl.minIndex);

	A_piecewise(address_arr,A_Padded_Map,max_num_pointer,svpar,sinoparams,recon_mask,imgparams);
	
	for (i = 0; i < n_y; i++)
	for (j = 0; j < n_x; j++) {
		free((void *)address_arr.addressA[i][j].countTheta);
		free((void *)address_arr.addressA[i][j].minIndex);
		free((void *)address_arr.addressB[i][j].val);
	}
    	multifree(address_arr.addressA,2);
    	multifree(address_arr.addressB,2);	    		
	return address_arr;
	
}



void readAmatrix(
	char *fname,
	struct AValues_char ** A_Padded_Map, 
	float * max_num_pointer,
	struct ImageParams3D *imgparams,
	struct SinoParams3DParallel *sinoparams,
	struct SVParams svpar)
{
	FILE *fp;
	int i, j, Nx, Ny;
	int M_nonzero;
	int SVLength=svpar.SVLength;
	int sum=svpar.Nsv;
	struct minStruct * bandMinMap = svpar.bandMinMap;
	struct maxStruct * bandMaxMap = svpar.bandMaxMap;
	int pieceLength = svpar.pieceLength;

	Nx = imgparams->Nx;
	Ny = imgparams->Ny;

	if ((fp = fopen(fname, "rb")) == NULL)
	{
		fprintf(stderr, "ERROR in readAmatrix: can't open file %s.\n", fname);
		exit(-1);
	}

	for (i =0; i < sum ; i ++ ){
		if(fread(bandMinMap[i].bandMin,sizeof(int),sinoparams->NViews,fp) < sinoparams->NViews) {
			fprintf(stderr, "ERROR in readAmatrix: %s terminated early.\n", fname);
			exit(-1);
		}
		if(fread(bandMaxMap[i].bandMax,sizeof(int),sinoparams->NViews,fp) < sinoparams->NViews) {
			fprintf(stderr, "ERROR in readAmatrix: %s terminated early.\n", fname);
			exit(-1);
		}

		for ( j = 0; j< (SVLength*2+1)*(SVLength*2+1) ; j ++)
		{
			if(fread(&M_nonzero, sizeof(int), 1, fp) < 1) {
				fprintf(stderr, "ERROR in readAmatrix: %s terminated early.\n", fname);
				exit(-1);
			}
			A_Padded_Map[i][j].length = M_nonzero;
			if(M_nonzero > 0){
				A_Padded_Map[i][j].val = (unsigned char *)get_spc(M_nonzero, sizeof(unsigned char));
				A_Padded_Map[i][j].pieceWiseWidth=(int *)get_spc(sinoparams->NViews/pieceLength,sizeof(int));
				A_Padded_Map[i][j].pieceWiseMin=(int *)get_spc(sinoparams->NViews/pieceLength,sizeof(int));
			}
			if(M_nonzero > 0){
				if(fread(A_Padded_Map[i][j].val, sizeof(unsigned char), M_nonzero, fp) < M_nonzero) {
					fprintf(stderr, "ERROR in readAmatrix: %s terminated early.\n", fname);
					exit(-1);
				}
				int num = (sinoparams->NViews)/pieceLength;
				if(fread(A_Padded_Map[i][j].pieceWiseMin,sizeof(int),num,fp) < num) {
					fprintf(stderr, "ERROR in readAmatrix: %s terminated early.\n", fname);
					exit(-1);
				}
				if(fread(A_Padded_Map[i][j].pieceWiseWidth,sizeof(int),num,fp) < num) {
					fprintf(stderr, "ERROR in readAmatrix: %s terminated early.\n", fname);
					exit(-1);
				}
			}	
		}
	}
	
	if(fread(&max_num_pointer[0],sizeof(float),Nx*Ny,fp) < Nx*Ny) {
		fprintf(stderr, "ERROR in readAmatrix: %s terminated early.\n", fname);
		exit(-1);
	}

	fclose(fp);
}

void writeAmatrix(
	char *fname,
	struct AValues_char ** A_Padded_Map, 
	float * max_num_pointer,
	struct ImageParams3D *imgparams,
	struct SinoParams3DParallel *sinoparams,
	struct SVParams svpar)
{
	FILE *fp;
	int i,j, Nx, Ny;
	int M_nonzero;
	int SVLength = svpar.SVLength;
	int sum = svpar.Nsv;
	struct minStruct * bandMinMap = svpar.bandMinMap;
	struct maxStruct * bandMaxMap = svpar.bandMaxMap;
	int pieceLength = svpar.pieceLength;

	Nx = imgparams->Nx;
	Ny = imgparams->Ny;

	if ((fp = fopen(fname, "wb")) == NULL)
	{
		fprintf(stderr, "ERROR in writeAmatrix: can't open file %s.\n", fname);
		exit(-1);
	}

	for (i =0; i < sum ; i ++ ){
		fwrite(bandMinMap[i].bandMin,sizeof(int),sinoparams->NViews,fp);
		fwrite(bandMaxMap[i].bandMax,sizeof(int),sinoparams->NViews,fp);
		for ( j = 0; j< (SVLength*2+1)*(SVLength*2+1) ; j ++){
			M_nonzero =  A_Padded_Map[i][j].length;
			fwrite(&M_nonzero, sizeof(int), 1, fp);
			if(M_nonzero > 0){
				fwrite(A_Padded_Map[i][j].val, sizeof(unsigned char), M_nonzero, fp);			
				fwrite(A_Padded_Map[i][j].pieceWiseMin,sizeof(int),(sinoparams->NViews)/pieceLength,fp);
				fwrite(A_Padded_Map[i][j].pieceWiseWidth,sizeof(int),(sinoparams->NViews)/pieceLength,fp);
			} 	    										
		}
	}
	fwrite(&max_num_pointer[0],sizeof(float),Nx*Ny,fp);
	fclose(fp);
}


