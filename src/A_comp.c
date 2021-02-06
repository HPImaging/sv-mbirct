
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "MBIRModularDefs.h"
#include "MBIRModularUtils.h"
#include "allocate.h"
#include "initialize.h"
#include "A_comp.h"
#ifndef MSVC    /* not included in MS Visual C++ */
#include <sys/time.h>
#endif

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
    float pi, ang, d1, d2, t, t_1, t_2, t_3, t_4, maxval, rc;
    float **pix_prof ; /* Detector-pixel profile, indexed by view angle and detector-pixel displacement */

    pix_prof = (float **)get_img(LEN_PIX, sinoparams->NViews, sizeof(float));

    pi = PI;
    rc = sin(pi/4.0);

    /* Compute 3 parameters of the profile function */
    /* Here the corresponding parameters are : maxval, d1 and d2 */

    for (i = 0; i < sinoparams->NViews; i++)
    {
        ang = sinoparams->ViewAngles[i];

        while(ang >= pi/2.0)
            ang -= pi/2.0;
        while(ang < 0.0)
            ang += pi/2.0;

        if (ang <= pi/4.0)
            maxval = imgparams->Deltaxy/cos(ang);
        else
            maxval = imgparams->Deltaxy/cos(pi/2.0-ang);

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
            if(t <= t_1 || t > t_4)
                pix_prof[i][j] = 0.0;
            else if(t <= t_2)
                pix_prof[i][j] = maxval*(t-t_1)/(t_2-t_1);
            else if(t <= t_3)
                pix_prof[i][j] = maxval;
            else
                pix_prof[i][j] = maxval*(t_4-t)/(t_4-t_3);
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
    static float dprof[LEN_DET];
    float t_0, x_0, y_0;
    int i, k, pr, ind_min, ind_max, pix_prof_ind, proj_count;
    float Aval, t_min, t_max, ang, x, y;
    float t, const1, const2, const3, const4;

    float Deltaxy = imgparams->Deltaxy;
    int NChannels = sinoparams->NChannels;
    float DeltaChannel = sinoparams->DeltaChannel;

    if (first_call == 1)
    {
        first_call = 0;

        /* delta profile */
        /*
        for(k=0;k<LEN_DET;k++)
        dprof[k] =0;
        dprof[(LEN_DET-1)/2]=1.0;
        */

        /* square profile */
        for (k = 0; k < LEN_DET; k++)
            dprof[k] = 1.0/(LEN_DET);

        /* triangular profile */
        /*
        float sum=0;
        for(k=0;k<LEN_DET;k++)
        {
            if(k<=(LEN_DET-1)/2)
                dprof[k]=2.0*k/(LEN_DET-1);
            else
                dprof[k]=2.0*(LEN_DET-1-k)/(LEN_DET-1);
            sum += dprof[k];
        }
        for(k=0;k<LEN_DET;k++)
            dprof[k] /= sum;
        */
    }

    /* WATCH THIS; ONLY FOR SQUARE PIXELS NOW   */

    t_0 = -(NChannels-1)*DeltaChannel/2.0 - sinoparams->CenterOffset * DeltaChannel;
    x_0 = -(imgparams->Nx-1)*Deltaxy/2.0;
    y_0 = -(imgparams->Ny-1)*Deltaxy/2.0;

    y = y_0 + im_row*Deltaxy;
    x = x_0 + im_col*Deltaxy;

    proj_count = 0;
    for (pr = 0; pr < sinoparams->NViews; pr++)
    {
        int countTemp=proj_count;
        int write=1;
        int minCount=0;

        ang = sinoparams->ViewAngles[pr];

        /* range for pixel profile.  Need profile to contain 2 pixel widths */
        t_min = y*cos(ang) - x*sin(ang) - Deltaxy;
        t_max = t_min + 2.0*Deltaxy;

        /* Relevant detector indices */
        ind_min = ceil((t_min-t_0)/DeltaChannel - 0.5);
        ind_max= floor((t_max-t_0)/DeltaChannel + 0.5);

        /* move on if voxel clearly out of range of detectors */
        if(ind_max<0 || ind_min>NChannels-1)
        {
            A_col->countTheta[pr]=0;
            A_col->minIndex[pr]=0;
            continue;
        }

        /* Fix this 4/91 to prevent over-reach at ends  */
        ind_min = (ind_min<0) ? 0 : ind_min;
        ind_max = (ind_max>=NChannels) ? NChannels-1 : ind_max;

        const1 = t_0 - DeltaChannel/2.0;
        const2 = DeltaChannel/(float)(LEN_DET-1);
        const3 = Deltaxy - (y*cos(ang) - x*sin(ang));
        const4 = (float)(LEN_PIX-1)/(2.0*Deltaxy);

        for (i = ind_min; i <= ind_max; i++)
        {
            #ifdef WIDE_BEAM
            /* step through values of detector profile, inner product with PIX prof */
            Aval = 0;
            for (k = 0; k < LEN_DET; k++)
            {
                t = const1 + (float)i*DeltaChannel + (float)k*const2;
                pix_prof_ind = (t+const3)*const4 +0.5;   /* +0.5 for rounding */
                if (pix_prof_ind >= 0 && pix_prof_ind < LEN_PIX)
                    Aval += dprof[k]*pix_prof[pr][pix_prof_ind];
            }
            #else
            /*** this block computes zero-beam-width projection model ****/
            int prof_ind = LEN_PIX*(t_0+i*DeltaChannel+const3)/(2.0*Deltaxy);

            if (prof_ind >= LEN_PIX || prof_ind < 0)
            {
                if (prof_ind == LEN_PIX)
                    prof_ind = LEN_PIX-1;
                else if (prof_ind == -1)
                    prof_ind = 0;
                else
                {
                    fprintf(stderr,"A_comp_ij() Error: input parameters inconsistant\n");
                    exit(-1);
                }
            }
            Aval = pix_prof[pr][prof_ind];
            #endif

            if (Aval > 0.0)
            {
                /*XW: record the starting position for each view. */
                if(write==1) {
                    minCount=i;
                    write=0;
                }
                A_Values[proj_count] = Aval;
                proj_count++;
            }
        }
        /* data type of ACol.countTheta is typically unsigned char --check for overflow */
        const int overflow_val = 1 << (8*sizeof(chanwidth_t));
        if(proj_count-countTemp >= overflow_val) {
            fprintf(stderr,"A_comp_ij() Error: overflow detected--check voxel/detector dimensions\n");
            exit(-1);
        }
        A_col->countTheta[pr] = proj_count-countTemp;
        A_col->minIndex[pr] = minCount;
    }

    A_col->n_index = proj_count;
}


void A_piecewise(
    struct ACol **ACol_ptr,
    struct AValues_char **AVal_ptr,
    struct AValues_char **A_Padded_Map,
    struct SVParams svpar,
    struct SinoParams3DParallel *sinoparams,
    struct ImageParams3D *imgparams)
{

    int i,j,jj,p,t;
    int Nx = imgparams->Nx;
    int Ny = imgparams->Ny;
    int NViews = sinoparams->NViews;
    int NChannels = sinoparams->NChannels;
    int SVLength = svpar.SVLength;
    int pieceLength = svpar.pieceLength;
    int NViewSets = NViews/svpar.pieceLength;
    struct minStruct * bandMinMap = svpar.bandMinMap;
    struct maxStruct * bandMaxMap = svpar.bandMaxMap;

    int *order = (int *) mget_spc(svpar.Nsv,sizeof(int));
    t=0;
    for(i=0; i<Ny; i+=(SVLength*2-svpar.overlap))
    for(j=0; j<Nx; j+=(SVLength*2-svpar.overlap)) {
        order[t]=i*Nx+j;  /* order is the first voxel coordinate, not the center */
        t++;
    }

    for(jj=0; jj<svpar.Nsv; jj++)
    for(i=0; i<(SVLength*2+1)*(SVLength*2+1); i++) {
        A_Padded_Map[jj][i].val=NULL;
        A_Padded_Map[jj][i].length=0;
    }

    for(i=0; i<Ny; i++)
    for(j=0; j<Nx; j++)
    if(ACol_ptr[i][j].n_index > 0)
    for(p=0; p<NViews; p++)
    {
        if(ACol_ptr[i][j].minIndex[p]==0 && ACol_ptr[i][j].countTheta[p]==0)
        {
            if(p!=0)
                ACol_ptr[i][j].minIndex[p] = ACol_ptr[i][j].minIndex[p-1];
            else
            {
                t=0;
                while(ACol_ptr[i][j].minIndex[t] == 0 && t<NViews-1)
                    t++;
                ACol_ptr[i][j].minIndex[p] = ACol_ptr[i][j].minIndex[t];
            }
        }
    }

    /* Note this gets *slower* if running more than a few threads */
    #pragma omp parallel private(i,j,p,t) num_threads(4)
    {
        int jx,jy,q,jx_new,jy_new;

        channel_t *bandMin = (channel_t *) mget_spc(NViews,sizeof(channel_t));
        channel_t *bandMax = (channel_t *) mget_spc(NViews,sizeof(channel_t));
        channel_t *bandWidth=(channel_t *) mget_spc(NViews,sizeof(channel_t));
        channel_t *bandWidthPW = (channel_t *) mget_spc(NViewSets,sizeof(channel_t));

        int SVSize = (2*SVLength+1)*(2*SVLength+1);
        channel_t **piecewiseMin = (channel_t **)multialloc(sizeof(channel_t),2,SVSize,NViewSets);
        channel_t **piecewiseMax = (channel_t **)multialloc(sizeof(channel_t),2,SVSize,NViewSets);
        channel_t **piecewiseWidth=(channel_t **)multialloc(sizeof(channel_t),2,SVSize,NViewSets);
        int *jx_list = (int *) mget_spc(SVSize,sizeof(int));
        int *jy_list = (int *) mget_spc(SVSize,sizeof(int));
        int *totalSum = (int *) mget_spc(SVSize,sizeof(int));

        #pragma omp for schedule(static)
        for(jj=0; jj<svpar.Nsv; jj++)
        {
            jy = order[jj] / Nx;
            jx = order[jj] % Nx;
            SVSize = 0;

            for(jy_new=jy; jy_new<=(jy+2*SVLength); jy_new++)
            for(jx_new=jx; jx_new<=(jx+2*SVLength); jx_new++)
            if(jy_new<Ny && jx_new<Nx) {
                if(ACol_ptr[jy_new][jx_new].n_index >0) {
                    jy_list[SVSize] = jy_new;
                    jx_list[SVSize] = jx_new;
                    SVSize++;
                }
            }

            //channel_t bandMin[NViews]__attribute__((aligned(64)));
            //channel_t bandMax[NViews]__attribute__((aligned(64)));
            for(p=0; p< NViews; p++)
                bandMin[p] = NChannels;

            for(i=0; i<SVSize; i++)
            {
                jy_new = jy_list[i];
                jx_new = jx_list[i];
                for(p=0; p< NViews; p++)
                {
                    if(ACol_ptr[jy_new][jx_new].minIndex[p] < bandMin[p])
                        bandMin[p] = ACol_ptr[jy_new][jx_new].minIndex[p];
                }
            }

            for(p=0; p< NViews; p++)
                bandMax[p]=bandMin[p];

            for(i=0; i<SVSize; i++)
            {
                jy_new = jy_list[i];
                jx_new = jx_list[i];
                for(p=0; p< NViews; p++) {
                    if((ACol_ptr[jy_new][jx_new].minIndex[p] + ACol_ptr[jy_new][jx_new].countTheta[p]) > bandMax[p])
                        bandMax[p] = ACol_ptr[jy_new][jx_new].minIndex[p] + ACol_ptr[jy_new][jx_new].countTheta[p];
                }
            }

            //channel_t bandWidth[NViews]__attribute__((aligned(64)));
            //channel_t bandWidthPW[NViewSets]__attribute__((aligned(64)));
            //#pragma vector aligned
            for(p=0; p< NViews; p++)
                bandWidth[p] = bandMax[p]-bandMin[p];

            for (p=0; p < NViewSets; p++)
            {
                int bandWidthMax = bandWidth[p*pieceLength];
                for(t=0; t<pieceLength; t++) {
                    if(bandWidth[p*pieceLength+t] > bandWidthMax) {
                        bandWidthMax = bandWidth[p*pieceLength+t];
                    }
                }
                bandWidthPW[p] = bandWidthMax;
            }

            //#pragma vector aligned
            for(p=0; p< NViews; p++) {
                if((bandMin[p]+bandWidthPW[p/pieceLength]) >= NChannels)
                    bandMin[p] = NChannels - bandWidthPW[p/pieceLength];
            }

            memcpy(&bandMinMap[jj].bandMin[0],&bandMin[0],sizeof(channel_t)*NViews);
            memcpy(&bandMaxMap[jj].bandMax[0],&bandMax[0],sizeof(channel_t)*NViews);

            //int totalSum[SVSize]__attribute__((aligned(64)));
            for(i=0; i<SVSize; i++)
            {
                jy_new = jy_list[i];
                jx_new = jx_list[i];
                for (p=0; p < NViewSets; p++)
                {
                    int pwMin = (int)ACol_ptr[jy_new][jx_new].minIndex[p*pieceLength]-(int)bandMin[p*pieceLength];
                    int pwMax = pwMin + ACol_ptr[jy_new][jx_new].countTheta[p*pieceLength];
                    for(t=0; t<pieceLength; t++)
                    {
                        int idx0 = (int)ACol_ptr[jy_new][jx_new].minIndex[p*pieceLength+t]-(int)bandMin[p*pieceLength+t];
                        int idx1 = idx0 + ACol_ptr[jy_new][jx_new].countTheta[p*pieceLength+t];
                        if(idx0 < pwMin)
                            pwMin = idx0;
                        if(pwMax < idx1)
                            pwMax = idx1;
                    }
                    piecewiseMin[i][p] = pwMin;
                    piecewiseMax[i][p] = pwMax;
                    piecewiseWidth[i][p] = (pwMax - pwMin);
                }
            }

            for(i=0; i<SVSize; i++)
            {
                totalSum[i]=0;
                //#pragma vector aligned
                for (p = 0; p < NViewSets; p++)
                    totalSum[i] += piecewiseWidth[i][p] * pieceLength;
            }

            unsigned char **AMatrixPadded= (unsigned char **) mget_spc(SVSize,sizeof(unsigned char *));
            unsigned char **AMatrixPaddedTranspose=(unsigned char **) mget_spc(SVSize,sizeof(unsigned char *));

            for(i=0;i<SVSize;i++) {
                AMatrixPadded[i] = (unsigned char *) mget_spc(totalSum[i],sizeof(unsigned char));
                AMatrixPaddedTranspose[i] = (unsigned char *) mget_spc(totalSum[i],sizeof(unsigned char));
            }

            for(i=0; i<SVSize; i++)
            {
                jy_new = jy_list[i];
                jx_new = jx_list[i];
                unsigned char * A_padded_pointer = &AMatrixPadded[i][0];
                unsigned char * newProjectionValueArrayPointer = &AVal_ptr[jy_new][jx_new].val[0];
                for (p=0; p < NViews; p++)
                {
                    int n_pad;
                    n_pad=(int)ACol_ptr[jy_new][jx_new].minIndex[p]-(int)piecewiseMin[i][p/pieceLength]-(int)bandMin[p];
                    #pragma vector aligned
                    for(t=0; t<n_pad; t++) {
                        *A_padded_pointer = 0;
                        A_padded_pointer++;
                    }
                    #pragma vector aligned
                    for(t=0; t<ACol_ptr[jy_new][jx_new].countTheta[p]; t++) {
                        *A_padded_pointer = *newProjectionValueArrayPointer;
                        A_padded_pointer++;
                        newProjectionValueArrayPointer++;
                    }
                    n_pad=(int)piecewiseMax[i][p/pieceLength]-(int)ACol_ptr[jy_new][jx_new].minIndex[p]-(int)ACol_ptr[jy_new][jx_new].countTheta[p]+(int)bandMin[p];
                    #pragma vector aligned
                    for(t=0; t<n_pad; t++) {
                        *A_padded_pointer = 0;
                        A_padded_pointer++;
                    }
                }
            }

            for(i=0; i<SVSize; i++)
            {
                unsigned char * A_padded_pointer = &AMatrixPadded[i][0];
                unsigned char * A_padd_Tranpose_pointer = &AMatrixPaddedTranspose[i][0];
                for (p=0; p < NViewSets; p++)
                {
                    for(q=0; q<piecewiseWidth[i][p]; q++) {
                        for(t=0; t<pieceLength; t++) {
                            A_padd_Tranpose_pointer[q*pieceLength+t] = A_padded_pointer[t*piecewiseWidth[i][p]+q];
                        }
                    }
                    A_padded_pointer += piecewiseWidth[i][p]*pieceLength;
                    A_padd_Tranpose_pointer += piecewiseWidth[i][p]*pieceLength;
                }
            }

            for(i=0;i<SVSize;i++)
            {
                jy_new = jy_list[i];
                jx_new = jx_list[i];
                int VoxelPosition = (jy_new-jy)*(2*SVLength+1)+(jx_new-jx);

                A_Padded_Map[jj][VoxelPosition].val = (unsigned char *)get_spc(totalSum[i], sizeof(unsigned char));
                A_Padded_Map[jj][VoxelPosition].pieceWiseMin = (channel_t *)get_spc(NViewSets,sizeof(channel_t));
                A_Padded_Map[jj][VoxelPosition].pieceWiseWidth = (channel_t *)get_spc(NViewSets,sizeof(channel_t));
                A_Padded_Map[jj][VoxelPosition].length = totalSum[i];
                memcpy(&A_Padded_Map[jj][VoxelPosition].val[0],&AMatrixPaddedTranspose[i][0],sizeof(unsigned char)*totalSum[i]);
                memcpy(&A_Padded_Map[jj][VoxelPosition].pieceWiseMin[0],&piecewiseMin[i][0],sizeof(channel_t)*NViewSets);
                memcpy(&A_Padded_Map[jj][VoxelPosition].pieceWiseWidth[0],&piecewiseWidth[i][0],sizeof(channel_t)*NViewSets);
            }

            for(i=0;i<SVSize;i++) {
                free((void *)AMatrixPadded[i]);
                free((void *)AMatrixPaddedTranspose[i]);
            }
            free((void *)AMatrixPadded);
            free((void *)AMatrixPaddedTranspose);
        } //omp for block

        multifree(piecewiseMin,2);
        multifree(piecewiseMax,2);
        multifree(piecewiseWidth,2);
        free((void *) bandMin);
        free((void *) bandMax);
        free((void *) bandWidth);
        free((void *) bandWidthPW);
        free((void *) jx_list);
        free((void *) jy_list);
        free((void *) totalSum);

    }   //omp parallel block

    free((void *) order);

}  /*** END A_piecewise() ***/



/* Compute Entire System Matrix */
/* The System matrix does not vary with slice for 3-D Parallel Geometry */
/* So, the method of compuatation is same as that of 2-D Parallel Geometry */

void A_comp(
    struct AValues_char **A_Padded_Map,
    float *Aval_max_ptr,
    struct SVParams svpar,
    struct SinoParams3DParallel *sinoparams,
    char *recon_mask,
    struct ImageParams3D *imgparams)
{
    int i,j,r;
    struct ACol **ACol_arr;
    struct AValues_char **AVal_arr;
    int NViews = sinoparams->NViews;
    int NChannels = sinoparams->NChannels;
    int Nx = imgparams->Nx;
    int Ny = imgparams->Ny;

    ACol_arr = (struct ACol **)multialloc(sizeof(struct ACol), 2, Ny, Nx);
    AVal_arr = (struct AValues_char **)multialloc(sizeof(struct AValues_char), 2, Ny, Nx);

    float **pix_prof = ComputePixelProfile3DParallel(sinoparams,imgparams);

    //struct timeval tm1,tm2;
    //unsigned long long tdiff;
    //gettimeofday(&tm1,NULL);

    #pragma omp parallel private(j,r)
    {
        struct ACol A_col_sgl;
        A_col_sgl.countTheta = (chanwidth_t *)get_spc(NViews,sizeof(chanwidth_t));
        A_col_sgl.minIndex = (channel_t *)get_spc(NViews,sizeof(channel_t));
        float *A_val_sgl = (float *)get_spc(NViews*NChannels, sizeof(float));

        #pragma omp for schedule(static)
        for (i=0; i<Ny; i++)
        for (j=0; j<Nx; j++)
        if(recon_mask[i*Nx+j])
        {
            A_comp_ij(i,j,sinoparams,imgparams,pix_prof,&A_col_sgl,A_val_sgl);
            ACol_arr[i][j].n_index = A_col_sgl.n_index;
            ACol_arr[i][j].countTheta = (chanwidth_t *) get_spc(NViews,sizeof(chanwidth_t));
            ACol_arr[i][j].minIndex = (channel_t *) get_spc(NViews,sizeof(channel_t));
            AVal_arr[i][j].val = (unsigned char *) get_spc(A_col_sgl.n_index, sizeof(unsigned char));

            float maxval = A_val_sgl[0];
            for (r = 0; r < A_col_sgl.n_index; r++) {
                if(A_val_sgl[r]>maxval)
                    maxval = A_val_sgl[r];
            }
            Aval_max_ptr[i*Nx+j] = maxval;

            for (r=0; r < A_col_sgl.n_index; r++)
                AVal_arr[i][j].val[r] = (unsigned char)((A_val_sgl[r])/maxval*255+0.5);

            for (r=0; r < NViews; r++) {
                ACol_arr[i][j].countTheta[r] = A_col_sgl.countTheta[r];
                ACol_arr[i][j].minIndex[r] = A_col_sgl.minIndex[r];
            }
        }
        else
        {
            ACol_arr[i][j].n_index = 0;
            Aval_max_ptr[i*Nx+j] = 0;
        }

        free((void *)A_val_sgl);
        free((void *)A_col_sgl.countTheta);
        free((void *)A_col_sgl.minIndex);
    }

    //gettimeofday(&tm2,NULL);
    //tdiff = 1000 * (tm2.tv_sec - tm1.tv_sec) + (tm2.tv_usec - tm1.tv_usec) / 1000;
    //fprintf(stdout,"matrix time 1 = %llu ms\n",tdiff);
    //gettimeofday(&tm1,NULL);

    A_piecewise(ACol_arr,AVal_arr,A_Padded_Map,svpar,sinoparams,imgparams);

    //gettimeofday(&tm2,NULL);
    //tdiff = 1000 * (tm2.tv_sec - tm1.tv_sec) + (tm2.tv_usec - tm1.tv_usec) / 1000;
    //fprintf(stdout,"matrix time 2 = %llu ms\n",tdiff);

    for (i=0; i<Ny; i++)
    for (j=0; j<Nx; j++)
    if(recon_mask[i*Nx+j]) {
        free((void *)ACol_arr[i][j].countTheta);
        free((void *)ACol_arr[i][j].minIndex);
        free((void *)AVal_arr[i][j].val);
    }
    multifree(ACol_arr,2);
    multifree(AVal_arr,2);

    free_img((void **)pix_prof);

}


void readAmatrix(
    char *fname,
    struct AValues_char **A_Padded_Map,
    float *Aval_max_ptr,
    struct ImageParams3D *imgparams,
    struct SinoParams3DParallel *sinoparams,
    struct SVParams svpar)
{
    FILE *fp;
    int i,j;
    int M_nonzero;

    int Nxy = imgparams->Nx * imgparams->Ny;
    int NViews = sinoparams->NViews;
    int NViewSets = sinoparams->NViews/svpar.pieceLength;

    if ((fp = fopen(fname, "rb")) == NULL) {
        fprintf(stderr, "ERROR in readAmatrix: can't open file %s.\n", fname);
        exit(-1);
    }

    for (i=0; i<svpar.Nsv ; i++)
    {
        if(fread(svpar.bandMinMap[i].bandMin,sizeof(channel_t),NViews,fp) < (size_t)NViews) {
            fprintf(stderr, "ERROR in readAmatrix: %s terminated early.\n", fname);
            exit(-1);
        }
        if(fread(svpar.bandMaxMap[i].bandMax,sizeof(channel_t),NViews,fp) < (size_t)NViews) {
            fprintf(stderr, "ERROR in readAmatrix: %s terminated early.\n", fname);
            exit(-1);
        }

        for (j=0; j< (svpar.SVLength*2+1)*(svpar.SVLength*2+1); j++)
        {
            if(fread(&M_nonzero, sizeof(int), 1, fp) < 1) {
                fprintf(stderr, "ERROR in readAmatrix: %s terminated early.\n", fname);
                exit(-1);
            }
            A_Padded_Map[i][j].length = M_nonzero;
            if(M_nonzero > 0)
            {
                A_Padded_Map[i][j].val = (unsigned char *)get_spc(M_nonzero, sizeof(unsigned char));
                A_Padded_Map[i][j].pieceWiseWidth = (channel_t *)get_spc(NViewSets,sizeof(channel_t));
                A_Padded_Map[i][j].pieceWiseMin = (channel_t *)get_spc(NViewSets,sizeof(channel_t));

                if(fread(A_Padded_Map[i][j].val, sizeof(unsigned char), M_nonzero, fp) < (size_t)M_nonzero) {
                    fprintf(stderr, "ERROR in readAmatrix: %s terminated early.\n", fname);
                    exit(-1);
                }
                if(fread(A_Padded_Map[i][j].pieceWiseMin,sizeof(channel_t),NViewSets,fp) < (size_t)NViewSets) {
                    fprintf(stderr, "ERROR in readAmatrix: %s terminated early.\n", fname);
                    exit(-1);
                }
                if(fread(A_Padded_Map[i][j].pieceWiseWidth,sizeof(channel_t),NViewSets,fp) < (size_t)NViewSets) {
                    fprintf(stderr, "ERROR in readAmatrix: %s terminated early.\n", fname);
                    exit(-1);
                }
            }
        }
    }

    if(fread(&Aval_max_ptr[0],sizeof(float),Nxy,fp) < (size_t)Nxy) {
        fprintf(stderr, "ERROR in readAmatrix: %s terminated early.\n", fname);
        exit(-1);
    }

    fclose(fp);
}


void writeAmatrix(
    char *fname,
    struct AValues_char **A_Padded_Map,
    float *Aval_max_ptr,
    struct ImageParams3D *imgparams,
    struct SinoParams3DParallel *sinoparams,
    struct SVParams svpar)
{
    FILE *fp;
    int i,j;
    int M_nonzero;
    int NViewSets = sinoparams->NViews/svpar.pieceLength;

    if ((fp = fopen(fname, "wb")) == NULL) {
        fprintf(stderr, "ERROR in writeAmatrix: can't open file %s.\n", fname);
        exit(-1);
    }

    for (i=0; i<svpar.Nsv; i++)
    {
        fwrite(svpar.bandMinMap[i].bandMin,sizeof(channel_t),sinoparams->NViews,fp);
        fwrite(svpar.bandMaxMap[i].bandMax,sizeof(channel_t),sinoparams->NViews,fp);
        for (j=0; j< (svpar.SVLength*2+1)*(svpar.SVLength*2+1); j++)
        {
            M_nonzero = A_Padded_Map[i][j].length;
            fwrite(&M_nonzero, sizeof(int), 1, fp);
            if(M_nonzero > 0) {
                fwrite(A_Padded_Map[i][j].val, sizeof(unsigned char), M_nonzero, fp);
                fwrite(A_Padded_Map[i][j].pieceWiseMin,sizeof(channel_t),NViewSets,fp);
                fwrite(A_Padded_Map[i][j].pieceWiseWidth,sizeof(channel_t),NViewSets,fp);
            }
        }
    }
    fwrite(&Aval_max_ptr[0],sizeof(float),imgparams->Nx*imgparams->Ny,fp);
    fclose(fp);
}


void AmatrixComputeToFile(
    struct ImageParams3D imgparams,
    struct SinoParams3DParallel sinoparams,
    char *Amatrix_fname,
    char verboseLevel)
{
    struct SVParams svpar;
    struct AValues_char **A_Padded_Map;
    float *Aval_max_ptr;
    char *ImageReconMask;   /* Image reconstruction mask (determined by ROI) */
    int i,j;
    #ifndef MSVC    /* not included in MS Visual C++ */
    struct timeval tm1,tm2;
    unsigned long long tdiff;
    #endif

    if(verboseLevel) {
        fprintf(stdout,"Computing system matrix...\n");
        #ifndef MSVC    /* not included in MS Visual C++ */
        gettimeofday(&tm1,NULL);
        #endif
    }

    initSVParams(&svpar,imgparams,sinoparams);  /* Initialize/allocate SV parameters */

    int Nx = imgparams.Nx;
    int Ny = imgparams.Ny;
    int SVLength = svpar.SVLength;
    int Nsv = svpar.Nsv;

    /* Allocate and generate recon mask based on ROIRadius */
    ImageReconMask = GenImageReconMask(&imgparams);

    A_Padded_Map = (struct AValues_char **)multialloc(sizeof(struct AValues_char),2,Nsv,(2*SVLength+1)*(2*SVLength+1));
    Aval_max_ptr = (float *) get_spc(Nx*Ny,sizeof(float));

    A_comp(A_Padded_Map,Aval_max_ptr,svpar,&sinoparams,ImageReconMask,&imgparams);

    if(verboseLevel>1) {
        #ifndef MSVC    /* not included in MS Visual C++ */
        gettimeofday(&tm2,NULL);
        tdiff = 1000 * (tm2.tv_sec - tm1.tv_sec) + (tm2.tv_usec - tm1.tv_usec) / 1000;
        fprintf(stdout,"\tmatrix time = %llu ms\n",tdiff);
        #endif
        fprintf(stdout,"Writing system matrix %s\n",Amatrix_fname);
    }
    else if(verboseLevel)
        fprintf(stdout,"Writing system matrix...\n");

    writeAmatrix(Amatrix_fname,A_Padded_Map,Aval_max_ptr,&imgparams,&sinoparams,svpar);

    /* Free memory */
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

}

