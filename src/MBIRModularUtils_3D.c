
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <libgen.h>  /* for dirname */
#include "allocate.h"
#include "MBIRModularUtils_2D.h"
#include "MBIRModularUtils_3D.h"


/*************************************************/
/*  Utilities for reading/printing .param files  */
/*************************************************/

/* Print Sinogram parameters */
void printSinoParams3DParallel(struct SinoParams3DParallel *sinoparams)
{
    fprintf(stdout, "SINOGRAM PARAMETERS:\n");
    fprintf(stdout, " - Number of sinogram views per slice    = %d\n", sinoparams->NViews);
    fprintf(stdout, " - Number of detector channels per slice = %d\n", sinoparams->NChannels);
    fprintf(stdout, " - Number of slices                      = %d\n", sinoparams->NSlices);
    fprintf(stdout, " - Spacing between Detector Channels     = %.7f (mm) \n", sinoparams->DeltaChannel);
    fprintf(stdout, " - Center of rotation offset             = %.7f (channels)\n", sinoparams->CenterOffset);
    fprintf(stdout, " - Spacing between slices                = %.7f (mm)\n", sinoparams->DeltaSlice);
    fprintf(stdout, " - First Slice Index                     = %d \n", sinoparams->FirstSliceNumber);
}

/* Utility for reading 3D parallel beam sinogram parameters */
/* Returns 0 if no error occurs */
int ReadSinoParams3DParallel(
	char *basename,					/* Input: Reads sinogram parameters from <fname>.sinoparams */
	struct SinoParams3DParallel *sinoparams)	/* Output: Reads sinogram parameters into data structure */
{
	FILE *fp;
	char fname[200];
	char tag[200], fieldname[200], fieldval_s[200], *ptr;
	char AngleListFname[200]=" ";
	int i, Nlines;
	char Geometry_flag=0;

	/* set defaults, also used for error checking below */
	sinoparams->NViews=0;		/* Number of view angles */
	sinoparams->NChannels=-1;	/* Number of channels in detector */
	sinoparams->DeltaChannel=0.0;	/* Detector spacing (mm) */
	sinoparams->NSlices=0;		/* Number of slices stored in Sino array */
	sinoparams->FirstSliceNumber=-1;	/* slice index coresponding to first slice in volume */
	sinoparams->DeltaSlice=0.0;	/* Spacing along slice direction (mm) */
	sinoparams->CenterOffset=0.0;	/* Offset of center-of-rotation ... */

	strcpy(fname,basename);
	strcat(fname,".sinoparams"); /* append file extension */
	if((fp=fopen(fname,"r")) == NULL) {
		fprintf(stderr,"ERROR in ReadSinoParams3DParallel: can't open file %s\n",fname);
		exit(-1);
	}

	Nlines=0;
	while(fgets(tag, 200, fp)!=NULL && Nlines<100)
		Nlines++;
	rewind(fp);
	//printf("Nlines=%d\n",Nlines);

	/* parse each line assuming a single ":" delimiter */
	for(i=0; i<Nlines; i++)
	{
		strcpy(fieldname," ");
		strcpy(fieldval_s," ");
		fgets(tag, 200, fp);
		ptr=strtok(tag,":\n\r");	// including the newline will keep it out of the last token
		if(ptr!=NULL) {
			//strcpy(fieldname,ptr);
			sscanf(ptr,"%s",fieldname);	// this won't include leading/trailing spaces, strcpy will
			ptr=strtok(NULL,":\n\r");
		}
		if(ptr!=NULL) sscanf(ptr,"%s",fieldval_s);
		//printf("|%s|%s|\n",fieldname,fieldval_s);

		if(strcmp(fieldname,"Geometry")==0)
		{
			Geometry_flag=1;
			if(strcmp(fieldval_s,"3DPARALLEL")!=0) {
				fprintf(stderr,"Error in %s: Geometry value \"%s\" unrecognized\n",fname,fieldval_s);
				exit(-1);
			}
		}
		else if(strcmp(fieldname,"NChannels")==0)
		{
			sscanf(fieldval_s,"%d",&(sinoparams->NChannels));
		}
		else if(strcmp(fieldname,"NViews")==0)
		{
			sscanf(fieldval_s,"%d",&(sinoparams->NViews));
		}
		else if(strcmp(fieldname,"NSlices")==0)
		{
			sscanf(fieldval_s,"%d",&(sinoparams->NSlices));
		}
		else if(strcmp(fieldname,"DeltaChannel")==0)
		{
			sscanf(fieldval_s,"%f",&(sinoparams->DeltaChannel));
		}
		else if(strcmp(fieldname,"CenterOffset")==0)
		{
			sscanf(fieldval_s,"%f",&(sinoparams->CenterOffset));
		}
		else if(strcmp(fieldname,"DeltaSlice")==0)
		{
			sscanf(fieldval_s,"%f",&(sinoparams->DeltaSlice));
		}
		else if(strcmp(fieldname,"FirstSliceNumber")==0)
		{
			sscanf(fieldval_s,"%d",&(sinoparams->FirstSliceNumber));
		}
		else if(strcmp(fieldname,"ViewAngleList")==0)
			sscanf(fieldval_s,"%s",AngleListFname);  // this won't include leading/trailing spaces, strcpy will
			//strcpy(AngleListFname,fieldval_s);
		else
			fprintf(stderr,"Warning: unrecognized field \"%s\" in %s, line %d\n",fieldname,fname,i+1);

	}  // done parsing

	fclose(fp);

	//printSinoParams3DParallel(sinoparams);

	/* do some error checking */
	if(Geometry_flag==0) {
		fprintf(stderr,"Error in %s: \"Geometry\" field unspecified\n",fname);
		exit(-1);
	}
	if(sinoparams->NViews<=0 || sinoparams->NChannels<=0 || sinoparams->NSlices<=0) {
		printSinoParams3DParallel(sinoparams);
		fprintf(stderr,"Error in %s: NViews, NChannels, NSlices must all be positive\n",fname);
		exit(-1);
	}
	if(sinoparams->DeltaChannel<=0 && sinoparams->NChannels>1) {
		printSinoParams3DParallel(sinoparams);
		fprintf(stderr,"Error in %s: DeltaChannel needs to be positive (mm)\n",fname);
		exit(-1);
	}
	if(sinoparams->DeltaSlice<=0 && sinoparams->NSlices>1) {
		printSinoParams3DParallel(sinoparams);
		fprintf(stderr,"Error in %s: DeltaSlice needs to be positive (mm)\n",fname);
		exit(-1);
	}
	if(sinoparams->FirstSliceNumber < 0) {
		printSinoParams3DParallel(sinoparams);
		fprintf(stderr,"Error in %s: FirstSliceNumber should be non-negative\n",fname);
		exit(-1);
	}

	/* form full pathname of ViewAngleList file; path relative to sinoparams directory */
	strcpy(fieldval_s,AngleListFname);  // tmp copy
	strcpy(tag,fname);
	ptr=dirname(tag);
	sprintf(AngleListFname,"%s/%s",ptr,fieldval_s);
	//printf("Views filename \"%s\"\n",AngleListFname);

	/* Read the view angles file */
	if((fp=fopen(AngleListFname,"r")) == NULL) {
		fprintf(stderr,"ERROR in ReadSinoParams3DParallel: can't open ViewAngle file %s\n",AngleListFname);
		exit(-1);
	}

	sinoparams->ViewAngles = (float *)get_spc(sinoparams->NViews, sizeof(float));

	for(i=0;i<sinoparams->NViews;i++)
	if(fscanf(fp,"%f\n",&(sinoparams->ViewAngles[i])) == 0) {
		fprintf(stderr, "ERROR in ReadSinoParams3DParallel: View angles file %s terminated early\n", AngleListFname);
		exit(-1);
	}

	fclose(fp);

	return(0);
}


/* Print Image parameters */
void printImageParams3D(struct ImageParams3D *imgparams)
{
    fprintf(stdout, "IMAGE PARAMETERS:\n");
    fprintf(stdout, " - Number of Pixels within a single slice in X direction = %d\n", imgparams->Nx);
    fprintf(stdout, " - Number of Pixels within a single slice in Y direction = %d\n", imgparams->Ny);
    fprintf(stdout, " - Number of Slices to reconstruct                       = %d \n", imgparams->Nz);
    fprintf(stdout, " - Pixel width  in XY plane                              = %.7f (mm)\n", imgparams->Deltaxy);
    fprintf(stdout, " - Spacing between slices                                = %.7f (mm)\n", imgparams->DeltaZ);
    fprintf(stdout, " - First Slice Index                                     = %d \n", imgparams->FirstSliceNumber);
    fprintf(stdout, " - ROIRadius                                             = %.7f (mm)\n", imgparams->ROIRadius);
}

/* Utility for reading 2D Image parameters */
/* Returns 0 if no error occurs */
int ReadImageParams3D(
	char *basename,				/* Input: Reads image type parameter from <fname>.imgparams */
	struct ImageParams3D *imgparams)	/* Output: Reads image parameters into data structure */
{
	FILE *fp;
	char fname[200];
	char tag[200], fieldname[200], fieldval_s[200], *ptr;
	int i, Nlines;

	/* set defaults, also used for error checking below */
	imgparams->Nx=0;
	imgparams->Ny=0;
	imgparams->Nz=0;
	imgparams->FirstSliceNumber=-1;
	imgparams->Deltaxy=0.0;
	imgparams->DeltaZ=0;
	imgparams->ROIRadius=0.0;

	strcpy(fname,basename);
	strcat(fname,".imgparams");
	if((fp=fopen(fname,"r")) == NULL) {
		fprintf(stderr,"ERROR in ReadImageParams3D: can't open file %s\n",fname);
		exit(-1);
	}

	Nlines=0;
	while(fgets(tag, 200, fp)!=NULL && Nlines<100)
		Nlines++;
	rewind(fp);
	//printf("Nlines=%d\n",Nlines);

	/* parse each line assuming a single ":" delimiter */
	for(i=0; i<Nlines; i++)
	{
		strcpy(fieldname," ");
		strcpy(fieldval_s," ");
		fgets(tag, 200, fp);
		ptr=strtok(tag,":\n\r");	// including the newline will keep it out of the last token
		if(ptr!=NULL) {
			//strcpy(fieldname,ptr);
			sscanf(ptr,"%s",fieldname);	// this won't include leading/trailing spaces, strcpy will
			ptr=strtok(NULL,":\n\r");
		}
		if(ptr!=NULL) sscanf(ptr,"%s",fieldval_s);
		//printf("|%s|%s|\n",fieldname,fieldval_s);

		if(strcmp(fieldname,"Nx")==0)
		{
			sscanf(fieldval_s,"%d",&(imgparams->Nx));
		}
		else if(strcmp(fieldname,"Ny")==0)
		{
			sscanf(fieldval_s,"%d",&(imgparams->Ny));
		}
		else if(strcmp(fieldname,"Nz")==0)
		{
			sscanf(fieldval_s,"%d",&(imgparams->Nz));
		}
		else if(strcmp(fieldname,"FirstSliceNumber")==0)
		{
			sscanf(fieldval_s,"%d",&(imgparams->FirstSliceNumber));
		}
		else if(strcmp(fieldname,"Deltaxy")==0)
		{
			sscanf(fieldval_s,"%f",&(imgparams->Deltaxy));
		}
		else if(strcmp(fieldname,"DeltaZ")==0)
		{
			sscanf(fieldval_s,"%f",&(imgparams->DeltaZ));
		}
		else if(strcmp(fieldname,"ROIRadius")==0)
		{
			sscanf(fieldval_s,"%f",&(imgparams->ROIRadius));
		}
		else
			fprintf(stderr,"Warning: unrecognized field \"%s\" in %s, line %d\n",fieldname,fname,i+1);

	}  // done parsing

	fclose(fp);

	//printImageParams3D(imgparams);

	/* do some error checking */
	if(imgparams->Nx<=0 || imgparams->Ny<=0 || imgparams->Nz<=0) {
		printImageParams3D(imgparams);
		fprintf(stderr,"Error in %s: Nx, Ny, Nz must all be positive\n",fname);
		exit(-1);
	}
	if(imgparams->Deltaxy<=0) {
		printImageParams3D(imgparams);
		fprintf(stderr,"Error in %s: Deltaxy needs to be positive (mm)\n",fname);
		exit(-1);
	}
	if(imgparams->DeltaZ<=0 && imgparams->Nz>1) {
		printImageParams3D(imgparams);
		fprintf(stderr,"Error in %s: DeltaZ needs to be positive (mm)\n",fname);
		exit(-1);
	}
	if(imgparams->FirstSliceNumber < 0) {
		printImageParams3D(imgparams);
		fprintf(stderr,"Error in %s: FirstSliceNumber should be non-negative\n",fname);
		exit(-1);
	}
	if(imgparams->ROIRadius<=0) {
		imgparams->ROIRadius = imgparams->Nx * imgparams->Deltaxy;
		fprintf(stderr,"Warning in %s: ROIRadius needs to be positive. Defaulting to %.4f (mm)\n",fname,imgparams->ROIRadius);
		printImageParams3D(imgparams);
	}

	return(0);
}


/* Print QGGMRF reconstruction parameters */
void printReconParamsQGGMRF3D(struct ReconParamsQGGMRF3D *reconparams)
{
    fprintf(stdout, "PRIOR PARAMETERS:\n");
    fprintf(stdout, " - Q-GGMRF Prior Parameter, q                            = %f\n", reconparams->p);
    fprintf(stdout, " - Q-GGMRF Prior Parameter, p                            = %f\n", reconparams->q);
    fprintf(stdout, " - Q-GGMRF Prior Parameter, T                            = %f\n", reconparams->T);
    fprintf(stdout, " - Prior Regularization parameter, SigmaX                = %.7f (mm^-1)\n", reconparams->SigmaX);
    fprintf(stdout, " - Scaling for weight matrix, SigmaY (W <- W/SigmaY^2)   = %.7f (mm^-1)\n", reconparams->SigmaY);
    fprintf(stdout, " - Prior weight for nearest neighbors within slice       = %.7f\n", reconparams->b_nearest);
    fprintf(stdout, " - Prior weight for diagonal neighbors within slice      = %.7f\n", reconparams->b_diag);
    fprintf(stdout, " - Prior weight for nearest neighbors in adjacent slices = %.7f\n", reconparams->b_interslice);
    fprintf(stdout, " - Inital image value                                    = %-10f (mm-1)\n", reconparams->InitImageValue);
    fprintf(stdout, " - Stop threshold for convergence                        = %.7f %%\n", reconparams->StopThreshold);
    fprintf(stdout, " - Maximum number of ICD iterations                      = %d\n", reconparams->MaxIterations);
    fprintf(stdout, " - Positivity constraint flag                            = %d\n", reconparams->Positivity);
}

/* Utility for reading QGGMRF reconstruction parameters */
/* Returns 0 if no error occurs */
int ReadReconParamsQGGMRF3D(
	char *basename,				/* base file name <fname>.reconparams */
	struct ReconParamsQGGMRF3D *reconparams) /* recon parameters data structure */
{
	FILE *fp;
	char fname[200];
	char tag[200], fieldname[200], fieldval_s[200], *ptr;
	double fieldval_f;
	int fieldval_d;
	int i, Nlines;
	char Prior_flag=0;

	/* set defaults, also used for error checking below */
	reconparams->InitImageValue=MUWATER;
	reconparams->p=1.2;
	reconparams->q=2.0;
	reconparams->T=0.1;
	reconparams->SigmaX=0.02;
	reconparams->SigmaY=1.0;
	reconparams->b_nearest=1.0;
	reconparams->b_diag=0.707;
	reconparams->b_interslice=1.0;
	reconparams->StopThreshold=1.0;
	reconparams->MaxIterations=20;
	reconparams->Positivity=1;

	strcpy(fname,basename);
	strcat(fname,".reconparams");
	if((fp=fopen(fname,"r")) == NULL) {
		fprintf(stderr,"ERROR in ReadReconParamsQGGMRF3D: can't open file %s\n",fname);
		exit(-1);
	}

	Nlines=0;
	while(fgets(tag, 200, fp)!=NULL && Nlines<100)
		Nlines++;
	rewind(fp);
	//printf("Nlines=%d\n",Nlines);

	/* parse each line assuming a single ":" delimiter */
	for(i=0; i<Nlines; i++)
	{
		strcpy(fieldname," ");
		strcpy(fieldval_s," ");
		fgets(tag, 200, fp);
		ptr=strtok(tag,":\n\r");	// including the newline will keep it out of the last token
		if(ptr!=NULL) {
			//strcpy(fieldname,ptr);
			sscanf(ptr,"%s",fieldname);	// this won't include leading/trailing spaces, strcpy will
			ptr=strtok(NULL,":\n\r");
		}
		if(ptr!=NULL) sscanf(ptr,"%s",fieldval_s);
		//printf("|%s|%s|\n",fieldname,fieldval_s);

		if(strcmp(fieldname,"PriorModel")==0)
		{
			Prior_flag=1;
			if(strcmp(fieldval_s,"QGGMRF")!=0) {
				fprintf(stderr,"Error in %s: PriorModel value \"%s\" unrecognized\n",fname,fieldval_s);
				exit(-1);
			}
		}
		else if(strcmp(fieldname,"InitImageValue")==0)
		{
			//sscanf(fieldval_s,"%lf",&(reconparams->InitImageValue));
			//Changed above to the following to retain default value if input doesn't make sense
			sscanf(fieldval_s,"%lf",&(fieldval_f));
			if(fieldval_f <= 0)
				fprintf(stderr,"Warning in %s: InitImageValue should be positive. Reverting to default.\n",fname);
			else
				reconparams->InitImageValue = fieldval_f;
		}
		else if(strcmp(fieldname,"p")==0)
		{
			sscanf(fieldval_s,"%lf",&(fieldval_f));
			if(fieldval_f < 1 || fieldval_f > 2)
				fprintf(stderr,"Warning in %s: p parameter should be in range [1,2]. Reverting to default.\n",fname);
			else
				reconparams->p = fieldval_f;
		}
		else if(strcmp(fieldname,"q")==0)
		{
			sscanf(fieldval_s,"%lf",&(fieldval_f));
			if(fieldval_f < 1 || fieldval_f > 2)
				fprintf(stderr,"Warning in %s: q parameter should be in range [1,2]. Reverting to default.\n",fname);
			else
				reconparams->q = fieldval_f;
		}
		else if(strcmp(fieldname,"T")==0)
		{
			sscanf(fieldval_s,"%lf",&(fieldval_f));
			if(fieldval_f <= 0)
				fprintf(stderr,"Warning in %s: T parameter should be positive. Reverting to default.\n",fname);
			else
				reconparams->T = fieldval_f;
		}
		else if(strcmp(fieldname,"SigmaX")==0)
		{
			sscanf(fieldval_s,"%lf",&(fieldval_f));
			if(fieldval_f <= 0)
				fprintf(stderr,"Warning in %s: SigmaX parameter should be positive. Reverting to default.\n",fname);
			else
				reconparams->SigmaX = fieldval_f;
		}
		else if(strcmp(fieldname,"SigmaY")==0)
		{
			sscanf(fieldval_s,"%lf",&(fieldval_f));
			if(fieldval_f <= 0)
				fprintf(stderr,"Warning in %s: SigmaY parameter should be positive. Reverting to default.\n",fname);
			else
				reconparams->SigmaY = fieldval_f;
		}
		else if(strcmp(fieldname,"b_nearest")==0)
		{
			sscanf(fieldval_s,"%lf",&(fieldval_f));
			if(fieldval_f <= 0)
				fprintf(stderr,"Warning in %s: b_nearest parameter should be positive. Reverting to default.\n",fname);
			else
				reconparams->b_nearest = fieldval_f;
		}
		else if(strcmp(fieldname,"b_diag")==0)
		{
			sscanf(fieldval_s,"%lf",&(fieldval_f));
			if(fieldval_f < 0)
				fprintf(stderr,"Warning in %s: b_diag parameter should be non-negative. Reverting to default.\n",fname);
			else
				reconparams->b_diag = fieldval_f;
		}
		else if(strcmp(fieldname,"b_interslice")==0)
		{
			sscanf(fieldval_s,"%lf",&(fieldval_f));
			if(fieldval_f < 0)
				fprintf(stderr,"Warning in %s: b_interslice parameter should be non-negative. Reverting to default.\n",fname);
			else
				reconparams->b_interslice = fieldval_f;
		}
		else if(strcmp(fieldname,"StopThreshold")==0)
		{
			sscanf(fieldval_s,"%lf",&(fieldval_f));
			if(fieldval_f < 0)
				fprintf(stderr,"Warning in %s: StopThreshold should be non-negative. Reverting to default.\n",fname);
			else
				reconparams->StopThreshold = fieldval_f;
		}
		else if(strcmp(fieldname,"MaxIterations")==0)
		{
			sscanf(fieldval_s,"%d",&(fieldval_d));
			if(fieldval_d < 1)
				fprintf(stderr,"Warning in %s: MaxIterations should be at least 1. Reverting to default.\n",fname);
			else
				reconparams->MaxIterations = fieldval_d;
		}
		else if(strcmp(fieldname,"Positivity")==0)
		{
			sscanf(fieldval_s,"%d",&(fieldval_d));
			if( strcmp(fieldval_s,"0") && strcmp(fieldval_s,"1") )
				fprintf(stderr,"Warning in %s: \"Positivity\" parameter options are 0/1. Reverting to default.\n",fname);
			else
				reconparams->Positivity = fieldval_d;
		}
		else
			fprintf(stderr,"Warning: unrecognized field \"%s\" in %s, line %d\n",fieldname,fname,i+1);

	}  // done parsing

	fclose(fp);

	//printReconParamsQGGMRF3D(reconparams);

	/* do some error checking */
	if(Prior_flag==0) {
		fprintf(stderr,"Error in %s: \"PriorModel\" field unspecified\n",fname);
		exit(-1);
	}
	if(reconparams->p > reconparams->q) {
		printReconParamsQGGMRF3D(reconparams);
		fprintf(stderr,"Error in %s: Need (p <= q) for convexity. (p<q for strict convexity)\n",fname);
		exit(-1);
	}

	return(0);
}


/**********************************************/
/*  Utilities for reading/writing 3D sinogram */
/**********************************************/

/* Utility for reading 3D parallel beam sinogram data */
/* Warning: Memory must be allocated before use */
/* Returns 0 if no error occurs */
int ReadSinoData3DParallel(
    char *basename,   /* Input: Reads sinogram data from <basename>_slice<index>.2Dsinodata for given index range */
    struct Sino3DParallel *sinogram)  /* Input/Output: Uses sinogram parameters and reads sinogram data into data structure */
{
    FILE *fp;
    char fname[200];
    int i,NSlices,FirstSliceNumber,M;
    
    NSlices = sinogram->sinoparams.NSlices;
    FirstSliceNumber = sinogram->sinoparams.FirstSliceNumber;
    M = sinogram->sinoparams.NViews * sinogram->sinoparams.NChannels;
    
    //printf("Reading 3-D Projection Data ... \n");
    for(i=0;i<NSlices;i++)
    {
        /* slice index currently formed from fixed number of digits */
        sprintf(fname,"%s_slice%.*d.2Dsinodata",basename,MBIR_MODULAR_MAX_NUMBER_OF_SLICE_DIGITS,i+FirstSliceNumber);
	//printf("filename: |%s|\n",fname);
        
        if((fp = fopen(fname,"r")) == NULL) {
            fprintf(stderr, "ERROR in ReadSinoData3DParallel: can't open file %s\n",fname);
            exit(-1);
        }
        if(fread(sinogram->sino[i],sizeof(float),M,fp)!=M) {
            fprintf(stderr, "ERROR in ReadSinoData3DParallel: file %s terminated early\n",fname);
            exit(-1);
        }
        fclose(fp);
    }
    return 0;
}


/* Utility for reading weights for 3D sinogram projections data */
/* Warning: Memory must be allocated before use */
/* Returns 0 if no error occurs */
int ReadWeights3D(
    char *basename,       /* Input: Reads sinogram data from <basename>_slice<index>.2Dweightdata for given index range */
    struct Sino3DParallel *sinogram) /* Input/Output: Uses sinogram parameters and reads sinogram data into data structure */
{
    FILE *fp;
    char fname[200];
    int i,NSlices,FirstSliceNumber,M;

    NSlices = sinogram->sinoparams.NSlices;
    FirstSliceNumber = sinogram->sinoparams.FirstSliceNumber;
    M = sinogram->sinoparams.NViews * sinogram->sinoparams.NChannels;

    //printf("Reading 3-D Sinogram Weights ... \n");
    for(i=0;i<NSlices;i++)
    {
        /* slice index currently formed from fixed number of digits */
        sprintf(fname,"%s_slice%.*d.2Dweightdata",basename,MBIR_MODULAR_MAX_NUMBER_OF_SLICE_DIGITS,i+FirstSliceNumber);
	//printf("filename: |%s|\n",fname);
        
        if((fp = fopen(fname,"r")) == NULL) {
            fprintf(stderr,"ERROR in ReadWeights3D: can't open file %s\n",fname);
            exit(-1);
        }
        if(fread(sinogram->weight[i],sizeof(float),M,fp)!=M) {
            fprintf(stderr,"ERROR in ReadWeights3D: file %s terminated early\n",fname);
            exit(-1);
        }
        fclose(fp);
    }
    return 0;
}


/* Utility for writing out 3D parallel beam sinogram parameters and data */
/* Returns 0 if no error occurs */
int WriteSino3DParallel(
    char *basename,	/* Input: Writes sinogram data to <basename>_slice<n>.2Dsinodata for given slice range */
    struct Sino3DParallel *sinogram)  /* Input: Sinogran parameters and data */
{
    FILE *fp;
    char fname[200];
    int i,NSlices,FirstSliceNumber,M;

    NSlices = sinogram->sinoparams.NSlices;
    FirstSliceNumber = sinogram->sinoparams.FirstSliceNumber;
    M = sinogram->sinoparams.NViews * sinogram->sinoparams.NChannels;

    //printf("Writing 3-D Sinogram Projection Data ... \n");
    for(i=0;i<NSlices;i++)
    {
        /* slice index currently formed from fixed number of digits */
        sprintf(fname,"%s_slice%.*d.2Dsinodata",basename,MBIR_MODULAR_MAX_NUMBER_OF_SLICE_DIGITS,i+FirstSliceNumber);
	//printf("filename: |%s|\n",fname);
        
        if((fp = fopen(fname,"w")) == NULL) {
            fprintf(stderr,"ERROR in WriteSino3DParallel: can't open file %s\n",fname);
            exit(-1);
        }
        if(fwrite(sinogram->sino[i],sizeof(float),M,fp)!=M) {
            fprintf(stderr,"ERROR in WriteSino3DParallel: can't write to file %s\n",fname);
            exit(-1);
        }
        fclose(fp);
    }
    return 0;
}


/* Utility for writing out weights for 3D parallel beam sinogram data */
/* Returns 0 if no error occurs */
int WriteWeights3D(
    char *basename,	/* Input: Writes sinogram weights to <basename>_slice<n>.2Dweightdata for given slice range */
    struct Sino3DParallel *sinogram) /* Input: Sinogram parameters and data */
{
    FILE *fp;
    char fname[200];
    int i,NSlices,FirstSliceNumber,M;

    NSlices = sinogram->sinoparams.NSlices;
    FirstSliceNumber = sinogram->sinoparams.FirstSliceNumber;
    M = sinogram->sinoparams.NViews * sinogram->sinoparams.NChannels;

    //printf("Writing 3-D Sinogram Weights ... \n");
    for(i=0;i<NSlices;i++)
    {
        /* slice index currently formed from fixed number of digits */
        sprintf(fname,"%s_slice%.*d.2Dweightdata",basename,MBIR_MODULAR_MAX_NUMBER_OF_SLICE_DIGITS,i+FirstSliceNumber);
	//printf("filename: |%s|\n",fname);
        
        if((fp = fopen(fname,"w")) == NULL) {
            fprintf(stderr,"ERROR in WriteWeights3D: can't open file %s\n",fname);
            exit(-1);
        }
        if(fwrite(sinogram->weight[i],sizeof(float),M,fp)!=M) {
            fprintf(stderr,"ERROR in WriteWeights3D: can't write to file %s\n",fname);
            exit(-1);
        }
        fclose(fp);
    }
    return 0;
}


/* Utility for allocating memory for Sino */
/* Returns 0 if no error occurs */
int AllocateSinoData3DParallel(struct Sino3DParallel *sinogram)  /* Input: Sinogram parameters data structure */
{
    //printf("Allocating Sinogram Memory ... \n");
    sinogram->sino   = (float **)multialloc(sizeof(float), 2, sinogram->sinoparams.NSlices,sinogram->sinoparams.NViews * sinogram->sinoparams.NChannels);
    sinogram->weight = (float **)multialloc(sizeof(float), 2, sinogram->sinoparams.NSlices,sinogram->sinoparams.NViews * sinogram->sinoparams.NChannels);
    return 0;
}


/* Utility for freeing memory allocated for ViewAngles and Sino */
/* Returns 0 if no error occurs */
int FreeSinoData3DParallel(struct Sino3DParallel *sinogram)  /* Input: Sinogram parameters data structure */
{
    multifree(sinogram->sino,2);
    multifree(sinogram->weight,2);
    return 0;
}


/*******************************************/
/* Utilities for reading/writing 3D images */
/*******************************************/

/* Utility for reading 3D image data */
/* Warning: Memory must be allocated before use */
/* Returns 0 if no error occurs */
int ReadImage3D(
    char *basename,	/* Input: Reads 2D image data from <basename>_slice<n>.2Dimgdata for given slice range */
    struct Image3D *Image)	/* Input/Output: Uses image parameters (dimensions) and reads images into structure */
{
    FILE *fp;
    char fname[200];
    int i,Nz,FirstSliceNumber,M;
    
    Nz = Image->imgparams.Nz;
    FirstSliceNumber = Image->imgparams.FirstSliceNumber;
    M = Image->imgparams.Nx * Image->imgparams.Ny;
    
    //printf("Reading 3-D Image ... \n");
    for(i=0;i<Nz;i++)
    {
        /* slice index currently formed from fixed number of digits */
        sprintf(fname,"%s_slice%.*d.2Dimgdata",basename,MBIR_MODULAR_MAX_NUMBER_OF_SLICE_DIGITS,i+FirstSliceNumber);
	//printf("filename: |%s|\n",fname);
        
        if((fp = fopen(fname,"r")) == NULL) {
            fprintf(stderr, "ERROR in ReadImage3D: can't open file %s\n",fname);
            exit(-1);
        }
        if(fread(Image->image[i],sizeof(float),M,fp)!=M) {
            fprintf(stderr, "ERROR in ReadImage3D: file %s terminated early\n",fname);
            exit(-1);
        }
        fclose(fp);
    }
    return 0;
}


/* Utility for writing 3D image parameters and data */
/* Returns 0 if no error occurs */
int WriteImage3D(
    char *basename,	/* Input: Writes image data to <basename>_slice<n>.2Dimgdata for given slice range */
    struct Image3D *Image)  /* Input: Image data structure (both data and params) */
{
    FILE *fp;
    char fname[200];
    int i,Nz,FirstSliceNumber,M;
    
    Nz = Image->imgparams.Nz;
    FirstSliceNumber = Image->imgparams.FirstSliceNumber;
    M = Image->imgparams.Nx * Image->imgparams.Ny;
    
    //printf("Writing 3-D Image ... \n");
    for(i=0;i<Nz;i++)
    {
        /* slice index currently formed from fixed number of digits */
        sprintf(fname,"%s_slice%.*d.2Dimgdata",basename,MBIR_MODULAR_MAX_NUMBER_OF_SLICE_DIGITS,i+FirstSliceNumber);
	//printf("filename: |%s|\n",fname);
        
        if((fp = fopen(fname,"w")) == NULL) {
            fprintf(stderr, "ERROR in WriteImage3D: can't open file %s\n",fname);
            exit(-1);
        }
        if(fwrite(Image->image[i],sizeof(float),M,fp)!=M) {
            fprintf(stderr, "ERROR in WriteImage3D: file %s terminated early\n",fname);
            exit(-1);
        }
        fclose(fp);
    }
    return 0;
}


/* Utility for allocating memory for Image */
/* Returns 0 if no error occurs */
int AllocateImageData3D(struct Image3D *Image)
{
    Image->image = (float **)multialloc(sizeof(float), 2, Image->imgparams.Nz, Image->imgparams.Nx * Image->imgparams.Ny);
    return 0;
}

/* Utility for freeing memory for Image */
/* Returns 0 if no error occurs */
int FreeImageData3D(struct Image3D *Image)
{
    multifree(Image->image,2);
    return 0;
}

