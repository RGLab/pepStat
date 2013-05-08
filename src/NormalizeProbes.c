#include <R.h>
#include <Rmath.h>
#include <float.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_sf_log.h>
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_sort_double.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_sort_vector.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

#include "config.h"

int convertSeq(char *seq, char **alphabet, int *lAlphabet);

void createPairMatrixCount(gsl_matrix *pairNumCount1, gsl_matrix *pairNumCount2, gsl_matrix *pairNumCount3,
  char **seq, char **alphabet, int *lAlphabet);
void createPairMatrixCountSimple(gsl_matrix *pairNumCount1, char **seq, char **alphabet, int *lAlphabet);

void createDesignMatrixPairBinnedSimple(
  gsl_matrix *pairNumCount1,
  int *lAlphabet,
  gsl_matrix *X,
  int* numZpep,
  double* Zpep);

void createDesignMatrixPairBinned(
  gsl_matrix *pairNumCount1,
  gsl_matrix *pairNumCount2,
  gsl_matrix *pairNumCount3,
  int *lAlphabet,
  gsl_matrix *X);

void createDesignMatrixPairBinnedRow(
  gsl_matrix *pairNumCount1,
  gsl_matrix *pairNumCount2,
  gsl_matrix *pairNumCount3,
  int *lAlphabet,
  gsl_vector *XRow, int i);

void createDesignMatrixZpep(gsl_matrix* X, int* numZpep, double* Zpep);

void NormalizeProbes(char **seq, double *y, double *yNormalized, int *nProbes, int *nArrays, int *method,
		int *robust, double *adjRSquare, double *RSquare, double *BIC, double *center, double *outBeta,
		int *betaLength, char **alphabet, int *lAlphabet, int *MATScaling, int *isVerbose, int *numZpep, double* Zpep);

void normArray(char **seq, double *y, double *yNormalized, int *nProbes, int *nArrays, int *method, int *robust,
		double *adjRSquare, double *RSquare, double *BIC, double *center, double *outBeta, int *betaLength,
		char **alphabet, int *lAlphabet, int *MATScaling, int *isVerbose, int *numZpep, double* Zpep,
		gsl_matrix *pairNumCount1, gsl_matrix *pairNumCount2, gsl_matrix *pairNumCount3,
		int nProbesTotal,
		int nVariables,
		int nBins,
		int nProbesPerBin,
		gsl_vector_view yVector,
		int j);

/** Main method to normalize probes**/
void NormalizeProbes(char **seq, double *y, double *yNormalized, int *nProbes, int *nArrays, int *method,
		int *robust, double *adjRSquare, double *RSquare, double *BIC, double *center, double *outBeta,
		int *betaLength, char **alphabet, int *lAlphabet, int *MATScaling, int *isVerbose, int *numZpep, double* Zpep)
{
	
	/** Declaring gsl variables **/
	gsl_vector_view yVector;

	gsl_matrix *pairNumCount1=NULL, *pairNumCount2=NULL, *pairNumCount3=NULL;
	
	int nProbesTotal=*nProbes;
	int j=0,nVariables;
	
	
	/** Number of bins used in the MAT normalization **/
	int nBins=10;
	/** Average number of probes per bins **/
	int nProbesPerBin=(nProbesTotal)/nBins;
		
	if(*method==1) /** This is the MAT model **/
	{
	  nVariables= *numZpep + *lAlphabet;
    // Rprintf("nLetter=%d\n",nVariables);
	  
	  pairNumCount1=gsl_matrix_calloc(nProbesTotal,*lAlphabet + *numZpep);
	  createPairMatrixCountSimple(pairNumCount1, seq, alphabet, lAlphabet);
	}
	
	else if(*method == 3)
	{
	  nVariables= *numZpep + 1;
	  pairNumCount1=gsl_matrix_calloc(nProbesTotal, nVariables);
	}
	
	else if(*method==2)
	{
    nVariables=0+3**lAlphabet;
    // Rprintf("nLetter=%d\n",nVariables);

    pairNumCount1=gsl_matrix_calloc(nProbesTotal,*lAlphabet);
    pairNumCount2=gsl_matrix_calloc(nProbesTotal,*lAlphabet);
    pairNumCount3=gsl_matrix_calloc(nProbesTotal,*lAlphabet);
    createPairMatrixCount(pairNumCount1, pairNumCount2, pairNumCount3, seq, alphabet, lAlphabet);
	}

	DO_NORMALIZE(normArray(seq, y, yNormalized, nProbes, 
						   nArrays, method, 
						   robust, adjRSquare, RSquare, 
						   BIC, center, outBeta, betaLength, alphabet, lAlphabet,
						   MATScaling, isVerbose, numZpep, Zpep,
						   pairNumCount1, pairNumCount2, 
						   pairNumCount3,
						   nProbesTotal,
						   nVariables, nBins, nProbesPerBin,
						   yVector, jj),
				 jj, *nArrays);
	
	if(*method==1 || *method == 3)
	{
	  gsl_matrix_free(pairNumCount1);
	}
	else if(*method==2)
	{
	  gsl_matrix_free(pairNumCount1);
	  gsl_matrix_free(pairNumCount2);
	  gsl_matrix_free(pairNumCount3);	  
	}
	
	if(*isVerbose)
	{
	  Rprintf("** End of NormalizeProbes procedure **\n");
	}
}
      
int convertSeq(char* seq, char **alphabet, int *lAlphabet)
{
	int i=0;
    for(i=0;i<*lAlphabet;i++)
    {
//      Rprintf("a=%s\n",alphabet[i]);
//      Rprintf("seq=%s\n",seq);

    	/* Compare the two strings */
		if(strcmp(seq, alphabet[i]) == 0)
		{
			return i;
		}
	}
	/* Check that we did find the amino acid */
	error("The amino acid is not in the alphabet!\n");
}

void createDesignMatrixPairBinnedSimple(
  gsl_matrix *pairNumCount1,
  int *lAlphabet,
  gsl_matrix *X,
  int* numZpep,
  double* Zpep)
{

  int nProbes=X->size1;
  int i,j;
  int length=*lAlphabet;
  
  for(i=0;i<nProbes;i++)
  {
  	/** Intercept **/
   // gsl_matrix_set(X,i,0,1);
    /** Nucleotide positional pair effects **/
    for(j=0;j<length;j++)
    {
      gsl_matrix_set(X,i,0+j,gsl_matrix_get(pairNumCount1,i,j));
    }
  }
  if(*numZpep > 0)
  {
	  for(i = 0;i < *numZpep; i++)
	  {
		  for(j = 0;j < nProbes; j++)
		  {
			  gsl_matrix_set(X, j, i + length, Zpep[i*nProbes + j]);
		  }
	  }
  }
}


void createDesignMatrixPairBinned(
  gsl_matrix *pairNumCount1,
  gsl_matrix *pairNumCount2,
  gsl_matrix *pairNumCount3,
  int *lAlphabet,
  gsl_matrix *X)
{

  int nProbes=X->size1;
  int i,j;
  int length=*lAlphabet;
  
  for(i=0;i<nProbes;i++)
  {
  	/** Intercept **/
//    gsl_matrix_set(X,i,0,1);
    /** Nucleotide positional pair effects **/
    for(j=0;j<length;j++)
    {
      gsl_matrix_set(X,i,0+j,gsl_matrix_get(pairNumCount1,i,j));
      gsl_matrix_set(X,i,0+length+j,gsl_matrix_get(pairNumCount2,i,j));
      gsl_matrix_set(X,i,0+2*length+j,gsl_matrix_get(pairNumCount3,i,j));
    }
  }
}

void createDesignMatrixPairBinnedRow(
  gsl_matrix *pairNumCount1,
  gsl_matrix *pairNumCount2,
  gsl_matrix *pairNumCount3,
  int *lAlphabet,
  gsl_vector *XRow, int i)
{
	int j, length=*lAlphabet;
	/** Intercept **/
  // gsl_vector_set(XRow,0,1);
	/** Nucleotide positional pair effects **/
	for(j=0;j<length;j++)
	{
		gsl_vector_set(XRow,0+j,gsl_matrix_get(pairNumCount1,i,j));
		gsl_vector_set(XRow,0+length+j,gsl_matrix_get(pairNumCount2,i,j));
		gsl_vector_set(XRow,0+2*length+j,gsl_matrix_get(pairNumCount3,i,j));
	}
}

void createPairMatrixCount(gsl_matrix *pairNumCount1, gsl_matrix *pairNumCount2, gsl_matrix *pairNumCount3,
  char **seq, char **alphabet, int *lAlphabet)
{
  int nProbes=pairNumCount1->size1;
  int i=0,j=0;
  int tempInt;
  int lBin;
  char strTmp[2];
  
  for(i=0;i<nProbes;i++)
  {
  	/* We divide the sequence in 3 bins */
  	lBin=(int)round(strlen(seq[i])/3.);
  	/* First bin */
    for(j=0;j<strlen(seq[i]);j++)
    {
	    strncpy(strTmp, seq[i]+strlen(seq[i])-1-j, 1);
	    strTmp[1]='\0';
		
    	if(j<lBin)
    	{
			tempInt = convertSeq(strTmp, alphabet, lAlphabet);
			gsl_matrix_set(pairNumCount1,i,tempInt,gsl_matrix_get(pairNumCount1,i,tempInt)+1);
    	}
    	else if((j>=lBin) & (j<2*lBin))
    	{
			tempInt = convertSeq(strTmp, alphabet, lAlphabet);
			gsl_matrix_set(pairNumCount2,i,tempInt,gsl_matrix_get(pairNumCount2,i,tempInt)+1);    	
    	}
    	else
    	{
			tempInt = convertSeq(strTmp, alphabet, lAlphabet);
			gsl_matrix_set(pairNumCount3,i,tempInt,gsl_matrix_get(pairNumCount3,i,tempInt)+1);
    	}
    }
  }
}

void createPairMatrixCountSimple(gsl_matrix *pairNumCount1, char **seq, char **alphabet, int *lAlphabet)
{
  int nProbes=pairNumCount1->size1;
  int i=0,j=0;
  int tempInt;
  int lBin;
  char strTmp[2];
  
  for(i=0;i<nProbes;i++)
  {
    for(j=0;j<strlen(seq[i]);j++)
    {
	    strncpy(strTmp, seq[i]+j, 1);
	    strTmp[1]='\0';
		
			tempInt = convertSeq(strTmp, alphabet, lAlphabet);
			gsl_matrix_set(pairNumCount1,i,tempInt,gsl_matrix_get(pairNumCount1,i,tempInt)+1);
    }
  }
}

void createDesignMatrixZpep(gsl_matrix* X, int* numZpep, double* Zpep)
{
	int nProbes = X->size1;
	int i, j;
	for(i = 0; i < nProbes; i++)
	{
		gsl_matrix_set(X,i,0,1);
		for(j = 0; j < *numZpep; j++)
		{
			gsl_matrix_set(X, i, j+1, Zpep[nProbes*j + i]);
		}
	}
	return;
}
	
void normArray(char **seq, double *y, double *yNormalized, int *nProbes, int *nArrays, int *method, 
               int *robust, double *adjRSquare, double *RSquare, double *BIC, double *center, double *outBeta,
               int *betaLength, char **alphabet, int *lAlphabet, int *MATScaling,int *isVerbose, int* numZpep, double* Zpep,
               gsl_matrix *pairNumCount1, gsl_matrix *pairNumCount2, gsl_matrix *pairNumCount3,
               int nProbesTotal,
               int nVariables,
               int nBins,
               int nProbesPerBin,
               gsl_vector_view yVector, int j)
{
  double RSS=0,TSS=0,meanY=0;
  int i=0,k=0,iterMax=10;
  gsl_vector_view xRow;
  gsl_vector *beta, *betaTmp, *weight, *fittedSorted;
  gsl_vector *sdBins=NULL, *yBin=NULL, *XR, *S, *work;
  gsl_permutation *indexSorted;
  gsl_matrix *X, *H, *V;
  FILE *file;
  double centerY=0;

  weight=gsl_vector_calloc(*nProbes);
  fittedSorted=gsl_vector_calloc(nProbesTotal);

  if((*MATScaling)==1)
  {
    sdBins=gsl_vector_calloc(nBins);
    yBin=gsl_vector_calloc(nProbesPerBin+nProbesTotal-nProbesPerBin*(nBins-1)+1);
  }

  XR=gsl_vector_calloc(nVariables);

  if(*method==1) /** This is the MAT model **/
  {
    nVariables = *lAlphabet + *numZpep;
    X=gsl_matrix_calloc(*nProbes,nVariables);
    createDesignMatrixPairBinnedSimple(pairNumCount1, lAlphabet, X, numZpep, Zpep);
  }
  else if(*method == 3)
  {
    nVariables = 1 + *numZpep;
    X=gsl_matrix_calloc(*nProbes,nVariables);
    createDesignMatrixZpep(X, numZpep, Zpep);
  }
  else if(*method==2)
  {
    nVariables=0+3**lAlphabet;
    X=gsl_matrix_calloc(*nProbes,nVariables);
    // Rprintf("*nProbes=%d\n",*nProbes);
    // Rprintf("*nVariables=%d\n",nVariables);
    createDesignMatrixPairBinned(pairNumCount1,
      pairNumCount2,
      pairNumCount3,
      lAlphabet,
      X);
      // file=fopen("matrixX1.txt","w");
      // gsl_matrix_fprintf(file, X, "%g");
      // fclose(file);
      
  }

  /*****************************************************************************************************/

  beta=gsl_vector_alloc(nVariables);
  betaTmp=gsl_vector_alloc(nVariables);
  H=gsl_matrix_calloc(nVariables,nVariables);
  V=gsl_matrix_calloc(nVariables,nVariables);
  S=gsl_vector_calloc(nVariables);
  work=gsl_vector_calloc(nVariables);
  
  indexSorted=gsl_permutation_calloc(nProbesTotal);

  // gsl_blas_dsyrk(CblasLower, CblasTrans, 1.0, X, 0, H);
  gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, X, X, 0.0, H);
  // file=fopen("matrixH1.txt","w");
  // gsl_matrix_fprintf(file, H, "%g");
  // fclose(file);
  
  
  /** Cholesky decomposition of H **/
  // gsl_linalg_cholesky_decomp(H);
  gsl_linalg_SV_decomp(H, V, S, work);
  
  if(*isVerbose)
  {
    Rprintf("** Start normalization of array %d **\n",j);
  }

  meanY=gsl_stats_mean(y+j*nProbesTotal,1,*nProbes);

  /** Select the array **/
  yVector=gsl_vector_view_array(y+j*nProbesTotal, *nProbes);
  
  /* Remove the overall mean */
  gsl_vector_add_constant (&yVector.vector, -meanY);
  
  /** Compute X'y **/
  gsl_blas_dgemv(CblasTrans, 1.0, X, &yVector.vector, 0.0, betaTmp);

  /** Compute (X'X)^{-1}X'y **/
  // gsl_linalg_cholesky_solve(H, betaTmp, beta);
  gsl_linalg_SV_solve (H, V, S, betaTmp, beta);


  for(i=0;i<nProbesTotal;i++)
  {
	  xRow=gsl_matrix_row(X,i);
	  /** Compute the fitted data X(X'X)^{-1}X'y **/
	  gsl_blas_ddot(&xRow.vector, beta, yNormalized+j*nProbesTotal+i);
	  /* Recenter everything */
    yNormalized[j*nProbesTotal+i]+=meanY;
    y[j*nProbesTotal+i]+=meanY;
	  RSS+=gsl_pow_2(y[j*nProbesTotal+i]-yNormalized[j*nProbesTotal+i]);
	  TSS+=gsl_pow_2(y[j*nProbesTotal+i]-meanY);
  }
  
  RSquare[j]=1.0-RSS/TSS;
  adjRSquare[j]=1.0-(nProbesTotal-1.0)/(nProbesTotal-nVariables-1.0)*RSS/TSS;
  BIC[j]=-*nProbes*log(RSS/(*nProbes))-nVariables*log(*nProbes);

  // gsl_vector_fprintf(stdout, beta, "%.5g");
  center[j]=meanY;
  /* If we specify to use the robust estimate */
  if(*robust==1)
  {
    if(*isVerbose)
    {
      Rprintf("** Start robust normalization of array %d **\n",j);      
    }
    
    for(k=0;k<=iterMax;k++)
    {
      /** Reinitialize everything **/
      gsl_matrix_set_zero(H);
      gsl_matrix_set_zero(V);
      gsl_vector_set_zero(S);
      gsl_vector_set_zero(work);
      
      
      
      for(i=0;i<nProbesTotal;i++)
      {
        /** Compute the weights (with 4 degrees of freedom) **/
        gsl_vector_set(weight,i,(4.0+1.0)/(4.0+*nProbes/RSS*gsl_pow_2(y[j*nProbesTotal+i]-yNormalized[j*nProbesTotal+i])));
        /** Multiply each row of X by the weight **/
        xRow=gsl_matrix_row(X,i);
        gsl_vector_scale(&xRow.vector, sqrt(gsl_vector_get(weight,i)));
      }

      /* compute sum u y/sum u */
      centerY=gsl_stats_wmean(weight->data, 1, y+j*nProbesTotal, 1,*nProbes);
      center[j]=centerY;
      
      for(i=0;i<nProbesTotal;i++)
      {
        /** Multiply each observation by its weight **/
        y[j*nProbesTotal+i]=sqrt(gsl_vector_get(weight,i))*(y[j*nProbesTotal+i]-centerY);
      }
      

      // file=fopen("matrixX2.txt","w");
      // gsl_matrix_fprintf(file, X, "%g");
      // fclose(file);

      // gsl_blas_dsyrk(CblasLower, CblasTrans, 1.0, X, 0, H);
      gsl_blas_dgemm (CblasTrans, CblasNoTrans, 1.0, X, X, 0.0, H);
      gsl_linalg_SV_decomp (H, V, S, work);
      /** Cholesky decomposition of H **/
      // gsl_linalg_cholesky_decomp(H);
      
      /** Compute X'y **/      
      gsl_blas_dgemv(CblasTrans, 1.0, X, &yVector.vector, 0.0, betaTmp);
            
      /** Compute (X'X)^{-1}X'y **/
      // gsl_linalg_cholesky_solve(H, betaTmp, beta);
      gsl_linalg_SV_solve (H, V, S, betaTmp, beta);
      

      RSS=0;TSS=0;
      for(i=0;i<nProbesTotal;i++)
      {
          xRow=gsl_matrix_row(X,i);
          /** Compute the fitted data X(X'X)^{-1}X'y **/
          gsl_blas_ddot(&xRow.vector, beta, yNormalized+j*nProbesTotal+i);
          /** reweight the data **/
          y[j*nProbesTotal+i]=y[j*nProbesTotal+i]/sqrt(gsl_vector_get(weight,i))+centerY;
          yNormalized[j*nProbesTotal+i]=yNormalized[j*nProbesTotal+i]/sqrt(gsl_vector_get(weight,i))+centerY;
          RSS+=gsl_pow_2(y[j*nProbesTotal+i]-yNormalized[j*nProbesTotal+i]);
          TSS+=gsl_pow_2(y[j*nProbesTotal+i]-meanY);
          /** reweight X **/
          gsl_vector_scale(&xRow.vector, 1./sqrt(gsl_vector_get(weight,i)));
      }
    }
    /** Compute the residual and total SS **/
    RSquare[j]=1.0-RSS/TSS;
    adjRSquare[j]=1.0-(nProbesTotal-1.0)/(nProbesTotal-nVariables-1.0)*RSS/TSS;
    BIC[j]=-*nProbes*log(RSS/(*nProbes))-nVariables*log(*nProbes);
  }

    // gsl_vector_fprintf(stdout, beta, "%.5g");

  /** Need to scale the values for MAT **/
  if((*MATScaling)==1)
  {
    if(*isVerbose)
    {
      Rprintf("** Start scaling of array %d **\n",j);
    }
    gsl_permutation_init(indexSorted);
    for(i=0;i<nProbesTotal;i++)
      gsl_vector_set(fittedSorted, i, yNormalized[j*nProbesTotal+i]);
    /** Sort the fitted values **/
    gsl_sort_vector_index(indexSorted, fittedSorted);
    /** Compute the standard deviation per bin and scale the values **/
    for(i=0;i<(nBins-1);i++)
    {
      for(k=0;k<nProbesPerBin;k++)
        gsl_vector_set(yBin,k,y[j*nProbesTotal+indexSorted->data[nProbesPerBin*i+k]]);
      gsl_vector_set(sdBins,i,gsl_stats_sd(yBin->data,1,nProbesPerBin));

      for(k=0;k<nProbesPerBin;k++)
        yNormalized[j*nProbesTotal+indexSorted->data[nProbesPerBin*i+k]]=(y[j*nProbesTotal+indexSorted->data[nProbesPerBin*i+k]]-yNormalized[j*nProbesTotal+indexSorted->data[nProbesPerBin*i+k]])/gsl_vector_get(sdBins,i);
    }
    /** For the last bin, we actually used a bit more probes **/
    for(k=0;k<(nProbesTotal-nProbesPerBin*(nBins-1));k++)
      gsl_vector_set(yBin,k,y[j*nProbesTotal+indexSorted->data[nProbesPerBin*(nBins-1)+k]]);
    gsl_vector_set(sdBins,nBins-1,gsl_stats_sd(yBin->data,1,nProbesTotal-nProbesPerBin*(nBins-1)));
    for(k=0;k<(nProbesTotal-nProbesPerBin*(nBins-1));k++)
      yNormalized[j*nProbesTotal+indexSorted->data[nProbesPerBin*(nBins-1)+k]]=(y[j*nProbesTotal+indexSorted->data[nProbesPerBin*(nBins-1)+k]]-yNormalized[j*nProbesTotal+indexSorted->data[nProbesPerBin*(nBins-1)+k]])/gsl_vector_get(sdBins,nBins-1);
  } 
  else /** if no MAT Scaling for MAT or not using MAT **/
  {
    for(i=0;i<nProbesTotal;i++)
    {
      yNormalized[j*nProbesTotal+i]=y[j*nProbesTotal+i]-yNormalized[j*nProbesTotal+i];
    }
  }
  /*****************************************************************************************************/ 

  *betaLength = nVariables;
  for(k=0; k<nVariables; k++)
  {
    outBeta[k+j*nVariables] = beta->data[k];
  }
  /*****************************************************************************************************/
  
  gsl_vector_free(beta);
  gsl_vector_free(betaTmp);
  gsl_matrix_free(X);
  gsl_matrix_free(H);
  gsl_matrix_free(V);
  gsl_vector_free(S);
  gsl_vector_free(work);
  gsl_vector_free(weight);
  gsl_vector_free(fittedSorted);
  if((*MATScaling==1))
  {
    gsl_vector_free(sdBins);
    gsl_vector_free(yBin);
  }
  gsl_permutation_free(indexSorted);
  gsl_vector_free(XR);
}


static const R_CMethodDef cMethods[] = {
  	{"NormalizeProbes", (DL_FUNC) &NormalizeProbes, 19},
		{NULL, NULL, 0}
};

void R_init_pepStat(DllInfo *info)
{
	R_registerRoutines(info, cMethods, NULL, NULL, NULL);
}



