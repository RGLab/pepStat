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


void getIndices(int *regions, int *nProbes, int *numRegions, int *StartRegion, int *EndRegion);
void callEnrichedRegions(double *MATScores, int *nProbes, int *position, double *dMerge, double *dMax, double *threshold, double *pValues, int *method, int *regions, int *verbose, int *seqNum, int *numRegions);
void MATNullDistribution(int *position, int *nProbes, double *dMax, double *MATScore, double *sigma0, double *mu0, int  *seqNum);
void MATpValue(int nProbes, double *MATScore, double sigma0, double mu0, double *pValues);
double MATcutoffFDR(int *position, int nProbes, double dMax, double *MATScores, double mu0, double FDR, int *regions, int *seqNum);
int mergeMATScores(int *position, int nProbes, double dMerge, double *MATScores, double m0, double cutoff, int sign, int *regions, int *Mum);
void MATScore(double *I, int *nProbes, int *nArraysI, int *position, double *dMax, double *MATScores, int *seqNum);

void callEnrichedRegions(double *MATScores, int *nProbes, int *position, double *dMerge, double *dMax, double *threshold, double *pValues, int *method, int *regions, int *verbose, int *seqNum, int *numRegions)
{

  double sigma0=0, mu0=0, cutoff=0;
  // double totalSeqsPosition=0;

  /** Compute the associated pValues **/
  MATNullDistribution(position, nProbes, dMax, MATScores, &sigma0, &mu0, seqNum);
  
  /** Compute the FDR cutoff **/
  if(*method==1) /** Based on MAT scores **/
  {
    if(*verbose)
    {
      printf("** Merging regions **\n");  
    }        
    numRegions[0] = mergeMATScores(position, *nProbes, *dMerge, MATScores, mu0, *threshold, 1, regions, seqNum);
  }
  else if(*method==2)
  {
    if(*verbose)
    {
      printf("** Calculating P-values **\n");  
    }    
    MATpValue(*nProbes, MATScores, sigma0, mu0, pValues);
    if(*verbose)
    {
      printf("** Merging regions **\n");  
    }        
    numRegions[0] = mergeMATScores(position, *nProbes, *dMerge, pValues, 0, *threshold, -1, regions, seqNum);
  }
  else if(*method==3)
  {
    if(*verbose)
    {
      printf("** Calculating FDR **\n");  
    }    
    
    cutoff=MATcutoffFDR(position, *nProbes, *dMerge, MATScores, mu0, *threshold, regions, seqNum);
    if(*verbose)
    {
      printf("** Merging regions **\n");  
    }    
    /** Form the regions using the threshold found **/
    numRegions[0] = mergeMATScores(position, *nProbes, *dMerge, MATScores, mu0, cutoff, 1, regions, seqNum);
  }
}

void getIndices(int *regions, int *nProbes, int *numRegions, int *StartRegion, int *EndRegion)
{
  int i=0;
  int count=0;
  
  for(i=1;i<=*numRegions;i++)
  {
    while(((regions[count]<i) | (regions[count]==0)) & (count<*nProbes))
    {
      count++;
    }
    StartRegion[i-1]=count+1;
    while((regions[count]==i) & (count<*nProbes))
    {
      count++;
    }
    EndRegion[i-1]=count;
  }
}


/*\ Assume data already sorted by sequence number, then by position */
void MATScore(double *I, int *nProbes, int *nArraysI, int *position, double *dMax, double *MATScores, int *seqNum)
{
  int pMin=0,pMax=0;
  double cMin=0,cMax=0;
  double MI=0;
  double *dataInRegion;
  gsl_vector *vectorInRegion;
  gsl_vector_view vectorTmp;
  int nProbesRegion=0, nProbesNotTrimmed=0;
  int i,j,k;
  double sigma0=0, mu0=0, cutoff=0;

  
  pMin=0;
  pMax=0;

  for(i=0;i<*nProbes;i++)
  {
    /* Current implementation is such that it is within 5 bp from pMin to i and another 5bp from i to pMax*/
    if((seqNum[pMin]!=seqNum[i]) && (seqNum[pMax]!=seqNum[i]))
    {
      pMin = i;
      pMax = i;
    }
    /* with this implementation, pMax and pMin must stay within the same sequence as i. The only time pMax and pMin are not the same as i is when moves out to the next sequence, and by then pMax and pMin should both be set as i to initialize the beginning of the new sequence. The else if clause checks the improbable state when either one of seqNum[pMin] or seqNum[pMax] is not the same as i*/
    else if((seqNum[pMin]!=seqNum[i])||(seqNum[pMax]!=seqNum[i]))
    {
      error("Check that your intensities are ordered by chromosome then by position \n");
    }
    /* It was sorted.. so no need abs(); pMin keeps moving WRT to i */
    while((pMin < i) &&  ((position[i] - position[pMin])>*dMax/2.) && (seqNum[pMin]==seqNum[i]))
    {
      pMin++;
    }
    /* should be position[pMin+1] so it is always within */
    while((pMax<*nProbes) && ((position[pMax+1]-position[i])<=*dMax/2.) && (seqNum[pMax+1]==seqNum[i]) && (pMax+1 < *nProbes))
    {
      pMax++;
    }
    /* If there is only one probe we skip it */
    
    
    if(pMax-pMin>0)
    {
      MI=0;
      nProbesNotTrimmed=0;

      /* For the Immunoprecipitated Array*/
      nProbesNotTrimmed=0;
      nProbesRegion=(pMax-pMin)**nArraysI;
      dataInRegion=I+pMin**nArraysI;

      /** Vector view of the window data **/
      vectorTmp=gsl_vector_view_array(dataInRegion, nProbesRegion);
      /** Allocate a vector to store the probe measurements within the window **/
      vectorInRegion=gsl_vector_alloc(nProbesRegion);
      /** Copy the probe measurements within the window **/
      gsl_vector_memcpy(vectorInRegion,&vectorTmp.vector);
      /** Sort the probe measurements within the window **/
      gsl_sort_vector(vectorInRegion);
      /** Compute the .1 and .9 quantiles **/
      cMin=gsl_stats_quantile_from_sorted_data(vectorInRegion->data,1, nProbesRegion, 0);
      cMax=gsl_stats_quantile_from_sorted_data(vectorInRegion->data,1, nProbesRegion, 1);

      /** Desallocate the memory **/
      gsl_vector_free(vectorInRegion);
      /** Computing the IP treated sample's trimmed mean **/
      for(k=pMin;k<pMax;k++)
      {
        for(j=0;j<*nArraysI;j++)
        {
          if((I[k**nArraysI+j]>=cMin) & (I[k**nArraysI+j]<=cMax))
          {
            MI+=I[k**nArraysI+j];
            nProbesNotTrimmed++;
          }
        }
      }

      /** Computing the MAT Score **/
      if(nProbesNotTrimmed>0)
      {
        MI=MI/nProbesNotTrimmed;
          MATScores[i]=MI;
      }
      else
      {
        MATScores[i]=0.0;
      }
    }
    else
    {
      MATScores[i]=0.0;
    }
  }
  
  // double totalSeqsPosition=0;

  /** Compute the associated pValues **/
  // MATNullDistribution(position, nProbes, dMax, MATScores, &sigma0, &mu0, seqNum);
  // MATpValue(*nProbes, MATScores, sigma0, mu0, pValues);
}

void MATNullDistribution(int *position, int *nProbes, double *dMax, double *MATScores, double *sigma0, double *mu0, int *seqNum)
{
  /** Number of regions used to estimate the null distribution **/
  /**Assume position vector sorted **/
  /**  int nRegionsMax=(int) (((position[nProbes-1]-position[0])/dMax)+5),nRegions=0; **/
  int nRegionsMax = 0, nRegions=0;
  double totalSeqsPosition=0;
  int i,p0,p1;
  int curSeqStartPos=-1;
  int curSeqNum=-1;
  
  for(i=0;i<*nProbes;i++)
  {
    if(curSeqNum!=seqNum[i])
    {
      curSeqNum=seqNum[i];
      curSeqStartPos=position[i];
    }
    if((i+1==*nProbes) | (curSeqNum!=seqNum[i+1]))
    {
      totalSeqsPosition+=(double)(position[i]-curSeqStartPos);
    }
    nRegionsMax = (int)(totalSeqsPosition/(*dMax))+5;
  }

  gsl_vector *NoOverLapMATScores=gsl_vector_calloc(nRegionsMax);
  gsl_vector *NoOverLapMATScoresMinus;

  p0=0;p1=0;
  while(p1<*nProbes)
  {
      /** If a score is zero, keep going **/
    while((p1<*nProbes) && (MATScores[p1]==0))
    {
      p1++;
    }
    gsl_vector_set(NoOverLapMATScores,nRegions,MATScores[p1]);
    nRegions++;
    
    while((p1<*nProbes) && (position[p1]-position[p0])<=(*dMax) && (seqNum[p1]==seqNum[p0]))
    {
      p1++;
    }
    p0=p1;	    
  }

  gsl_sort(NoOverLapMATScores->data, 1, nRegions);
  *mu0=gsl_stats_median_from_sorted_data(NoOverLapMATScores->data, 1, nRegions);
  /*\ Half the Data? */
  /** Multiply by sqrt(2) since we estimated the variance from half of the data **/
  
  NoOverLapMATScoresMinus=gsl_vector_calloc(nRegions);
  
  for(i=0;i<(int)nRegions/2;i++)
  {
    gsl_vector_set(NoOverLapMATScoresMinus,i,gsl_vector_get(NoOverLapMATScores,i));
  }
  for(i=(int)nRegions/2;i<nRegions;i++)
  {
    gsl_vector_set(NoOverLapMATScoresMinus,i,2**mu0-gsl_vector_get(NoOverLapMATScores,i-(int)nRegions/2));
  }
  // gsl_sort(NoOverLapMATScoresMinus->data, 1, nRegions);
  *sigma0=gsl_stats_sd(NoOverLapMATScoresMinus->data, 1, nRegions);
  // *sigma0=sqrt(2)*gsl_stats_sd(NoOverLapMATScores->data, 1, nRegions/2);
  
  gsl_vector_free(NoOverLapMATScores);
  gsl_vector_free(NoOverLapMATScoresMinus);
}

void MATpValue(int nProbes, double *MATScores, double sigma0, double mu0, double *pValues)
{
  /** Number of regions used to estimate the null distribution **/
  int i;
  for(i=0;i<nProbes;i++)
  {
    pValues[i]=1-gsl_cdf_gaussian_P(MATScores[i]-mu0, sigma0);
  }
}

double MATcutoffFDR(int *position, int nProbes, double dMerge, double *MATScores, double mu0, double FDR, int *regions, int *seqNum)
{

  /** Number of regions used to estimate the null distribution **/
  double step=0.05, estimatedFDR=1;
  int nPositive=0,nNegative=0;

  // int INITIALSTEP=1;
  // int notInitial=0;
  double proposeCutoff=0.1;

  while((proposeCutoff<50) && (estimatedFDR > FDR))
  {
    nPositive=mergeMATScores(position, nProbes, dMerge, MATScores, mu0, proposeCutoff, 1, regions, seqNum);
    nNegative=mergeMATScores(position, nProbes, dMerge, MATScores, mu0, -proposeCutoff, -1, regions, seqNum);
    if(nPositive>0)
    {
      estimatedFDR=fmin2(nNegative/((double)nPositive),1.0);
    }
    else
    {
      estimatedFDR=0;
    }
    /* Increase the cutoff */
    proposeCutoff=proposeCutoff+step;
 
  }
  return(proposeCutoff);  
}

int mergeMATScores(int *position, int nProbes, double dMerge, double *MATScores, double m0, double cutoff, int sign, int *regions, int *seqNum)
{
  int i=0,p=0;
  int w=0,ww=0,max=0;
  int nRegions=0;

  while(p<nProbes)
  {

    /*First detection of a MAT Region */
    if(((MATScores[p]-m0>cutoff) && (sign==1)) | ((MATScores[p]-m0<cutoff) && (sign==-1)))
    {
      nRegions++;
      regions[p]=nRegions;
      w=p+1;
      max=p;
      ww=p;
      /*Scan through the whole MAT Region */
      /** Here we make sure we are on the same chromosome **/
      while((w<nProbes) && ((position[w]-position[ww])<=dMerge) && (seqNum[w]==seqNum[ww]))
      {
        if(((MATScores[w]-m0>cutoff) && (sign==1)) || ((MATScores[w]-m0<cutoff) && (sign==-1)))
        {
          max=w;
          ww=max;
        }
        w++;
      }
      for(i=p;i<=max;i++)
      {
        regions[i]=nRegions;
      }
      p=w;
    }
    else
    {
      regions[p]=0;
      p++;
    }
  }
  return(nRegions);
}
