#include <stdio.h>
#include <string>
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <set>
#include <algorithm>
#include <cstdlib>
#include <vector>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <limits>

using namespace std;

int DEBUG=1; // used for experimental and debugging purposes, e.g. testing of the VT feature


// Begin: some assistant functions
int compareint (const void * a, const void * b)
{
    return ( *(int*)a - *(int*)b );
}

// float comparator
int comparefloat(const void *a, const void  *b)
{
    const float *da = (const float *) a;
    const float *db = (const float *) b;
    return (*da > *db) - (*da < *db);
}

// double comparator
int comparedouble(const void *a, const void  *b)
{
    const double *da = (const double *) a;
    const double *db = (const double *) b;
    return (*da > *db) - (*da < *db);
}

double wilson_lower(int ncount, int pcount)
{
    //  float z = 3.7190165; // for alpha=10^-4 one-sided
    double z = 1.96;
    double phat = pcount / double(ncount);
    return max((double)0,((phat + z*z/(2*ncount) - z * sqrt((phat*(1-phat)+z*z/(4*ncount))/ncount))/(1+z*z/ncount)));
}
double wilson_upper(int ncount, int pcount)
{
    double z = 1.96;
    //  float z = 3.7190165;
    double phat = pcount / double(ncount);
    return min((double)1,((phat + z*z/(2*ncount) + z * sqrt((phat*(1-phat)+z*z/(4*ncount))/ncount))/(1+z*z/ncount)));
}

double orconf_lower(double OR, int a, int b, int c, int d)
{
    double z = 1.96;
    return exp(log(OR)-z*sqrt(1/float(a)+1/float(b)+1/float(c)+1/float(d)));
}

double orconf_upper(double OR, int a, int b, int c, int d)
{
    double z = 1.96;
    return exp(log(OR)+z*sqrt(1/float(a)+1/float(b)+1/float(c)+1/float(d)));
}


// End: some assistant functions
// #
// Begin: Main SMA function

void calc_singlemarker(struct MAP *map, uint64_t*** BinSNPsCCFlagsMC, uint32_t nwords, int ncases, int ncontrols, struct WINDOW *window, int nwindows, int nlinestfam, int nlinestped, double pthresh, int optimalrare, int NCT, uint64_t ***BinCarriers, int nsim, string outputname, int minindiv, bool odds, bool oddsconf)
{

    //* Initiate and print out the table's header of the results
    fstream singlemarkerout;
    string singlemarkerfile=outputname+"_singlemarker.txt";
    singlemarkerout.open(singlemarkerfile.c_str(), ios::out);
    if(!singlemarkerout) die("singlemarkerout can not be opened!");
    
    singlemarkerout << "#INFO:";
    singlemarkerout<<"SIMULATIONS="<<nsim<<";";
    singlemarkerout << "NCT="<<NCT<<";";
    singlemarkerout<<"\n";
    singlemarkerout<<"chr\tpos_bp\tpos_var_nr\tcarriers\tp";
    if(odds) singlemarkerout<<"\tOR";
    if(odds && oddsconf) singlemarkerout<<"[95%CI]";
    singlemarkerout<<"\n";

    double *pvaluesMC=new double[nsim]();
    fill_n(pvaluesMC,nsim,1);
    
    double *statsMC=new double[nsim]();
    fill_n(statsMC,nsim,0);
    
    // Loop on all levels (here there is one level)
    for(int l=0; l<nwindows; l++) 
    {	
	// Loop on all window sizes (here one window size, namely 1 is available)
        for(int m=0; m<window[l].n_level; m++) 
	{
	    uint64_t** dummy = new uint64_t*[1]();
            if(!dummy) die("Memory allocation error in dummy!");

            dummy[0]=new uint64_t[nwords]();
            if(!dummy) die("Memory allocation error in dummy!");
	
	    
	    int nInd=ncases+ncontrols;	// # all samples 
	    double ncnc=(double)ncases*(double)ncontrols;	// # all possible values of the tuple space (case,control)
	    
	    // Loop on all variants
            for(int n1=0; n1<window[l].n; n1++) 
	    {
                for (uint32_t p=0; p<nwords; p++) 
		{
                    dummy[0][p]=BinCarriers[l][n1][p]; // all variants
                }
                
		int sY=0;
                int sX=0;
                double stat;
		double OR_COLL;
                double pvalue;
		//
		for (uint32_t p=0; p<nwords; p++)
		{
		    sY += __builtin_popcountll(dummy[0][p] & BinSNPsCCFlagsMC[0][p][1]);
		    sX += __builtin_popcountll(dummy[0][p] & BinSNPsCCFlagsMC[0][p][2]);
                }
                
		int sXY=sX+sY;
		int diff=sX*ncontrols-sY*ncases;

                if((nInd-sX-sY)!=0) 
		{
            	    stat=nInd*(double)diff*(double)diff/(((double)(ncnc)*(double)(sXY)*(double)(nInd-sXY)));
                }
		    else
		{
            	    stat=0;
                }
                
                pvalue=alglib::chisquarecdistribution(1.0, stat);  // P-value calculator
		
		//  Calculate for once the odds ratios (if odds==1)
                if (odds)
		{
		    if(sX == 0 || sY==0 || ncases==sX) OR_COLL=-9999;
		    else OR_COLL=(double)sX*((double)ncontrols-(double)sY)/((double)sY*((double)ncases-(double)sX));
		}
			                                         
		if(pvalue<=pthresh) 
		{
			singlemarkerout<<map[window[l].index[n1]].chr<<"\t"<<map[window[l].index[n1]].pos<<"\t"<<n1+1<<"\t"<<window[l].Ind[n1]<<"\t"<<pvalue;
			if (odds) singlemarkerout<<"\t"<<OR_COLL;
			if(odds && oddsconf){
			  if(OR_COLL!=-9999 && sX!=0 && ncases-sX!=0 && sY!=0 && ncontrols-sY!=0){
			    singlemarkerout<<"["<<orconf_lower(OR_COLL, sX, ncases-sX, sY, ncontrols-sY)<<","<<orconf_upper(OR_COLL, sX, ncases-sX, sY, ncontrols-sY)<<"]";
			  }
			}
			singlemarkerout<<"\n";
		}
                
                for(int n=1; n<nsim+1; n++)
		{
		    int sY=0;
                    int sX=0;
                    double stat;
		    //
		    for (uint32_t p=0; p<nwords; p++)
		    {
			sY += __builtin_popcountll(dummy[0][p] & BinSNPsCCFlagsMC[n][p][1]);
			sX += __builtin_popcountll(dummy[0][p] & BinSNPsCCFlagsMC[n][p][2]);
			
                    }
                    
		    int sXY=sX+sY;
		    int diff=sX*ncontrols-sY*ncases;

                    if((nInd-sX-sY)!=0) 
		    {
                	stat=nInd*(double)diff*(double)diff/(((double)(ncnc)*(double)(sXY)*(double)(nInd-sXY)));
                    }
		    else
		    {
                	stat=0;
                    }

		    if(stat>statsMC[n-1]) statsMC[n-1]=stat;
		    
		}
    	    }
    	
        delete[] dummy[0];
	delete[] dummy;
	
	}
    }
    
    singlemarkerout.close();
    
    
    for(int s=0; s<nsim; s++)
    {
           pvaluesMC[s]=alglib::chisquarecdistribution(1.0, statsMC[s]);
    }
    
    //
                                                                       
    if(nsim>0) // Print out the P-values from the permutations 
    {
	fstream pvalsMC;
	string pvalsMCfile=outputname+"_singlemarker.pvals";
	pvalsMC.open(pvalsMCfile.c_str(), ios::out); // Open a file of list of P-values from the permutations
	if(!pvalsMC) 
	{
	    die("pvalsMC can not be opened!");
	}
	for(int n=0; n<nsim; n++) 
	{
	    pvalsMC<<pvaluesMC[n]<<endl;
	}
	pvalsMC.close(); //close the *.pvales file

        sort(pvaluesMC, pvaluesMC + nsim); // sort the P-values in ascending order
        pvalsMCfile=outputname+"_singlemarker.srt.pvals";
        pvalsMC.open(pvalsMCfile.c_str(), ios::out); // open a file of ordered list of P-values from the permutations
        if(!pvalsMC) 
	{
	    die("pvalsMC can not be opened!");
	}
	for(int n=0; n<nsim; n++) 
	{
	    pvalsMC<<pvaluesMC[n]<<endl;
	}
	pvalsMC.close(); // Close the *.srt.pvals file
    }
    
    // output: alpha at 0.05 in the *.log file for the correction of multiple testing.
    if((nsim+1)>=100) 
    {
	double alpha=pvaluesMC[(int)floor((nsim+1)*0.05)-1];
	ostringstream alphas;
	alphas<<alpha;
	logg(alphas.str()+"\tnominal alpha corresponding to global alpha=0.05");
    }

    delete[] pvaluesMC;

    cout<<"\ndone!"<<endl;

}

// End: Main SMA function

// Begin: Main GECS with variable binning for all possible bins (without reductions)
// NOTICE: This function will be usefull for plotting aims in small genomic regions
 
void vb_ft_allbins(struct MAP *map, uint64_t*** BinSNPsCCFlagsMC, uint32_t nwords, int ncases, int ncontrols,  struct WINDOW *window, int nwindows, int nlinestfam, int nlinestped, double pthresh, int optimalrare, int NCT, uint64_t ***BinCarriers, int nsim, string outputname, int minindiv, bool odds, bool oddsconf)
{
 
    int nindiv=ncases+ncontrols;
    int differentbin=0;
    long unsigned int nbinsnaive=0;

    // initiate and print out the table's header of the results
    fstream rareVB;
    string rarefileVB=outputname+"_gecs_allbins.txt";
    rareVB.open(rarefileVB.c_str(), ios::out);
    rareVB << "#INFO:";
    rareVB<<"NSIM="<<nsim<<";";
    rareVB << "NCT="<<NCT<<";";
    rareVB << "PTHRESHOLD="<<pthresh<<";";
    rareVB << "ALL_BINS=1;";
    rareVB<<"\n";
    rareVB<<"chr\tstart_bp\tend_bp\tstart_var_nr\tend_var_nr\tcarriers\tp";
    if (odds) rareVB<<"\tOR";
    if(odds && oddsconf) rareVB<<"[95%CI]";

    rareVB<<"\n";


    for(int l=0; l<nwindows; l++) 
    {
	nbinsnaive+=(long unsigned int)window[l].n*(window[l].n-1)/2;
	if(optimalrare==0) 
	{
    	    window[l].n_level=0;  //  for levels=1
	}
	for(int m=0; m<window[l].n_level+1; m++) 
	{
    	    uint64_t** dummy1 = new uint64_t*[window[l].n_level+1]();
    	    if(!dummy1) die("Memory allocation error in dummy1!");
    	    uint64_t** dummy2 = new uint64_t*[window[l].n_level+1]();
	    if(!dummy2) die("Memory allocation error in dummy1!");
	    uint64_t** dummy = new uint64_t*[window[l].n_level+1]();
    	    if(!dummy) die("Memory allocation error in dummy1!");
	    dummy1[m]=new uint64_t[nwords]();
    	    if(!dummy1) die("Memory allocation error in dummy1!");
    	    dummy2[m]=new uint64_t[nwords]();
    	    if(!dummy2) die("Memory allocation error in dummy2!");
    	    dummy[m]=new uint64_t[nwords]();
    	    if(!dummy) die("Memory allocation error in dummy2!");

	    int* Ind1 = new int[window[l].n_level+1]();
    	    if(!Ind1) die("Memory allocation error in Ind1!");
      
    	    for(int n1=0; n1<window[l].n; n1++) 
	    {
		memset(dummy[m], 0, nwords*sizeof(*dummy));
		memset(dummy1[m], 0, nwords*sizeof(*dummy2));
		memset(dummy2[m], 0, nwords*sizeof(*dummy2));
		Ind1[m]=0;
		for (uint32_t p=0; p<nwords; p++) 
		{
		    dummy1[m][p]=BinCarriers[l][n1][p]; // all but last variant
		    dummy2[m][p]=0; // all but first variant
		    dummy[m][p]=dummy1[m][p]; // all variants
		    Ind1[m]+=__builtin_popcountll(dummy[m][p]);
		}
		for(int n2=n1; n2<window[l].n; n2++) 
		{
		    differentbin=0; // if new variants are contributed by last SNP in bin, set to 1
		    if(n2!=n1) 
		    {
			for (uint32_t p=0; p<nwords; p++) 
			{
	    		    dummy[m][p] |= BinCarriers[l][n2][p];
	    		    dummy2[m][p] |= BinCarriers[l][n2][p];
	    		    Ind1[m]+=__builtin_popcountll(~dummy1[m][p] & BinCarriers[l][n2][p]);
	    		    if(differentbin==0 && dummy[m][p] != dummy1[m][p]) 
			    {
				differentbin=1;
	    		    }
			    if(differentbin==1) 
			    {
				dummy1[m][p]=dummy[m][p]; // dummy1 becomes dummy for next iteration
	    		    }
			}
			if(Ind1[m]>=(nindiv-minindiv)) break;
		    }
		    int sY=0;
		    int sX=0;
		    double stat;
		    double OR_COLL;
		    double pvalue;
		    for (uint32_t p=0; p<nwords; p++) 
		    {
			sY += __builtin_popcountll(dummy[m][p] & BinSNPsCCFlagsMC[0][p][1]);
			sX += __builtin_popcountll(dummy[m][p] & BinSNPsCCFlagsMC[0][p][2]);
		    }
		    if((ncases+ncontrols-sX-sY)!=0) 
		    {
			stat=((double)ncases + (double)ncontrols)*((double)sX*(double)ncontrols - (double)sY*(double)ncases)*((double)sX*(double)ncontrols-(double)sY*(double)ncases)/((double)ncases*(double)ncontrols*((double)sX+(double)sY)*((double)ncases+(double)ncontrols-(double)sX-(double)sY));
		    }
		    else 
		    {
			stat=0;
		    }
		    if(odds) 
		    {    
			if(sX == 0 || sY==0 || ncases==sX)
			{
			    OR_COLL=-9999;
			}
			else
			{
			    OR_COLL=(double)sX*((double)ncontrols-(double)sY)/((double)sY*((double)ncases-(double)sX));
			}
		    }
		    
		    pvalue=alglib::chisquarecdistribution(1.0, stat);
		    if(pvalue<=pthresh) 
		    {
	    		rareVB<<map[window[l].index[n1]].chr<<"\t"<<map[window[l].index[n1]].pos<<"\t"<<map[window[l].index[n2]].pos<<"\t"<<n1+1<<"\t"<<n2+1<<"\t"<<Ind1[0]<<"\t"<<pvalue;

	    		if (odds) rareVB<<"\t"<<OR_COLL;
			if(odds && oddsconf){
			  if(OR_COLL!=-9999 && sX!=0 && ncases-sX!=0 && sY!=0 && ncontrols-sY!=0){
			    rareVB<<"["<<orconf_lower(OR_COLL, sX, ncases-sX, sY, ncontrols-sY)<<","<<orconf_upper(OR_COLL, sX, ncases-sX, sY, ncontrols-sY)<<"]";
			  }
			}
			rareVB<<endl;
		    }
		}
	    }
      
	    for(int m=0; m<window[l].n_level+1; m++) 
	    {
		delete[] dummy1[m];
		delete[] dummy2[m];
		delete[] dummy[m];
    	    }
    	    delete[] dummy1;
    	    delete[] dummy2;
    	    delete[] dummy;
    	    delete[] Ind1;
	}

    }
    rareVB.close();
    cout<<"\ndone!"<<endl;
}

void vb_ft(struct MAP *map, uint64_t*** BinSNPsCCFlagsMC, uint32_t nwords, int ncases, int ncontrols, struct WINDOW *window, int nwindows, int nlinestfam, int nlinestped, double pthresh, int optimalrare, int NCT, uint64_t ***BinCarriers, int nsim, string outputname, int minindiv, int verbose, bool odds, bool oddsconf)
{
    long unsigned int nbinsdistinct=0;
    long unsigned int nbinsnaive=0;
    int maxindiv=nlinestfam;
    if(minindiv!=0)
    {
	maxindiv=nlinestfam-minindiv+1;
    }

    fstream rareVB;
    ostringstream fileVB;
    fileVB<<outputname<<"_gecs_"<<NCT<<".txt";
    string rarefileVB = fileVB.str();
    rareVB.open(rarefileVB.c_str(), ios::out);
    rareVB << "#INFO:";
    
    rareVB<<"NSIM="<<nsim<<";";
    rareVB << "NCT="<<NCT<<";";
    rareVB << "PTHRESHOLD="<<pthresh<<";";
    rareVB<<"\n";
    rareVB<<"chr\tstart_bp\tend_bp\tstart_var_nr\tend_var_nr\tcarriers\tp";
    if(verbose) rareVB<<"\tstat";
    if(odds) rareVB<<"\tOR";
    if(odds && oddsconf) rareVB<<"[95%CI]";

    rareVB<<"\n";

    double *pvaluesMC=new double[nsim]();
    fill_n(pvaluesMC,nsim,1);
    double *statsMC=new double[nsim]();
    fill_n(statsMC,nsim,0);
    uint64_t *dummy1 = new uint64_t[nwords]();
    if(!dummy1) die("Memory allocation error in dummy1!");
    uint64_t *dummy2 = new uint64_t[nwords]();
    if(!dummy2) die("Memory allocation error in dummy2!");
    uint64_t *dummy = new uint64_t[nwords]();
    if(!dummy) die("Memory allocation error in dummy!");

    int nInd=ncases+ncontrols;
    double ncnc=(double)ncases*(double)ncontrols;

    for(int l=0; l<nwindows; ++l) 
    {
	nbinsnaive+=(long unsigned int)window[l].n*(window[l].n-1)/2;
	for(int n1=0; n1<window[l].n; n1++) 
	{
	    fill_n(dummy,nwords,0ULL);
    	    fill_n(dummy1,nwords,0ULL);
    	    fill_n(dummy2,nwords,0ULL);
    	    int differentbin=0;
    	    if(n1<window[l].n-1) // check if the start variant is non-distinct from the next variant
	    { 
		for (uint32_t p=0; p<nwords; p++) 
		{
		    if(BinCarriers[l][n1][p]!=BinCarriers[l][n1+1][p])
		    {
			differentbin=1;
			break;
		    }
		}
    	    }
	    else if(n1==window[l].n-1) differentbin=1;
    	    if(differentbin==0) continue;
    	    for (uint32_t p=0; p<nwords; p++) 
	    {
		dummy1[p]=BinCarriers[l][n1][p]; // all but last variant
		dummy[p]=dummy1[p]; // all variants
    	    }
	    for(int n2=n1; n2<window[l].n; n2++) 
	    {
		int Ind1=0;
		differentbin=0; // if new variants are contributed by last SNP in bin, set to 1
		int skiprow=1; // if new variants are contributed by first SNP in bin, set to 0
		if(n2!=n1) 
		{
		    for (uint32_t p=0; p<nwords; p++) 
		    {
			dummy2[p] |= BinCarriers[l][n2][p];
		    }
		    for (uint32_t p=0; p<nwords; p++) 
		    {
			if((BinCarriers[l][n1][p] & ~dummy2[p])!=0) 
			{
	    		    skiprow=0;
	    		    break;
			}
		    }
		    if(skiprow==1) break; // if first variant at n1 does not contribute new carriers anymore: continue with next n1
		    //check if first variant still contributes carriers
		    for (uint32_t p=0; p<nwords; p++) 
		    {
			dummy[p] |= BinCarriers[l][n2][p];
		    }
		    for (uint32_t p=0; p<nwords; p++) 
		    {
			if(dummy[p] != dummy1[p]) 
			{
	    		    differentbin=1;
	    		    break;
			}
		    }
		    if(differentbin==0) continue; // if no new carriers: continue with next n2
		    else { // dummy1 becomes dummy for next iteration
			for (uint32_t p=0; p<nwords; p++) 
			{
	    		    dummy1[p]=dummy[p];
			}
		    }
		}
 		int sY=0;
		int sX=0;
		double stat;
	        double OR_COLL;
		double pvalue;

		for (uint32_t p=0; p<nwords; p++) {
		    Ind1+=__builtin_popcountll(dummy[p]);
		}
	    
		if(Ind1>=maxindiv) break; // maxindiv=nlinestfam, unless MININDIV!=0 or 1
		else if(Ind1<minindiv) continue;

		for (uint32_t p=0; p<nwords; p++) 
		{
		    sY+=__builtin_popcountll(dummy[p] & BinSNPsCCFlagsMC[0][p][1]);
		}
		sX=Ind1-sY;
		int sXY=sX+sY;
		int diff=sX*ncontrols-sY*ncases;
	
		if((ncases+ncontrols-sX-sY)!=0) 
		{
		    stat=nInd*(double)diff*(double)diff/(((double)(ncnc)*(double)(sXY)*(double)(nInd-sXY)));
		}
		else
		{
		    stat=0;
		}
		
		pvalue=alglib::chisquarecdistribution(1.0, stat);
		
		if(n1!=n2) 
		{
		    nbinsdistinct+=1;
		}
		if(pvalue<=pthresh) 
		{
		    rareVB<<map[window[l].index[n1]].chr<<"\t"<<map[window[l].index[n1]].pos<<"\t"<<map[window[l].index[n2]].pos<<"\t"<<n1+1<<"\t"<<n2+1<<"\t"<<Ind1;
		    if(verbose) rareVB<<"\t"<<stat;
		    rareVB<<"\t"<<pvalue;
		
		if(odds)
                {
                    if(sX == 0 || sY==0 || ncases==sX)
                    {
                        OR_COLL=-9999;
                    }
                    else
                    {
                        OR_COLL=(double)sX*((double)ncontrols-(double)sY)/((double)sY*((double)ncases-(double)sX));
                    }
                    rareVB<<"\t"<<OR_COLL;
		    if(odds && oddsconf){
		      if(OR_COLL!=-9999 && sX!=0 && ncases-sX!=0 && sY!=0 && ncontrols-sY!=0){
			rareVB<<"["<<orconf_lower(OR_COLL, sX, ncases-sX, sY, ncontrols-sY)<<","<<orconf_upper(OR_COLL, sX, ncases-sX, sY, ncontrols-sY)<<"]";
		      }
		    }
                }
		if(DEBUG==1){
		  rareVB<<"\t";
		  for (uint32_t p=0; p<nwords; p++) 
		    {
		      rareVB<<dummy[p];
		      if(p!=nwords-1) rareVB<< " ";
		    }
		}

		rareVB<<"\n";
		}

		for(int n=1; n<nsim+1; n++) 
		{
		    sY=0;
		    for (uint32_t p=0; p<nwords; p++) 
		    {
			sY+=__builtin_popcountll(dummy[p] & BinSNPsCCFlagsMC[n][p][1]);
		    }
		    sX=Ind1-sY;
		    int sXY=sX+sY;
		    int diff=sX*ncontrols-sY*ncases;
		    if((ncases+ncontrols-sX-sY)!=0) 
		    {
			stat=nInd*(double)diff*(double)diff/(((ncnc)*(double)(sXY)*(double)(nInd-sXY)));
		    }
		    else
		    {
			stat=0;
		    }

		    if(stat>statsMC[n-1]) statsMC[n-1]=stat;
		}
	    }
	}
    }
    delete[] dummy1;
    delete[] dummy2;
    delete[] dummy;

    rareVB.close();

    for(int s=0; s<nsim; s++)
    {
	//cout<<statsMC[s]<<endl;
	pvaluesMC[s]=alglib::chisquarecdistribution(1.0, statsMC[s]);
    }

    if(nsim>0) 
    {
	fstream pvalsMC;
	ostringstream pvalsMCf;
        pvalsMCf<<outputname<<"_gecs_"<<NCT<<".pvals";
        string pvalsMCfile = pvalsMCf.str();
        pvalsMC.open(pvalsMCfile.c_str(), ios::out);

	if(!pvalsMC) 
	{
    	    die("pvalsMC can not be opened!");
	}
	for(int n=0; n<nsim; n++) 
	{
    	    pvalsMC<<pvaluesMC[n]<<endl;
	}
	pvalsMC.close();
	
	sort(pvaluesMC, pvaluesMC + nsim);

	ostringstream srtpvalsMCf;
        srtpvalsMCf<<outputname<<"_gecs_"<<NCT<<".srt.pvals";
        string srtpvalsMCfile = srtpvalsMCf.str();
        pvalsMC.open(srtpvalsMCfile.c_str(), ios::out);
	if(!pvalsMC) 
	{
    	    die("pvalsMC can not be opened!");
	}
	for(int n=0; n<nsim; n++) 
	{
    	    pvalsMC<<pvaluesMC[n]<<endl;
	}
	pvalsMC.close();
    }
    double reduction=(1-nbinsdistinct/(double)nbinsnaive)*100;
    ostringstream reductions;
    ostringstream nbinsnaives;
    ostringstream nbinsdistincts;

    reductions<<reduction;
    nbinsnaives<<nbinsnaive;
    nbinsdistincts<<nbinsdistinct;
    
    logg(nbinsnaives.str()+"\tbins (naive computation, not counting single variants)");
    logg(nbinsdistincts.str()+"\tdistinct bins (not counting single variants)");
    logg(reductions.str()+"%\treduction");

    // output: alpha at 0.05
    if((nsim+1)>=100) 
    {
	double alpha=pvaluesMC[(int)floor((nsim+1)*0.05)-1];
	ostringstream alphas;
	alphas<<alpha;
	logg(alphas.str()+"\tnominal alpha corresponding to global alpha=0.05");
    }

    delete[] statsMC;
    delete[] pvaluesMC;

    cout<<"\ndone!"<<endl;

}

void vb_vt(struct MAP *map, uint64_t*** BinSNPsCCFlagsMC, uint32_t nwords, int ncases, int ncontrols, struct WINDOW *window, int nwindows, int nlinestfam, int nlinestped, double pthresh, int optimalrare, int NCT, uint64_t ***BinCarriers, int nsim, string outputname, int minindiv, int *Ind)
    {

      int nindiv=ncases+ncontrols;
      fstream rareVB;
      string rarefileVB=outputname+"VBVT.txt";
      rareVB.open(rarefileVB.c_str(), ios::out);
      rareVB << "#INFO:";
      rareVB<<"NSIM="<<nsim<<";";
      rareVB << "NCT="<<NCT<<";";
      rareVB << "VT="<<optimalrare<<";";
      rareVB << "PTHRESHOLD="<<pthresh<<";";

      rareVB<<"\n";
      rareVB<<"chr\tstart_bp\tend_bp\tstart_var_nr\tend_var_nr\tcarriers\tp";
      rareVB<<"\n";


      double *pvaluesMC=new double[nsim]();
      fill_n(pvaluesMC,nsim,1);
      double *statsMC=new double[nsim]();
      fill_n(statsMC,nsim,0);


      for(int l=0; l<nwindows; l++) {
	for(int n1=0; n1<window[l].n; n1++) { // startpos: go thourgh all variants that are present at highest level
	  //cout<<l<<" "<<n1<<" "<<window[l].Ind[n1]<<endl;
	  int startlevel=window[l].level_at_NCT[window[l].Ind[n1]-1];
	  //	  cout<<l<<"  "<<window[l].level_at_NCT[window[l].Ind[n1]-1]<<" " <<window[l].n_level<<endl;
	  int levels_at_startpos=window[l].n_level-startlevel; 

	  int differentbin=0;
	  if(n1<window[l].n-1){
	    for (uint32_t p=0; p<nwords; p++) {	      
	      if(BinCarriers[l][n1][p]!=BinCarriers[l][n1+1][p]){
		differentbin=1;
		break;
		}
	    }
	  }
	  if(n1==window[l].n-1) differentbin=1;

	  if(differentbin==0) continue; // if adjacent non distinct bins

	  int* Ind1 = new int[levels_at_startpos](); // keep track of number of individuals per level
	  if(!Ind1) die("Memory allocation error in Ind1!");
	  int* Ind2 = new int[levels_at_startpos](); // number of individuals per level without the first variant
	  if(!Ind2) die("Memory allocation error in Ind2!");
	  uint64_t** dummy1 = new uint64_t*[levels_at_startpos]();
	  if(!dummy1) die("Memory allocation error in dummy1!");
	  uint64_t** dummy2 = new uint64_t*[levels_at_startpos]();
	  if(!dummy2) die("Memory allocation error in dummy1!");
	  uint64_t** dummy = new uint64_t*[levels_at_startpos]();
	  if(!dummy) die("Memory allocation error in dummy1!");

	  for(int m=0; m<levels_at_startpos; m++) {
	    dummy1[m]=new uint64_t[nwords]();
	    if(!dummy1) die("Memory allocation error in dummy1!");
	    dummy2[m]=new uint64_t[nwords]();
	    if(!dummy2) die("Memory allocation error in dummy2!");
	    dummy[m]=new uint64_t[nwords]();
	    if(!dummy) die("Memory allocation error in dummy2!");
	  }

	  Ind1[0]=Ind[n1]; // lowest level=0, startlevel rescaled

	  for (uint32_t p=0; p<nwords; p++) {
	    dummy1[0][p]=BinCarriers[l][n1][p]; // all but last variant
	    dummy2[0][p]=0; // all but first variant
	    dummy[0][p]=dummy1[0][p]; // all variants
	    Ind1[0]+=__builtin_popcountll(dummy[0][p]);
	  }	  
	  int* validlevel=new int[levels_at_startpos](); // vector indicating whether the level needs to be analyzed. encoding: 0=never seen; 1=active; -9: invalid for this level and above
	  validlevel[0]=1; // level at startposition is always active, until no level is active
	  int skipstartpos=1;
	  for(int n2=n1; n2<window[l].n; n2++) { // variation of end position
	    int skipendpos=1;
	    int currentlevel=window[l].level_at_NCT[window[l].Ind[n2]-1];
	    int lastvalidlevel=0; // if the level has not been encounted before, we need to compute this bin from the bin at a lower level, with position=n2-1xs
	    if(currentlevel<=startlevel) {
	      currentlevel=startlevel; // level has to be at least the level of startposition
	    }
	    // go from currentlevel to maxlevel
	    //cout<<n1<<" "<<n2<<" "<<startlevel<<" "<<currentlevel<<" "<<validlevel[currentlevel]<<" 0"<<endl;
	    for(int m=currentlevel-startlevel; m<levels_at_startpos; m++) {		  
	      int skiplevels=1;
	      if(validlevel[m]==1 || (m==currentlevel-startlevel && validlevel[m]!=-9)) {
		if(n1!=n2) {
		  Ind2[m]=0; 
		  differentbin=0; // if new variants are contributed by last SNP in bin, set to 1
		  if(m==currentlevel-startlevel){ // lowest level for this variant
		    if(currentlevel>startlevel){ // if higher than startpos, find lastvalidlevel
		      for(int h=currentlevel-startlevel; h>=0; h--){
			if(validlevel[h]==1){
			  lastvalidlevel=h;
			  break;
			}
		      }
		    }
		    for (uint32_t p=0; p<nwords; p++) {
		      dummy[m][p]=dummy[lastvalidlevel][p] | BinCarriers[l][n2][p]; // all variants
		      dummy2[m][p]=dummy2[lastvalidlevel][p] | BinCarriers[l][n2][p]; // all but first variant
		      Ind2[m]+=__builtin_popcountll(dummy[m][p]);
		      if(differentbin==0 && dummy[m][p]!=dummy1[lastvalidlevel][p]) {
			differentbin=1;
			validlevel[m]=1;
		      }
		      if(differentbin==1) {
			dummy1[m][p]=dummy[m][p]; // dummy1 becomes dummy for next iteration
		      }
		      if((skiplevels==1 || skipstartpos==1) && dummy[m][p]!=dummy2[m][p]) {
			skiplevels=0;
		      }
		    }
		  } // next bin is computed
		  else{ // next bins are at higher levels
		    for (uint32_t p=0; p<nwords; p++) {
		      dummy[m][p] |= BinCarriers[l][n2][p]; // all variants
		      dummy2[m][p] |= BinCarriers[l][n2][p]; // all but first variant
		      Ind2[m]+=__builtin_popcountll(dummy[m][p]);
		      if(differentbin==0 && dummy[m][p]!=dummy1[m][p]) { // check if current end position is redundant
			differentbin=1;
			validlevel[m]=1;
		      }
		      if(differentbin==1) {
			dummy1[m][p]=dummy[m][p]; // dummy1 becomes dummy for next iteration
		      }
		      if((skiplevels==1 || skipstartpos==1) && dummy[m][p]!=dummy2[m][p]) {
			skiplevels=0;
		      }
		    }
		  }
		  if(skiplevels==1 || Ind2[m]>=(nindiv-minindiv)) {
		    if(m==0) {
		      skipstartpos=1;
		    } else if(m>0) {
		      for(int h=m; h<levels_at_startpos; h++) {
			validlevel[h]=-9;
		      }
		    }
		    break;
		  } else if(differentbin==0 || nindiv<minindiv) {
		    continue;
		  }
		} else if(n1==n2) {
		  skipstartpos=0; skiplevels=0;
		}
		if(differentbin==1) lastvalidlevel=m;		
		if(differentbin==0 || skiplevels==1) break;
		for(int n=0; n<nsim+1; n++) {
		  int sY=0;
		  int sX=0;
		  double stat;
		  //	    double OR_COLL;
		  double pvalue;
		  for (uint32_t p=0; p<nwords; p++) {
		    sY += __builtin_popcountll(dummy[m][p] & BinSNPsCCFlagsMC[n][p][1]);
		    sX += __builtin_popcountll(dummy[m][p] & BinSNPsCCFlagsMC[n][p][2]);
		  }
		  //	  cout<<ncases<<" "<<ncontrols<<" "<<sX<<" "<<sY<<endl;
		  if((ncases+ncontrols-sX-sY)!=0) {
		    stat=((double)ncases + (double)ncontrols)*((double)sX*(double)ncontrols - (double)sY*(double)ncases)*((double)sX*(double)ncontrols-(double)sY*(double)ncases)/((double)ncases*(double)ncontrols*((double)sX+(double)sY)*((double)ncases+(double)ncontrols-(double)sX-(double)sY));
		  } else {
		    stat=0;
		  }

		  if(n==0 ) {
		    pvalue=alglib::chisquarecdistribution(1.0, stat);
		    if(pvalue<=pthresh){
		      if(n1!=n2) {
			Ind1[m]=Ind2[m];
		      }
		      //				cout<<map[window[l].levelpos[window[l].n_level-1][n1]].chr<<"\t"<<map[window[l].levelpos[window[l].n_level-1][n1]].pos<<"\t"<<map[window[l].levelpos[window[l].n_level-1][n2]].pos<<"\t"<<n1+1<<"\t"<<n2+1<<"\t"<<Ind1[m]<<"\t"<<pvalue<<endl;
		      rareVB<<map[window[l].index[n1]].chr<<"\t"<<map[window[l].index[n1]].pos<<"\t"<<map[window[l].index[n2]].pos<<"\t"<<n1+1<<"\t"<<n2+1<<"\t"<<Ind1[m]<<"\t"<<pvalue;

		      if(DEBUG==1){
			rareVB<<"\t";
			for (uint32_t p=0; p<nwords; p++) 
			  {
			    rareVB<<" "<<dummy[m][p];
			    if(p!=nwords-1) rareVB<< " ";
			  }
		      }
		      rareVB<<endl;
		    }
		  } else if(n>0){
		    if(stat>statsMC[n-1]){
		      statsMC[n-1]=stat;
		    }
		    // if(pvalue<pvaluesMC[n-1]) {
		    //   pvaluesMC[n-1]=pvalue;
		    // }
		  }
		} // MC-simulations
	      } // validlevel check
	      //	  cout<<validlevel[m]<<endl;
	      else if(validlevel[m]==-9) break; // -9 indicates that all higher levels are non-distinct
	      if(!differentbin) break;

	      if(n1==n2) break; // special case: at bin=variant at startpos, do not go through higher levels - they are all non-distinct
	    } // levels
	    //	    if(skipstartpos!=0)cout<<skipstartpos<<endl;
	    if(skipstartpos==1) break;
	    //	if(skiprow==1) break;
	  } // endpos
	  delete[] validlevel;
	  //	  cout<<n1<<" "<<levels_at_startpos<<endl;
	  for(int m=0; m<levels_at_startpos; m++) {
	    delete[] dummy1[m];
	    delete[] dummy2[m];
	    delete[] dummy[m];
	  }
	  delete[] dummy1;
	  delete[] dummy2;
	  delete[] dummy;
	  delete[] Ind1;
	  delete[] Ind2;
	} // startpos
      } // chromosomes
      rareVB.close();

      for(int s=0; s<nsim; s++){
	pvaluesMC[s]=alglib::chisquarecdistribution(1.0, statsMC[s]);
      }

      if(nsim>0) {

	sort(pvaluesMC, pvaluesMC + nsim);

	fstream pvalsMC;
	string pvalsMCfile=outputname+".pvals";
	pvalsMC.open(pvalsMCfile.c_str(), ios::out);
	if(!pvalsMC) {
	  die("pvalsMC can not be opened!");
	}
	for(int n=0; n<nsim; n++) {
	  pvalsMC<<pvaluesMC[n]<<endl;
	}
	pvalsMC.close();

	ifstream rareVBin;
	rareVBin.open(rarefileVB.c_str(), ios::in);
	if(!rareVBin) {
	  die("rareVBin can not be opened!");
	}

	string rarefileVBout=outputname+"corrVB.txt";
	fstream rareVBout;
	rareVBout.open(rarefileVBout.c_str(), ios::out);
	if(!rareVBout) {
	  die("rareVBout can not be opened!");
	}

	string line;
	int chr,posbp1,posbp2,pos1,pos2,carr;
	double pvaluein;
	getline(rareVBin,line);
	rareVBout<<line<<endl;
	getline(rareVBin,line);
	rareVBout<<line<<"\tp_corr"<<endl;

	while(true) {
	  rareVBin>>chr>>posbp1>>posbp2>>pos1>>pos2>>carr>>pvaluein;
	  if(rareVBin.eof()) {
	    break;
	  }
	  int pcounter=0;
	  double pvaluecorr;
	  for(int n=0; n<nsim; n++) {
	    if(pvaluesMC[n]<pvaluein) {
	      pcounter++;
	    } else {
	      break;
	    }
	  }
	  pvaluecorr=(double)(pcounter+1)/(double)(nsim+1);
	  rareVBout<<chr<<"\t"<<posbp1<<"\t"<<posbp2<<"\t"<<pos1<<"\t"<<pos2<<"\t"<<carr<<"\t"<<pvaluein<<"\t"<<pvaluecorr<<endl;
	}
	rareVBin.close();
	rareVBout.close();

      }

      delete[] pvaluesMC;

      cout<<"\ndone!"<<endl;

    }
