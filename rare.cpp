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


void calc_singlemarker(struct MAP *map, uint64_t*** BinSNPsCCFlagsMC, uint32_t nwords, int ncases, int ncontrols, struct WINDOW *window, int nwindows, int nlinestfam, int nlinestped, double pthresh, int optimalrare, int NCT, uint64_t ***BinCarriers, int nsim, string outputname,  bool vb_bmp, int minindiv)
{

    fstream singlemarkerout;
    string singlemarkerfile=outputname+"Singlemarker.txt";
    singlemarkerout.open(singlemarkerfile.c_str(), ios::out);
    if(!singlemarkerout) die("singlemarkerout can not be opened!");

    singlemarkerout << "#INFO:";
    singlemarkerout<<"NSIM="<<nsim<<";";
    singlemarkerout << "NCT="<<NCT<<";";
    singlemarkerout<<"\n";
    singlemarkerout<<"chr\tpos_bp\tpos_var_nr\tcarriers\tp";
    singlemarkerout<<"\n";

    double *pvaluesMC=new double[nsim]();
    fill_n(pvaluesMC,nsim,1);


    for(int l=0; l<nwindows; l++) {
        for(int m=0; m<window[l].n_level; m++) {
            uint64_t** dummy = new uint64_t*[1]();
            if(!dummy) die("Memory allocation error in dummy!");
            dummy[0]=new uint64_t[nwords]();
            if(!dummy) die("Memory allocation error in dummy!");
            //      cout<<"window["<<l<<"].n_at_level["<<m<<"] "<<window[l].n_at_level[m]<<endl;
            for(int n1=0; n1<window[l].n; n1++) {
                for (uint32_t p=0; p<nwords; p++) {
                    dummy[0][p]=BinCarriers[l][n1][p]; // all variants
                }
                for(int n=0; n<nsim+1; n++) {
                    int sY=0;
                    int sX=0;
                    double stat;
                    //	    double OR_COLL;
                    double pvalue;
                    for (uint32_t p=0; p<nwords; p++) {
		      sY += __builtin_popcountll(dummy[0][p] & BinSNPsCCFlagsMC[n][p][1]);
		      sX += __builtin_popcountll(dummy[0][p] & BinSNPsCCFlagsMC[n][p][2]);
                    }
                    //	  cout<<ncases<<" "<<ncontrols<<" "<<sX<<" "<<sY<<endl;
                    if((ncases+ncontrols-sX-sY)!=0) {
                        stat=((double)ncases + (double)ncontrols)*((double)sX*(double)ncontrols - (double)sY*(double)ncases)*((double)sX*(double)ncontrols-(double)sY*(double)ncases)/((double)ncases*(double)ncontrols*((double)sX+(double)sY)*((double)ncases+(double)ncontrols-(double)sX-(double)sY));
                    } else {
                        stat=0;
                    }
                    /*if(sX == 0 || sY==0 || ncases==sX){
                      OR_COLL=-9999;
                    }
                    else{
                      OR_COLL=(double)sX*((double)ncontrols-(double)sY)/((double)sY*((double)ncases-(double)sX));
                    }
                    if(vb_bmp){
                      if(OR_COLL>0 && stat>statmax_pos){
                    statmax_pos=stat;
                      }
                      else if(OR_COLL<0 && stat>statmax_neg){
                    statmax_neg=stat;
                      }
                      }*/
							  pvalue=alglib::chisquarecdistribution(1.0, stat);
                    if(n==0) {
		      singlemarkerout<<map[window[l].index[n1]].chr<<"\t"<<map[window[l].index[n1]].pos<<"\t"<<n1+1<<"\t"<<window[l].Ind[n1]<<"\t"<<pvalue<<endl;
                    } else {
                        if(pvalue<pvaluesMC[n-1]) {
                            pvaluesMC[n-1]=pvalue;
                        }
                    }
                }
            }
            for(int m=0; m<1; m++) {
                delete[] dummy[m];
            }
            delete[] dummy;
        }
    }
    singlemarkerout.close();


    if(nsim>0) {
      fstream pvalsMC;
      string pvalsMCfile=outputname+"_singlemarker.pvals";
      pvalsMC.open(pvalsMCfile.c_str(), ios::out);
      if(!pvalsMC) {
	die("pvalsMC can not be opened!");
      }
      for(int n=0; n<nsim; n++) {
	pvalsMC<<pvaluesMC[n]<<endl;
      }
      pvalsMC.close();

      sort(pvaluesMC, pvaluesMC + nsim);
      pvalsMCfile=outputname+"_singlemarker_srt.pvals";
      pvalsMC.open(pvalsMCfile.c_str(), ios::out);
      if(!pvalsMC) {
	die("pvalsMC can not be opened!");
      }
      for(int n=0; n<nsim; n++) {
	pvalsMC<<pvaluesMC[n]<<endl;
      }
      pvalsMC.close();
    }
// output: alpha at 0.05
    if((nsim+1)>=100) {
        double alpha=pvaluesMC[(int)floor((nsim+1)*0.05)-1];
        ostringstream alphas;
        alphas<<alpha;
        logg(alphas.str()+"\tnominal alpha corresponding to global alpha=0.05");
    }



    delete[] pvaluesMC;

    cout<<"\ndone!"<<endl;

}

void vb_ft_allbins(struct MAP *map, uint64_t*** BinSNPsCCFlagsMC, uint32_t nwords, int ncases, int ncontrols,  struct WINDOW *window, int nwindows, int nlinestfam, int nlinestped, double pthresh, int optimalrare, int NCT, uint64_t ***BinCarriers, int nsim, string outputname,  bool vb_bmp, int minindiv)
{

  //  int equivstartsnps=0;
  //  nvbstart=0;
  //  int startcounter=0;
  //  int shift=0; // shift is the number of variants on previous chromosomes

 
  int nindiv=ncases+ncontrols;

  long unsigned int nbinsdistinct=0;
  long unsigned int nbinsnaive=0;

  fstream rareVB;
  string rarefileVB=outputname+"VB_allbins.txt";
  rareVB.open(rarefileVB.c_str(), ios::out);
  rareVB << "#INFO:";

  rareVB<<"NSIM="<<nsim<<";";
  rareVB << "NCT="<<NCT<<";";
  rareVB << "VT="<<optimalrare<<";";
  rareVB << "PTHRESHOLD="<<pthresh<<";";
  rareVB << "ALL_BINS=1;";
  rareVB<<"\n";
  rareVB<<"chr\tstart_bp\tend_bp\tstart_var_nr\tend_var_nr\tcarriers\tp";
  /*
    if(verbose==2){
    rareVB<<"\tdistinct_rare_start\tdistinct_rare_end";
    }
    rareVB<<"\tnRV";
    rareVB<<"\tCarrierLevel";

    if(vb_print_perm_stat==0){
    rareVB<<"\tp_asymp";
    }
    else if(vb_print_perm_stat==1){
    rareVB<<"\tteststat_asymp";
    }

    if(nsim!=0 && vb_binwise_corr==true){
    rareVB<<"\tp_bw_corr";
    }
    if(nsim!=0 && vb_print_perm_stat==0){
    rareVB<<"\tp_gw_corr";
    }
  */
  rareVB<<"\n";




  // new approach
  //  nRareLimits=new int[nwindows]();


  
  int differentbin=0;
  
  for(int l=0; l<nwindows; l++) {
    nbinsnaive+=(long unsigned int)window[l].n*(window[l].n-1)/2;
    //    cout<<window[l].n_level<<endl;
    if(optimalrare==0) {
      window[l].n_level=0;  //  for levels=1
    }
    for(int m=0; m<window[l].n_level+1; m++) {
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
      
      for(int n1=0; n1<window[l].n; n1++) {
	memset(dummy[m], 0, nwords*sizeof(*dummy));
	memset(dummy1[m], 0, nwords*sizeof(*dummy2));
	memset(dummy2[m], 0, nwords*sizeof(*dummy2));
	//	cout<<n1<<endl;
	Ind1[m]=0;
	for (uint32_t p=0; p<nwords; p++) {
	  dummy1[m][p]=BinCarriers[l][n1][p]; // all but last variant
	  dummy2[m][p]=0; // all but first variant
	  dummy[m][p]=dummy1[m][p]; // all variants
	  Ind1[m]+=__builtin_popcountll(dummy[m][p]);
	}
	for(int n2=n1; n2<window[l].n; n2++) {
	  //	  	  cout<<n1<<" "<<n2<<" "<<Ind1[m]<<" "<<minindiv<<endl;
	  differentbin=0; // if new variants are contributed by last SNP in bin, set to 1
	  int skiprow=1; // if new variants are contributed by first SNP in bin, set to 0
	  //	  cout<<Ind1[m]<<" "<<minindiv<<endl;
	  if(n2!=n1) {
	    //                        Ind1[m]=0;
	    for (uint32_t p=0; p<nwords; p++) {
	      dummy[m][p] |= BinCarriers[l][n2][p];
	      dummy2[m][p] |= BinCarriers[l][n2][p];
	      Ind1[m]+=__builtin_popcountll(~dummy1[m][p] & BinCarriers[l][n2][p]);
	      if(differentbin==0 && dummy[m][p] != dummy1[m][p]) {
		differentbin=1;
	      }
	      if(differentbin==1) {
		dummy1[m][p]=dummy[m][p]; // dummy1 becomes dummy for next iteration
	      }
	    }
	    if(Ind1[m]>=(nindiv-minindiv)) break;
	  }
	  int sY=0;
	  int sX=0;
	  double stat;
	  //	    double OR_COLL;
	  double pvalue;
	  for (uint32_t p=0; p<nwords; p++) {
	    sY += __builtin_popcountll(dummy[m][p] & BinSNPsCCFlagsMC[0][p][1]);
	    sX += __builtin_popcountll(dummy[m][p] & BinSNPsCCFlagsMC[0][p][2]);
	  }
	  //			sX=Ind1[m]-sY;
	  //	  cout<<ncases<<" "<<ncontrols<<" "<<sX<<" "<<sY<<endl;
	  if((ncases+ncontrols-sX-sY)!=0) {
	    stat=((double)ncases + (double)ncontrols)*((double)sX*(double)ncontrols - (double)sY*(double)ncases)*((double)sX*(double)ncontrols-(double)sY*(double)ncases)/((double)ncases*(double)ncontrols*((double)sX+(double)sY)*((double)ncases+(double)ncontrols-(double)sX-(double)sY));
	  } else {
	    stat=0;
	  }
	  /*if(sX == 0 || sY==0 || ncases==sX){
	    OR_COLL=-9999;
	    }
	    else{
	    OR_COLL=(double)sX*((double)ncontrols-(double)sY)/((double)sY*((double)ncases-(double)sX));
	    }
	    if(vb_bmp){
	    if(OR_COLL>0 && stat>statmax_pos){
	    statmax_pos=stat;
	    }
	    else if(OR_COLL<0 && stat>statmax_neg){
	    statmax_neg=stat;
	    }
	    }*/
	  pvalue=alglib::chisquarecdistribution(1.0, stat);
	    if(pvalue<=pthresh) {
	      rareVB<<map[window[l].index[n1]].chr<<"\t"<<map[window[l].index[n1]].pos<<"\t"<<map[window[l].index[n2]].pos<<"\t"<<n1+1<<"\t"<<n2+1<<"\t"<<Ind1[0]<<"\t"<<pvalue<<endl;
	    }
	}
      }
      
      for(int m=0; m<window[l].n_level+1; m++) {
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

void vb_ft(struct MAP *map, uint64_t*** BinSNPsCCFlagsMC, uint32_t nwords, int ncases, int ncontrols, struct WINDOW *window, int nwindows, int nlinestfam, int nlinestped, double pthresh, int optimalrare, int NCT, uint64_t ***BinCarriers, int nsim, string outputname,  bool vb_bmp, int minindiv, int verbose)
{
  long unsigned int nbinsdistinct=0;
  long unsigned int nbinsnaive=0;
  int maxindiv=nlinestfam;
  if(minindiv!=0){
    maxindiv=nlinestfam-minindiv+1;
  }
  
  fstream rareVB;
  string rarefileVB=outputname+"VB.txt";
  rareVB.open(rarefileVB.c_str(), ios::out);
  rareVB << "#INFO:";

  rareVB<<"NSIM="<<nsim<<";";
  rareVB << "NCT="<<NCT<<";";
  rareVB << "VT="<<optimalrare<<";";
  rareVB << "PTHRESHOLD="<<pthresh<<";";

  rareVB<<"\n";
  rareVB<<"chr\tstart_bp\tend_bp\tstart_var_nr\tend_var_nr\tcarriers";
  if(verbose) rareVB<<"\tstat";
  rareVB<<"\tp";
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

  
  for(int l=0; l<nwindows; ++l) {
    nbinsnaive+=(long unsigned int)window[l].n*(window[l].n-1)/2;
    for(int n1=0; n1<window[l].n; n1++) {
      fill_n(dummy,nwords,0ULL);
      fill_n(dummy1,nwords,0ULL);
      fill_n(dummy2,nwords,0ULL);
      int differentbin=0;
      if(n1<window[l].n-1) { // check if the start variant is non-distinct from the next variant
	for (uint32_t p=0; p<nwords; p++) {
	  if(BinCarriers[l][n1][p]!=BinCarriers[l][n1+1][p]){
	    differentbin=1;
	    break;
	  }
	}
      } else if(n1==window[l].n-1) differentbin=1;
      if(differentbin==0) continue;
      for (uint32_t p=0; p<nwords; p++) {
	dummy1[p]=BinCarriers[l][n1][p]; // all but last variant
	//	dummy2[p]=0; // all but first variant
	dummy[p]=dummy1[p]; // all variants
	//	Ind1=window[l].Ind[n1];
      }
      for(int n2=n1; n2<window[l].n; n2++) {
	int Ind1=0;
	differentbin=0; // if new variants are contributed by last SNP in bin, set to 1
	int skiprow=1; // if new variants are contributed by first SNP in bin, set to 0
	// int newInd=0;
	if(n2!=n1) {
	  for (uint32_t p=0; p<nwords; p++) {
	    dummy2[p] |= BinCarriers[l][n2][p];
	  }
	  for (uint32_t p=0; p<nwords; p++) {
	    if((BinCarriers[l][n1][p] & ~dummy2[p])!=0) {
	      skiprow=0;
	      break;
	    }
	  }
	  if(skiprow==1) break; // if first variant at n1 does not contribute new carriers anymore: continue with next n1
	  //check if first variant still contributes carriers
	  for (uint32_t p=0; p<nwords; p++) {
	    dummy[p] |= BinCarriers[l][n2][p];
	  }
	  for (uint32_t p=0; p<nwords; p++) {
	    if(dummy[p] != dummy1[p]) {
	      differentbin=1;
	      break;
	    }
	  }
	  if(differentbin==0) continue; // if no new carriers: continue with next n2
	  else { // dummy1 becomes dummy for next iteration
	    for (uint32_t p=0; p<nwords; p++) {
	      dummy1[p]=dummy[p];
	    }
	  }
	}
 	int sY=0;
	int sX=0;
	double stat;
	//	    double OR_COLL;
	double pvalue;

	for (uint32_t p=0; p<nwords; p++) {
	  Ind1+=__builtin_popcountll(dummy[p]);
	}

	if(Ind1>=maxindiv) break; // maxindiv=nlinestfam, unless MININDIV!=0 or 1
	else if(Ind1<minindiv) continue;

	   for (uint32_t p=0; p<nwords; p++) {
	  sY+=__builtin_popcountll(dummy[p] & BinSNPsCCFlagsMC[0][p][1]);
	}
	sX=Ind1-sY;
	int sXY=sX+sY;
	int diff=sX*ncontrols-sY*ncases;
	
	if((ncases+ncontrols-sX-sY)!=0) {
	  stat=nInd*(double)diff*(double)diff/(((double)(ncnc)*(double)(sXY)*(double)(nInd-sXY)));
	  //	  cout<<stat<<endl;
	} else{
	  stat=0;
	}

	pvalue=alglib::chisquarecdistribution(1.0, stat);
	if(n1!=n2) {
	  nbinsdistinct+=1;
	}

	if(pvalue<=pthresh) {
	  rareVB<<map[window[l].index[n1]].chr<<"\t"<<map[window[l].index[n1]].pos<<"\t"<<map[window[l].index[n2]].pos<<"\t"<<n1+1<<"\t"<<n2+1<<"\t"<<Ind1;
	  if(verbose) rareVB<<"\t"<<stat;
	  rareVB<<"\t"<<pvalue<<endl;
	}

	for(int n=1; n<nsim+1; n++) {
	  sY=0;
	  for (uint32_t p=0; p<nwords; p++) {
	    sY+=__builtin_popcountll(dummy[p] & BinSNPsCCFlagsMC[n][p][1]);
	  }
	  sX=Ind1-sY;
	  if((ncases+ncontrols-sX-sY)!=0) {
	    stat=nInd*(double)diff*(double)diff/(((ncnc)*(double)(sXY)*(double)(nInd-sXY)));
	    //	    stat=((double)(ncases + ncontrols)*(double)(sX*ncontrols - sY*ncases)*(double)(sX*ncontrols-sY*ncases))/(((double)(ncases*ncontrols)*(double)(sX+sY)*(double)(ncases+ncontrols-sX-sY)));
	  } else{
	    stat=0;
	  }
	  /*if(sX == 0 || sY==0 || ncases==sX){
	    OR_COLL=-9999;
	    }
	    else{
	    OR_COLL=(double)sX*((double)ncontrols-(double)sY)/((double)sY*((double)ncases-(double)sX));
	    }*/
	  if(stat>statsMC[n-1]){
	    statsMC[n-1]=stat;
	  }
	}
      }
    }
  }
  delete[] dummy1;
  delete[] dummy2;
  delete[] dummy;

  rareVB.close();

  for(int s=0; s<nsim; s++){
    pvaluesMC[s]=alglib::chisquarecdistribution(1.0, statsMC[s]);
  }


  if(nsim>0) {
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

    sort(pvaluesMC, pvaluesMC + nsim);
    pvalsMCfile=outputname+"_srt.pvals";
    pvalsMC.open(pvalsMCfile.c_str(), ios::out);
    if(!pvalsMC) {
      die("pvalsMC can not be opened!");
    }
    for(int n=0; n<nsim; n++) {
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
  if((nsim+1)>=100) {
    double alpha=pvaluesMC[(int)floor((nsim+1)*0.05)-1];
    ostringstream alphas;
    alphas<<alpha;
    logg(alphas.str()+"\tnominal alpha corresponding to global alpha=0.05");
  }

  delete[] statsMC;
  delete[] pvaluesMC;


  cout<<"\ndone!"<<endl;

}

    void vb_vt(struct MAP *map, uint64_t*** BinSNPsCCFlagsMC, uint32_t nwords, int ncases, int ncontrols, struct WINDOW *window, int nwindows, int nlinestfam, int nlinestped, double pthresh, int optimalrare, int NCT, uint64_t ***BinCarriers, int nsim, string outputname,  bool vb_bmp, int minindiv, int *Ind)
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
	    //	  cout<<n1<<" "<<n2<<" "<<startlevel<<" "<<currentlevel<<" "<<validlevel[currentlevel]<<" 0"<<endl;
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
		      if(differentbin==0 && dummy[m][p]!=dummy1[lastvalidlevel][p]) { // check if current end position is redundant
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
		    //		  skiprow=0;
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
		  /*if(sX == 0 || sY==0 || ncases==sX){
		    OR_COLL=-9999;
		    }
		    else{
		    OR_COLL=(double)sX*((double)ncontrols-(double)sY)/((double)sY*((double)ncases-(double)sX));
		    }
		    if(vb_bmp){
		    if(OR_COLL>0 && stat>statmax_pos){
		    statmax_pos=stat;
		    }
		    else if(OR_COLL<0 && stat>statmax_neg){
		    statmax_neg=stat;
		    }
		    }*/
		  //	    cout<<OR_COLL<<" " <<statmax_pos<<" "<<statmax_neg<<" "<<stat<<" "<<endl;

		  if(n==0 ) {
		    pvalue=alglib::chisquarecdistribution(1.0, stat);
		    if(pvalue<=pthresh){
		      if(n1!=n2) {
			Ind1[m]=Ind2[m];
		      }
		      //				cout<<map[window[l].levelpos[window[l].n_level-1][n1]].chr<<"\t"<<map[window[l].levelpos[window[l].n_level-1][n1]].pos<<"\t"<<map[window[l].levelpos[window[l].n_level-1][n2]].pos<<"\t"<<n1+1<<"\t"<<n2+1<<"\t"<<Ind1[m]<<"\t"<<pvalue<<endl;
		      rareVB<<map[window[l].index[n1]].chr<<"\t"<<map[window[l].index[n1]].pos<<"\t"<<map[window[l].index[n2]].pos<<"\t"<<n1+1<<"\t"<<n2+1<<"\t"<<Ind1[m]<<"\t"<<pvalue<<endl;
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


void outRareVB_BMP(struct MAP *map, bool optimalrare, int nsim, string outputname, fstream &logfile, double raref, struct WINDOW *window, int nlinestped, int verbose, int ncasesqc, int ncontrolsqc, uint64_t*** BinSNPsCCFlagsMC, int nwords, string bmpname, uint64_t ***BinCarriers)
{

  cout << "\nBMP output of rare variant VB analysis is being generated..." << endl;
  bmpname=bmpname+".bmp";
  unsigned char bmpfileheader[14] = {'B','M', 0,0,0,0, 0,0, 0,0, 54,0,0,0};
  unsigned char bmpinfoheader[40] = {40,0,0,0, 0,0,0,0, 0,0,0,0, 1,0, 24,0};
  unsigned char bmppad[3] = {0,0,0};

  int w,h;

  int snps2plot=window[0].n;

  w=snps2plot;
  h=snps2plot;


  int filesize = 54 + 3*snps2plot*snps2plot;
  bmpfileheader[ 2] = (unsigned char)(filesize    );
  bmpfileheader[ 3] = (unsigned char)(filesize>> 8);
  bmpfileheader[ 4] = (unsigned char)(filesize>>16);
  bmpfileheader[ 5] = (unsigned char)(filesize>>24);
  bmpinfoheader[ 4] = (unsigned char)(       w    );
  bmpinfoheader[ 5] = (unsigned char)(       w>> 8);
  bmpinfoheader[ 6] = (unsigned char)(       w>>16);
  bmpinfoheader[ 7] = (unsigned char)(       w>>24);
  bmpinfoheader[ 8] = (unsigned char)(       h    );
  bmpinfoheader[ 9] = (unsigned char)(       h>> 8);
  bmpinfoheader[10] = (unsigned char)(       h>>16);
  bmpinfoheader[11] = (unsigned char)(       h>>24);

  //    fstream rareVB;


  FILE *rareVB;
  rareVB=fopen(bmpname.c_str(), "wb");

  fwrite(bmpfileheader,1,14,rareVB);
  fwrite(bmpinfoheader,1,40,rareVB);


  int imgstat=0;
  if(imgstat==1) {


    unsigned char *img=NULL;
    img = (unsigned char *)malloc(3*w*h);
    if(!img) {
      die("Memory allocation error in img!\n");
    }
    memset(img,255,3*w*h*sizeof(unsigned char));



    float **stats=NULL;
    stats = new float *[w];
    if(!stats) {
      die("Memory allocation error in stats!\n");
    }
    for(int i=0; i<h; i++) {
      stats[i] = new float [w]();
      if(!stats[i]) {
	die("Memory allocation error in stats[i]!\n");
      }

    }

    uint64_t* dummy1 = new uint64_t[nwords];
    if(!dummy1) die("Memory allocation error in dummy1!");
    float statmax_pos=0;
    float statmax_neg=0;

    for(int j1=0; j1<snps2plot; j1++) {
      dummy1=(uint64_t *)calloc(nwords, sizeof(uint64_t));
      if(!dummy1) die("Memory allocation error in dummy1!");
      for(int i1=j1; i1<snps2plot; i1++) {
	int sX=0,sY=0;
	double stat=0;

	for (int p=0; p<nwords; p++) {
	  dummy1[p] |= BinCarriers[0][i1][p];
	}
	for(int p=0; p<nwords; p++) {
	  sY += __builtin_popcountll(dummy1[p] & BinSNPsCCFlagsMC[0][p][1]);
	  sX += __builtin_popcountll(dummy1[p] & BinSNPsCCFlagsMC[0][p][2]);
	}
	stat=((double)ncasesqc + (double)ncontrolsqc)*((double)sX*(double)ncontrolsqc - (double)sY*(double)ncasesqc)*((double)sX*(double)ncontrolsqc-(double)sY*(double)ncasesqc)/((double)ncasesqc*(double)ncontrolsqc*((double)sX+(double)sY)*((double)ncasesqc+(double)ncontrolsqc-(double)sX-(double)sY));
	if((ncasesqc+ncontrolsqc-sY-sX)==0 || (sX == 0 && sY==0)) stat=0;


	float propNUC=float(sY)/float(ncontrolsqc);
	float propNAC=float(sX)/float(ncasesqc);
	if((propNUC-propNAC)>EPS) { // protective
	  stats[j1][i1]=-(float)stat;
	  if((stat-statmax_neg)>EPS) {
	    statmax_neg=stat;
	  }
	} else if((propNAC-propNUC)>EPS) { // damaging
	  stats[j1][i1]=(float)stat;
	  if((stat-statmax_pos)>EPS) {
	    statmax_pos=stat;
	  }
	}
	//      cout<< stats[j1][i1]<<endl;
      }
      delete[] dummy1;
    }
    cout<<"statmax_neg "<<statmax_neg<<endl;
    cout<<"statmax_pos "<<statmax_pos<<endl;
    float statmax=statmax_pos;
    if(statmax_neg>statmax_pos) statmax=statmax_neg;

    for(int i=0; i<w; i++) {
      for(int j=i; j<h; j++) {
	int y=i;
	int x=j;
	float r=0;
	float g=0;
	float b=0;
	if(stats[i][j]>0) {
	  r=255;
	  g=255-(stats[i][j]/statmax)*255;
	  b=255-(stats[i][j]/statmax)*255;
	} else if(stats[i][j]<0) {
	  b=255-(-stats[i][j]/statmax)*255;
	  r=255-(-stats[i][j]/statmax)*255;
	  g=255;
	}

	if ((r-255)>EPS) r=0;
	if ((g-255)>EPS) g=0;
	if ((b-255)>EPS) b=255;

	img[(x+y*w)*3+2] = (unsigned char)(r);
	img[(x+y*w)*3+1] = (unsigned char)(g);
	img[(x+y*w)*3+0] = (unsigned char)(b);

	//

	//               if (r<EPS) r=0;
	// 		            if (g<EPS) g=0;
	// 			            if (b<EPS) b=0;

	// 				          if(fabs(stats[i][j])<EPS){
	// 					              break;
	// 						      }
      }
    }

    for(int i=0; i<h; i++) {
      fwrite(img+(w*(h-i-1)*3),3,w,rareVB);
      fwrite(bmppad,1,(4-(w*3)%4)%4,rareVB);
    }

    free(stats);

  }

  else if(imgstat==0) {
    
    for(int j1=snps2plot-1; j1>=0; j1--) {
      unsigned char *img=NULL;
      img = (unsigned char *)calloc(3*w, sizeof(unsigned char));
      //      if(!img) die("Memory allocation error in img!");
      memset(img,255,3*w*sizeof(unsigned char));

      uint64_t* dummy1 = new uint64_t[nwords]();

      // for AND only
      
      //      for (int p=0; p<nwords; p++) {
      //	dummy1[p]=~dummy1[p];
      //	}
       
      if(!dummy1) die("Memory allocation error in dummy1!");
      
      for(int i1=j1; i1<snps2plot; i1++) {


	uint64_t dummy2=0;

	//	int sX=0,sY=0;
	//	double stat=0;
	

	int y=j1;
	int x=i1;

	int Ind=0;
	for (int p=0; p<nwords; p++) {
	  dummy1[p] |= BinCarriers[0][i1][p];
	  //	  dummy1[p] &=  (BinCarriers[0][i1][p]);
	  Ind+=__builtin_popcountll(dummy1[p]);
	  //	  cout<<Ind<<endl;
	}

	int nlinestfam=ncasesqc+ncontrolsqc;
	if(Ind==nlinestfam){
	//		if(Ind==0){
	  break;
	}

	uint64_t* dummycolor = new uint64_t[3]();
	if(!dummycolor) die("Memory allocation error in dummycolor!");

	for(int l=0; l<nwords/3; l++) {
	  dummycolor[2]^=dummy1[l];
	}
	dummycolor[2]=dummycolor[2]%255;

	for(int m=nwords/3; m<2*(nwords/3); m++) {
	  dummycolor[1]^=dummy1[m];
	}
	dummycolor[1]=dummycolor[1]%255;

	for(int n=2*(nwords/3); n<nwords; n++) {
	  dummycolor[0]^=dummy1[n];
	}
	dummycolor[0]=dummycolor[0]%255;

	// for(int n=0; n<nwords-1; n++) {
	//   dummy2^=dummy1[n];
	// }
	// dummycolor[0]=(dummy2 & 0xFFFFF00000000000)%255;
	// dummycolor[1]=(dummy2 & 0x00000FFFFF000000)%255;
	// dummycolor[2]=(dummy2 & 0x0000000000FFFFFF)%255;


	/*
	  int tmpcounter=0;
	  int indthird=nlinestfam/3;
	  int indrest=nlinestfam%3;
	  int c=1;
	  int word=0;
	  if(Ind!=nlinestfam){
	  for(int n=0; n<nlinestfam; n++) {
	  if(n>(indthird*c-1)){
	  c++;
	  }
	  if(n!=0 && n%64==0){
	  word++;
	  }
	  //		    dummycolor[n%3]+=getbit64(dummy1[word],n%64);
	  dummycolor[n%3]+=((dummy1[n/64] >> (n%64)) & 1);
	  }
	  dummycolor[0]-=255*Ind/nlinestfam;
	  dummycolor[1]-=255*Ind/nlinestfam;
	  dummycolor[2]-=255*Ind/nlinestfam;
		  
	  //		  dummycolor[0]^=dummycolor[1];
	  //		  dummycolor[1]^=dummycolor[2];
	  //		  dummycolor[2]^=dummycolor[0];
		  
	  } else{
	  for(int n=0; n<3; n++){
	  dummycolor[n]=0;
	  }
	  }
	*/
	img[(x)*3+2] = (unsigned char)(dummycolor[0]);
	img[(x)*3+1] = (unsigned char)(dummycolor[1]);
	img[(x)*3+0] = (unsigned char)(dummycolor[2]);
	delete[] dummycolor;
      }
      fwrite(img,3,w,rareVB);
      fwrite(bmppad,1,(4-(w*3)%4)%4,rareVB);

      delete[] dummy1;

      free(img);

    }


  //   for(int i=0; i<w; i++)
  //     {
  // 	for(int j=i; j<h; j++)
  // 	  {
  // 	    int y=i;
  // 	    int x=j;
  // 	    float r=0;
  // 	    float g=0;
  // 	    float b=0;
  // 	    if(stats[i][j]>0){
  // 	      r=255;
  // 	      g=255-(stats[i][j]/statmax_pos)*255;
  // 	      b=255-(stats[i][j]/statmax_pos)*255;
  // 	    }
  // 	    else if(stats[i][j]<0){
  // 	      b=255;
  // 	      r=255-(-stats[i][j]/statmax_neg)*255;
  // 	      g=255-(-stats[i][j]/statmax_neg)*255;
  // 	    }

  // 	    if ((r-255)>EPS) r=255;
  // 	    if ((g-255)>EPS) g=255;
  // 	    if ((b-255)>EPS) b=255;

  // 	    img[(x+y*w)*3+2] = (unsigned char)(r);
  // 	    img[(x+y*w)*3+1] = (unsigned char)(g);
  // 	    img[(x+y*w)*3+0] = (unsigned char)(b);


  // 	    //      if (r<EPS) r=0;
  // 	    //      if (g<EPS) g=0;
  // 	    //      if (b<EPS) b=0;

  // 	    //      if(fabs(stats[i][j])<EPS){
  // 	    //      break;
  // 	  }
  //     }
  //}


}

  fclose(rareVB);
  cout << "Rare variant VB BMP output is finished!"<<endl;
  logfile << "Rare variant VB BMP output is finished!"<<endl;
  exit(0);
}


void vctable_write(struct MAP *map, bool optimalrare, int nsim, string outputname, fstream &logfile, double raref, struct WINDOW *window, int nlinestped, int verbose, int ncasesqc, int ncontrolsqc, uint64_t*** BinSNPsCCFlagsMC, int nwords, uint64_t ***BinCarriers, int nwindows, int nlinestfam)
{



  fstream rareVBcas;
  fstream rareVBcon;
  string outnamecas=outputname+"casvctable.txt";
  string outnamecon=outputname+"convctable.txt";
  rareVBcas.open(outnamecas.c_str(), ios::out);
  rareVBcon.open(outnamecon.c_str(), ios::out);

  int counter=0;
  for(int i=0; i<nwindows; i++){
    for(int j=1448; j<2187; j++){
      for(int p=0; p<nwords; p++){
	for(int pos=0; pos<64; pos++){
	  if((BinSNPsCCFlagsMC[0][p][1] & mask64bit[pos])==0){ // case
	    //	    counter++;
	    if(counter>nlinestfam) break;
	    if((BinCarriers[i][j][p] & mask64bit[pos])==0) rareVBcas<<"0 ";
	    else rareVBcas<<"1 ";
	  }
	  if((BinSNPsCCFlagsMC[0][p][2] & mask64bit[pos])==0){ // control
	    if((BinCarriers[i][j][p] & mask64bit[pos])==0) rareVBcon<<"0 ";
	    //	    counter++;
	    else rareVBcon<<"1 ";
	    if(counter>nlinestfam) break;
	  }
	}
	if(counter>nlinestfam) break;
      }
      rareVBcas<<"\n";
      rareVBcon<<"\n";
      if(counter>nlinestfam) break;
    }
  }
  rareVBcas.close();
  rareVBcon.close();
  cout << "Table output is finished!"<<endl;
  logfile << "Table variant VB BMP output is finished!"<<endl;
  exit(0);
}
