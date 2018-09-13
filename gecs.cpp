#include <cmath>
#include <stdint.h>
#include <ctime>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>
#include <algorithm>
#include <iterator>

int __builtin_popcountll(long long unsigned int x);

#define ITMAX 200
#define FPMIN 1.0e-30
#define EPS 1.19209290e-07F

#include "alglib/specialfunctions.h" 
#include "alglib/specialfunctions.cpp" 
#include "alglib/ap.h"
#include "alglib/ap.cpp"
#include "alglib/alglibinternal.cpp"
#include "alglib/alglibinternal.h"

#include "gecs.hpp"
#include "rare.cpp"


using namespace std;

string inputname = "";
string outputname = "gecs_";
fstream logfile;
string logname = ".log";
string bmpname = "";

string famfile  = "";
string bpedfile = "";
string bmapfile = "";
uint8_t filemode = 0;

int nlinestfam = 0;
int ncases = 0;
int ncontrols = 0;
int nlinestped = 0;

stringstream sstm;
uint32_t nwordsSNPs = 0;
int nwordsalleles = 0;
uint64_t*** BinSNPsCCFlagsMC   = NULL;

struct MAP *map = NULL;

// variables for rare variant analysis
uint64_t ***BinCarriers=NULL;
uint64_t **BinCarriersPermuted=NULL;
uint64_t **vi=NULL; // binary encoded variants (1=carrier; 0=non-carrier)

// parameters
double raref = 0;
int verbose=0;
bool optimalrare = false;
bool singlemarker = false;
bool vb_bmp = false;
int NCT = 0;
int nwindows=0;
int vctable=0;
int allbins=0;
int nsim = 0;
int minindiv=0;
double pthresh=1;

// random seeds
int ix = 2505;
int iy = 11831;
int iz = 23492;



int main(int argc, char *argv[]) {
  srand((int)time(NULL));

  if (argc<2){
    die("Parameter file required!");
  }

  cout<<"-- GECS: genomic exhaustive collapsing scan --\n"<<endl;
  
  // read paramters
  ifstream paramfile;
  paramfile.open(argv[1], ios::in);  
  if(!paramfile) die("Parameter file can not be opened!");
  cout<<"Parameters read from "<<argv[1]<<":\n"<<endl;
  for(string line; getline(paramfile, line);){
    istringstream iss(line);
    string id, val;
    bool error = false;
    if (!(iss >> id)) error = true;	 
    else if (id[0] == '#' || id[0] == '/') continue;
    //    else if (!(iss >> val >> ws) || iss.get() != EOF) error = true;
    if (!(iss >> val)) error = true;	 
    if(!error){
      if(id=="BFILE"){
	cout<<id<<" "<<val<<endl;
	famfile=val+".fam";
	bmapfile=val+".bim";
	bpedfile=val+".bed";
	inputname="BFILE";
      }
      else if(id=="OUTPUTNAME" || id=="OUTPUT" || id=="OUT"){
	cout<<id<<" "<<val<<endl;
	outputname=val;
      }
      else if(id=="OUTBMP" || id=="BMP"){
	cout<<id<<" "<<val<<endl;
	vb_bmp=true;
	bmpname=val;
      }
      else if(id=="SINGLEMARKER" || id=="SM"){
	cout<<id<<" "<<val<<endl;
	singlemarker=atoi(val.c_str());
	if(singlemarker!=0 && singlemarker!=1) die("SINGLEMARKER must be 0 or 1!");
      }
      else if(id=="SIMULATION" || id=="SIMULATIONS" || id=="PERMUTATIONS" || id=="PERMUTATION"){
	cout<<id<<" "<<val<<endl;
	nsim=atoi(val.c_str());
      }
      else if(id=="NCT"){
	cout<<id<<" "<<val<<endl;
	NCT=atoi(val.c_str());
	if(NCT==0) die("NCT must be > 0!");
      }
      else if(id=="VT"){
	cout<<id<<" "<<val<<endl;
	optimalrare=atoi(val.c_str());
	if(optimalrare!=0 && optimalrare!=1) die("VT must be 0 or 1!");
      }
      else if(id=="VERBOSE"){
	cout<<id<<" "<<val<<endl;
	verbose=atoi(val.c_str());
	if(verbose!=0 && verbose!=1) die("VERBOSE must be 0 or 1!");
      }
      else if(id=="ALL_BINS" || id=="ALLBINS"){
	cout<<id<<" "<<val<<endl;
	allbins=atoi(val.c_str());
	if(allbins!=0 && allbins!=1) die("ALL_BINS must be 0 or 1!");
      }
      else if(id=="VCTABLE" || id=="CVTABLE"){
	cout<<id<<" "<<val<<endl;
	vctable=atoi(val.c_str());
	if(vctable!=0 && vctable!=1) die("VCTABLE must be 0 or 1!");
      }
      else if(id=="PTHRESHOLD"){
	cout<<id<<" "<<val<<endl;
	pthresh=atof(val.c_str());
      }
      else if(id=="MININDIV" || id=="MINDINDIV"){
	cout<<id<<" "<<val<<endl;
	minindiv=atof(val.c_str());
      }
    }
    else die("Problem with parameter file!"); 
  }
  paramfile.close();
  
  // output WARNINGS: parameter sanity check
  if(inputname!="BFILE")die("This is not a valid parameter file! Keyword BFILE is required!");
  if(FILE *file = fopen(famfile.c_str(), "r")) fclose(file);
  else die(famfile+" does not exist!");
  if(FILE *file = fopen(bmapfile.c_str(), "r")) fclose(file);
  else die(bmapfile+" does not exist!");
  if(FILE *file = fopen(bpedfile.c_str(), "r")) fclose(file);
  else die(bpedfile+" does not exist!");
  if(allbins==1){
    logg("WARNING: ALL_BINS was selected. SIMULATIONS was set to 0.");
    singlemarker=0;
    nsim=0;
    optimalrare=0;
  }
  if(vctable==1){
    logg("WARNING: VCTABLE was selected. SIMULATIONS was set to 0.");
    singlemarker=0;
    nsim=0;
    optimalrare=0;
  }
  if(pthresh==0) logg("WARNING: PTHRESHOLD is 0. No bins will be written to output file!");
  if(nsim>1 && nsim%10==0) logg("WARNING: the number of permutations is divisible by 10. We strongly advise to use 99..9 permutations instead");
  if(nsim==0) logg("WARNING: analysis is conducted with 0 permutations. No correction for multiple testing will be performed");
  else if(nsim<99) logg("WARNING: analysis is conducted with fewer than 99 permutations. This is not recommended for significance testing at alpha=0.05"); 

  logname = outputname + ".log";
  logfile.open(logname.c_str(), ios::out);
  logfile.clear();

  // count lines in fam file
  nlinestped = 0;
  string line;
  int number_of_lines = 0;
  string famline;
  ifstream ffile;
  ffile.open(famfile.c_str(), ios::in);  
  while (getline(ffile, famline)) nlinestfam++;
  cout << "\nNumber of lines in famfile: " << nlinestfam<<endl;
  ffile.close();
  nwordsSNPs=nlinestfam/64;
  if(nlinestfam%64!=0) nwordsSNPs++;
  nwordsalleles=(nlinestfam*2)/8;
  if((nlinestfam*2)%8!=0) nwordsalleles++;

  if(singlemarker){
    optimalrare=0;
    NCT=nlinestfam;
  }
  
  // allocate permuted affection status
  BinSNPsCCFlagsMC=new uint64_t**[nsim+1];
  if(!BinSNPsCCFlagsMC)die("Memory allocation error in BinSNPsCCFlagsMC!");
  for(int n = 0; n < nsim+1; n++) { // MC simulations
    BinSNPsCCFlagsMC[n]=new uint64_t*[nwordsSNPs];
    if(!BinSNPsCCFlagsMC[n])die("Memory allocation error in BinSNPsCCFlagsMC[n]!");
    for(uint32_t p=0; p<nwordsSNPs; p++) {
      BinSNPsCCFlagsMC[n][p]=new uint64_t[3]();
      if(!BinSNPsCCFlagsMC[n][p])die("Memory allocation error in BinSNPsCCFlagsMC[n][p]!");
    }
  }
  

  
  // read fam file
  stringstream strm;

  ffile.open(famfile.c_str(), ios::in);  
  if(!ffile) die("fam file can not be opened!");
  int counter=0;
  int wordcount=0;
  for(string line; getline(ffile, line);){
    istringstream iss(line);
    string fid,iid,pid,mid,sex,aff;
    if(!(iss>>fid>>iid>>pid>>mid>>sex>>aff)) die("famfile is malformed!");
    else{
      if(pid!="0"){
	strm<<fid<<" "<<iid<<" has a father: "<<pid<<". Currently, analysis is only possible for unrelated individuals!";
	die(strm.str());
      }
      else if(mid!="0"){
	strm<<fid<<" "<<iid<<" has a mother: "<<mid<<". Currently, analysis is only possible for unrelated individuals!";
	die(strm.str());
      }
      else if(aff!="1" && aff!="2" ){
	strm<<fid<<" "<<iid<<" has affection status: "<<aff<<". Analysis is only possible for binary phenotypes (1=control; 2=case)!";
	die(strm.str());
      }
      if(aff=="1") {
	BinSNPsCCFlagsMC[0][wordcount][1]|=mask64bit[counter];
	ncontrols++;
      }
      else if(aff=="2") {
	BinSNPsCCFlagsMC[0][wordcount][2]|=mask64bit[counter];
	ncases++;
      }
    }
    counter++;
    if(counter==64) {
      wordcount++;
      counter=0;
    }
  }
  ffile.close();
  nlinestfam=ncases+ncontrols;
  strm<<"famfile ok, "<<ncases<<" cases and "<<ncontrols<<" controls"<<endl;
  logg(strm.str());
  strm.str(string());

  // read bmap file
  ifstream bifile;
  bifile.open(bmapfile.c_str(), ios::in);  
  if(!bifile) die("bim file can not be opened!");
  for(string line; getline(bifile, line);){
    istringstream iss(line);
    string chr,rs,cm,bp,a1,a2;
    if(!(iss>>chr>>rs>>cm>>bp>>a1>>a2)) die("bimfile is malformed!");
    else{
      if(chr=="0"){
	strm<<"Variant "<<rs<<" is on chromosome 0!";
	die(strm.str());
      }
      else if(bp=="0"){
	strm<<"Variant "<<rs<<" is on pb position 0!";
	die(strm.str());
      }
      else if(a1=="0" || a2=="0"){
	strm<<"Variant "<<rs<<" has alleles "<<a1<<" and "<<a2<<"! The variant appears to be monomorphic!"<<endl;
	logg(strm.str());
	strm.str(string());
      }
    }
    nlinestped++;
  }
  bifile.close();
  nlinestfam=ncases+ncontrols;
  strm<<"bimfile ok, "<<nlinestped<<" variants"<<endl;
  logg(strm.str());
  strm.str(string());

  // read bed file
  ifstream befile(bpedfile.c_str(), ios::binary);
  char *magicnumber=new char[2];
  befile.read(magicnumber, 2);
  if((int)magicnumber[0]!=108 || (int)magicnumber[1]!=27) die("Magic number of bed file is wrong!");
  delete[] magicnumber;
  char *bmode=new char[1];
  befile.read(bmode, 1);
  if((int)bmode[0]==0)die("bed file is in individual-major-mode!");
  delete[] bmode;
  // allocate memory
  vi=new uint64_t*[nlinestped];
  if(!vi)die("Memory allocation error in vi!");

  for(int i=0; i<nlinestped; i++){
    vi[i]=new uint64_t[nwordsSNPs]();
    if(!vi[i])die("Memory allocation error in vi[i]!");
  }
  

  uint64_t x=1;
  char *buffer=new char[nwordsalleles];
  int rarevariants=0;
  int *ncarriers=new int[nlinestped]();

  for(int i=0; i<nlinestped; i++){
    uint32_t countword=0;
    uint32_t counter=0;
    befile.read(buffer,nwordsalleles);
    //    cout<<buffer[0]<<endl;
    for(int j=0; j<nwordsalleles; j++){
      for(uint8_t k=0; k<4; k++){
	counter++;
	if((counter+64*countword)>nlinestfam) break;
	//	uint8_t k2=(k*2+1);
	//      cout<<(buffer[j]&mask8bit[k]>>k2)<<endl;
	if(((uint8_t)buffer[j]&mask8bit[k])==0){
	  //	  cout<<counter+64*countword<<" "<<nwordsSNPs<<" "<<countword<<endl;
	  vi[i][countword] |=mask64bit[counter-1];
	  // cout<<i<<" "<<counter+64*countword<<endl;

	}	 
      }
      if((counter+64*countword)>nlinestfam) break;
      if(counter==64){
	//	cout<<" "<<vi[i][countword]<<" "<<endl;
	countword++;
	counter=0;
	}
    }
    for(int l=0; l<nwordsSNPs; l++){
      ncarriers[i]+=__builtin_popcountll(vi[i][l]);
    }
    if(ncarriers[i]<=NCT && ncarriers[i]>0){
      rarevariants++;
    }
  }
  if(rarevariants==0){
    cout<<"No rare variants left! Choose a higher NCT!"<<endl;
    exit(1);
  }
  else if(rarevariants==1){
    cout<<"1 rare variant left! Choose a higher NCT!"<<endl;
    exit(1);
  }
  strm<<rarevariants<<" rare variants"<<endl;
  logg(strm);
  strm.str(string());
  delete[] buffer;
  befile.close();
  

  // read in mapfile again, this time without common variants: fill struct rv[nrv] with attributes string chrom, string pos 
  //  struct MAP *map=NULL;
  map=new struct MAP[rarevariants];
  int *rv_to_allv=NULL;
  rv_to_allv=new int[rarevariants]();
  //  ifstream bifile;
  bifile.open(bmapfile.c_str(), ios::in);  
  if(!bifile) die("bim file can not be opened!");
  int countrv=0;
  int countv=0;
  string chromtmp="";
  for(string line; getline(bifile, line);){
    istringstream iss(line);
    string chr,rs,cm,a1,a2;
    int bp;
    if(!(iss>>chr>>rs>>cm>>bp>>a1>>a2)) die("bimfile is malformed!");
    else{
      if(chr=="0"){
	strm<<"Variant "<<rs<<" is on chromosome 0!";
	die(strm.str());
      }
      else if(bp==0){
	strm<<"Variant "<<rs<<" is on pb position 0!";
	die(strm.str());
      }
      else if(a1=="0" || a2=="0"){
	strm<<"Variant "<<rs<<" has alleles "<<a1<<" and "<<a2<<"! The variant appears to be monomorphic!";
	logg(strm.str());
      }
      if(chromtmp!=chr){
	chromtmp=chr;
	nwindows++;
      }
      if(ncarriers[countv]<=NCT && ncarriers[countv]>0){
	map[countrv].chr=chr;
	map[countrv].rs=rs;
	map[countrv].pos=bp;
	rv_to_allv[countrv]=countv;
	countrv++;
      }
      countv++;
    }
  }
  bifile.close();
  if(nwindows==1) strm<<nwindows<<" chromosome\n";
  else strm<<nwindows<<" chromosomes\n";
  logg(strm);
  strm.str(string());
  
  window=new struct WINDOW[nwindows];
  if(!window)die("Memory allocation error in window struct!");

  chromtmp="";
  int contigcount=0;
  for(int i=0; i<rarevariants; i++){
    if(chromtmp!=map[i].chr) {
      chromtmp=map[i].chr;
      contigcount++;      
      window[contigcount-1].n=0;
    }
    window[contigcount-1].n++;
  }

  BinCarriers=new uint64_t**[nwindows];
  for(int i=0; i<nwindows; i++){
    BinCarriers[i]=new uint64_t*[window[i].n];
    for(int j=0; j<window[i].n; j++){
      BinCarriers[i][j]=new uint64_t[nwordsSNPs]();
    }
    window[i].index=new int[window[i].n];
  }

  int rvcount=0;
  for(int i=0; i<nwindows; i++){
    window[i].Ind=new int[window[i].n]();
    for(int j=0; j<window[i].n; j++){
      window[i].index[j]=rvcount;
      for(int k=0; k<nwordsSNPs; k++) {
	BinCarriers[i][j][k]=vi[rv_to_allv[rvcount]][k];
	window[i].Ind[j]+=__builtin_popcountll(BinCarriers[i][j][k]);
      }
      rvcount++;
    }
  }

  delete[] ncarriers;
  delete[] rv_to_allv;
    
  for(int i=0; i<nlinestped; i++){
    delete[] vi[i];
  }
  delete[] vi;
  



  if(nwindows>1 && allbins==1) die("ALL_BINS can be run on only one contig!");
  else if(nwindows>1 && vb_bmp==1) die("OUTBMP can be run on only one contig!");


  int *variantsPerChromosome=new int[nwindows]();
  int *Ind=new int[nlinestped]();

  // int *chromosomeCount=new int[25]();  // chromosome count value; in case that not all chromsomes are in the dataset, count those that are present; assume no chromsome has value 0
  // int currentchr=0;
  // for(int m1=0; m1<nlinestped-1; m1++){ 
  //   if(map[m1].chr!=map[m1+1].chr){
  //     currentchr++;
  //     chromosomeCount[atoi(map[m1+1].chr)-1]=currentchr;
  //   }
  // }
  // if(currentchr+1!=nwindows) die("Bug! currentchr+1!=nwindows for some reason.");
 

  int **levelcounter=NULL;

  
  // for FT (only one level)
  if(!optimalrare){
    for(int l=0; l<nwindows; l++){
      window[l].n_level=1;
    }
  }
  // for VT: count occupied levels
  else if(optimalrare){
    levelcounter=new int*[nwindows];
    for(int l=0; l<nwindows; l++){
      levelcounter[l]=new int[NCT+1]();
    }
  }
  
  // count n_at_level

   // for VT: determine n_level
  if(optimalrare){
   for(int l=0; l<nwindows; l++){ 
     for(int m=0; m<window[l].n; m++){ 
       levelcounter[l][window[l].Ind[m]]++;
     }
   }
    for(int l=0; l<nwindows; l++){
      int tmplevelcount=0;
      for(int m=1; m<=NCT; m++){
	if(levelcounter[l][m]!=0){
	  tmplevelcount++;
	  //	  cout<<l<<" "<<m<<" "<<levelcounter[l][m]<<" "<<tmplevelcount<<endl;

	}
      }
      window[l].n_level=tmplevelcount;
      window[l].NCT_at_level=new int[tmplevelcount];
      window[l].level_at_NCT=new int[NCT];
      fill_n(window[l].level_at_NCT,NCT,-9);
      int tmplevelcount2=0;
      for(int m=0; m<NCT; m++){
	if(levelcounter[l][m+1]!=0){
	  window[l].level_at_NCT[m]=tmplevelcount2;
	  window[l].NCT_at_level[tmplevelcount2]=m+1;
	  tmplevelcount2++;
	}	
      }
    }
    for(int l=0; l<nwindows; l++){
      delete[] levelcounter[l];
    }
    delete[] levelcounter;
  }
  /*
    if(optimalrare){
    }
  */
  // allocate space for levelpos
  // if(!optimalrare){
  //   for(int l=0; l<nwindows; l++){
  //     if(!singlemarker) cout<<"Rare SNPs on contig "<<l+1<<": "<<window[l].n_at_level[0]<<endl;
  //     else if(singlemarker) cout<<"SNPs on contig "<<l+1<<": "<<window[l].n_at_level[0]<<endl;
  //     window[l].levelpos[0]=new int[window[l].n_at_level[0]]();
  //   }
  //   // fill levelpos
  //   cout<<"nwindows "<<nwindows<<endl;
  //   int tmpcount=0;
  //   //    cout<<"l "<<l<<endl;
  //   for(int m1=0; m1<nlinestped; m1++){
  //     if(Ind[m1]>0 && Ind[m1]<=NCT){
  // 	window[chromosomeCount[atoi(map[m1].chr)-1]].levelpos[0][tmpcount]=m1; // for FT (only one level)
  // 	//	cout<<"here: window["<<chromosomeCount[atoi(map[m1].chr)-1]<<"].levelpos[0]["<<tmpcount<<"] "<<window[chromosomeCount[atoi(map[m1].chr)-1]].levelpos[0][tmpcount]<<endl;
  // 	tmpcount++;
  // 	if(tmpcount==variantsPerChromosome[chromosomeCount[atoi(map[m1].chr)-1]]){
  // 	  tmpcount=0;
  // 	}
  //     }
  //   }
  // }
  // else if(optimalrare){
  //   for(int m1=0; m1<nlinestped; m1++){
  //     if(Ind[m1]>0 && Ind[m1]<=NCT){
  // 	window[chromosomeCount[atoi(map[m1].chr)-1]].n_at_level[window[chromosomeCount[atoi(map[m1].chr)-1]].level_at_NCT[Ind[m1]-1]]++; // increment n_at_level at locally maximal level
  // 	for (int m2=window[chromosomeCount[atoi(map[m1].chr)-1]].level_at_NCT[Ind[m1]-1]+1; m2<window[chromosomeCount[atoi(map[m1].chr)-1]].n_level; m2++){ // .. and all levels above
  // 	  window[chromosomeCount[atoi(map[m1].chr)-1]].n_at_level[m2]++;
  // 	}
  //     }
  //   }
  //   for(int l=0; l<nwindows; l++){
  //     for(int m1=0; m1<window[l].n_level; m1++){
  // 	//	cout<<"window["<<l<<"].n_at_level["<<m1<<"] "<<window[l].n_at_level[m1]<<endl;
  // 	window[l].levelpos[m1]=new int[window[l].n_at_level[m1]]();
  //     }
  //   }

  //   for(int l=0; l<nwindows; l++){
  //     for(int m=0; m<window[l].n_level; m++){
  // 	//	cout<<"window["<<l<<"].NCT_at_level["<<m<<"] "<<window[l].NCT_at_level[m]<<endl;
  //     }
  //   }
  //   for(int l=0; l<nwindows; l++){
  //     for(int m=0; m<NCT; m++){
  // 	//	cout<<"window["<<l<<"].level_at_NCT["<<m<<"] "<<window[l].level_at_NCT[m]<<endl;
  //     }
  //   }
  //   //    int tmpcount=0;
  //   int **tmpcount_on_level=new int*[nwindows];
  //   for(int l=0; l<nwindows; l++){
  //     tmpcount_on_level[l]=new int[window[l].n_level]();
  //   }


  //   //    cout<<"l "<<l<<endl;
  //   for(int m1=0; m1<nlinestped; m1++){
  //     if(Ind[m1]>0 && Ind[m1]<=NCT){
  // 	//	cout<<"m1 "<<m1<<" chromosomeCount["<<atoi(map[m1].chr)-1<<"] "<<chromosomeCount[atoi(map[m1].chr)-1]<<" window[chromosomeCount["<<atoi(map[m1].chr)-1<<"].n_level "<<window[chromosomeCount[atoi(map[m1].chr)-1]].n_level<<" window["<<chromosomeCount[atoi(map[m1].chr)-1]<<"].level_at_NCT["<<Ind[m1]-1<<"] "<<window[chromosomeCount[atoi(map[m1].chr)-1]].level_at_NCT[Ind[m1]-1]<<" Ind["<<m1<<"] "<<Ind[m1]<<" window["<<chromosomeCount[atoi(map[m1].chr)-1]<<"].n_level "<<window[chromosomeCount[atoi(map[m1].chr)-1]].n_level<<" window["<<chromosomeCount[atoi(map[m1].chr)-1]<<"].n_at_level["<<window[chromosomeCount[atoi(map[m1].chr)-1]].level_at_NCT[Ind[m1]-1]<<"] "<<window[chromosomeCount[atoi(map[m1].chr)-1]].n_at_level[window[chromosomeCount[atoi(map[m1].chr)-1]].level_at_NCT[Ind[m1]-1]] <<" "<<endl;
  // 	window[chromosomeCount[atoi(map[m1].chr)-1]].levelpos[window[chromosomeCount[atoi(map[m1].chr)-1]].level_at_NCT[Ind[m1]-1]][tmpcount_on_level[chromosomeCount[atoi(map[m1].chr)-1]][window[chromosomeCount[atoi(map[m1].chr)-1]].level_at_NCT[Ind[m1]-1]]]=m1;
  // 	tmpcount_on_level[chromosomeCount[atoi(map[m1].chr)-1]][window[chromosomeCount[atoi(map[m1].chr)-1]].level_at_NCT[Ind[m1]-1]]++;
  // 	for (int m2=window[chromosomeCount[atoi(map[m1].chr)-1]].level_at_NCT[Ind[m1]-1]+1; m2<window[chromosomeCount[atoi(map[m1].chr)-1]].n_level; m2++){ // .. and all levels below
  // 	  //	  cout<<m2<<" "<<window[chromosomeCount[atoi(map[m1].chr)-1]].n_at_level[m2]<<" "<<m1<<" "<<tmpcount<<" "<<endl;
  // 	  window[chromosomeCount[atoi(map[m1].chr)-1]].levelpos[m2][tmpcount_on_level[chromosomeCount[atoi(map[m1].chr)-1]][m2]]=m1;
  // 	  tmpcount_on_level[chromosomeCount[atoi(map[m1].chr)-1]][m2]++;
  // 	}
  // 	// tmpcount++;
  // 	// if(tmpcount==variantsPerChromosome[chromosomeCount[atoi(map[m1].chr)-1]]-1){
  // 	//   tmpcount=0;
  // 	// }
  //     }
  //   }
  //   for(int l=0; l<nwindows; l++){
  //     delete[] tmpcount_on_level[l];
  //   }
  //   delete[] tmpcount_on_level;

  // }
  // permute affection status
  for (int n = 1; n <= nsim; n++) {
    srand(int(time(NULL)));
    uint64_t** BinCarriersPermuted=new uint64_t*[nwordsSNPs];
    for(int j=0; j<nwordsSNPs; j++){
      BinCarriersPermuted[j]=new uint64_t[3]();
      BinCarriersPermuted[j][1]=0xFFFFFFFFFFFFFFFFull;
    }
    permute(BinCarriersPermuted, nlinestfam, ncases, &ix, &iy, &iz);
    int newcases=0;
    int newcontrols=0;
    for (uint32_t p=0; p<nwordsSNPs; p++) {
      for(int i=0; i<3; i++) {
	BinSNPsCCFlagsMC[n][p][i]=BinCarriersPermuted[p][i];
      }
      newcontrols+=__builtin_popcountll(BinCarriersPermuted[p][1]);
      newcases+=__builtin_popcountll(BinCarriersPermuted[p][2]);
    }
    for(int j=0; j<nwordsSNPs; j++){
      delete[] BinCarriersPermuted[j];
    }
    delete[] BinCarriersPermuted;
  }

  
  if(singlemarker){
    calc_singlemarker(map, BinSNPsCCFlagsMC, nwordsSNPs, ncases, ncontrols, window,  nwindows, nlinestfam, nlinestped, pthresh, optimalrare, NCT, BinCarriers, nsim, outputname, vb_bmp, minindiv);
  }
  else if(!optimalrare && !vctable){
    if(allbins) vb_ft_allbins(map, BinSNPsCCFlagsMC, nwordsSNPs, ncases, ncontrols, window, nwindows, nlinestfam, nlinestped, pthresh, optimalrare, NCT, BinCarriers, nsim, outputname, vb_bmp, minindiv);
    else vb_ft(map, BinSNPsCCFlagsMC, nwordsSNPs, ncases, ncontrols, window, nwindows, nlinestfam, nlinestped, pthresh, optimalrare, NCT, BinCarriers, nsim, outputname, vb_bmp, minindiv, verbose);
  }
  else if(optimalrare){
    vb_vt(map, BinSNPsCCFlagsMC, nwordsSNPs, ncases, ncontrols, window, nwindows, nlinestfam, nlinestped, pthresh, optimalrare, NCT, BinCarriers, nsim, outputname, vb_bmp, minindiv, Ind);
  }
  
  if(vctable) {
    vctable_write(map, optimalrare, nsim, outputname, logfile, raref, window, nlinestped, verbose, ncases, ncontrols, BinSNPsCCFlagsMC, nwordsSNPs, BinCarriers, nwindows, nlinestfam);
  }

  if(vb_bmp) {
    outRareVB_BMP(map, optimalrare, nsim, outputname, logfile, raref, window, nlinestped, verbose, ncases, ncontrols, BinSNPsCCFlagsMC, nwordsSNPs, bmpname, BinCarriers);
  }


  for(int l=0; l<nwindows; l++){
    for(int j=0; j<window[l].n; j++){
      delete[] BinCarriers[l][j];
    }
    delete[] BinCarriers[l];
  }
  delete[] BinCarriers;

  
  for(int l=0; l<nwindows; l++){
    // delete[] window[l].n_at_level;
    // delete[] window[l].levelpos;
    if(optimalrare){
      delete[] window[l].NCT_at_level;
      delete[] window[l].level_at_NCT;
      delete[] window[l].index;
      delete[] window[l].Ind;      
    }
  }
 
  delete[] window;
  delete[] variantsPerChromosome;
  //  delete[] chromosomeCount;
  delete[] Ind;
  

  for(int n = 0; n < nsim+1; n++) { // MC simulations
    for(uint32_t p=0; p<nwordsSNPs; p++) {
      delete[] BinSNPsCCFlagsMC[n][p];
    }
  }
  for(int n = 0; n < nsim+1; n++) { // MC simulations
    delete[] BinSNPsCCFlagsMC[n];
  }
  delete[] BinSNPsCCFlagsMC;

 
  delete[] map;
  
  logfile.close();

}
