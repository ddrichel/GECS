const uint8_t mask8bit[4]={1,4,16,64};
const uint64_t mask64bit[64] = {0x8000000000000000ull,0x4000000000000000ull,0x2000000000000000ull,0x1000000000000000ull,0x800000000000000ull,0x400000000000000ull,0x200000000000000ull,0x100000000000000ull,0x80000000000000ull,0x40000000000000ull,0x20000000000000ull,0x10000000000000ull,0x8000000000000ull,0x4000000000000ull,0x2000000000000ull,0x1000000000000ull,0x800000000000ull,0x400000000000ull,0x200000000000ull,0x100000000000ull,0x80000000000ull,0x40000000000ull,0x20000000000ull,0x10000000000ull,0x8000000000ull,0x4000000000ull,0x2000000000ull,0x1000000000ull,0x800000000ull,0x400000000ull,0x200000000ull,0x100000000ull,0x80000000ull,0x40000000ull,0x20000000ull,0x10000000ull,0x8000000ull,0x4000000ull,0x2000000ull,0x1000000ull,0x800000ull,0x400000ull,0x200000ull,0x100000ull,0x80000ull,0x40000ull,0x20000ull,0x10000ull,0x8000ull,0x4000ull,0x2000ull,0x1000ull,0x800ull,0x400ull,0x200ull,0x100ull,0x80ull,0x40ull,0x20ull,0x10ull,0x8ull,0x4ull,0x2ull,0x1ull};

static const char *USAGE_MESSAGE =
"GECS - Genomic Exhaustive Collapsing Scan for association analysis of rare variants\n"
"Copyright George Kanoungi, Dmitriy Drichel (2019)\n"
"Version 1.1\n\n"
  
"All parameters to GECS can be passed only by the parameter file ~/path_to/*.param\n"
"Usage: ./gecs ~/path_to/*.param\n\n"
"Content of param file:\n\n"
"BFILE\t<string>\t// Prefix of the plink binary file\n"
"PERMUTATIONS\t<integer>\t// Number of permutations, default=0\n"
"SINGLEMARKER\t<bool>\t// Perform pointwise single marker analysis (SMA) instead of the region-based exhaustive scan, default=0\n"
"NCT\t<integer>\t// Threshold of the number of carriers, (equivalent to MAFT)\n"
"MAFT\t<double>\t// Minor allele frequency threshold, alternative to the NCT keyword\n"
"OR\t<bool>\t// Whether the odds ratios should be calculated, default=0\n"
"CORRECTED_P\t<bool>\t// Whether corrected p-values with Wilson score interval of CI 95% shall be computed, default=0\n"
"PTHRESHOLD\t<double>\t// Max. nominal p-value for bins to be written to output files, default=1\n"
"ALLBINS\t<bool>\t// Whether locally not-distinct bins should be written to output (useful for plotting of small regions), default=0\n"
"OUTPUTNAME\t<string>\t// Output prefix\n\n\n"
"For more information please check the README file\n\n"
"--ENJOY!!--\n\n";


void logg(string logstring){
  extern fstream logfile;
  logfile.clear();
  cout<<logstring<<endl;
  logfile<<logstring<<endl;
}
void logg(stringstream& ss){
  logg(ss.str());
}

void die(string logstring){
  extern fstream logfile;
  logg(logstring);
  exit(1);
}


double HillRandom(int *ix, int *iy, int *iz) {
    double random = 0;
    double v;
    *ix = 171 * ((*ix) % 177) - 2 * (*ix / 177);
    *iy = 172 * ((*iy) % 176) - 35 * (*iy / 176);
    *iz = 170 * ((*iz) % 178) - 63 * (*iz / 178);
    if (*ix < 0) *ix = *ix + 30269;
    if (*iy < 0) *iy = *iy + 30307;
    if (*iz < 0) *iz = *iz + 30323;
    random = modf(((float) (*ix)) / 30269.0 + ((float) (*iy)) / 30307.0 + ((float) (*iz)) / 30323.0, &v);
    return random;
}


void permute(uint64_t** BinCarriersPermuted, int nlinestfam, int ncases, int *ix, int *iy, int *iz) {
  double random1;
  int i, newcases = 0;
  int word, pos=0;
  while (newcases < ncases) {
    random1 = HillRandom(ix, iy, iz);
    i = (int)(nlinestfam * random1);
    word=i/64;
    pos=i%64;
    if ((BinCarriersPermuted[word][1] & mask64bit[pos]) != 0) {
      BinCarriersPermuted[word][1] &= ~(mask64bit[pos]);
      BinCarriersPermuted[word][2] |= mask64bit[pos];
      newcases++;
    }
  }
}


struct WINDOW {
  int n; // rare variant count per contig
  int *Ind;
  int *index; // rare variant count per contig
  int **levelpos; // relative positions of SNPs in bim file,  levelwise
  int n_level;
    //  int start;
    //  int end;
    //  int start1,start2,end1,end2; // For RAREINTER

  int *n_at_level; // how many SNPs per level

  int *NCT_at_level;
  int *level_at_NCT;
  int **wordswithcarriers;
  int *nwordswithcarriers;
  //  double *maf_at_level; // How high is MAF at level
};

struct WINDOW *window;



struct MAP
{
  //  char chr[3];        // Chromosome
  string chr;
  int pos;            // physical position on chromosome
  string rs;           // ID
  /**/double rsNum;       // ggf. Probleme mit accuracy und Länge der rs-Nummern, da [[alpha]] auch durch Ziffern ersetzt wird -> ersetze durch long
  unsigned int dist;  // genetic distance
  unsigned int line;  // line in input file
  string minorA;
  string majorA;
  bool analysis_in;   // preselection for analysis
  bool matching_in;   // preselection for matching
  bool done;
  bool qcin;          // QC
  bool in;
  int missing;        // Missings insgesamt
  float missingrate;
  
  double maf;          // MAF for qc (case+control)
  double controlmaf;   // MAF for showing results (controls)
  double casemaf;      // MAF for showing results (cases)

  double mafa;          // MAF for qc (case+control) MAF_ADJUST
  double controlmafa;   // MAF for showing results (controls) MAF_ADJUST
  double casemafa;      // MAF for showing results (cases) MAF_ADJUST
  
  double mafr;          // MAF for qc (case+control) RARE
  double controlmafr;   // MAF for showing results (controls) RARE
  double casemafr;      // MAF for showing results (cases) RARE
  
  char *gene;
  int location;
  int locationToGene;
  int codingStatus;
  double p;
  double pmod;
};


