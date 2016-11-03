#ifndef _NB_EXPERT_HH_
#define _NB_EXPERT_HH_

#include "nb_param.hh" 
#include "nb_cpp_utils.h"

#include <stdio.h>
#include <iostream>
#include <vector>
class Expert
{
public:
  /// types of tasks
  enum ACTION {
    INCLDENSITY,   //0 -> INCLDENSITY
    TMIN,     //1 *
    TMAX,     //2 *
    MEDIAN,   //3 *
    LERROR,   //4 *
    RERROR,   //5 *
    MEAN,     //6 *
    TRIM,     //7 *
    QUANTILE, //8 *
    INVQUANT, //9 *
    INVQINCL, //10 *
    CONDDENSITY, //11 -> CONDDENSITY
    MAXLIKELI, //12 -> MAXLIKELI
    PLOT,     //13 *
    RNDCOND,  //14 ?
    RNDINCL,  //15 ?
    BINCLASS, //16 *
    REGR,      //17 ?
    NOOP = -1
  };

  Expert(const char *a,int debug = -2, bool writeout = false, ec_t** ec=NULL,
      log_func_t log_f = NULL,void* log_enclosed=NULL,
      delete_enclosed_func_t log_delete_enclosed=NULL);

  Expert(const std::string a,int debug = -2, bool writeout = false, ec_t** ec=NULL,
      log_func_t log_f = NULL,void* log_enclosed=NULL,
      delete_enclosed_func_t log_delete_enclosed=NULL);

  bool check_num_inputs(int num_inputs);

  Expert(float* myExpertise, int debug = -2, bool writeout = false, ec_t** ec=NULL,
      log_func_t log_f = NULL,void* log_enclosed=NULL,
      delete_enclosed_func_t log_delete_enclosed=NULL);
  virtual ~Expert();

  float nb_expert(ACTION ACT,double* X,
      float ARGUMENT=0.0);
  float nb_expert(ACTION ACT,float* X,float ARGUMENT=0.0);
  float nb_expert(double* X,ACTION ACT=Expert::BINCLASS,
      float ARGUMENT=0.0)
  {
    return nb_expert(ACT,X,ARGUMENT);
  }
  float nb_expert(float* X,ACTION ACT=Expert::BINCLASS,float ARGUMENT=0.0)
  {
    return nb_expert(ACT,X,ARGUMENT);
  }

  // Expert entries
  void     NB_EXPERT_FILLTABXS(float* STABXS);
  void     NB_EXPERT_DEFGSPLINE(int ModeIn, int NIn, float RegIn);
  void     NB_EXPERT_GETPINP(float* XPrePro);
  float    NB_EXPERT_FTMEAN(double (*f)(float*));
  static float NB_RNDM2(int& a, int dbg);


  // convert Expertise files (ASCII) to a file which can be included
  // in C++ code holding the Expertise as an array.
  void convertToArray(const char *ExpertiseFileName);
			  
  void SetHistFile(const char* thisName){HistoFileName = thisName;};
  void SetRootFile(const char* thisName){ SetHistFile(thisName); }

  float* getFtArray(){return plotArray;};

  void Print();
  char* getErrorMessage(){return errorMessage;};

  // Use these methods only after calling nb_expert with the same event, e.g.
  // with the action NOOP.
  float  NB_INVQUANT(float value);
  float  NB_QUANTILE(float value);
  float  NB_BINCLASS();
  float  NB_TMEAN(float value);
  float  TABX[NB_MAXNODE][NB_NVALUE];
  float  NB_INVQINCL(float argument);
  float  NB_RNDINCL(float argument);
  float  NB_CONDDENSITY(float& arg);
  void   NB_DERIVATIVE(float t, float* der);
  
  bool is_valid;
  int* cg_of_nb_cols;
  char* filename;
  int n_nb_cols;


  Expert* boost;
  Expert* next_expert;
  char* descriptor;
  int boost_sample_no;
  bool is_classify_density;
  bool tell_colindices_called;
  
  bool getWriteOutData();
  void setWriteOutData(bool value);
  int tellinputs_mode;
  float c_interface_v[NB_MAXDIM];
  common_t* com;         //common struct for logging and such  
  double NETOUT[NB_LEVMAX];
  
  int NODE1;
  int NODE2;
  int NODE3;
  int NODE4;


private:

  //
  // functions
  //

  // default constructor is private since
  // constructor needs some arguments
  Expert();
  void init_expert(const std::string a,int debug = -2,
      bool writeout = false, ec_t** ec=NULL,
      log_func_t log_f = NULL,void* log_enclosed=NULL,
      delete_enclosed_func_t log_delete_enclosed=NULL);
  void initialise_e();

  void   NB_DefPolynomials_e();

  void   NBBOOK1(int a,const char* b,int c,float d,float e,float f,int i);
  void   NBFILL_e(int a,float& b,float c,float& d);

  void   NB_DEFEXPERTISE_e();

  void   NB_DEFFT_e();
  void   NBPAK_e(int a, float* b);

  void   NB_PREPRO2_e();
  void   NB_CHOUTH(int& a,float* b,float* c,float* d,float* e,
		   int& f,float* g,int& h, int dbg);
  void   NB_CHOUTC(float* b,float* c,float* d,
		   int& e,float* f,int& g);
  void   NB_FORWE_e(int* a,float* b, double ( *f)(common_t**,float *),float* d,int& e,int& g);
  void   NB_SPLINEF2_e(int& histoID);
  void   NB_SPLINEF2_WHEN_NEEDED_e();
  void   NB_PERFPLOT2(int histoID);
  void   NB_FTXDEF(int& mmax);
  void   NB_FTXDEF_e(float* array,int nBins);
  void   NB_TRANSBACKTAIL(float& randomN,float& result);
  void   NB_TRANSBACK(float& eq,float& result);

  void   NB_READEXPERTISEC_e(const char* ExName, float* Ex);

  float  NB_TRANSGLE(float& a,float* aa);
  float  NB_FTMEAN(double (*f)(float*));
  float  NB_TOMARGINAL_e(float* X);
  float  NB_BSKFUN_e(float& t,float& Der);
			
  double NB_F3_e(float a);
  
  void PrintArrays(const char*);
  void InfoAndLicence();

  void writeErrorMessage(int varIndex, float limit, float current,
			 int status,int prepro);

  // computes and returns TABL (conditional PDF) for nBins 
  void getPDF(float* inputArray,int nBins,float* array);
  // checks the input array for INF, NAN and the range seen
  // in the training for each variable
  void checkInputRange_e(float* X);
  // checks if the input array is the same as the last one or not
  int isNewEvent(float* X);
  // initializes the arrays with the new event
  void calcNetOutput_e(float *X);

  void nb_prepare_boostdiag();
  void nb_prepro2_boost_e();
  void nb_diag_boostnet_e();
private:
  
  friend void fillPdfFunction(Expert*,float *x,int nBins,float* array);

  bool writeoutdata;
  //  int Debug;
  int LDEFSPLINE;

  float XSAVE[NB_MAXNODE];
  int NumPreproVar;      //variable selection, 0 if not used

  int NewEx;             //assume new expertise
  int MODEGS; 
  int NGS; 
  float REGGS;

  float T[NB_MAXNODE][NB_LEVMAX+1];
  float Weights[NB_MAXNODE][NB_MAXNODE][NB_MAXLAYER];
  int NODES[NB_MAXLAYER];

  int NLEVEL; //OUTPUT levels
  float OUTLEVEL[NB_LEVMAX];
  float XMEAN[NB_LEVMAX+1];
  float XSHAPE[NB_LEVMAX]; 

  float SigFrac;              //used for classification: Ratio Signal/BG

  // for training delta+density
  int nlLevel;
  float tDeltaVal[NB_MAXTDELTA];
  float tDeltaFrac[NB_MAXTDELTA];

  int LSHAPE,LLOG;
  int IFIXSHAPE,IPRUNE,IFIXORDER;
  int PREPROC;           //needed for PreProc2

  int AutoVarSelect;
  float TABG[NB_NVALUE];  
  float TABD[NB_NVALUE];  // g(s)
  float TABF[NB_NVALUE];  // f(t)
  float TABL[NB_NVALUE];
  float TABXS[NB_NVALUE+20];

  float RsfTable[NB_MAXNODE][NB_nRsfBins];
  int   nMapKey[NB_MAXNODE];
  float MapKeyValue[NB_MAXNODE][NB_MaxMapKey];  
  float MapKeyTrans[NB_MAXNODE][NB_MaxMapKey]; 

  float ITABY[NB_NIVALUE]; 
  float AA[(NB_MAXNODE-1)*(NB_MAXNODE-1)];
  float DIAG[NB_MAXNODE-1];
  float SCRATCH[NB_MAXNODE];

  float THETA[NB_MAXNODE*(NB_MAXNODE-1)/2];
  float CHEBY[NB_MAXNODE],CSHAPE[NB_LEVMAX];
  double CTH[NB_MAXNODE-2][NB_MAXNODE-2];
  double STH[NB_MAXNODE-2][NB_MAXNODE-2];

  static const int IVERSION=20070129;

  float NB_EXPERT;

  int eventCounter;   //number of current call to nb_expert

  int NLAYER;
  int MXNODE;
  int NEWEVT;           
  int NTABL;

  int ITER;              //number of iterations, from expertise
  
  //array holding number of variable sorted by significance
  int ISigSort[NB_MAXNODE-1];
  // inverse of ISigSort array: for an input variable with index i,
  // stores the corresponding index in the TABX and PreproFlag arrays
  int inverseISigSort[NB_MAXNODE-1];

  // number of input nodes (after pruning, signifcance cuts, etc.)
  int NVar;

  int PreproFlag[NB_MAXNODE]; //preprocessing flag for each variable

  int nMargDim;
  int MargVarid[NB_MaxMargDim];

  float PreproPar[NB_MaxPreproPar][NB_MAXNODE]; //store additional parameters for individual variable prepro 
  float A[NB_MAXNODE][NB_MAXLAYER];     
  float IN[1][NB_MAXDIM];        //input field for preproc2

  float XPRE[NB_MAXDIM];

  float TABX2[NB_MAXNODE][NB_NVALUE];
  float RsfTable2[NB_MAXNODE][NB_nRsfBins];

  // static definition is too big for CDF, instanttiate with new in the constructor
  float *EXPERTISE;

  std::string HistoFileName;      // ROOT file holding histos
  
  // array containing the PDF for an event 
  //(result of NB_EXPERT called with the option PLOT)
  float plotArray[100]; 

  float MargCoeff0;
  float MargCoeff[NB_MaxMargBins][NB_MaxMargDim];

  char errorMessage[500];
  // knot positions from NB_DEFFT
  float TT[33];
  // number of knots from NB_DEFFT
  int NKNOT;
  // spline coefficients from NB_DEFFT
  float CP[100];
  //SPLINEF2 performed for current event 
  bool spline_performed;
  //variables for boost training
  float TABX3[NB_NVALUE][NB_MAXNODE];
  float RsfTable3[NB_nRsfBins][NB_MAXNODE];     
  float a_boost[(NB_MAXNODE-1)*(NB_MAXNODE-1)];
  float diag_boost[NB_MAXNODE-1];
  float inv_sigma_boostvar[NB_MAXNODE];               
  bool is_boost;        
  bool is_net_plus_diag;
  bool is_preboost;
  int ioffset;
  int NumPreproVar_save;
  int preproc_save;
  int node1_save;
  float mean_boostvar[NB_MAXNODE];
  int   ipos_add_vars[NB_MAXNODE];
  int iversiont;
  int iterations;
  int final_diagfit;
  std::stringstream ss;
};

#endif // Expert.hh
