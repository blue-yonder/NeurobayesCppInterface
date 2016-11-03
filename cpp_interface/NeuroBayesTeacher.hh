/*
 * NeuroBayes Teacher
 *
 * 
 * Revision history:
 * 2004-June-05 Ulrich Kerzel
 *    initial revision
 */


#ifndef NEUROBAYESTEACHER_HH_
#define NEUROBAYESTEACHER_HH_

#include "nb_param.hh"
#include "nb_cpp_utils.h"

#include <string> 
#include <vector>

class NeuroBayesTeacher{

public:
  // Get the singleton instance.
  static NeuroBayesTeacher* Instance(ec_t** ec=NULL, int debug=-1, log_func_t log_f=NULL,
      void* log_enclosed=NULL, delete_enclosed_func_t log_delete_enclosed=NULL);

  // NeuroBayes steering functions

  void NB_DEF (bool resetInput=true);
  void NB_DEF_TASK (const char* thisTask);
  void NB_DEF_TASK (std::string & thisTask);
  void NB_DEF_DEBUG(int   thisDebug);
  void NB_DEF_PRE  (int   thisPre);
  void NB_DEF_INITIALPRUNE (int   thisIprune);
  void NB_DEF_NODE1(int   thisNode1);
  void NB_DEF_NODE2(int   thisNode2);
  void NB_DEF_NODE3(int   thisNode3);
  void NB_DEF_REG  (const char* thisReg);
  void NB_DEF_LEARNDIAG    (int   thisValue);
  void NB_DEF_LOSS  (const char* thisLoss);
  void NB_DEF_SHAPE (const char* thisShape);
  void NB_DEF_METHOD(const char* thisMethod);
  void NB_DEF_MOM   (float thisMom);
  void NB_DEF_EPOCH (int   thisEpoch);
  void NB_DEF_ITER  (int   thisIter);
  void NB_DEF_RTRAIN(float thisRtrain);
  void NB_DEF_SPEED (float thisSpeed);
  void NB_DEF_MAXLEARN (float thisMaxlearn);
  void NB_DEF_RELIMPORTANCE(float thisRelimportance);
  void NB_DEF_SURRO (float thisSurro);
  void NB_DEF_PRUNEMIN (float thisPrunemin);
  void NB_DEF_PRUNEMAX (float thisPrunemax);
  void NB_DEF_PRUNERESULT (float thisPruneresult);
  void NB_DEF_QUANTILE (float thisQuantile);
  // parameter thisDbg in NB_RANVIN(thisJseed, thisJwarm, thisDbg) is deprecated
  // and will be ignored
  void NB_RANVIN (int   thisJseed = 4711, 
		  int   thisJwarm = 10, 
		  int   thisDbg   = -2);  
  void NB_DEF_LOSSWGT (float thisWeight);
  void NB_DEF_TDELTA (float thisWeight);
  void NB_DEF_WEIGHT_MODE(int mode);
  void NB_DEF_SPLOT_MODE(int mode);

  // has to be static: can be called without a teacher-instance is created
  // needed by useTargetDistribution before callTeacher is done
  static void NB_TABDEF1(float* targetDist, float* targetWeight,
      int numTarget, float* targetTab, int numTab,
      common_t* com1=NULL);

  void SetIndividualPreproFlag (int thisIvar, int thisFlag,const char* varname = "");
  void SetIndividualPreproParameter(int thisIvar, int thisParNr, 
				    float thisValue);

  // define training target, weight:
  void SetTarget(float thisTarget);
  void SetWeight(float thisWeight,float thisWeight2 = 1);
  void SetNextInput(int thisNvar,float* thisVars);

  // hint for efficient memory allocation
  void SetNEvents(int nevt);

  // start network training
  void TrainNet(bool write_output_files=true);

  // misc. setups
  void SetOutputFile(const char* thisName);
  void SetHistosFile(const char* thisName);
  void SetCArrayFile(const char* thisName);

  float* nb_get_expertise();

  //create correl_signi.txt file
  void nb_infoout(float* weightsum,float* total,int* keep,int* rank,
      float* single,float* added,float* global,float* loss,
      int* nvar,int* index);
  void nb_correl_signi(const char filename_txt[],
      const char filename_html[]);
  void nb_correl_signi(char** varnames,const char filename_txt[],
      const char filename_html[]);

  char** nb_get_varnames(int* n_var_all);
  int* nb_get_individual_prepro_flags(int* n_var_all);

  // just to check if things work
  void SayHello();

  // Destructor private because it is a singleton
  virtual ~NeuroBayesTeacher();
  //
  
  // Logging and error handling
  common_t* com;
  
  // singleton
  static NeuroBayesTeacher* instance;
  
private:

  void NB_DEF_WEIGHT_FACTOR();
  
  // status
  static unsigned int instanceCounter;

  // Constructor private because it is a singleton
  NeuroBayesTeacher(ec_t** ec=NULL, int debug=-1, log_func_t log_f=NULL,
      void* log_enclosed=NULL, delete_enclosed_func_t log_delete_enclosed=NULL);

  // NeuroBayes variables
  std::vector<std::string> varnames;
  std::vector<int> prepros;
  int maxEvent;
  int storedEvents;
  int weight_mode;
  float trainingTarget1;
  float trainingTarget2;
  float eventWeight;
  float eventWeight2;
  double wsum;
  double w2sum;
  
  //net input/output arrays 
  std::vector<float> inarray;
  
  //NeuroBayes Expertise
  float Expertise[NB_NEXPERTISE];

  std::string ExpertiseFileName;
  std::string HistosFileName;      // file holding histos, etc
  std::string CArrayFileName;

  bool writeCArray;
  std::stringstream ss;

}; //class NeuroBayesTeacher

#endif
