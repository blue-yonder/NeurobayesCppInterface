#include <cstring>
#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include "NeuroBayesTeacher.hh"

// declare nb's FORTRAN functions
extern "C"
{
  void nb_def_ (common_t ** com);
  void nb_def_task_e_ (common_t ** com, const char *, int);
  void nb_def_node1_e_ (common_t ** com, const int &);
  void nb_def_node2_e_ (common_t ** com, const int &);
  void nb_def_node3_e_ (common_t ** com, const int &);
  void nb_def_debug_ (common_t ** com, const int &);
  void nb_def_learndiag_e_ (common_t ** com, const int &);

  void nb_def_reg_e_ (common_t ** com, const char *, int);
  void nb_def_pre_e_ (common_t ** com, const int &);
  void nb_def_initialprune_e_ (common_t ** com, const int &);
  void nb_def_loss_e_ (common_t ** com, const char *, int);
  void nb_def_shape_e_ (common_t ** com, const char *, int);
  void nb_def_method_e_ (common_t ** com, const char *, int);
  void nb_def_mom_e_ (common_t ** com, const float &);
  void nb_def_epoch_ (common_t ** com, const int &);
  void nb_def_rtrain_ (common_t ** com, const float &);
  void nb_def_iter_ (common_t ** com, const int &);
  void nb_def_maxlearn_ (common_t ** com, const float &);
  void nb_def_speed_ (common_t ** com, const float &);
  void nb_def_relimportance_ (common_t ** com, const float &);
  void nb_def_surro_ (common_t ** com, const float &);
  void nb_def_prunemin_ (common_t ** com, const float &);
  void nb_def_prunemax_ (common_t ** com, const float &);
  void nb_def_pruneresult_ (common_t ** com, const float &);
  void nb_def_losswgt_ (common_t **, const float &);
  void nb_def_tdelta_ (common_t **, const float &);
  void nb_tabdef1_e_ (common_t **, float *, float *, int *, float *, int *,
		      int *, double *);
  void nb_def_weight_mode_ (common_t **, const int &);
  void nb_def_splot_mode_ (common_t **, const int &);
  void nb_def_weight_factor_ (common_t **, const double &);
  void nb_def_preproflag_ (common_t **, const int &node, const int &flag);
  void nb_def_prepropar_ (common_t **, const int &node, const int &par,
			  const float &value);

  void nb_ranvin_ (common_t ** com, int *jseed, int *jwarm);

  void nb_prepro_and_network_e_ (common_t ** com, const int &, float *,
				 float[]);

  void nb_saveexpertise_ (common_t **, const char *, float[], int);
  void nb_saveascarray_ (common_t **, const char *, float[], int);

  //fill correl_signi arrays
  void nb_infoout_e_ (common_t ** com, float *weightsum, float *total,
		      int *keep, int *rank, float *single, float *added,
		      float *global, float *loss, int *nvar, int *index);

  // functions from histogram interface
  void nb_histo_init_ ();
  void nb_histo_save_ (common_t ** com, const char *, int);

}				//extern C

// status
unsigned int
  NeuroBayesTeacher::instanceCounter = 0;
NeuroBayesTeacher *
  NeuroBayesTeacher::instance = 0;

void
_NeuroBayesTeacher_destructor () __attribute__ ((destructor));
void
_NeuroBayesTeacher_destructor ()
{
  if (NeuroBayesTeacher::instance != 0)
    {
      delete
	NeuroBayesTeacher::instance;
      NeuroBayesTeacher::instance = 0;
    }
}


// Constructors
NeuroBayesTeacher::NeuroBayesTeacher (ec_t ** ec, int debug, log_func_t log_f,
				      void *log_enclosed,
				      delete_enclosed_func_t
				      log_delete_enclosed)
{
  // life counter
  instanceCounter++;

  // initialisations
  //
  // 1. perform a licence check and set the debug flag (quiet)
  // 2.

  com =
    nb_init_common (debug, ec, rand_double1, NULL,
		    1);
  if (nb_get_debug (&com) > -3 && log_f)
    {
      nb_register_logging (&com, log_f, log_enclosed, log_delete_enclosed);
    }

  if (nb_cpp_handle_error (com))
    return;

  // by default, do not write out C++ array file
  writeCArray = false;
  // determine max. number of events
  // 5 for safety margin
  int maxint = 2147483647;
  maxEvent = maxint - 1;	//NB_MAXPATTERN-NB_MaxPreproPar-5

  // reset 
  storedEvents = 0;
  trainingTarget1 = 0;
  trainingTarget2 = 0;
  eventWeight = 1.0;
  eventWeight2 = 1.0;
  wsum = 0;
  w2sum = 0;
  weight_mode = 0;
  inarray.clear ();
  nb_def_ (&com);
}				// constructor

NeuroBayesTeacher *
NeuroBayesTeacher::Instance (ec_t ** ec, int debug, log_func_t log_f,
			     void *log_enclosed,
			     delete_enclosed_func_t log_delete_enclosed)
{
  if (!instance)
    {
      // KCC refused to compile: instance = new NeuroBayesTeacher::NeuroBayesTeacher();
      instance = new NeuroBayesTeacher (ec, debug, log_f,
					log_enclosed, log_delete_enclosed);
    }
  else
    {
      nb_delete_common (&(instance->com));
      instance->com =
	nb_init_common (debug, ec, rand_double1,
			NULL, 1);
      if (nb_get_debug (&(instance->com)) > -3 && log_f)
	{
	  nb_register_logging (&(instance->com), log_f, log_enclosed,
			       log_delete_enclosed);
	}
    }
  return instance;
}				// Instance()

// Destructors
NeuroBayesTeacher::~NeuroBayesTeacher ()
{
  instanceCounter--;
  nb_delete_common (&com);
}				// destructor

// general NeuroBayes settings
void
NeuroBayesTeacher::NB_DEF (bool resetInput)
{
  if (nb_get_debug (&com) >= -1)
    {
      ss << "*** reset NeuroBayes Teacher ***" << std::endl;
      nb_cpp_log (ss, -1, &com);
    }

  // reset variables
  storedEvents = 0;
  trainingTarget1 = 0;
  trainingTarget2 = 0;
  eventWeight = 1.0;
  eventWeight2 = 1.0;
  wsum = 0;
  w2sum = 0;
  weight_mode = 0;
  varnames.clear ();
  prepros.clear ();

  if (resetInput)
    {
      if (nb_get_debug (&com) >= -1)
	{
	  ss << "NeuroBayesTeacher: reset input values ...";
	  nb_cpp_log (ss, -1, &com);
	}
      inarray.clear ();
      if (nb_get_debug (&com) >= -1)
	{
	  ss << "... done " << std::endl;
	  nb_cpp_log (ss, -1, &com);
	}
    }
  if (nb_get_debug (&com) >= -1)
    {
      ss << "NeuroBayesTeacher: reset Expertise ...";
      nb_cpp_log (ss, -1, &com);
    }
  memset (Expertise, 0, NB_NEXPERTISE * sizeof (float));
  if (nb_get_debug (&com) >= -1)
    {
      ss << "... done " << std::endl;
      nb_cpp_log (ss, -1, &com);
      ss << "NeuroBayesTeacher: reset internal variables...";
      nb_cpp_log (ss, -1, &com);
    }
  nb_def_ (&com);
  if (nb_get_debug (&com) >= -1)
    {
      ss << "... done " << std::endl;
      nb_cpp_log (ss, -1, &com);
    }
}

void
NeuroBayesTeacher::NB_DEF_TASK (const char *thisTask)
{
  std::string myTask = thisTask;
  NB_DEF_TASK (myTask);
}

void
NeuroBayesTeacher::NB_DEF_TASK (std::string & myTask)
{
  if (myTask.size () < 4)
    myTask += "    ";
  nb_def_task_e_ (&com, myTask.c_str (), myTask.size ());
  nb_cpp_handle_error (com);
}

void
NeuroBayesTeacher::NB_DEF_DEBUG (int thisDebug)
{
  nb_def_debug_ (&com, thisDebug);
}

void
NeuroBayesTeacher::NB_DEF_LEARNDIAG (int thisValue)
{
  nb_def_learndiag_e_ (&com, thisValue);
  nb_cpp_handle_error (com);
}

void
NeuroBayesTeacher::NB_DEF_PRE (int thisPre)
{
  nb_def_pre_e_ (&com, thisPre);
  nb_cpp_handle_error (com);
}

void
NeuroBayesTeacher::NB_DEF_INITIALPRUNE (int thisIprune)
{
  nb_def_initialprune_e_ (&com, thisIprune);
  nb_cpp_handle_error (com);
}

void
NeuroBayesTeacher::NB_DEF_NODE1 (int thisNode1)
{
  nb_def_node1_e_ (&com, thisNode1);
  nb_cpp_handle_error (com);
}

void
NeuroBayesTeacher::NB_DEF_NODE2 (int thisNode2)
{
  nb_def_node2_e_ (&com, thisNode2);
  nb_cpp_handle_error (com);
}

void
NeuroBayesTeacher::NB_DEF_NODE3 (int thisNode3)
{
  nb_def_node3_e_ (&com, thisNode3);
  nb_cpp_handle_error (com);
}

void
NeuroBayesTeacher::NB_DEF_REG (const char *thisReg)
{
  std::string myReg = thisReg;
  if (myReg.size () < 4)
    myReg += "    ";
  nb_def_reg_e_ (&com, myReg.c_str (), myReg.size ());
  nb_cpp_handle_error (com);
}

void
NeuroBayesTeacher::NB_DEF_LOSSWGT (float thisWeight)
{
  nb_def_losswgt_ (&com, thisWeight);
}

void
NeuroBayesTeacher::NB_DEF_LOSS (const char *thisLoss)
{
  std::string myLoss = thisLoss;
  if (myLoss.size () < 4)
    myLoss += "    ";
  nb_def_loss_e_ (&com, myLoss.c_str (), myLoss.size ());
  nb_cpp_handle_error (com);
}

void
NeuroBayesTeacher::NB_DEF_SHAPE (const char *thisShape)
{
  std::string myShape = thisShape;
  // hang some spaces to avoid crashes
  if (myShape.size () < 4)
    myShape += "    ";
  nb_def_shape_e_ (&com, myShape.c_str (), myShape.size ());
  nb_cpp_handle_error (com);
}

void
NeuroBayesTeacher::NB_DEF_METHOD (const char *thisMethod)
{
  std::string myMethod = thisMethod;
  // hang some spaces to avoid crashes
  if (myMethod.size () < 4)
    myMethod += "    ";
  nb_def_method_e_ (&com, myMethod.c_str (), myMethod.size ());
  nb_cpp_handle_error (com);
}

void
NeuroBayesTeacher::NB_DEF_MOM (float thisMom)
{
  nb_def_mom_e_ (&com, thisMom);
  nb_cpp_handle_error (com);
}

void
NeuroBayesTeacher::NB_DEF_EPOCH (int thisEpoch)
{
  nb_def_epoch_ (&com, thisEpoch);
}

void
NeuroBayesTeacher::NB_DEF_ITER (int thisIter)
{
  nb_def_iter_ (&com, thisIter);
}

void
NeuroBayesTeacher::NB_DEF_RTRAIN (float thisRtrain)
{
  nb_def_rtrain_ (&com, thisRtrain);
}

void
NeuroBayesTeacher::NB_DEF_SPEED (float thisSpeed)
{
  nb_def_speed_ (&com, thisSpeed);
}

void
NeuroBayesTeacher::NB_DEF_MAXLEARN (float thisMaxlearn)
{
  nb_def_maxlearn_ (&com, thisMaxlearn);
}

void
NeuroBayesTeacher::NB_DEF_RELIMPORTANCE (float thisRelimportance)
{
  nb_def_relimportance_ (&com, thisRelimportance);
}

void
NeuroBayesTeacher::NB_DEF_SURRO (float thisSurro)
{
  nb_def_surro_ (&com, thisSurro);
}

void
NeuroBayesTeacher::NB_DEF_PRUNEMIN (float thisPrunemin)
{
  nb_def_prunemin_ (&com, thisPrunemin);
}

void
NeuroBayesTeacher::NB_DEF_PRUNEMAX (float thisPrunemax)
{
  nb_def_prunemax_ (&com, thisPrunemax);
}

void
NeuroBayesTeacher::NB_DEF_PRUNERESULT (float thisPruneresult)
{
  nb_def_pruneresult_ (&com, thisPruneresult);
}

void
NeuroBayesTeacher::NB_DEF_QUANTILE (float thisQuantile)
{
  nb_def_pruneresult_ (&com, thisQuantile);
}

void
NeuroBayesTeacher::NB_DEF_TDELTA (float delta)
{
  nb_def_tdelta_ (&com, delta);
}

void
NeuroBayesTeacher::NB_DEF_WEIGHT_MODE (int mode)
{
  weight_mode = mode;
  nb_def_weight_mode_ (&com, mode);
}

void
NeuroBayesTeacher::NB_DEF_SPLOT_MODE (int mode)
{
  nb_def_splot_mode_ (&com, mode);
}

void
NeuroBayesTeacher::NB_DEF_WEIGHT_FACTOR ()
{
  if (weight_mode == 2 && w2sum != 0)
    {
      double weight_factor = wsum / w2sum;
      nb_def_weight_factor_ (&com, weight_factor);
    }
}

// wrapper for fortran func tabdef1
// field 'targetDist' of dim 'numTarget' divided into 'numTa' quantiles and 
// the 'targetTab' filed  contains the borders of the qauntiles
void
NeuroBayesTeacher::NB_TABDEF1 (float *targetDist, float *targetWeight,
			       int numTarget, float *targetTab, int numTab,
			       common_t * com1)
{
  const int size = numTarget;
  int tmpArray[size];
  double tmpWSum[(size + 1)];
  nb_tabdef1_e_ (&com1, &targetDist[0], &targetWeight[0],
		 &numTarget, &targetTab[0], &numTab, &tmpArray[0],
		 &tmpWSum[0]);
  nb_cpp_handle_error (com1);
}

void
NeuroBayesTeacher::NB_RANVIN (int thisJseed, int thisJwarm, int thisDbg)
{
  init_random_number_generator (thisJseed);
}				//NB_RANVIN

void
NeuroBayesTeacher::SetIndividualPreproFlag (int thisIvar, int thisFlag,
					    const char *varname)
{
  varnames.push_back (varname);
  prepros.push_back (thisFlag);
  // pass-through
  nb_def_preproflag_ (&com, thisIvar + 2, thisFlag);
}				// SetIndividualPreproFlag

void
NeuroBayesTeacher::SetIndividualPreproParameter (int thisIvar, int thisParNr,
						 float thisValue)
{

  // FORTRAN counts from 1, and first node is bias node
  int node = thisIvar + 2;

  // pass-through
  nb_def_prepropar_ (&com, node, thisParNr + 1, thisValue);
}				// SetIndividualPreproParameter

// misc. settings
void
NeuroBayesTeacher::SetTarget (float thisTarget)
{
  if (nb_cpp_handle_error (com))
    return;
  if (isnan (thisTarget) || isinf (thisTarget))
    {
      ss << "NeuroBayesTeacher::SetTarget  ERROR: "
	<< "NAN or INF passed as target value" << std::endl;
      char *msg = stream_c_str (ss);
      if (nb_get_debug (&com) >= -2)
	nb_c_log (msg, -2, &com);
      nb_set_error (&com, invalid_data_exc, msg);
      free (msg);
      nb_cpp_handle_error (com);
      return;
    }
  trainingTarget1 = thisTarget;
}				// SetTarget

void
NeuroBayesTeacher::SetWeight (float thisWeight, float thisWeight2)
{
  if (nb_cpp_handle_error (com))
    return;
  eventWeight = thisWeight;
  eventWeight2 = thisWeight2;
  if (isnan (eventWeight))
    {
      ss << "NeuroBayesTeacher::SetWeight  ERROR: "
	<< "NAN passed as weight1 in event " << storedEvents << std::endl;
      ss << "\n Please, correct the error" << std::endl;
      char *msg = stream_c_str (ss);
      if (nb_get_debug (&com) >= -2)
	nb_c_log (msg, -2, &com);
      nb_set_error (&com, invalid_data_exc, msg);
      free (msg);
      nb_cpp_handle_error (com);
      return;
    }

  if (isnan (eventWeight2))
    {
      ss << "NeuroBayesTeacher::SetWeight  ERROR: "
	<< "NAN passed as weight2 in event " << storedEvents << std::endl;
      ss << "\n Please, correct the error" << std::endl;
      char *msg = stream_c_str (ss);
      if (nb_get_debug (&com) >= -2)
	nb_c_log (msg, -2, &com);
      nb_set_error (&com, invalid_data_exc, msg);
      free (msg);
      nb_cpp_handle_error (com);
      return;
    }
  wsum += eventWeight2 * eventWeight;
  w2sum += eventWeight * eventWeight * eventWeight2;
}				//SetWeight

void
NeuroBayesTeacher::SetNextInput (int numVariables, float *thisVars)
{
  if (nb_cpp_handle_error (com))
    return;
  if (storedEvents < maxEvent)
    {

      // input variable sanity checks
      for (int i = 0; i != numVariables; ++i)
	{
	  if (isnan (thisVars[i]))
	    {
	      ss << "NeuroBayesTeacher::SetNextInput  ERROR: "
		<< "NAN passed as input in event " << storedEvents
		<< ". column: " << i << std::endl;
	      ss << "\n Please, correct the error" << std::endl;
	      char *msg = stream_c_str (ss);
	      if (nb_get_debug (&com) >= -2)
		nb_c_log (msg, -2, &com);
	      nb_set_error (&com, invalid_data_exc, msg);
	      free (msg);
	      nb_cpp_handle_error (com);
	      return;
	    }

	  if (isinf (thisVars[i]))
	    {
	      ss << "NeuroBayesTeacher::SetNextInput  ERROR: "
		<< "INF passed as input in event " << storedEvents
		<< ". column: " << i << std::endl;
	      ss << "\n Please, correct the error" << std::endl;
	      char *msg = stream_c_str (ss);
	      if (nb_get_debug (&com) >= -2)
		nb_c_log (msg, -2, &com);
	      nb_set_error (&com, invalid_data_exc, msg);
	      free (msg);
	      nb_cpp_handle_error (com);
	      return;
	    }
	}

      /*
         inarray index structure for one event:

         0    : target
         1    : first input variable
         2    : second input variable
         .
         .
         .
         numVariables    : last input variable
         numVariables + 1: empty space
         .
         .
         .
         NB_MAXNODE - 1: empty space
         NB_MAXNODE + 0: event weight
         NB_MAXNODE + 1: target 1
         NB_MAXNODE + 2: target 2
         NB_MAXNODE + 3: empty space for internal boost
         .
         .
         .
         NB_MAXDIM - 1: empty space for internal boost

         from nb_param.hh:
         NB_MAXDIM = NB_MAXNODE + 3 + NB_MAXNODE

       */

      // target
      inarray.push_back (trainingTarget1);

      // input variables
      inarray.insert (inarray.end (), &thisVars[0], &thisVars[numVariables]);

      // empty space
      inarray.insert (inarray.end (), NB_MAXNODE - numVariables - 1, 0.);

      // event weight and targets
      inarray.push_back (eventWeight * eventWeight2);
      inarray.push_back (trainingTarget1);
      inarray.push_back (trainingTarget2);

      // empty space
      inarray.insert (inarray.end (), NB_MAXNODE, 0.);

      // increase event counter
      storedEvents++;

    }
  else
    {
      ss << "NeuroBayesTeacher::SetNextInput "
	<< "Number of events too high for your version. Abort " << std::endl;
      char *msg = stream_c_str (ss);
      if (nb_get_debug (&com) >= -2)
	nb_c_log (msg, -2, &com);
      nb_set_error (&com, invalid_data_exc, msg);
      free (msg);
      nb_cpp_handle_error (com);
      return;
    }

  return;
}				// SetNextInput

void
NeuroBayesTeacher::SetNEvents (int nevt)
{
  // This method tells the inarray to allocate enough space for the given
  // number of events. Until the given number is reached, no reallocation
  // is necessary when filling input events, speeding up the process and
  // avoiding copying large chunks of memory around.

  // check upper limit
  if (nevt > maxEvent)
    {
      nevt = maxEvent;

      if (nb_get_debug (&com) >= -2)
	{
	  // tell the user
	  ss << "NeuroBayesTeacher::SetNEvents warning" << std::endl;
	  ss << "  You requested space for " << nevt << "events," << std::
	    endl;
	  ss << "  but " << maxEvent << " is still the hard limit" << std::
	    endl;
	  ss << "  for your version." << std::endl;
	  ss << "  I'll allocate memory for " << nevt << " events." << std::
	    endl;
	  nb_cpp_log (ss, -2, &com);
	}
    }

  inarray.reserve (nevt * NB_MAXDIM);

}				// SetNEvents

void
NeuroBayesTeacher::SetOutputFile (const char *thisName)
{
  ExpertiseFileName = thisName;
  return;
}				// SetOutputFile

void
NeuroBayesTeacher::SetHistosFile (const char *thisName)
{
  HistosFileName = thisName;
  return;
}				//SetHistosFilename

void
NeuroBayesTeacher::SetCArrayFile (const char *thisName)
{

  writeCArray = true;
  CArrayFileName = thisName;

}				// SetCArrayFile

// train network
void
NeuroBayesTeacher::TrainNet (bool write_output_files)
{
  if (nb_cpp_handle_error (com))
    return;
  NB_DEF_WEIGHT_FACTOR ();
  nb_histo_init_ ();
  nb_prepro_and_network_e_ (&com, storedEvents, &inarray[0], Expertise);
  if (nb_cpp_handle_error (com))
    return;
  if (write_output_files)
    {
      nb_saveexpertise_ (&com, ExpertiseFileName.c_str (),
			 Expertise, ExpertiseFileName.size ());
      // save Expertise as C++ array?
      if (writeCArray)
	{
	  nb_saveascarray_ (&com, CArrayFileName.c_str (),
			    Expertise, CArrayFileName.size ());
	}			//save as C++ array
    }
  if (nb_get_debug (&com) > quiet_dbg)
    {
      if ((int) HistosFileName.size () < 1)
	{
	  HistosFileName = "ahist.txt";
	}
      nb_histo_save_ (&com, HistosFileName.c_str (),
		      (int) HistosFileName.size ());
    }
}				//TrainNet

void
NeuroBayesTeacher::SayHello ()
{
  if (nb_get_debug (&com) >= -1)
    {
      ss << "NeuroBayes Teacher(R) " << std::endl;
      nb_cpp_log (ss, -1, &com);
    }
}

float *
NeuroBayesTeacher::nb_get_expertise ()
{
  return &Expertise[0];
}

void
NeuroBayesTeacher::nb_infoout (float *weightsum, float *total, int *keep,
			       int *rank, float *single, float *added,
			       float *global, float *loss, int *nvar,
			       int *index)
{
  nb_infoout_e_ (&com, weightsum, total, keep, rank, single, added, global,
		 loss, nvar, index);
  nb_cpp_handle_error (com);
  return;
}

char **
NeuroBayesTeacher::nb_get_varnames (int *n_var_all)
{
  int size = varnames.size ();
  *(n_var_all) = size;
  if (size == 0)
    return NULL;
  else
    {
      char **varnames_c = (char **) malloc (size * sizeof (char *));
      for (int i = 0; i < size; i++)
	varnames_c[i] = strdup (varnames[i].c_str ());
      return varnames_c;
    }
}

int *
NeuroBayesTeacher::nb_get_individual_prepro_flags (int *n_var_all)
{
  int size = prepros.size ();
  *(n_var_all) = size;
  if (size == 0)
    return NULL;
  else
    {
      int *prepros_c = (int *) malloc (size * sizeof (int));
      for (int i = 0; i < size; i++)
	prepros_c[i] = prepros[i];
      return prepros_c;
    }
}

void
NeuroBayesTeacher::nb_correl_signi (char **varnames_,
				    const char filename_txt[],
				    const char filename_html[])
{
  if (nb_get_debug (&com) >= -2)
    {
      nb_c_log ((char *) ("*** Deprecated warning\n"
			  "Use NeuroBayesTeacher::nb_correl_signi("
			  "const char filename_txt[],const char filename_html[])"
			  "and set the varnames directly with: \nNeuroBayesTeacher::"
			  "SetIndividualPreproFlag("
			  "int thisIvar, int thisFlag,const char* varname)\n"),
		-2, &com);
    }
  float weightsum, total, single[NB_MAXNODE], global[NB_MAXNODE],
    added[NB_MAXNODE], loss[NB_MAXNODE];
  int keep, nvar;
  int j;
  int rank[NB_MAXNODE], index[NB_MAXNODE];

  nb_infoout (&weightsum, &total, &keep, rank, single, added, global, loss,
	      &nvar, index);
  if (nb_cpp_handle_error (com))
    return;
  float signi = sqrt (weightsum);

  //now open file and write out arrays
  FILE *file = fopen (filename_txt, "w");

  fprintf (file, "\n total correlation to target:  %3.3f%%\n", total * 100.);
  fprintf (file, " total significance:  %3.3f sigma\n", total * signi);
  fprintf (file,
	   " (additional signif. , only this var , loss when removed , global corr. to others)\n\n");

  for (int i = 0; i < nvar; i++)
    {
      j = index[i] - 1;
      fprintf (file, "%3d%s%4d:%6.2f %6.2f %6.2f %5.1f%%   %s\n", i + 1,
	       ".:  variable", j + 2, fabs (added[j] * signi),
	       fabs (single[j] * signi), fabs (loss[j] * signi),
	       100. * global[j], varnames_[j]);
    }
  fprintf (file, "\n Keep only %d most significant input variables\n", keep);
  fclose (file);
  //Now write out html file
  FILE *html_file = fopen (filename_html, "w");
  //head of html file
  fprintf (html_file,
	   "<html>\n<head>\n<title>%s</title>\n</head>\n<body>\n\n<h1>%s</h1>\n<table border=\"1\">\n",
	   filename_html, filename_html);
  //headline table
  fprintf (html_file,
	   "<tr>\n<th>nrank</th>\n<th>nvar</th>\n<th>additional signif</th>\n<th>only this var</th>\n<th>loss when removed</th>\n<th>global corr. to others [%%]</th>\n");
  fprintf (html_file, "<h3>total correlation to target:  %3.3f%%</h3>\n",
	   total * 100.);
  fprintf (html_file, "<h3>total significance:  %3.3f sigma</h3>\n",
	   total * signi);
  fprintf (html_file,
	   "<h6>(use firefox AddOn TableTools to sort the columns)</h6>\n");
  char tr_options[50];
  for (int i = 0; i < nvar; i++)
    {
      j = index[i] - 1;
      if (i + 1 > keep)
	strcpy (tr_options, "tr bgcolor=\"#C0C0C0\"");
      else
	strcpy (tr_options, "tr");
      fprintf (html_file,
	       "<%s>\n<td>%3d</td><td>%4d</td><td>%6.2f</td><td>%6.2f</td><td>%6.2f</td><td>%5.1f</td><td>%s</td></tr>\n",
	       tr_options, i + 1, j + 2, fabs (added[j] * signi),
	       fabs (single[j] * signi), fabs (loss[j] * signi),
	       100. * global[j], varnames_[j]);
    }

  fprintf (html_file, "</table>\n</body>\n</html>\n");
  fclose (html_file);
}

void
NeuroBayesTeacher::nb_correl_signi (const char filename_txt[],
				    const char filename_html[])
{
  float weightsum, total, single[NB_MAXNODE], global[NB_MAXNODE],
    added[NB_MAXNODE], loss[NB_MAXNODE];
  int keep, nvar;
  int j;
  int rank[NB_MAXNODE], index[NB_MAXNODE];

  nb_infoout (&weightsum, &total, &keep, rank, single, added, global, loss,
	      &nvar, index);
  if (nb_cpp_handle_error (com))
    return;

  float signi = sqrt (weightsum);

  //now open file and write out arrays
  FILE *file = fopen (filename_txt, "w");

  fprintf (file, "\n total correlation to target:  %3.3f%%\n", total * 100.);
  fprintf (file, " total significance:  %3.3f sigma\n", total * signi);
  fprintf (file,
	   " (additional signif. , only this var , loss when removed , global corr. to others)\n\n");
  for (int i = 0; i < nvar; i++)
    {
      j = index[i] - 1;
      fprintf (file, "%3d%s%4d:%6.2f %6.2f %6.2f %5.1f%%   %s %d #%d\n",
	       i + 1, ".:  variable", j + 2, fabs (added[j] * signi),
	       fabs (single[j] * signi), fabs (loss[j] * signi),
	       100. * global[j], varnames[j].c_str (), prepros[j], j + 2);
    }
  fprintf (file, "\n Keep only %d most significant input variables\n", keep);
  fclose (file);
  //Now write out html file
  FILE *html_file = fopen (filename_html, "w");
  //head of html file
  fprintf (html_file,
	   "<html>\n<head>\n<title>%s</title>\n</head>\n<body>\n\n<h1>%s</h1>\n<table border=\"1\">\n",
	   filename_html, filename_html);
  //headline table
  fprintf (html_file,
	   "<tr>\n<th>nrank</th>\n<th>nvar</th>\n<th>additional signif</th>\n<th>only this var</th>\n<th>loss when removed</th>\n<th>global corr. to others [%%]</th>\n");
  fprintf (html_file, "<h3>total correlation to target:  %3.3f%%</h3>\n",
	   total * 100.);
  fprintf (html_file, "<h3>total significance:  %3.3f sigma</h3>\n",
	   total * signi);
  fprintf (html_file,
	   "<h6>(use firefox AddOn TableTools to sort the columns)</h6>\n");
  char tr_options[50];
  for (int i = 0; i < nvar; i++)
    {
      j = index[i] - 1;
      if (i + 1 > keep)
	strcpy (tr_options, "tr bgcolor=\"#C0C0C0\"");
      else
	strcpy (tr_options, "tr");
      fprintf (html_file,
	       "<%s>\n<td>%3d</td><td>%4d</td><td>%6.2f</td><td>%6.2f</td><td>%6.2f</td><td>%5.1f</td><td>%s %d #%d</td></tr>\n",
	       tr_options, i + 1, j + 2, fabs (added[j] * signi),
	       fabs (single[j] * signi), fabs (loss[j] * signi),
	       100. * global[j], varnames[j].c_str (), prepros[j], j + 2);
    }

  fprintf (html_file, "</table>\n</body>\n</html>\n");
  fclose (html_file);
}
