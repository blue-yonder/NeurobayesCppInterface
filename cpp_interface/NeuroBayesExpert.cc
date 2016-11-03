#include <cstdlib>
#include <cmath>
#include "NeuroBayesExpert.hh"
#include <iostream>
#include <string.h>
#include <stdio.h>

// these are the Fortran NeuroBayes functions
// since Fortran passes all arguments by
// reference, each argument has to be <var type>*
extern "C"
{

  double nb_f3_e_ (common_t **, float *);
  void nb_wrap_f3_e_ (common_t **, float *, double *result);
  void nb_wrap_invqincl_ (common_t **, float *, float *tabxs,
			  int *nllevel, float *tdeltaval, float *tdeltafrac,
			  float *result);
  void nb_wrap_rndincl_e_ (common_t **, float *, int *debug, float *tabxs,
			   int *nllevel, float *tdeltaval, float *tdeltafrac,
			   float *result);
  void nb_invquant2_ (common_t **, float *t, float *tabg, float *tabxs,
		      float *result);
  void nb_wrap_quantile_ (common_t **, float *, float *tabg, float *tabxs,
			  int *nllevel, float *tdeltaval, double *netout,
			  int *dbg, float *result);
  void nb_wrap_tmean_e_ (common_t **, float *TABG, float *TABXS, float *TRIM,
			 int *nllevel, float *tdeltaval, double *netout,
			 int *dbg, float *result);
  void nb_wrap_eval_bspline_e_ (common_t **, float *x, int *k, float *t,
				int *nkn, float *c, float *der, float *der2,
				float *der3, float *result);
  void nb_wrap_transgle_ (common_t **, float *, float *, float *result);
  void nb_wrap_rndm2_ (common_t **, int *, int *dbg, float *result);

  void nb_wrap_conddensity_ (common_t **, float *, int *dbg, float *tabxs,
			     float *tabd, float *TT, int *NKNOT, float *CP,
			     int *nllevel, float *tdeltaval,
			     float *tdeltafrac, double *netout,
			     float *result);

  void nb_defpolynomials_e_ (common_t **, int *, float *, int *, int *,
			     float *, float *, float *);

  void nb_defexpertise_e_ (common_t **, int *NLAYER, int *NODES, int *ITER,
			   int *IPRUNE, float *W, float *TABX, float *ITABY,
			   float *AE, float *DIAG, float *CHEBY, float *THETA,
			   float *EXPERTISE, int *NumPreproVar, int *PreProc,
			   int *AutoVarSelect, int *ISigSort, float *SigFrac,
			   int *PreproFlag, float *PreproPar, int *Debug,
			   float *RsfTable, int *nMapKey, float *MapKeyValue,
			   float *MapKeyTrans, float *TABX2, float *RsfTable2,
			   int *nlevel, int *node1, int *lshape, int *llog,
			   int *ifixorder, int *ifixshape, int *nvar,
			   int *n_marg_dim, int *marg_varid,
			   float *marg_coeff, float *marg_coeff0,
			   int *nllevel, float *tdeltaval, float *tdeltafrac,
			   float *TABX3, float *RsfTable3, float *a_boost,
			   float *inv_sigma_boostvar, float *diag_boost,
			   float *mean_boostvar, int *ipos_add_vars);

  void nb_defft_e_ (common_t **, float *, float *, float *, float *, int *,
		    float *TT, int *NKNOT, float *CP, int *iversiont);
  void nb_prepro2_e_ (common_t **, float *IN, int *NEVTS1, int *NVARM1,
		      float *TABX, float *ITABY, int *MXNODE, float *SCRATCH,
		      double *CTH, double *STH, float *A, float *DIAG,
		      float *CHEBY, float *THETA, int *IFIXORDER,
		      int *PREPROC, int *Node1, int *LSHAPE, int *NEWEX,
		      int *ISigSort, int *IPRUNE, int *NumPreproVar,
		      int *PreproFlag, float *PreproPar, int *Debug,
		      float *RsfTable, int *nmapkey, float *mapkeyvalue,
		      float *mapkeytrans, int *doHisto, float *OUTLEVEL,
		      float *T, float *TABX2, float *RsfTable2,
		      int *IFixShape, int *NLEVEL, int *nllevel,
		      float *tdeltaval, float *tdeltafrac,
		      bool * is_preboost);


  void nb_chouth_ (common_t **, int *, float *, float *, float *, float *,
		   int *, float *, int *, int *dbg);
  void nb_choutc_ (common_t **, float *, float *, float *,
		   int *, float *, int *, int *nllevel, float *tdeltafrac);
  void nb_forwe_e_ (common_t **, int *, float *,
		    double (*f) (common_t **, float *), float *, int *,
		    int *);
  void nb_splinef2_e_ (common_t **, int *, double *, float *, float *, int *,
		       int *, float *, int *, int *, int *);
  void nb_splinef2_when_needed_e_ (common_t **, double *, float *, int *,
				   int *, float *, int *, int *, int *,
				   bool *);
  void nb_perfplot2_ (common_t **, int *, double *, int *, float *, int *);
  void nb_ftxdef_ (common_t **, float *, float *, float *, float *, int *,
		   int *nllevel, double *netout, int *dbg);
  void nb_transbacktail_ (common_t **, float *, float *, float *, int *dbg);
  void nb_transback_ (common_t **, float *, float *, float *,
		      int *nlLevel, double *netout, int *dbg);

  void nb_converttocarray_e_ (common_t **, const char *ExName, int);

  void nb_readexpertise_e_ (common_t **, const char *ExName, float *Ex, int);

  void nb_derivative_ (common_t **, float *t, float *deri, float *TABXS,
		       float *TT, int *NKNOT, float *CP, int *nllevel,
		       float *tdeltaval, float *tdeltafrac);
  void nb_wrap_ftmean_ (common_t **, float *TABG, float *TABXS,
			double (*f) (float *), int *nllevel, float *tdeltaval,
			double *netout, int *dbg, float *result);

  void nb_wrap_tomarginal_e_ (common_t **, float *X, float *TABX,
			      int *nMapKey, float *MapKeyValue,
			      int *PreproFlag, float *PreproPar, int *Debug,
			      int *n_marg_dim, int *marg_varid,
			      float *marg_coeff, float *marg_coeff0,
			      float *result);

  void nb_wrap_decomposepreproflag_ (common_t **, int *PreproFlag,
				     int *PreproFlag1, int *PreproFlag10,
				     int *PreproFlag100, int *Debug,
				     int *result);

  void nb_prepare_boostdiag_ (bool * is_preboost, bool * is_boost,
			      bool * is_net_plus_diag, int *ifixshape,
			      int *preproc, int *iterations, int *ioffset,
			      int *iversiont);

  void nb_prepro2_boost_e_ (common_t **, float *IN, int *inum, float *scratch,
			    bool * is_boost, bool * is_preboost,
			    int *NumPreproVar_save, int *NumPreproVar,
			    int *preproc_save, int *preproc, int *NODE1,
			    int *node1_save, int *IFIXSHAPE, int *nlevel,
			    int *nllevel, float *outlevel,
			    float *inv_sigma_boostvar, float *a_boost,
			    float *diag_boost, int *NVar, int *nodes,
			    float *mean_boostvar, int *ISigSort, int *nMapKey,
			    float *MapKeyValue, int *ipos_add_vars,
			    int *iversiont, int *final_diagfit);

  void nb_diag_boostnet_e_ (common_t **, float *IN, int *NEVTS1, float *A,
			    int *ioffset, int *nlevel, int *nllevel,
			    int *lshape, float *tabx3, float *rsftable3,
			    int *PREPROC, int *preproc_save, bool * is_boost,
			    int *node1, int *node1_save, int *NumPreproVar,
			    int *NumPreproVar_save, int *ifixshape,
			    int *debug, bool * is_preboost, double *NETOUT,
			    int *NVar, int *nodes, int *iversiont,
			    float *cshape, int *final_diagfit);
  void nb_regul_inputs_ (float *IN, int *NVar, int *inum);
  // Functions called from NeuroBayes libs

  void nbbook1_ (int *, const char *, int *, float *, float *, float *, int);
  void nbfill_e_ (common_t **, int *, float *, float *, float *);

  int nb_def_debugexpert_ (common_t **, const int &);

  void nbpak_e_ (common_t **, int *, float *);
  void nb_histo_init_ ();
  void nb_histo_save_ (const char *, int);

  int nb_check_netout_ (double *netout, int *nlevel);
  void nb_interpolate_tabg_ (double *netout, int *nlevel, float *tabg);

}

/// simple wrapper that uses default arguments for rand function, etc
static common_t *
nb_init_common_simple (int debug_lvl, ec_t ** ec1)
{
  return nb_init_common (debug_lvl, ec1, rand_double1,
			 NULL, 0);
}

double
Expert::NB_F3_e (float a)
{
  double result = 0.;
  nb_wrap_f3_e_ (&com, &a, &result);
  nb_cpp_handle_error (com);
  return result;
}

float
Expert::NB_TRANSGLE (float &a, float *aa)
{
  float result = 0.;
  nb_wrap_transgle_ (&com, &a, aa, &result);
  return result;
}

float
Expert::NB_RNDM2 (int &a, int dbg)
{
  return rand_double ();
}

float
Expert::NB_CONDDENSITY (float &arg)
{
  float result = 0.;
  NB_SPLINEF2_WHEN_NEEDED_e ();
  if (nb_cpp_handle_error (com))
    return result;
  nb_wrap_conddensity_ (&com, &arg, &(com->log->debug_lvl), TABXS, TABD, TT,
			&NKNOT, CP, &nlLevel, &tDeltaVal[0], &tDeltaFrac[0],
			&NETOUT[0], &result);
  return result;
}

float
Expert::NB_BSKFUN_e (float &t, float &Der)
{
  float result = 0.;
  int order = 4;
  float Der2, Der3;
  nb_wrap_eval_bspline_e_ (&com, &t, &order, TT, &NKNOT, CP, &Der, &Der2,
			   &Der3, &result);
  nb_cpp_handle_error (com);
  return result;
}

void
Expert::NB_DefPolynomials_e ()
{
  float XLEVEL = 0;
  nb_cpp_handle_error (com);
  nb_defpolynomials_e_ (&com, &NVar, &XLEVEL, &(com->log->debug_lvl), &NLEVEL,
			OUTLEVEL, XMEAN, &T[0][0]);
}

void
Expert::NBBOOK1 (int a, const char *b, int c, float d, float e, float f,
		 int i)
{
  nbbook1_ (&a, b, &c, &d, &e, &f, i);
}

void
Expert::NBFILL_e (int a, float &b, float c, float &d)
{
  nbfill_e_ (&com, &a, &b, &c, &d);
  nb_cpp_handle_error (com);
}

void
Expert::NB_DEFEXPERTISE_e ()
{
  nb_defexpertise_e_ (&com, &NLAYER, NODES, &ITER, &IPRUNE, &Weights[0][0][0],
		      &TABX[0][0], ITABY, AA, DIAG, CHEBY, THETA, EXPERTISE,
		      &NumPreproVar, &PREPROC, &AutoVarSelect,
		      ISigSort, &SigFrac, PreproFlag, &PreproPar[0][0],
		      &(com->log->debug_lvl), &RsfTable[0][0], nMapKey,
		      &MapKeyValue[0][0], &MapKeyTrans[0][0], &TABX2[0][0],
		      &RsfTable2[0][0], &NLEVEL, &NODE1, &LSHAPE, &LLOG,
		      &IFIXORDER, &IFIXSHAPE, &NVar, &nMargDim, &MargVarid[0],
		      &MargCoeff[0][0], &MargCoeff0, &nlLevel, &tDeltaVal[0],
		      &tDeltaFrac[0], &TABX3[0][0], &RsfTable3[0][0],
		      &a_boost[0], &inv_sigma_boostvar[0], &diag_boost[0],
		      &mean_boostvar[0], &ipos_add_vars[0]);
  nb_cpp_handle_error (com);
}

void
Expert::NB_DEFFT_e ()
{
  nb_defft_e_ (&com, &TABX[0][0], TABXS, ITABY, TABF, &LSHAPE, TT, &NKNOT, CP,
	       &iversiont);
  nb_cpp_handle_error (com);
}

void
Expert::NBPAK_e (int a, float *b)
{
  nbpak_e_ (&com, &a, b);
  nb_cpp_handle_error (com);
}

void
Expert::NB_PREPRO2_e ()
{
  int nVar = NVar - 1;
  int nEv = 1;
  int newExpert = 1;
  int doHisto = 0;
  nb_prepro2_e_ (&com, &IN[0][0], &nEv, &nVar, &TABX[0][0], ITABY, &MXNODE,
		 SCRATCH, &CTH[0][0], &STH[0][0], AA, DIAG, CHEBY, THETA,
		 &IFIXORDER, &PREPROC, &NODE1, &LSHAPE, &newExpert, ISigSort,
		 &IPRUNE, &NumPreproVar, PreproFlag, &PreproPar[0][0],
		 &(com->log->debug_lvl), &RsfTable[0][0], nMapKey,
		 &MapKeyValue[0][0], &MapKeyTrans[0][0], &doHisto, OUTLEVEL,
		 &T[0][0], &TABX2[0][0], &RsfTable2[0][0], &IFIXSHAPE,
		 &NLEVEL, &nlLevel, &tDeltaVal[0], &tDeltaFrac[0],
		 &is_preboost);
  nb_cpp_handle_error (com);
}

void
Expert::NB_CHOUTH (int &a, float *b, float *c, float *d, float *e,
		   int &f, float *g, int &h, int dbg)
{
  nb_chouth_ (&com, &a, b, c, d, e, &f, g, &h, &dbg);
}

void
Expert::NB_CHOUTC (float *b, float *c, float *d, int &e, float *f, int &g)
{
  nb_choutc_ (&com, b, c, d, &e, f, &g, &nlLevel, &tDeltaFrac[0]);
}

void
Expert::NB_FORWE_e (int *a, float *b, double (*f) (common_t **, float *),
		    float *d, int &e, int &g)
{
  nb_forwe_e_ (&com, a, b, f, d, &e, &g);
  nb_cpp_handle_error (com);
}

void
Expert::NB_SPLINEF2_e (int &histoID)
{
  nb_splinef2_e_ (&com, &histoID, NETOUT, TABG, TABD, &MODEGS, &NGS, &REGGS,
		  &NLEVEL, &(com->log->debug_lvl), &nlLevel);
  nb_cpp_handle_error (com);
}

void
Expert::NB_SPLINEF2_WHEN_NEEDED_e ()
{
  nb_splinef2_when_needed_e_ (&com, NETOUT, TABD, &MODEGS, &NGS, &REGGS,
			      &NLEVEL, &(com->log->debug_lvl), &nlLevel,
			      &spline_performed);
  nb_cpp_handle_error (com);
}

void
Expert::NB_PERFPLOT2 (int histoID)
{
  nb_perfplot2_ (&com, &histoID, NETOUT, &NLEVEL, OUTLEVEL, &nlLevel);
}

void
Expert::NB_FTXDEF (int &mmax)
{
  nb_ftxdef_ (&com, TABXS, TABF, TABD, TABL, &mmax,
	      &nlLevel, &NETOUT[0], &(com->log->debug_lvl));
}

void
Expert::NB_FTXDEF_e (float *array, int nBins)
{
  float t, s, der;
  float gs = 0;
  float step = TABXS[NB_NVALUE - 1] - TABXS[0];
  float Der = 0;

  for (int ii = 0; ii < nBins; ++ii)
    {
      t = TABXS[0] + step * (ii + 0.5) / (float) nBins;
      s = NB_TRANSGLE (t, TABXS);
      // computation of TABF[ii]
      NB_BSKFUN_e (t, Der);
      if (nb_cpp_handle_error (com))
	return;
      der = -Der / 2.;
      NB_TRANSBACK (s, gs);
      array[ii] = gs * der;
      if (nb_get_debug (&com) >= 2)
	{
	  ss << "t = " << t << "\ts = " << s << "\tgs = " << gs << "\tder = "
	    << der << "\tarray[" << ii << "] = " << array[ii] << std::endl;
	  nb_cpp_log (ss, 2, &com);
	}
    }
}

void
Expert::NB_TRANSBACKTAIL (float &randomN, float &result)
{
  nb_transbacktail_ (&com, &randomN, TABXS, &result, &(com->log->debug_lvl));
}

void
Expert::NB_TRANSBACK (float &eq, float &result)
{
  nb_transback_ (&com, &eq, TABD, &result,
		 &nlLevel, &NETOUT[0], &(com->log->debug_lvl));
}

float
Expert::NB_QUANTILE (float value)
{
  float result = 0.;
  nb_wrap_quantile_ (&com, &value, TABG, TABXS,
		     &nlLevel, &tDeltaVal[0], &NETOUT[0],
		     &(com->log->debug_lvl), &result);
  return result;
}

float
Expert::NB_BINCLASS ()
{
  return (float) NETOUT[0];
}

float
Expert::NB_INVQUANT (float value)
{
  float result = 0.;
  nb_invquant2_ (&com, &value, &TABG[0], &TABXS[0], &result);
  return result;
}

float
Expert::NB_TMEAN (float value)
{
  float result = -999;
  nb_wrap_tmean_e_ (&com, TABG, TABXS, &value,
		    &nlLevel, &tDeltaVal[0], &NETOUT[0],
		    &(com->log->debug_lvl), &result);
  nb_cpp_handle_error (com);
  return result;
}

float
Expert::NB_INVQINCL (float argument)
{
  float result = 0;
  nb_wrap_invqincl_ (&com, &argument, TABXS,
		     &nlLevel, &tDeltaVal[0], &tDeltaFrac[0], &result);
  return result;
}

float
Expert::NB_RNDINCL (float argument)
{
  float result = 0;
  nb_wrap_rndincl_e_ (&com, &argument, &(com->log->debug_lvl), TABXS,
		      &nlLevel, &tDeltaVal[0], &tDeltaFrac[0], &result);
  nb_cpp_handle_error (com);
  return result;
}

void
Expert::NB_DERIVATIVE (float t, float *der)
{
  nb_derivative_ (&com, &t, der, TABXS, TT, &NKNOT, CP, &nlLevel,
		  &tDeltaVal[0], &tDeltaFrac[0]);
}

void
Expert::NB_READEXPERTISEC_e (const char *ExName, float *Ex)
{
  std::string exName = ExName;
  nb_readexpertise_e_ (&com, exName.c_str (), Ex, exName.size ());
  nb_cpp_handle_error (com);
}

float
Expert::NB_FTMEAN (double (*f) (float *))
{
  float result = 0.;
  nb_wrap_ftmean_ (&com, TABG, TABXS, f,
		   &nlLevel, &tDeltaVal[0], &NETOUT[0],
		   &(com->log->debug_lvl), &result);
  return result;
}

float
Expert::NB_TOMARGINAL_e (float *X)
{
  float result = 0.;
  nb_wrap_tomarginal_e_ (&com, X, TABX[0], nMapKey, MapKeyValue[0],
			 PreproFlag, PreproPar[0], &(com->log->debug_lvl),
			 &nMargDim, MargVarid, MargCoeff[0], &MargCoeff0,
			 &result);
  nb_cpp_handle_error (com);
  return result;
}

void
Expert::nb_prepare_boostdiag ()
{
  nb_prepare_boostdiag_ (&is_preboost, &is_boost, &is_net_plus_diag,
			 &IFIXSHAPE, &PREPROC, &iterations, &ioffset,
			 &iversiont);
}

void
Expert::nb_prepro2_boost_e ()
{
  int n_events = 1;
  nb_prepro2_boost_e_ (&com, &IN[0][0], &n_events, SCRATCH, &is_boost,
		       &is_preboost, &NumPreproVar_save, &NumPreproVar,
		       &preproc_save, &PREPROC, &NODE1, &node1_save,
		       &IFIXSHAPE, &NLEVEL, &nlLevel, OUTLEVEL,
		       &inv_sigma_boostvar[0], a_boost, diag_boost, &NVar,
		       NODES, &mean_boostvar[0], ISigSort, nMapKey,
		       MapKeyValue[0], &ipos_add_vars[0], &iversiont,
		       &final_diagfit);
  nb_cpp_handle_error (com);
}

void
Expert::nb_diag_boostnet_e ()
{
  int n_events = 1;
  nb_diag_boostnet_e_ (&com, &IN[0][0], &n_events, &A[0][0], &ioffset,
		       &NLEVEL, &nlLevel, &LSHAPE, &TABX3[0][0],
		       &RsfTable3[0][0], &PREPROC, &preproc_save, &is_boost,
		       &NODE1, &node1_save, &NumPreproVar, &NumPreproVar_save,
		       &IFIXSHAPE, &(com->log->debug_lvl), &is_preboost,
		       &NETOUT[0], &NVar, NODES, &iversiont, &CSHAPE[0],
		       &final_diagfit);
  nb_cpp_handle_error (com);
}

Expert::Expert ()
{
  std::cout << "NeuroBayes(R): wrong constructor called, abort" << std::endl;
}

bool
Expert::check_num_inputs (int num_inputs)
{
  if (num_inputs > NODE1 - 1)
    {
      if (nb_get_debug (&com) >= -2)
	{
	  ss << "ERROR: num_inputs = " << num_inputs
	    << ", which is bigger than NODE1 - 1 = " << NODE1 -
	    1 << std::endl;
	  nb_cpp_log (ss, -2, &com);
	}
      return false;
    }
  return true;
}

Expert::Expert (float *myExpertise, int myDebug, bool writeout,
		ec_t ** ec, log_func_t log_f, void *log_enclosed,
		delete_enclosed_func_t log_delete_enclosed)
{
  com = nb_init_common_simple (myDebug, ec);
  if (nb_get_debug (&com) > -3 && log_f)
    {
      nb_register_logging (&com, log_f, log_enclosed, log_delete_enclosed);
    }
  cg_of_nb_cols = NULL;
  n_nb_cols = -1;
  tellinputs_mode = 0;
  filename = NULL;
  if (nb_get_debug (&com) >= 1)
    {
      ss << "NeuroBayes Expert  -- constructor" << std::endl;
      nb_cpp_log (ss, 1, &com);
    }
  writeoutdata = writeout;

  EXPERTISE = new float[NB_NEXPERTISE];
  if (myExpertise == NULL)
    {
      ss << "NeuroBayes Expert  -- Expertise array is NULL" << std::endl;
      char *msg = stream_c_str (ss);
      if (nb_get_debug (&com) >= -2)
	nb_c_log (msg, -2, &com);
      nb_set_error (&com, invalid_arg_exc, msg);
      free (msg);
      nb_cpp_handle_error (com);
      return;
    }

  int expSize = (int) myExpertise[19];
  if (nb_get_debug (&com) >= 1)
    {
      ss << "read size " << expSize << " and max is " << NB_NEXPERTISE <<
	std::endl;
      nb_cpp_log (ss, 1, &com);
    }
  if (expSize < 1)
    {
      ss << "Expert constructor error: the array passed might be corrupted" <<
	std::endl;
      char *msg = stream_c_str (ss);
      if (nb_get_debug (&com) >= -2)
	nb_c_log (msg, -2, &com);
      nb_set_error (&com, invalid_arg_exc, msg);
      free (msg);
      nb_cpp_handle_error (com);
      return;
    }

  // fill the expertise array by copying the content of the passed array
  // (this is to avoid problem with the destructor, which deletes EXPERTISE
  // -> if the EXPERTISE is not initialized with new we end up trying to delete
  // an object we do not own)
  for (int k = 0; k < NB_NEXPERTISE; ++k)
    if (k < expSize)
      EXPERTISE[k] = myExpertise[k];
    else
      EXPERTISE[k] = 0;

  // now call initialisation routine
  if (nb_get_debug (&com) >= 1)
    {
      ss << " initialise NeuroBayes(R) Expert " << std::endl;
      nb_cpp_log (ss, 1, &com);
    }
  initialise_e ();
  if (nb_cpp_handle_error (com))
    return;


  if (nb_get_debug (&com) >= 1)
    {
      ss << "NeuroBayes Expert - constructor ends" << std::endl;
      nb_cpp_log (ss, 1, &com);
    }
}

Expert::Expert (const std::string ExpertiseName, int myDebug,
		bool writeout, ec_t ** ec, log_func_t log_f,
		void *log_enclosed,
		delete_enclosed_func_t log_delete_enclosed)
{
  init_expert (ExpertiseName, myDebug, writeout, ec,
	       log_f, log_enclosed, log_delete_enclosed);
}

Expert::Expert (const char *ExpertiseName, int myDebug,
		bool writeout, ec_t ** ec, log_func_t log_f,
		void *log_enclosed,
		delete_enclosed_func_t log_delete_enclosed)
{
  std::string Name = "";
  if (ExpertiseName)
    {
      Name = std::string (ExpertiseName);
    }
  init_expert (Name, myDebug, writeout, ec,
	       log_f, log_enclosed, log_delete_enclosed);
}

void
Expert::init_expert (const std::string ExpertiseName, int myDebug,
		     bool writeout, ec_t ** ec, log_func_t log_f,
		     void *log_enclosed,
		     delete_enclosed_func_t log_delete_enclosed)
{
  com = nb_init_common_simple (myDebug, ec);
  if (nb_get_debug (&com) > -3 && log_f)
    {
      nb_register_logging (&com, log_f, log_enclosed, log_delete_enclosed);
    }
  cg_of_nb_cols = NULL;
  n_nb_cols = -1;
  filename = NULL;
  is_valid = true;
  tellinputs_mode = 0;
  // create arrays
  EXPERTISE = new float[NB_NEXPERTISE];
  if (nb_cpp_handle_error (com))
    return;

  if (nb_get_debug (&com) >= 0)
    {
      ss << "NeuroBayes Expert  -- constructor" << std::endl;
      nb_cpp_log (ss, 0, &com);
    }

  writeoutdata = writeout;

  // reading in Expertise
  if (nb_get_debug (&com) >= 1)
    {
      ss << "NeuroBayes Expert: reading in Expertise" << std::endl;
      nb_cpp_log (ss, 1, &com);
    }
  if (nb_get_debug (&com) >= 1)
    {
      ss << " now read expertise " << std::endl;
      ss << " size of expertise  " << NB_NEXPERTISE << std::endl;
      nb_cpp_log (ss, 1, &com);
    }				// if debug

  NB_READEXPERTISEC_e (ExpertiseName.c_str (), EXPERTISE);
  if (nb_cpp_handle_error (com))
    return;
  // now call initialisation routine
  if (nb_get_debug (&com) >= 1)
    {
      ss << " initialise NeuroBayes(R) Expert " << std::endl;
      nb_cpp_log (ss, 1, &com);
    }

  initialise_e ();
  if (nb_cpp_handle_error (com))
    return;

  if (nb_get_debug (&com) >= 1)
    {
      ss << "NeuroBayes Expert - constructor ends" << std::endl;
      nb_cpp_log (ss, 1, &com);
    }
}

void
Expert::initialise_e ()
{

  if (nb_get_debug (&com) >= 1)
    {
      ss << "NeuroBayes Expert - initialise begins" << std::endl;
      nb_cpp_log (ss, 1, &com);
    }

  InfoAndLicence ();

  LDEFSPLINE = 0;
  memset (XSAVE, 0, NB_MAXNODE * sizeof (float));
  memset (&T[0][0], 0, NB_MAXNODE * (NB_LEVMAX + 1) * sizeof (float));

  NumPreproVar = 0;		//variable selection, 0 if not used

  // there is only one Expertise per class instance of 
  // NeuroBayes Expert, i.e. the Expertise is always new
  // in the constructor
  NewEx = 1;			//assume new expertise
  NEWEVT = 0;
  eventCounter = 0;

  MODEGS = 0;
  NGS = 0;
  REGGS = 0.;

  memset (&TABX[0][0], 0, NB_MAXNODE * (NB_NVALUE) * sizeof (float));

  NLAYER = NB_MAXLAYER;
  MXNODE = NB_MAXNODE;

  // need this call here also if writeout is false
  // because a histogram is booked inside NB_DEFFT
  nb_histo_init_ ();

  // decode expertise
  if (nb_get_debug (&com) >= 2)
    {
      ss << "NeuroBayes Expert: now decode Expertise" << std::endl;
      nb_cpp_log (ss, 2, &com);
    }
  NB_DEFEXPERTISE_e ();
  if (nb_cpp_handle_error (com))
    return;

  iversiont = int (EXPERTISE[0]) + 20000000;
  iterations = int (EXPERTISE[14]);
  nb_prepare_boostdiag ();	// check version and set flags

  if (nb_get_debug (&com) >= 2)
    {
      PrintArrays ("After DefExpertise");
      ss << "NeuroBayes Expert (C++): initialise orthogonal polynomials" <<
	std::endl;
      nb_cpp_log (ss, -2, &com);
    }

  NB_DefPolynomials_e ();
  if (nb_cpp_handle_error (com))
    return;

  if (nb_get_debug (&com) >= 2)
    {
      ss << "NeuroBayes Expert (C++): orth. pol. initialised " << std::endl;
      nb_cpp_log (ss, 2, &com);
    }

  if (nb_get_debug (&com) >= 1)
    {
      for (int ii = 0; ii < NB_MAXNODE; ii++)
	{
	  ss << "Expert (C++): T(" << ii << ",ii):" << std::endl;
	  for (int j = 0; j < NB_LEVMAX + 1; j++)
	    ss << T[j][ii] << "  ";
	  ss << std::endl;
	}
      nb_cpp_log (ss, 1, &com);
    }

  //   adjustments for target delta functions
  if (nlLevel > 0)
    {
      SigFrac = tDeltaFrac[nlLevel - 1];
      for (int k = NLEVEL - 1; k >= 0; --k)
	OUTLEVEL[nlLevel + k] = OUTLEVEL[k];
      for (int k = 0; k < nlLevel - 1; ++k)
	OUTLEVEL[k] = 0.5 * (tDeltaVal[k] + tDeltaVal[k + 1]);
      OUTLEVEL[nlLevel - 1] = 0.5 * tDeltaVal[nlLevel - 1];
    }

  // define levels (based on number of nodes in output layer)
  // for classification
  if (IFIXSHAPE == 1)
    {
      switch (LSHAPE)
	{
	case 0:
	  XSHAPE[0] = -log (1.0 / SigFrac - 1.0);
	  break;
	case 1:
	  for (int k = 0; k < NLEVEL; ++k)
	    XSHAPE[k] = log (1. / OUTLEVEL[k] - 1.0);
	  SigFrac = 1.;
	  break;
	case 2:
	  SigFrac = tDeltaFrac[nlLevel - 1];
	  for (int k = 0; k < nlLevel; ++k)
	    XSHAPE[k] = -log (1. / tDeltaFrac[k] - 1.0);
	  for (int k = nlLevel; k < NLEVEL + nlLevel; ++k)
	    XSHAPE[k] = -log (1.0 / (SigFrac * (1 - OUTLEVEL[k])) - 1.0);
	  break;
	}
      if (nb_get_debug (&com) >= 1)
	{
	  ss << "Expert: Define Levels:" << std::endl;
	  for (int i = 0; i < NLEVEL + nlLevel; ++i)
	    ss << "SigFrac = " << SigFrac << "\tOUTLEVEL[" << i << "] = " <<
	      OUTLEVEL[i] << "\tXSHAPE[" << i << "] = " << XSHAPE[i] <<
	      "\tF3(XSHAPE[" << i << "]) = " << NB_F3_e (XSHAPE[i]) << std::
	      endl;
	  if (nb_cpp_handle_error (com))
	    return;
	  nb_cpp_log (ss, 1, &com);
	}
    }

  //     new expertise, calc everthing
  if (nb_get_debug (&com) >= 2)
    {
      ss << "NeuroBayes Expert: calc. incl. density" << std::endl;
      nb_cpp_log (ss, 2, &com);
    }

  if (LSHAPE != 0)
    {
      NB_DEFFT_e ();
      if (nb_cpp_handle_error (com))
	return;
    }

  if (nb_get_debug (&com) >= 2)
    {
      ss << "NeuroBayes Expert: incl. density calculated" << std::endl;
      nb_cpp_log (ss, 2, &com);
    }

  // initialize and fill the inverse array of ISigSort 
  for (int i = 0; i < NB_MAXNODE - 1; ++i)
    inverseISigSort[i] = 0;
  if (NumPreproVar > 0)
    {
      for (int i = 0; i < NODE1 - 1; ++i)
	inverseISigSort[ISigSort[NODE1 - 2 - i] - 2] = i + 1;
    }
  else
    {
      // case in which there is no significance cut
      for (int i = 0; i < NB_MAXNODE - 1; ++i)
	inverseISigSort[i] = i + 1;
    }

  if (nb_get_debug (&com) >= 1)
    {
      for (int i = 0; i < NODE1 - 1; ++i)
	{
	  if (inverseISigSort[i] != 0)
	    {
	      ss << "inverseISigSort[" << i << "] = " << inverseISigSort[i]
		<< "\tprepro = " << PreproFlag[inverseISigSort[i]]
		<< "\tTABX extremes = " << TABX[inverseISigSort[i]][0]
		<< "\t" << TABX[inverseISigSort[i]][NB_NVALUE -
						    1] << std::endl;
	      nb_cpp_log (ss, 1, &com);
	    }
	}
    }
  if (writeoutdata)
    {
      NBBOOK1 (498, "TABX   ", 101, 0., 101., 0., 8);
      NBPAK_e (498, &TABX[0][0]);
      if (nb_cpp_handle_error (com))
	return;
      NBBOOK1 (499, "TABF ", 100, TABX[0][0], TABX[0][100], 0., 5);
      NBPAK_e (499, TABF);
      if (nb_cpp_handle_error (com))
	return;
      NBBOOK1 (987, "E likelihood sum  ", 100, TABX[0][0], TABX[0][100], 0.0,
	       18);
      NBBOOK1 (986, "s likelihood sum  ", 101, 0., 1.00001, 0.0, 18);
      //book histogrammes for purity/efficiency plots
      for (int jj = 0; jj < NLEVEL + nlLevel; jj++)
	{
	  NBBOOK1 (100 + jj + 1, "F(NETOUT) FOR BACKGROUND  ", 100, -1., 1.,
		   0., 26);
	  NBBOOK1 (200 + jj + 1, "F(NETOUT) FOR SIGNAL  ", 100, -1., 1., 0.,
		   22);
	}
    }

  if (nb_get_debug (&com) >= 1)
    {
      ss << "-------------------------------------" << std::endl;
      ss << "Expert: actual network topology:" << std::endl;
      ss << "NODES[1]     = " << NODES[0] << std::endl;
      ss << "NODES[2]     = " << NODES[1] << std::endl;
      ss << "NODES[3]     = " << NODES[2] << std::endl;
      ss << "nominal node1= " << NODE1 << std::endl;
      ss << "-------------------------------------" << std::endl;
      ss << "NeuroBayes Expert: initialise ends" << std::endl;
      nb_cpp_log (ss, 1, &com);
    }

}

Expert::~Expert ()
{
  // save histograms
  if (writeoutdata)
    {
      if (HistoFileName.size () < 1)
	HistoFileName = "expertHistos.hist";
      nb_histo_save_ (HistoFileName.c_str (), (int) HistoFileName.size ());
    }
  delete[]EXPERTISE;
  nb_delete_common (&com);
}

void
Expert::Print ()
{
  if (nb_get_debug (&com) >= -1)
    {
      ss << "Greetings" << std::endl;
      nb_cpp_log (ss, -1, &com);
    }
}

float
Expert::nb_expert (ACTION key, double *input, float ARGUMENT)
{

  float X[NB_MAXNODE];

  memset (X, 0, NB_MAXNODE * sizeof (float));	//reset

  for (int i = 0; i < NODE1 - 1; ++i)
    X[i] = (float) input[i];

  return nb_expert (key, X, ARGUMENT);
}

void
Expert::writeErrorMessage (int varIndex, float limit,
			   float current, int status, int prepro)
{
  char help[100] = "";
  sprintf (errorMessage,
	   "****************************************************\n");
  strcat (errorMessage,
	  "** NeuroBayes Expert Warning (only last one shown)**\n");
  strcat (errorMessage,
	  "****************************************************\n");
  sprintf (help, "Variable: %3i with flag %i value", varIndex, prepro);
  strcat (errorMessage, help);
  if (status == 0)
    strcat (errorMessage, " too low\nLowest");
  else if (status == 1)
    strcat (errorMessage, " too high\nHighest");
  else
    strcat (errorMessage, " unknown");
  sprintf (help, " known input: %5.4f\n", limit);
  strcat (errorMessage, help);
  sprintf (help, "Current input     : %5.4f\n", current);
  strcat (errorMessage, help);
  strcat (errorMessage,
	  "****************************************************\n");
}

void
Expert::checkInputRange_e (float *X)
{
  // loop to check the values passed (if INF or NAN issue an error message and abort)
  for (int ii = 0; ii < NODE1 - 1; ++ii)
    {
      if (isnan (X[ii]))
	{
	  ss << ii << "\t" << NODE1 << std::endl;
	  ss << " Expert::nb_expert  ERROR: NAN passed as input" << std::endl;
	  for (int kk = 0; kk < NODE1; kk++)
	    {
	      ss << " input variable " << kk << " = " << X[kk] << std::endl;
	    }
	  ss << "\n Please, correct the error" << std::endl;
	  char *msg = stream_c_str (ss);
	  if (nb_get_debug (&com) >= -2)
	    nb_c_log (msg, -2, &com);
	  nb_set_error (&com, invalid_data_exc, msg);
	  free (msg);
	  nb_cpp_handle_error (com);
	  return;
	}
      if (isinf (X[ii]))
	{
	  ss << "Expert::nb_expert  ERROR: INF passed as input" << std::endl;
	  for (int kk = 0; kk < NODE1; kk++)
	    {
	      ss << " input variable " << kk << " = " << X[kk] << std::endl;
	    }
	  ss << "\n Please, correct the error" << std::endl;
	  char *msg = stream_c_str (ss);
	  if (nb_get_debug (&com) >= -2)
	    nb_c_log (msg, -2, &com);
	  nb_set_error (&com, invalid_data_exc, msg);
	  free (msg);
	  nb_cpp_handle_error (com);
	  return;
	}

      // Check here if the variables passed are within the range given during the training
      // -- since the variable indices in the PreproFlag and in the TABX array are sorted
      //    according to the variable ranking, we have to check the original index with
      //    the array ISigSort
      int indexToCheck = inverseISigSort[ii];
      // array for the decomposed prepro-flags
      // 0 -> PreproFlag1
      // 1 -> PreproFlag10
      // 2 -> PreproFlag100
      int preproFlags[3] = { 0, 0, 0 };
      int errorCode = 0;
      nb_wrap_decomposepreproflag_ (&com, &(PreproFlag[indexToCheck]),
				    &(preproFlags[0]), &(preproFlags[1]),
				    &(preproFlags[2]), &(com->log->debug_lvl),
				    &errorCode);

      if ((errorCode == 0 || errorCode == 3) &&
	  PreproFlag[indexToCheck] != -999. &&
	  ((preproFlags[1] != 3 && preproFlags[1] != 4
	    && preproFlags[1] != 7 && preproFlags[1] != 8
	    && preproFlags[1] != 9) || int (X[ii]) != -999.))
	{

	  if (X[ii] < TABX[indexToCheck][0])
	    writeErrorMessage (ii + 1, TABX[indexToCheck][0], X[ii], 0,
			       PreproFlag[indexToCheck]);
	  if (X[ii] > TABX[indexToCheck][NB_NVALUE - 1])
	    writeErrorMessage (ii + 1, TABX[indexToCheck][NB_NVALUE - 1],
			       X[ii], 1, PreproFlag[indexToCheck]);
	}
    }
  return;
}

int
Expert::isNewEvent (float *X)
{
  int newEvent = 0;
  for (int n = 0; n < NODE1 - 1; ++n)
    {
      if (X[n] != XSAVE[n])
	{
	  newEvent = 1;
	  XSAVE[n] = X[n];
	}
    }
  return newEvent;
}

// The x array contains the input variables for NeuroBayes, not the target.
// Please store your target somewhere else or at indices >= NODE1 - 1 in the
// x array.
void
Expert::calcNetOutput_e (float *X)
{
  //  new event, calc network output:
  //  perform preprocessing
  //  obsolete: one event at a time (i.e. absn(argument no 2)=1)
  //  obsolete: argument no 2 (target) is negative here
  for (int j = 0; j < NODE1 - 1; ++j)
    IN[0][j + 1] = X[j];

  if (nb_get_debug (&com) >= 2)
    {
      ss << "Expert (C++) after filling array IN" << std::endl;
      ss << "       array IN" << std::endl;
      ss << "---------------------------" << std::endl;
      for (int j = 0; j < NB_MAXNODE; ++j)
	ss << "IN(" << j << ",1) = " << IN[0][j] << std::endl;
      nb_cpp_log (ss, 2, &com);
      PrintArrays ("before PREPRO2");
    }

  NB_PREPRO2_e ();
  if (nb_cpp_handle_error (com))
    return;

  if (is_preboost)
    {
      nb_prepro2_boost_e ();
      if (nb_cpp_handle_error (com))
	return;
    }
  int one_event = 1;
  if (iterations > 0 && iversiont >= NB_DATE_BOOST2)
    nb_regul_inputs_ (&IN[0][0], &NVar, &one_event);
  if (nb_get_debug (&com) >= 2)
    PrintArrays ("after PREPRO2");

  A[0][0] = 1.0;		//bias node
  for (int n = 1; n < NVar; ++n)
    {
      A[n][0] = IN[0][n];
      XPRE[n] = IN[0][n];
    }

  if (nb_get_debug (&com) >= 2)
    {
      ss << " Expert (C++) after filling array A" << std::endl;
      ss << " Expert (C++) NODE1        = " << NODE1 << std::endl;
      ss << " Expert (C++) NumPreproVar = " << NumPreproVar << std::endl;
      ss << " Expert (C++) NVar       = " << NVar << std::endl;
      ss << " Expert (C++) array IN" << std::endl;
      ss << "---------------------------" << std::endl;

      for (int m = 0; m < NB_MAXNODE; ++m)
	ss << "IN(" << m << ",1) = " << IN[0][m] << std::endl;

      ss << "array A" << std::endl;
      ss << "---------------------------" << std::endl;
      for (int n = 0; n < NB_MAXLAYER; ++n)
	{
	  ss << "A for layer " << n + 1 << std::endl;
	  for (int j = 0; j < NB_MAXNODE; ++j)
	    ss << A[j][n] << "  ";
	  ss << std::endl;
	}
      nb_cpp_log (ss, 2, &com);
    }
  if (LSHAPE == 0 && IFIXSHAPE == 3)
    {
      if (iversiont >= NB_DATE_BOOST2 && is_boost)
	{
	  CSHAPE[0] = IN[0][NB_MAXNODE + 2];
	}
      else
	{
	  CSHAPE[0] = A[1][0];
	}
    }
  else if (LSHAPE >= 1)
    {
      switch (IFIXSHAPE)
	{
	case 2:
	  if (PREPROC / 10 == 4)
	    NB_CHOUTH (NVar, &A[0][0], CHEBY, CSHAPE, OUTLEVEL,
		       NLEVEL, &T[0][0], IFIXORDER, nb_get_debug (&com));
	  else
	    NB_CHOUTC (&A[0][0], CHEBY, CSHAPE, NLEVEL, &T[0][0], IFIXORDER);
	  break;
	case 3:
	  for (int i = 0; i < NLEVEL + nlLevel; i++)
	    {
	      if (iversiont >= NB_DATE_BOOST2 && is_boost)
		{
		  CSHAPE[i] = IN[0][NB_MAXNODE + NLEVEL + nlLevel + i + 1];
		}
	      else
		{
		  CSHAPE[i] = IN[0][IPRUNE + i];
		}
	    }
	  break;
	default:
	  break;
	}
    }
  if (nb_get_debug (&com) >= 2)
    PrintArrays ("CHOUTC/CHOUTH");
  NB_FORWE_e (NODES, &Weights[0][0][0], nb_f3_e_, &A[0][0], MXNODE, NLAYER);
  if (nb_cpp_handle_error (com))
    return;

  if (is_net_plus_diag || is_boost)
    {
      if (is_net_plus_diag)
	final_diagfit = 3;
      nb_diag_boostnet_e ();
      if (nb_cpp_handle_error (com))
	return;
    }
  else
    {
      for (int nn = 0; nn < NLEVEL + nlLevel; nn++)
	{
	  if (LLOG == 3)
	    NETOUT[nn] = A[nn][NB_MAXLAYER - 1];
	  else
	    {
	      switch (IFIXSHAPE)
		{
		case 1:
		  NETOUT[nn] = NB_F3_e (A[nn][NB_MAXLAYER - 1] + XSHAPE[nn]);
		  break;
		case 2:
		case 3:
		  NETOUT[nn] = NB_F3_e (A[nn][NB_MAXLAYER - 1] + CSHAPE[nn]);
		  break;
		default:
		  NETOUT[nn] = NB_F3_e (A[nn][NB_MAXLAYER - 1]);
		}
	      if (nb_cpp_handle_error (com))
		return;
	    }
	}
    }
  if (LSHAPE >= 1)
    {
      if (nb_check_netout_ (&NETOUT[0], &NLEVEL) == 0)
	{
	  nb_interpolate_tabg_ (&NETOUT[0], &NLEVEL, &TABG[0]);
	  spline_performed = false;
	}
      else
	{
	  ++eventCounter;
	  int histoID = -eventCounter;
	  if (eventCounter < 100.)
	    histoID = eventCounter;
	  NB_SPLINEF2_e (histoID);
	  if (nb_cpp_handle_error (com))
	    return;
	  spline_performed = true;
	  if (histoID > 0)
	    NB_PERFPLOT2 (eventCounter);
	}
    }
  if (IFIXSHAPE == 4)
    {				// marginal sums
      NETOUT[0] = NB_TOMARGINAL_e (X);
      if (nb_cpp_handle_error (com))
	return;
      NETOUT[0] = 2 * NETOUT[0] - 1;
    }
  return;
}


float
Expert::nb_expert (ACTION key, float *input, float ARGUMENT)
{				//nb_expert
  //used histograms
  //498                         TABX
  //499                         TABF
  //497                         inclusive cumulative distribution 
  //496                         inclusive distribution 
  //495                         inclusive distribution derivative 
  //494                         inclusive distribution curvature 
  //493                         input for 497 fit (compare//////) 
  //4                         smoothed TABX array 
  //987                         E likelihood sum
  //986                         s (0:1) likelihood sum
  //50000+i (1=1,10)            f(t)
  //10300+id                    spline
  //10400+id                    diff spline
  //10600+id                    d2 spline
  //10700+id                    d3 spline      
  float nb_expert_delta = -999.99;
  if (nb_cpp_handle_error (com))
    return nb_expert_delta;

  // check what parameter we got
  if (isnan (ARGUMENT) || isinf (ARGUMENT))
    {
      ss << " Expert::nb_expert  ERROR: " << ARGUMENT <<
	" passed as parameter for ACTION " << key << " !" << std::endl;
      char *msg = stream_c_str (ss);
      if (nb_get_debug (&com) >= -2)
	nb_c_log (msg, -2, &com);
      nb_set_error (&com, invalid_data_exc, msg);
      free (msg);
      nb_cpp_handle_error (com);
      return nb_expert_delta;
    }

  memset (&A[0][0], 0, NB_MAXNODE * NB_MAXLAYER * sizeof (float));
  memset (&IN[0][0], 0, NB_MAXDIM * sizeof (float));
  errorMessage[0] = '\0';
  //1: use only one event at a time

  // copy input field
  float X[NB_MAXNODE];

  memset (X, 0, NB_MAXNODE * sizeof (float));	//reset 

  for (int jj = 0; jj < NODE1 - 1; ++jj)
    X[jj] = input[jj];

  //  check if this is a new event
  NEWEVT = isNewEvent (X);

  if (NEWEVT == 1)
    {
      NTABL = 0;
      checkInputRange_e (X);
      if (nb_cpp_handle_error (com))
	return nb_expert_delta;
    }

  if (NewEx != 1 && NEWEVT != 1 && LDEFSPLINE == 1)
    {
      if (nb_get_debug (&com) >= -2)
	{
	  ss <<
	    "WARNING, NB_DEFSPLINE HAD PROBABLY BEEN CALLED \nAT WRONG PLACE"
	    << std::
	    endl <<
	    "NB_DEFSPLINE has to be called BEFORE the first\ncall to NB_EXPERT"
	    << std::
	    endl << "with either -- a new EXPERTISE \n\t or -- a new EVENT" <<
	    std::endl;
	  nb_cpp_log (ss, -2, &com);
	}
    }
  LDEFSPLINE = 0;
  if (NewEx == 1 || NEWEVT == 1)
    {
      calcNetOutput_e (X);
      if (nb_cpp_handle_error (com))
	return nb_expert_delta;
    }
  int mMax = 0;
  float randomN = 0.;
  int intArg = 0;
  float parameter = 0;
  int nBins = 100;
  if (LSHAPE == 0 && (key == BINCLASS || key == NOOP))
    {
      NB_EXPERT = (float) NETOUT[0];
    }
  else if (LSHAPE != 0)
    {
      switch (key)
	{			// switch key
	case INCLDENSITY:
	  parameter = -999.;
	  NB_DERIVATIVE (ARGUMENT, &parameter);
	  NB_EXPERT = parameter;
	  break;
	case TMIN:
	  NB_EXPERT = TABX[0][0];
	  break;
	case TMAX:
	  NB_EXPERT = TABX[0][NB_NVALUE - 1];
	  break;
	case MEDIAN:
	  NB_EXPERT = NB_QUANTILE (0.5);
	  break;
	case LERROR:
	  NB_EXPERT = NB_QUANTILE (0.8413);
	  break;
	case RERROR:
	  NB_EXPERT = NB_QUANTILE (0.1587);
	  break;
	case MEAN:
	  NB_EXPERT = NB_TMEAN (0.5);
	  break;
	case TRIM:
	  NB_EXPERT = NB_TMEAN (ARGUMENT);
	  break;
	case QUANTILE:
	  NB_EXPERT = NB_QUANTILE (ARGUMENT);
	  break;
	case INVQUANT:
	  NB_EXPERT = NB_INVQUANT (ARGUMENT);
	  break;
	case INVQINCL:
	  NB_EXPERT = NB_INVQINCL (ARGUMENT);
	  break;
	case CONDDENSITY:
	  NB_EXPERT = NB_CONDDENSITY (ARGUMENT);
	  if (nb_cpp_handle_error (com))
	    return nb_expert_delta;
	  break;
	case PLOT:
	  if (eventCounter <= 100)
	    {
	      NB_SPLINEF2_WHEN_NEEDED_e ();
	      if (nb_cpp_handle_error (com))
		return nb_expert_delta;
	      for (int nBin = 0; nBin < nBins; ++nBin)
		plotArray[nBin] = 0.;
	      if (NTABL == 0)
		{
		  NB_FTXDEF (mMax);
		  NTABL = 1;
		}
	      NB_EXPERT = 1.;
	      if (writeoutdata)
		{
		  ARGUMENT = 50000;
		  int iPlot = (int) (ARGUMENT) + eventCounter;
		  NBBOOK1 (iPlot, "f(t)  ", 100, TABX[0][0], TABX[0][100], 0.,
			   7);
		  NBPAK_e (iPlot, TABL);
		}
	      for (int nBin = 0; nBin < 100; ++nBin)
		plotArray[nBin] = TABL[nBin];
	    }
	  break;
	case RNDCOND:
	  intArg = int (ARGUMENT + 0.5);
	  randomN = NB_RNDM2 (intArg, nb_get_debug (&com));
	  NB_EXPERT = NB_QUANTILE (randomN);
	  break;
	case RNDINCL:
	  intArg = int (ARGUMENT + 0.5);
	  NB_EXPERT = NB_RNDINCL (intArg);
	  if (nb_cpp_handle_error (com))
	    return nb_expert_delta;
	  break;
	case BINCLASS:
	  NB_EXPERT = (float) NETOUT[0];
	  break;
	case REGR:
	  if (nb_get_debug (&com) >= -2)
	    {
	      ss <<
		"nb_expert called with action REGR: this action is obsolete!"
		<< std::endl;
	      nb_cpp_log (ss, -2, &com);
	    }
	  NB_EXPERT = nb_expert_delta;
	  break;
	case NOOP:
	  break;
	default:
	  if (nb_get_debug (&com) >= -2)
	    {
	      ss << "NB_EXPERT: UNKNOWN ACTION " << std::endl;
	      ss << "please select one of the possible ACTIONS:" << std::endl;
	      ss << " " << std::endl;
	      ss << "incl. quantities independent of input vector X:" << std::
		endl;
	      ss << " " << std::endl;
	      ss << "TMIN  -- minimum value of inclusive distribution" <<
		std::endl;
	      ss << "TMAX  -- maximum value of inclusive distribution" <<
		std::endl;
	      ss << "INCLDENSITY  -- PDF of inclusive distribution" << std::
		endl;
	      ss << "INVQINCL -- incl. prob. for t larger than ARGUMENT" <<
		std::endl;
	      ss << "RNDINCL -- random number acc. to incl. distri." << std::
		endl;
	      ss << " " << std::endl;
	      ss << "quantities depending on input vector X:" << std::endl;
	      ss << " " << std::endl;
	      ss << "MEDIAN -- median " << std::endl;
	      ss << "LERROR -- -1 sigma " << std::endl;
	      ss << "RERROR -- +1 sigma " << std::endl;
	      ss << "QUANTILE -- quantile, ARGUMENT=Quantile 0...1.  " <<
		std::endl;
	      ss << "MEAN    -- mean value " << std::endl;
	      ss << "TRIM    -- trimmed mean value, with ARGUMENT  " << std::
		endl;
	      ss << "CONDDENSITY --conditional density at t. ARGUMENT=t" <<
		std::endl;
	      ss << "PLOT -- plot density in histogram (int)(ARGUMENT)" <<
		std::endl;
	      ss << "RNDCOND -- random number acc. to cond. prob dens." <<
		std::endl;
	      ss << "INVQUANT -- cond. prob. for t larger than ARGUMENT" <<
		std::endl;
	      ss << "BINCLASSification -- binary classification" << std::endl;
	      ss <<
		"PPCLASSification -- binary classification PrePro only for 22FIXSHAPE TOTAL"
		<< std::endl;
	      ss << "REGRession -- regression" << std::endl;
	      nb_cpp_log (ss, -2, &com);
	    }
	  NB_EXPERT = nb_expert_delta;
	}
    }
  else
    {
      if (nb_get_debug (&com) >= -2)
	{
	  ss << "Illegal action " << key << " required for a classification"
	    << std::endl;
	  nb_cpp_log (ss, -2, &com);
	}
      NB_EXPERT = nb_expert_delta;
    }

  if (nb_get_debug (&com) >= 2)
    {
      ss << "At end of EXPERT:" << std::endl;
      ss << "      returning:" << NB_EXPERT << std::endl;
      for (int m = 0; m < NB_NVALUE; ++m)
	{
	  ss << "Expert: TabL(" << m << ")= " << TABL[m] << std::endl;
	  ss << "        TabD(" << m << ")= " << TABD[m] << std::endl;
	}
      nb_cpp_log (ss, 2, &com);
    }

  if (writeoutdata)
    {
      for (int n = 0; n < 101; ++n)
	{
	  if (n < 100)
	    {
	      float TT =
		TABX[0][0] + (TABX[0][100] - TABX[0][0]) * (n - 0.5) / 100.;
	      NBFILL_e (987, TT, 0., TABL[n]);
	      if (nb_cpp_handle_error (com))
		return nb_expert_delta;
	    }
	  float parameter = (float) (n - 1) / 100.;
	  NBFILL_e (986, parameter, 0., TABD[n]);
	  if (nb_cpp_handle_error (com))
	    return nb_expert_delta;
	}
    }

  NewEx = 0;
  return NB_EXPERT;
}


void
Expert::PrintArrays (const char *topic)
{
  ss << "---------------------------" << std::endl;
  ss << "---------------------------" << std::endl;
  ss << "Expert (C++): after        " << topic << std::endl;
  ss << "---------------------------" << std::endl;
  ss << "---------------------------" << std::endl;
  ss << "        NLayer   = " << NLAYER << std::endl;
  ss << "        NODE1    = " << NODE1 << std::endl;
  ss << "        NODES[1] = " << NODES[0] << std::endl;
  ss << "        NODES[2] = " << NODES[1] << std::endl;
  ss << "        NODES[3] = " << NODES[2] << std::endl;
  ss << "        Iter     = " << ITER << std::endl;
  ss << "        IPRUNE   = " << IPRUNE << std::endl;
  ss << "        LShape   = " << LSHAPE << std::endl;
  ss << "        NewEx    = " << NewEx << std::endl;
  ss << "        NumPreproVar= " << NumPreproVar << std::endl;
  ss << "        PREPROC  = " << PREPROC << std::endl;
  ss << "        AutoVarSelect= " << AutoVarSelect << std::endl;

  ss << " Expert(C++)       IN: " << std::endl;
  ss << "---------------------------" << std::endl;
  for (int jj = 0; jj < NB_MAXNODE; ++jj)
    ss << "In[" << jj << "][1] = " << IN[0][jj] << std::endl;

  ss << "  Expert(C++) array A" << std::endl;
  ss << "---------------------------" << std::endl;
  for (int n = 0; n < NB_MAXLAYER; ++n)
    for (int j = 0; j < NB_MAXNODE; ++j)
      ss << "A(" << n << "," << j << ")= " << A[j][n] << std::endl;

  ss << "  Expert(C++)   TABX:" << std::endl;
  ss << "---------------------------" << std::endl;
  for (int j = 0; j < NB_MAXNODE; ++j)
    {
      ss << "TABX for var " << j + 1 << std::endl;
      for (int n = 0; n < NB_NVALUE; ++n)
	ss << TABX[j][n] << "  ";
      ss << std::endl;
    }

  ss << "  Expert(C++)     ITabY:" << std::endl;
  ss << "---------------------------" << std::endl;
  for (int n = 0; n < NB_NIVALUE; ++n)
    ss << ITABY[n] << "  ";
  ss << std::endl;

  ss << "  Expert(C++)     AA:" << std::endl;
  ss << "---------------------------" << std::endl;
  for (int n = 0; n < (NB_MAXNODE - 1) * (NB_MAXNODE - 1); ++n)
    ss << AA[n] << "  ";
  ss << std::endl;

  ss << "  Expert(C++)     DIag" << std::endl;
  ss << "---------------------------" << std::endl;
  for (int n = 0; n < NB_MAXNODE - 1; ++n)
    ss << DIAG[n] << "  ";
  ss << std::endl;

  ss << "  Expert(C++)      CHEBY:" << std::endl;
  ss << "---------------------------" << std::endl;
  for (int n = 0; n < NB_MAXNODE; ++n)
    ss << CHEBY[n] << "  ";
  ss << std::endl;

  ss << "  Expert(C++)  array CSHAPE:" << std::endl;
  ss << "---------------------------" << std::endl;
  for (int n = 0; n < NB_LEVMAX; ++n)
    ss << CSHAPE[n] << "  ";
  ss << std::endl;

  ss << "  Expert(C++)  array outlevel" << std::endl;
  ss << "---------------------------" << std::endl;
  for (int n = 0; n < NB_LEVMAX; ++n)
    ss << OUTLEVEL[n] << "  ";
  ss << std::endl;

  ss << "  Expert(C++)         THETA" << std::endl;
  ss << "---------------------------" << std::endl;
  for (int n = 0; n < NB_MAXNODE * (NB_MAXNODE - 1) / 2; ++n)
    ss << THETA[n] << "  ";
  ss << std::endl;

  ss << "Expert (C++): nMapKey" << std::endl;
  ss << "---------------------------" << std::endl;
  for (int i = 0; i < NB_MAXNODE; ++i)
    ss << nMapKey[i] << "  ";
  ss << std::endl;

  ss << "Expert (C++): MapKeyValue" << std::endl;
  ss << "---------------------------" << std::endl;
  for (int i = 0; i < NB_MAXNODE; ++i)
    {
      ss << "Expert (C++): MapKeyValue for var " << i + 1 << std::endl;
      for (int ii = 0; ii < NB_MaxMapKey; ++ii)
	ss << MapKeyValue[i][ii] << "  ";
      ss << std::endl;
    }

  ss << "Expert (C++): MapKeyTrans" << std::endl;
  ss << "---------------------------" << std::endl;
  for (int i = 0; i < NB_MAXNODE; ++i)
    {
      ss << "Expert (C++): MapKeyTrans for var " << i + 1 << std::endl;
      for (int ii = 0; ii < NB_MaxMapKey; ++ii)
	ss << MapKeyValue[i][ii] << "  ";
      ss << std::endl;
    }
  nb_cpp_log (ss, 2, &com);
}

void
Expert::InfoAndLicence ()
{

  if (nb_get_debug (&com) >= -1)
    {
      ss << std::endl;
	  ss << "  Blue Yonder NeuroBayes(R) Expert" << std::endl;
      ss << "  Algorithms by Michael Feindt" << std::endl;
      ss << "  Implementation by Blue Yonder GmbH" << std::endl;
      ss << "  http://www.blue-yonder.com" << std::endl;
      ss << std::endl;
      nb_cpp_log (ss, -1, &com);
    }


  if (nb_get_debug (&com) >= 0)
    {
      ss << "Version " << IVERSION << std::endl;
      ss << "Library compiled with:" << std::endl;
      ss << "NB_MAXNODE   = " << NB_MAXNODE << std::endl;
      ss << "-----------------------------------" << std::endl;
      nb_cpp_log (ss, 0, &com);
    }
}


void
Expert::NB_EXPERT_FILLTABXS (float *STABXS)
{
  // no call to NB_DEFEXPERTISE since this is done
  // in constructor of Expert class
  NB_DEFFT_e ();
}

void
Expert::NB_EXPERT_DEFGSPLINE (int ModeIn, int NIn, float RegIn)
{
  LDEFSPLINE = 1;
  MODEGS = ModeIn;
  NGS = NIn;
  REGGS = RegIn;
}


void
Expert::NB_EXPERT_GETPINP (float *XPrePro)
{
  for (int i = 0; i < NB_MAXDIM; i++)
    XPrePro[i] = IN[0][i];
}


float
Expert::NB_EXPERT_FTMEAN (double (*f) (float *))
{
  float returnValue = -999.9;
  if (nb_cpp_handle_error (com))
    return returnValue;

  if (TABXS[0] == TABXS[NB_NVALUE - 1])
    {
      ss << "NB_EXPERT must be called before NB_EXPERT_FTMEAN" << std::endl;
      ss << "-> STOP" << std::endl;
      char *msg = stream_c_str (ss);
      if (nb_get_debug (&com) >= -2)
	nb_c_log (msg, -2, &com);
      nb_set_error (&com, invalid_arg_exc, msg);
      free (msg);
      nb_cpp_handle_error (com);
      return returnValue;
    }
  returnValue = NB_FTMEAN (f);

  return returnValue;
}

void
Expert::convertToArray (const char *ExpertiseFileName)
{
  std::string name = ExpertiseFileName;
  nb_converttocarray_e_ (&com, name.c_str (), name.size ());
}

void
Expert::getPDF (float *inputArray, int nBins, float *array)
{
  if (nb_cpp_handle_error (com))
    return;
  // copy input field
  float X[(const int) NODE1 - 1];
  memset (X, 0, (NODE1 - 1) * sizeof (float));	//reset 

  memset (array, 0, nBins * sizeof (float));	//reset 

  for (int jj = 0; jj < NODE1 - 1; ++jj)
    {
      X[jj] = inputArray[jj];
    }

  if (isNewEvent (X) || NewEx == 1)
    {
      NTABL = 0;
      // check the input range
      checkInputRange_e (X);
      if (nb_cpp_handle_error (com))
	return;
      // compute the distribution
      calcNetOutput_e (X);
      if (nb_cpp_handle_error (com))
	return;
    }
  if (NTABL == 1 && nBins == 100)
    {
      // the PDF has already been computed
      array = TABL;
    }
  else
    {
      // compute the PDF
      NB_FTXDEF_e (array, nBins);
      if (nb_cpp_handle_error (com))
	return;
      NTABL = 1;
    }
  NewEx = 0;
}

bool
Expert::getWriteOutData ()
{
  return writeoutdata;
}

void
Expert::setWriteOutData (bool value)
{
  writeoutdata = value;
}

