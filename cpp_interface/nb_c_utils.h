/* == Neurobayes c utils == */
/* by Martin Hahn, neurobayes@blue-yonder.com */

#ifndef NB_C_UTILS_H
#define NB_C_UTILS_H

#include <stdlib.h>
#include <math.h>
#ifndef __cplusplus
#include <stdbool.h>
#endif

#ifdef __cplusplus
extern "C"
{
#endif

  typedef struct ec_t
  {
    int kind;
    char *reason;
    char *stacktrace;
  } ec_t;

  typedef double (*rand_func_t) (void *rand_enc);
  typedef void (*log_func_t) (char *msg, int debug_lvl, void *enc);
  typedef void (*delete_enclosed_func_t) (void *enc);
  typedef struct log_t
  {
    int debug_lvl;
    log_func_t f;
    void *enclosed;
    delete_enclosed_func_t delete_enclosed;
  } log_t;

  typedef struct common_t
  {
    ec_t **ec;
    ec_t *ecf;
    log_t *log;
    rand_func_t rand;
    void *rand_enc;
    int is_teacher;
  } common_t;

/* == Common  == */

  common_t *nb_init_common (int debug_lvl, ec_t ** ec1,
			    rand_func_t rand_func, void *rand_enc,
			    int is_teacher);

/* Fortran functions */
  void nb_delete_common (common_t ** com);
  void nb_delete_common_and_ecf (common_t ** com);

/* == Error Codes == */
  int nb_get_error (common_t ** com);
  void nb_set_error (common_t ** com, int input_error_code, char *msg);
  void nb_set_error1 (common_t * com, int input_error_code, char *msg);
  int nb_c_handle_error (common_t ** com);
  int nb_get_is_teacher (common_t ** com);
  void free_ec1 (ec_t * ec);

/* == Random Generator == */
  float nb_rand_double_f (common_t ** com1);

  enum
  { invalid_arg_exc = 1, file_not_found_exc, file_open_exc,
    null_pointer_exc, not_found_exc, duplication_exc, errno_exc,
    assertion_exc, division_by_zero_exc, test_skipped_exc, bounds_check_exc,
    not_implemented_exc, unix_command_exc, invalid_data_exc, dynv_exc,
    gen_exc, sql_open_exc, sql_query_failed_exc, algo_exc,
    out_of_memory_exc
  };

/* == Logging == */
  typedef enum
  { quiet_dbg = -3, error_dbg = -2, normal_dbg = -1,
    verbose_dbg = 0, vverbose_dbg = 1
  } debug_t;

  void nb_c_log (char *msg, int debug_lvl, common_t ** com);
  void nb_c_logging_delete_enclosed (common_t ** com);
  void nb_register_logging (common_t ** com, log_func_t log,
			    void *enc,
			    delete_enclosed_func_t delete_enclosed);
  void nb_c_logging_set_default_logging (common_t ** com);

/* Fortran functions */
  int nb_get_debug (common_t ** com);


/*
 * For a description of __thread see
 * https://gcc.gnu.org/onlinedocs/gcc-4.2.4/gcc/Thread_002dLocal.html Usage
 * here is probably to avoid implementing a mutex. When several threads
 * initialize the random number generator at the same time, they could
 * overwrite an initialization from another thread.
 */
extern __thread struct random_data my_random_data;
extern __thread bool my_random_state_initialized;
extern __thread double my_random_state;

#define DEFAULT_RANDOM_SEED 314159

static inline void
init_random_number_generator (unsigned int seed)
{
  initstate_r (seed, (char *) &my_random_state, sizeof (my_random_state),
	       &my_random_data);
  my_random_state_initialized = true;
}

static inline long
random_threadsafe ()
{
  int32_t result;
  if (! my_random_state_initialized)
    {
      init_random_number_generator (DEFAULT_RANDOM_SEED);
    }
  random_r (&my_random_data, &result);
  return result;
}

static inline double
rand_double ()
{
  return random_threadsafe () / ((double) RAND_MAX + 1);
}

// Dummy argument <enc> for neurobayes random generator
static inline double
rand_double1 (void *enc)
{
  return random_threadsafe () / ((double) RAND_MAX + 1);
}

#ifdef __cplusplus
}				// extern "C"
#endif

#endif
