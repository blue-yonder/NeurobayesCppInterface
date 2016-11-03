/* == Neurobayes cpp utils == */
/* by Martin Hahn, neurobayes@blue-yonder.com */
#ifndef NB_CPP_UTILS_H
#define NB_CPP_UTILS_H

#include "nb_c_utils.h"
#include <iostream>
#include <sstream>
#include <string>

/* == Logging == */
void nb_cpp_log (std::stringstream & ss, int debug_lvl, common_t ** com);
char *stream_c_str (std::stringstream & ss);
/* == Error Codes == */
bool nb_cpp_handle_error (common_t * com);

#endif
