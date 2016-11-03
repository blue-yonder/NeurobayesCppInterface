/* == Neurobayes cpp utils == */
/* by Martin Hahn, neurobayes@blue-yonder.com */

#include "nb_cpp_utils.h"
#include <string.h>
#include <stdlib.h>

/* == Logging == */
void
nb_cpp_log (std::stringstream & ss, int debug_lvl, common_t ** com)
{
  nb_c_log ((char *) (ss.str ().c_str ()), debug_lvl, com);
  ss.str ("");
}

char *
stream_c_str (std::stringstream & ss)
{
  char *msg = strdup (ss.str ().c_str ());
  ss.str ("");
  return msg;
}

/* == Error Codes == */
bool
nb_cpp_handle_error (common_t * com)
{
  if (com->ec && *(com->ec))
    {				// error occured
      if (com->ecf)		// no error code wanted
	exit ((*(com->ec))->kind);
      return true;
    }
  else
    return false;
}
