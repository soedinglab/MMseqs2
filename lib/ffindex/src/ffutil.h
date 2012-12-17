/*
 * Ffindex
 * written by Andy Hauser <hauser@genzentrum.lmu.de>.
 * Please add your name here if you distribute modified versions.
 * 
 * Ffindex is provided under the Create Commons license "Attribution-ShareAlike
 * 3.0", which basically captures the spirit of the Gnu Public License (GPL).
 * 
 * See:
 * http://creativecommons.org/licenses/by-sa/3.0/
 * 
 * Ffindex is a very simple database for small files. The files are stored
 * concatenated in one big data file, seperated by '\0'. A second file
 * contains a plain text index, giving name, offset and length of of the small
 * files.
 */

#include <err.h>
#include <errno.h>
#include <string.h>


int fferror_print(char *sourcecode_filename, int line, const char *function_name, const char *message);

char* ffnchomp(char *s, size_t len);

/* vim: ts=2 sw=2 et
*/
