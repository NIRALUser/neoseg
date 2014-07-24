
////////////////////////////////////////////////////////////////////////////////
//
// Functions for extracting information from a given file name string
//
////////////////////////////////////////////////////////////////////////////////

#ifndef _muFile_h
#define _muFile_h

#include <string>

#if defined(WIN32) && !defined(__CYGWIN__)
// FLTK on Windows  uses forward slash
//#define MU_DIR_SEPARATOR '\\'
#define MU_DIR_SEPARATOR '/'
#else
#define MU_DIR_SEPARATOR '/'
#endif

namespace mu
{

std::string get_path(const char* fn);
std::string get_name(const char* fn);
std::string get_ext(const char* fn);

bool create_dir(const char* dir);
char** list_dir(const char* dir);

}

#endif
