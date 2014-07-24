
#include "muFile.h"

#include <iostream>

#include <errno.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>

// For directory creation
#if defined(WIN32) && ! defined (__CYGWIN__)
#include <direct.h>
#include <io.h>
#else
#include <unistd.h>
#include <pwd.h>
#endif

namespace mu
{

std::string
get_path(const char* fn)
{

  int n = strlen(fn);

  int pathlen = 0;

  // Get length of path string
  for (int i = (n-1); i >= 0; i--)
  {
    if (fn[i] == MU_DIR_SEPARATOR)
    {
      pathlen = i;
      break;
    }
  }

  if (pathlen <= 0)
    return std::string(".");

  char* path = new char[pathlen+2];

  for (int j = 0; j < pathlen; j++)
    path[j] = fn[j];

  path[pathlen] = MU_DIR_SEPARATOR;
  path[pathlen+1] = 0;

  std::string retstr(path);
  delete [] path;

  return retstr;

}

std::string
get_name(const char* fn)
{

  int n = strlen(fn);

  // Skip file extension
  int nameend = n-1;
  for (int k = (n-1); k >= 0; k--)
  {
    if (fn[k] == '.')
    {
      nameend = k-1;
      break;
    }
    if (fn[k] == MU_DIR_SEPARATOR)
      break;
  }

  int offt = 0;

  // Get offset
  for (int i = nameend; i >= 0; i--)
  {
    if (fn[i] == MU_DIR_SEPARATOR)
    {
      offt = i+1;
      break;
    }
  }

  int namelen = (nameend - offt) + 1;

  if (namelen <= 0)
    return std::string("");

  char* name = new char[namelen+1];

  for (int j = 0; j < namelen; j++)
    name[j] = fn[j+offt];

  name[namelen] = 0;

  std::string retstr(name);
  delete [] name;

  return retstr;

}

std::string
get_ext(const char* fn)
{

  int n = strlen(fn);

  int offt = n+1;

  for (int i = (n-1); i >= 0; i--)
  {
    if (fn[i] == '.')
    {
      offt = i+1;
      break;
    }
    if (fn[i] == MU_DIR_SEPARATOR)
      break;
  }

  int extlen = (n - offt);

  if (extlen <= 0)
    return std::string("");

  char* ext = new char[extlen+1];

  for (int j = 0; j < extlen; j++)
    ext[j] = fn[j+offt];

  ext[extlen] = 0;

  std::string retstr(ext);
  delete [] ext;

  return retstr;

}

bool
create_dir(const char* dir)
{

  // Create directory, ignore EEXIST errors...
#if defined(WIN32) && ! defined (__CYGWIN__)
  if (mkdir(dir))
#else
  if (mkdir(dir, 0777))
#endif
    if (errno != EEXIST)
    {
      std::cerr << "Cannot create directory: " << dir << std::endl;
      return false;
    }

  return true;

}


char**
list_dir(const char* dir)
{

//TODO

  return 0;

}

} // namespace mu

// Test, for UNIX style file names
#if 0

#include <assert.h>

int main(int argc, char** argv)
{

  char* s1 = "/abc/def/ghi";
  char* s2 = "/abc/def/ghi.jkl";
  char* s3 = "abc";
  char* s4 = "abc.def";
  char* s5 = "/abc/def/";

  std::cout << "s2 = " << s2 << std::endl;
  std::cout << "get_path(s2) = " << mu::get_path(s2) << std::endl;
  std::cout << "get_name(s2) = " << mu::get_name(s2) << std::endl;
  std::cout << "get_ext(s2) = " << mu::get_ext(s2) << std::endl;

  std::cout << "s4 = " << s4 << std::endl;
  std::cout << "get_path(s4) = " << mu::get_path(s4) << std::endl;
  std::cout << "get_name(s4) = " << mu::get_name(s4) << std::endl;
  std::cout << "get_ext(s4) = " << mu::get_ext(s4) << std::endl;

  assert(strcmp(mu::get_path(s1).c_str(),"/abc/def/") == 0);
  assert(strcmp(mu::get_name(s1).c_str(),"ghi") == 0);
  assert(strcmp(mu::get_ext(s1).c_str(),"") == 0);

  assert(strcmp(mu::get_path(s2).c_str(),"/abc/def/") == 0);
  assert(strcmp(mu::get_name(s2).c_str(),"ghi") == 0);
  assert(strcmp(mu::get_ext(s2).c_str(),"jkl") == 0);

  assert(strcmp(mu::get_path(s3).c_str(),".") == 0);
  assert(strcmp(mu::get_name(s3).c_str(),"abc") == 0);
  assert(strcmp(mu::get_ext(s3).c_str(),"") == 0);

  assert(strcmp(mu::get_path(s4).c_str(),".") == 0);
  assert(strcmp(mu::get_name(s4).c_str(),"abc") == 0);
  assert(strcmp(mu::get_ext(s4).c_str(),"def") == 0);

  assert(strcmp(mu::get_path(s5).c_str(),"/abc/def/") == 0);
  assert(strcmp(mu::get_name(s5).c_str(),"") == 0);
  assert(strcmp(mu::get_ext(s5).c_str(),"") == 0);

  return 0;

}

#endif // Test
